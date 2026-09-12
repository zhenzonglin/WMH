"""Read-only SAS/CSV census before image access or patient modelling."""
from __future__ import annotations

import html
from collections import Counter

import numpy as np
import pandas as pd
import pyreadstat

from .common import DataError, dump_json, outdir, read_csv, read_json, resolve
from .fields import FIELDS, UNKNOWN_98
from .harmonize import sas_date
from .sas_extract import inventory, iter_sas_chunks, resolve_owners

BASE = ["code_n", "AGE", "GENDER", "D_DIAG", "BSL_HCY", "BSL_B12", "BSL_B9", "BSL_CYSC",
        "H_SMK", "H_DRINK", "H_HYPT", "H_DIAB", "H_STROKE", "ONSET_D", "I_BLDSAMP_DT", "IMG_ONSET_TO_MRI_D"]
EVENT = ["y1_is", "y1_is_dd"]
M3 = ["F3_BLDSAMP_D", "M03_HCY", "M03_B12", "M03_B9", "M03_CYSC"]
FUNCTION = ["F12_MRS", "H_MRS", "A_NIHSS", "IMG_C_TOAST"]
PROFILES = {"H1": BASE, "H2": BASE + EVENT, "H3": BASE + EVENT + M3, "H4": BASE + FUNCTION}


def inspect_values(series: pd.Series, source: str, fmt: str, overrides: dict) -> tuple[pd.Series, pd.Series]:
    """Return usable and invalid masks, without changing or imputing source values."""
    raw = series.astype("string").fillna("").str.strip()
    name, kind, codes = FIELDS[source]
    missing = raw.eq("") | raw.str.match(r"^\.?[A-Z_]$|^\.$", na=False)
    if kind == "id":
        good = raw.ne("")
        return good, ~good
    if source in UNKNOWN_98:
        missing |= raw.isin(["98", "98.0"])
    if kind in {"date", "datetime"}:
        clean = raw.mask(missing, "")
        try:
            parsed = sas_date(clean.rename(source), fmt, overrides.get(source))
            return parsed.notna(), pd.Series(False, index=series.index)
        except DataError:
            # The report requests an explicit format/coding correction before prepare.
            return pd.Series(False, index=series.index), ~missing
    numeric = pd.to_numeric(raw.mask(missing), errors="coerce")
    good = numeric.notna() & np.isfinite(numeric)
    if codes:
        good &= numeric.isin(codes)
    if name in {"hcy", "b12", "folate", "cysc", "creatinine", "hcy3", "b123", "folate3", "cysc3"}:
        good &= numeric.gt(0)
    if name in {"pre_mrs", "mrs12", "nihss"}:
        good &= numeric.between(0, 42 if name == "nihss" else 5) & numeric.mod(1).eq(0)
    return good.fillna(False), (~missing & ~good).fillna(False)


def clinical_audit(cfg: dict) -> dict:
    folder = outdir(cfg) / "audit"
    folder.mkdir(parents=True, exist_ok=True)
    read_errors: list[dict] = []
    inv = inventory(cfg, errors=read_errors)
    owner, source_issues = resolve_owners(inv, cfg)
    csv_source = cfg["inputs"].get("clinical_csv")
    if csv_source:
        path = resolve(cfg, csv_source)
        raw = read_csv(path)
        spelling = {c.lower(): c for c in raw.columns}
        mapping = {c: spelling[c.lower()] for c in FIELDS if c.lower() in spelling}
        fp = cfg["inputs"].get("clinical_formats_json")
        formats = read_json(resolve(cfg, fp)) if fp else {}
        inv = [{"path": str(path), "basename": path.name, "rows": len(raw), "matched": list(mapping),
                "field_columns": mapping, "columns": list(raw), "formats": formats, "csv": True}]
        owner = {c: str(path) for c in mapping if c != "code_n"}
        source_issues, read_errors = [], []
    files, variables, frames, id_sets = [], [], [], {}
    global_invalid = Counter()
    overrides = cfg["sas"].get("date_formats", {})
    merged_formats = {}
    for rec in inv:
        selected = rec["matched"]
        if not selected:
            files.append({"file": rec["path"], "rows_metadata": rec["rows"], "rows_scanned": None,
                          "whitelist_fields": 0, "status": "NO_WHITELIST_FIELDS"})
            continue
        tally = {c: Counter() for c in selected}
        parts, ids, rows = [], [], 0
        try:
            if rec.get("csv"):
                pairs = [(read_csv(resolve(cfg, csv_source)).rename(columns={v: k for k, v in rec["field_columns"].items()}), None)]
            else:
                pairs = iter_sas_chunks(rec["path"], selected, cfg["sas"], rec["field_columns"])
            for chunk, _ in pairs:
                rows += len(chunk)
                for col in selected:
                    fmt = rec["formats"].get(rec["field_columns"][col], rec["formats"].get(col, ""))
                    good, bad = inspect_values(chunk[col], col, fmt, overrides)
                    tally[col].update(rows=len(chunk), usable=int(good.sum()), invalid=int(bad.sum()))
                    if owner.get(col) == rec["path"]:
                        merged_formats[col] = fmt
                        global_invalid[col] += int(bad.sum())
                if "code_n" in chunk:
                    ids.extend(chunk.code_n.astype("string").fillna("").tolist())
                    owned = [c for c in selected if owner.get(c) == rec["path"]]
                    if owned:
                        parts.append(chunk[["code_n", *owned]].copy())
        except (DataError, OSError, ValueError, pyreadstat.ReadstatError) as exc:
            read_errors.append({"path": rec["path"], "reason": str(exc)})
            files.append({"file": rec["path"], "rows_scanned": rows, "status": "READ_ERROR"})
            continue
        nonempty = [x for x in ids if str(x).strip()]
        counts = Counter(nonempty)
        duplicate_rows = len(nonempty)-len(counts)
        missing_ids = len(ids)-len(nonempty)
        ids_valid = "code_n" in selected and duplicate_rows == 0 and missing_ids == 0
        id_sets[rec["path"]] = set(nonempty)
        files.append({"file": rec["path"], "rows_metadata": rec["rows"], "rows_scanned": rows,
                      "unique_nonempty_ids": len(counts), "duplicate_rows": duplicate_rows,
                      "duplicate_id_count": sum(v > 1 for v in counts.values()), "empty_ids": missing_ids,
                      "has_code_n": "code_n" in selected, "whitelist_fields": len(selected),
                      "status": "ID_OK" if ids_valid else "ID_REVIEW_REQUIRED"})
        if parts and ids_valid:
            frames.append(pd.concat(parts, ignore_index=True))
        elif parts or any(v == rec["path"] for v in owner.values()):
            source_issues.append({"file": rec["path"], "reason": "Selected source lacks unique nonempty code_n"})
        for c, t in tally.items():
            variables.append({"file": rec["path"], "variable": c, "actual_column": rec["field_columns"][c],
                              "selected_source": owner.get(c) == rec["path"] or c == "code_n",
                              "rows": t["rows"], "usable": t["usable"], "invalid": t["invalid"],
                              "missing": t["rows"]-t["usable"]-t["invalid"],
                              "format": rec["formats"].get(rec["field_columns"][c], "")})
    merged = pd.DataFrame(columns=["code_n"])
    for part in frames:
        part["code_n"] = part.code_n.astype(str)
        merged = part if merged.empty else merged.merge(part, on="code_n", how="outer", validate="one_to_one")
    usable = {}
    for c in FIELDS:
        usable[c] = inspect_values(merged[c], c, merged_formats.get(c, ""), overrides)[0] if c in merged \
            else pd.Series(False, index=merged.index)
    present = {c for rec in inv for c in rec["matched"]}
    coverage = [{"variable": c, "analysis_name": FIELDS[c][0], "present_anywhere": c in present,
                 "selected_source": owner.get(c, "ID in each selected file" if c == "code_n" else ""),
                 "n_usable_selected_union": int(usable[c].sum()),
                 "expected_unit_or_type": FIELDS[c][1]} for c in FIELDS]
    profiles = []
    for h, columns in PROFILES.items():
        required = [c for c in columns if c != "y1_is_dd"]
        complete = pd.concat([usable[c] for c in required], axis=1).all(axis=1)
        if "y1_is" in columns and "y1_is" in merged:
            ev = pd.to_numeric(merged.y1_is, errors="coerce")
            times = pd.to_numeric(merged.get("y1_is_dd", pd.Series(np.nan, index=merged.index)), errors="coerce")
            complete &= ev.eq(0) | (ev.eq(1) & times.between(0, 365))
        profiles.append({"hypothesis": h, "missing_fields": ";".join(sorted(set(columns)-present)),
                         "n_observed_all_listed_values": int(complete.sum()),
                         "interpretation": "Observed clinical overlap before MRI/QC/time eligibility; not final sample size"})
    comparisons = []
    for i, (a, aid) in enumerate(id_sets.items()):
        for b, bid in list(id_sets.items())[i+1:]:
            comparisons.append({"file_a": a, "file_b": b, "shared_unique_ids": len(aid & bid)})
    required_missing = sorted(set(PROFILES["H2"])-present)
    # Event times can be entirely empty when no patient has an event; other core inputs cannot.
    required_empty = [c for c in BASE + ["y1_is"] if c in present and not usable[c].any()]
    tracking_pairs = [(f"F{v}_DATE", f"F{v}_DEATH") for v in [3, 6, 12]]
    tracking_present = any(set(pair) <= present for pair in tracking_pairs)
    followup_issue = [] if tracking_present else ["No visit date + survival-status field pair supplied"]
    invalid_core = {c: n for c, n in global_invalid.items() if n and c in PROFILES["H2"]}
    blocking = bool(required_missing or required_empty or read_errors or source_issues or followup_issue or invalid_core or merged.empty)
    result = {"status": "REVIEW_REQUIRED" if blocking else "READY_FOR_EXTRACTION", "mode": cfg["mode"],
              "sas_or_csv_files": len(inv), "patients_in_selected_source_union": len(merged),
              "unique_ids_in_any_scanned_file": len(set().union(*id_sets.values())) if id_sets else 0,
              "whitelist_fields_present": len(present), "whitelist_fields_expected": len(FIELDS),
              "required_missing": required_missing, "required_all_missing_or_invalid": required_empty,
              "invalid_core_values": invalid_core, "source_issues": source_issues, "read_errors": read_errors,
              "followup_issues": followup_issue, "profiles": profiles,
              "definition": "Clinical field/value census. Cohort eligibility is calculated by prepare after imaging joins."}
    dump_json(folder / "summary.json", result)
    dump_json(folder / "sas_inventory.json", inv)
    dump_json(folder / "variable_sources.json", owner)
    for name, rows in [("files", files), ("variables_by_file", variables), ("field_coverage", coverage),
                       ("hypothesis_counts", profiles), ("file_intersections", comparisons)]:
        pd.DataFrame(rows).to_csv(folder / f"{name}.csv", index=False, encoding="utf-8-sig")
    text = f"# 临床数据审计\n\n状态：{result['status']}\n\n"
    text += f"- 来源文件：{len(inv)}\n- 合并患者数（已选来源并集）：{len(merged)}\n- 字段存在：{len(present)}/{len(FIELDS)}\n"
    text += "\n## 按假设的原始临床完整值交集\n\n"
    text += "\n".join(f"- {r['hypothesis']}：{r['n_observed_all_listed_values']}；缺字段：{r['missing_fields'] or '无'}" for r in profiles)
    text += "\n\n这些是影像连接前的观察值交集，不是最终样本量；协变量散在缺失仍按方案处理。\n\n详细问题见summary.json；逐文件行数、重复ID和缺失见files.csv、variables_by_file.csv。\n"
    (folder / "SUMMARY.md").write_text(text, encoding="utf-8")
    tables = "".join(f"<h2>{title}</h2>{pd.DataFrame(rows).to_html(index=False, escape=True)}" for title, rows in
                     [("文件与ID", files), ("字段覆盖与有效值", coverage), ("假设交集", profiles)])
    (folder / "report.html").write_text('<meta charset="utf-8"><title>WMH clinical audit</title>'
        '<style>body{font:16px sans-serif;margin:30px}table{border-collapse:collapse}td,th{padding:7px;border:1px solid #ddd}</style>'
        f"<h1>临床数据审计：{html.escape(result['status'])}</h1><pre>{html.escape(text)}</pre>"+tables, encoding="utf-8")
    return result
