"""Read-only screenshot pages: saved field provenance and aggregate value checks."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..common import DataError, read_json, resolve, sha256, unique_ids
from .registry import FOLDERS

FOCUS = {"bp": ("IMG_ICAS", "H_CHD", "BSL_CYSC"),
         "ceramide": ("BSL_Cer_16_0", "BSL_Cer_24_0", "BSL_Cer_24_1"),
         "cec": ("CEC",), "kidney": ("M03_CYSC",)}
PATTERNS = {"ICAS": r"icas|颅内.*(?:狭窄|动脉)|intracranial",
            "CHD": r"chd|冠心病|冠状动脉|coronary",
            "CERAMIDE": r"ceramide|神经酰胺|cer.*(?:16|24)|c(?:16|24)[_: .]*[01]",
            "CEC": r"(?:^|_)cec(?:$|_)|efflux|外排|胆固醇流出"}


def saved_runs(cfg):
    runs = {}
    for study, folder in FOLDERS.items():
        pointer = read_json(Path(cfg["_out"]) / "studies" / folder / "latest_prepared.json")
        if pointer:
            runs[study] = Path(pointer["path"])
    if not runs:
        raise DataError("No saved study audits. Run wmh-study audit --study all first.")
    return runs


def field_page(runs):
    lines = ["[1/3] SAVED FIELD PRESENCE (counts refer to all extracted clinical rows, not eligible N)"]
    for study, sources in FOCUS.items():
        root = runs.get(study)
        if root is None or not (root / "field_audit.csv").is_file():
            lines.append(f"{study}: field audit unavailable")
            continue
        fields = pd.read_csv(root / "field_audit.csv", keep_default_na=False)
        state = read_json(root / "status.json")
        a = state.get("audit", {})
        lines.append(f"{study}: run={root.name}; status={state.get('status')}; eligible={a.get('eligible_n')}")
        for source in sources:
            match = fields.loc[fields.source.eq(source)]
            if match.empty:
                lines.append(f"  {source}: NOT_IN_SAVED_AUDIT")
                continue
            r = match.iloc[0]
            present = str(r.present).lower() == "true"
            flag = "ABSENT_COLUMN" if not present else "PRESENT_NO_VALID_VALUES" if int(r.observed) == 0 else "PRESENT"
            file = Path(str(r.get("file", ""))).name if r.get("file", "") else "-"
            lines.append(f"  {source}: {flag}; valid={r.observed}/{r.patients}; invalid={r.invalid}; file={file}")
        if a.get("eligible_n") == 0:
            lines.append("  Covariate missingness NOT ASSESSED: empty cohort does not mean all source covariates are missing.")
    lines.append("ABSENT_COLUMN means absent from this extraction; candidate metadata names are on page 2.")
    return lines


def metadata_records(runs):
    records = {}
    for root in runs.values():
        inventory = read_json(root / "source_inventory.json")
        if isinstance(inventory, list):
            for rec in inventory:
                # Preserve different saved file versions; do not silently select an owner.
                key = (rec["path"], rec.get("mtime_ns"), rec.get("bytes"))
                records[key] = rec
    return sorted(records.values(), key=lambda r: (r["path"], r.get("mtime_ns", 0)))


def candidate_page(runs):
    records = metadata_records(runs)
    lines = ["[2/3] SAVED SAS METADATA SEARCH (names/labels only; candidates are NOT mapped)",
             f"Saved file versions={len(records)}; current/new SAS files are not rescanned."]
    for rec in records:
        lines.append(f"  FILE {Path(rec['path']).name}: columns={len(rec.get('columns', []))}")
    for topic, pattern in PATTERNS.items():
        found = []
        for rec in records:
            for name in rec.get("columns", []):
                label = str((rec.get("labels") or {}).get(name) or "")
                if re.search(pattern, name, re.IGNORECASE) or re.search(pattern, label, re.IGNORECASE):
                    found.append((Path(rec["path"]).name, name, " ".join(label.split())))
        lines.append(f"{topic}: candidates={len(found)} (first 8)")
        lines.extend(f"  {file} | {name} | {label[:90]}" for file, name, label in found[:8])
        if not found:
            lines.append("  NO MATCH in saved metadata; this does not prove the assay was never measured.")
        elif len(found) > 8:
            lines.append(f"  {len(found)-8} further matches retained in local source_inventory.json.")
    if not records:
        lines.append("No saved SAS inventory (for example CSV input); no fuzzy replacement is attempted.")
    return lines


def raw_columns(root, names):
    path = root / "extracted/clinical_raw.csv"
    if not path.is_file():
        cfg = read_json(root / "config_snapshot.json")
        source = cfg.get("inputs", {}).get("clinical_csv")
        if not source:
            raise DataError("Saved extracted CSV unavailable")
        path = resolve(cfg, source)
        manifest = read_json(root / "source_files.json")
        if not isinstance(manifest, list) or not any(r.get("sha256") == sha256(path) for r in manifest):
            raise DataError("Clinical CSV changed since saved audit; no values inspected")
    header = pd.read_csv(path, nrows=0)
    spellings = {c.casefold(): c for c in header}
    if len(spellings) != len(header.columns):
        raise DataError("Case-insensitive duplicate columns")
    mapping = {spellings[n.casefold()]: n for n in ["code_n", *names] if n.casefold() in spellings}
    if "code_n" not in mapping.values():
        raise DataError("Saved raw CSV has no exact code_n")
    data = pd.read_csv(path, usecols=list(mapping), dtype=str, keep_default_na=False).rename(columns=mapping)
    unique_ids(data, "code_n", "diagnostic input")
    return data


def numeric_parts(values):
    text = values.astype("string").str.strip()
    missing = text.str.match(r"^\.?[A-Z_]$|^\.$|^$", na=False)
    numeric = pd.to_numeric(text.mask(missing), errors="coerce").astype(float)
    finite = numeric.notna() & np.isfinite(numeric)
    masks = {"sas_or_blank_missing": missing, "unparseable": ~missing & numeric.isna(),
             "nonfinite": numeric.notna() & ~np.isfinite(numeric),
             "zero": finite & numeric.eq(0), "negative": finite & numeric.lt(0), "positive": finite & numeric.gt(0)}
    return numeric, masks


def value_page(runs):
    lines = ["[3/3] LOCAL VALUE CHECK (aggregate counts only; no IDs or individual laboratory values)"]
    if "bp" in runs:
        root = runs["bp"]
        lines.append(f"bp run={root.name}")
        try:
            raw = raw_columns(root, ["H_CHD", "IMG_ICAS"])
            eligible = pd.read_csv(root / "eligible.csv", usecols=["patient_id"], dtype=str)
            selected = raw.loc[raw.code_n.isin(eligible.patient_id)]
            for source, codes in (("H_CHD", [0, 1]), ("IMG_ICAS", [1, 2, 3])):
                if source not in raw:
                    lines.append(f"  {source}: ABSENT_COLUMN")
                    continue
                for scope, d in (("all_extracted", raw), ("eligible", selected)):
                    value, parts = numeric_parts(d[source])
                    counts = "; ".join(f"code{c}={int(value.eq(c).sum())}" for c in codes)
                    other = int((~parts['sas_or_blank_missing'] & ~value.isin(codes)).sum())
                    lines.append(f"  {source} {scope} N={len(d)}; {counts}; missing={int(parts['sas_or_blank_missing'].sum())}; other={other}")
            lines.append("  H_CHD blanks are not assigned 0; IMG_ICAS code3 is dictionary unknown, not normal.")
        except (DataError, OSError, ValueError) as exc:
            lines.append(f"  BP value check unavailable: {type(exc).__name__} (check local saved inputs)")
    if "kidney" in runs:
        root = runs["kidney"]
        lines.append(f"kidney run={root.name}")
        try:
            raw = raw_columns(root, ["M03_CYSC"])
            if "M03_CYSC" not in raw:
                raise DataError("M03_CYSC absent")
            eligible = pd.read_csv(root / "eligible.csv", usecols=["patient_id", "state60"], dtype={"patient_id": str})
            _, parts = numeric_parts(raw.M03_CYSC)
            lines.append("  M03_CYSC all_extracted: "+"; ".join(f"{k}={int(v.sum())}" for k, v in parts.items()))
            invalid = parts["unparseable"] | parts["nonfinite"] | parts["zero"] | parts["negative"]
            for label, ids in (("eligible_before_outcome", eligible.patient_id),
                               ("observed_five_year_outcome", eligible.loc[eligible.state60.notna(), "patient_id"])):
                selected = raw.code_n.isin(ids)
                lines.append(f"  {label}: raw_rows={int(selected.sum())}; invalid_M03_CYSC={int((invalid & selected).sum())}")
            lines.append("  Invalid values were masked by existing preparation; counts locate them, without approving imputation.")
        except (DataError, OSError, ValueError) as exc:
            lines.append(f"  Kidney value check unavailable: {type(exc).__name__} (check local saved inputs)")
    lines.append("Read-only: no recoding, cohort exclusion, fitting, source scan or result-pointer changes.")
    return lines


def diagnose(cfg, page=1):
    runs = saved_runs(cfg)
    pages = {1: field_page, 2: candidate_page, 3: value_page}
    return [f"WMH STUDY DIAGNOSTICS | mode={cfg['mode']} | page={page}",
            "READ ONLY. Saved audit snapshots; candidate names require definition/units/coding confirmation.",
            *pages[page](runs)]
