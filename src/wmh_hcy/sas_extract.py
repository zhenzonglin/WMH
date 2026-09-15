"""Read-only SAS discovery and chunked whitelist extraction."""
from __future__ import annotations

import glob
from pathlib import Path

import pandas as pd
import pyreadstat

from .common import DataError, dump_json, outdir, read_csv, record_run, resolve, unique_ids
from .fields import FIELDS


def iter_sas_chunks(path, selected, options, column_map=None):
    """SAS names are case-insensitive; record actual spellings and preserve ID strings."""
    column_map = column_map or {c: c for c in selected}
    actual = [column_map[c] for c in selected]
    for chunk, meta in pyreadstat.read_file_in_chunks(
        pyreadstat.read_sas7bdat, str(path), chunksize=options["chunksize"],
        usecols=actual, disable_datetime_conversion=True, user_missing=True,
        encoding=options.get("encoding"),
    ):
        chunk = chunk.rename(columns={v: k for k, v in column_map.items()})
        if "code_n" in chunk:
            ids = chunk["code_n"]
            stored_type = (getattr(meta, "readstat_variable_types", None) or {}).get(column_map["code_n"])
            if pd.api.types.is_numeric_dtype(ids) or stored_type == "double":
                ids = pd.to_numeric(ids, errors="coerce")
                if ((ids.dropna() % 1) != 0).any():
                    raise DataError("Numeric code_n contains noninteger IDs")
                chunk["code_n"] = ids.astype("Int64").astype("string").fillna("")
                if options.get("id_width"):
                    good = chunk.code_n.ne("")
                    chunk.loc[good, "code_n"] = chunk.loc[good, "code_n"].str.zfill(options["id_width"])
        yield chunk, meta


def sas_to_csv(path, target, selected, options, column_map=None):
    """Preserve raw dates, SAS special missing codes, and text in UTF-8 CSV."""
    count, missing_meta = 0, {}
    for chunk, meta in iter_sas_chunks(path, selected, options, column_map):
        chunk.to_csv(target, index=False, mode="w" if count == 0 else "a",
                     header=count == 0, encoding="utf-8")
        count += len(chunk)
        missing_meta.update(meta.missing_user_values or {})
    if count == 0:
        pd.DataFrame(columns=selected).to_csv(target, index=False)
    return count, missing_meta


def inventory(cfg: dict, errors: list | None = None) -> list[dict]:
    paths = sorted({Path(p).resolve() for pattern in cfg["inputs"]["sas_globs"]
                    for p in glob.glob(str(resolve(cfg, pattern)), recursive=True)
                    if Path(p).is_file() and Path(p).suffix.lower() == ".sas7bdat"})
    records = []
    for path in paths:
        try:
            _, meta = pyreadstat.read_sas7bdat(str(path), metadataonly=True,
                                         encoding=cfg["sas"].get("encoding"))
        except (OSError, ValueError, pyreadstat.ReadstatError) as exc:
            if errors is None:
                raise DataError(f"{path.name}: SAS metadata read failed: {exc}") from exc
            errors.append({"path": str(path), "reason": str(exc)})
            continue
        spelling = {c.lower(): c for c in meta.column_names}
        field_columns = {c: spelling[c.lower()] for c in FIELDS if c.lower() in spelling}
        records.append({
            "path": str(path), "basename": path.name, "bytes": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns, "rows": meta.number_rows,
            "columns": meta.column_names, "labels": meta.column_names_to_labels,
            "formats": meta.original_variable_types, "encoding": meta.file_encoding,
            "matched": list(field_columns), "field_columns": field_columns,
        })
    return records


def resolve_owners(inv: list[dict], cfg: dict) -> tuple[dict, list[dict]]:
    owner, issues = {}, []
    for source in FIELDS:
        if source == "code_n":
            continue
        choices = [r for r in inv if source in r["matched"]]
        preference = cfg.get("variable_sources", {}).get(source)
        if preference:
            choices = [r for r in choices if preference in [r["basename"], r["path"]]]
        if len(choices) > 1 or (preference and len(choices) != 1):
            issues.append({"variable": source, "reason": "Select one source in variable_sources",
                           "candidates": [r["path"] for r in choices]})
        elif choices:
            owner[source] = choices[0]["path"]
    return owner, issues


def doctor(cfg: dict) -> dict:
    inv = inventory(cfg)
    supplied = {c for row in inv for c in row["matched"]}
    source_csv = cfg["inputs"].get("clinical_csv")
    if source_csv and resolve(cfg, source_csv).is_file():
        supplied.update(read_csv(resolve(cfg, source_csv)).columns)
    required = ["code_n", "AGE", "D_DIAG", "BSL_HCY", "y1_is", "y1_is_dd",
                "ONSET_D", "I_BLDSAMP_DT", "GENDER", "BSL_B12",
                "BSL_B9", "BSL_CYSC", "H_SMK", "H_DRINK", "H_HYPT", "H_DIAB", "H_STROKE"]
    image_paths = [cfg["inputs"].get(k) for k in ["derivatives_root", "imaging_csv"]]
    image_ready = any(p and resolve(cfg, p).exists() for p in image_paths)
    result = {
        "status": "INPUTS_PRESENT" if set(required) <= supplied and image_ready else "INPUTS_REQUIRED",
        "mode": cfg["mode"], "sas": inv,
        "missing_required_sources": sorted(set(required) - supplied),
        "missing_core_sources": [s for s, (a, _, _) in FIELDS.items()
                                 if a in ["sex", "b12", "folate", "cysc", "smoking", "drinking",
                                          "hypertension", "diabetes"] and s not in supplied],
        "imaging_input_present": image_ready,
        "message": "Presence is not patient-level eligibility or model validation.",
    }
    dump_json(outdir(cfg) / "doctor.json", result)
    return result


def extract(cfg: dict) -> Path:
    inv = inventory(cfg)
    if not inv:
        raise DataError("No SAS files found. Set inputs.sas_globs; patient data are not bundled.")
    out = outdir(cfg) / "extracted"
    out.mkdir(parents=True, exist_ok=True)
    owner, issues = resolve_owners(inv, cfg)
    if issues:
        raise DataError(f"Ambiguous SAS sources: {issues}")
    merged = None
    formats = {}
    extracts = []
    for index, rec in enumerate(inv):
        selected = [s for s, p in owner.items() if p == rec["path"]]
        if not selected:
            continue
        if "code_n" not in rec["matched"]:
            raise DataError(f"{rec['basename']}: code_n absent; supply an explicitly mapped source")
        selected = ["code_n"] + selected
        target = out / f"{index:03d}_{Path(rec['basename']).stem}.csv"
        count, missing_meta = sas_to_csv(rec["path"], target, selected, cfg["sas"], rec["field_columns"])
        part = read_csv(target)
        unique_ids(part, "code_n", rec["basename"])
        merged = part if merged is None else merged.merge(part, on="code_n", how="outer", validate="one_to_one")
        formats.update({s: rec["formats"].get(rec["field_columns"][s], "") for s in selected})
        dump_json(target.with_suffix(".metadata.json"), {**rec, "selected": selected,
                  "extracted_rows": count, "sas_special_missing": missing_meta})
        extracts.append({"csv": str(target), "rows": count, "columns": selected})
    if merged is None:
        raise DataError("No whitelist fields matched")
    target = out / "clinical_raw.csv"
    merged.to_csv(target, index=False, encoding="utf-8")
    dump_json(out / "formats.json", formats)
    dump_json(out / "manifest.json", {"files": extracts, "variable_owner": owner})
    record_run(cfg, "extract", {"patients": len(merged)})
    return target
