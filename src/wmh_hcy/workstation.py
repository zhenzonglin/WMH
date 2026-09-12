"""Portable workstation setup and explicit audit → image join → analysis stages."""
from __future__ import annotations

import shutil
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .common import DataError, dump_json, load_config, outdir, read_csv, resolve, unique_ids


def discover_sustain(directory: str | Path) -> dict:
    root = Path(directory).expanduser().resolve()
    if not root.is_dir():
        raise DataError(f"SuStaIn path is not a directory: {root}")

    def contains_subjects(path):
        return path.is_dir() and any(s.is_dir() and ((s / "wmh").is_dir() or (s / "t1").is_dir())
                                     for s in path.glob("sub-*"))

    if contains_subjects(root):
        selected = root
    else:
        candidates = [root / "derivatives", root / "outputs/derivatives", root / "output/derivatives"]
        if (root / "derivatives").is_dir():
            candidates.extend(p for p in (root / "derivatives").iterdir() if p.is_dir())
        found = sorted({p for p in candidates if contains_subjects(p)})
        if len(found) != 1:
            raise DataError(f"Expected one directory containing sub-*/wmh or sub-*/t1; found {len(found)}. "
                            "Pass the exact derivatives directory with --sustain-dir.")
        selected = found[0]
    participants = next((p for p in [root / "participants.tsv", selected / "participants.tsv",
                                     selected.parent / "participants.tsv"] if p.is_file()), None)
    return {"derivatives_root": str(selected), "participants_tsv": str(participants) if participants else "",
            "subject_directories": sum(p.is_dir() for p in selected.glob("sub-*")),
            "existing_qc": (selected / "tables/qc_reviews.tsv").is_file()}


def configure(config_path: str, sas_dir: str | None = None, sustain_dir: str | None = None,
              id_map: str | None = None, qc_csv: str | None = None) -> dict:
    target = Path(config_path).resolve()
    root = target.parent.parent
    if target.parent.name != "config" or not (root / "pyproject.toml").is_file():
        raise DataError("Keep the workstation configuration in this repository's config/ directory")
    if not target.name.endswith(".local.yml"):
        raise DataError("Use a *.local.yml configuration; the versioned template stays portable")
    template = target if target.is_file() else root / "config/analysis.yml"
    cfg = yaml.safe_load(template.read_text(encoding="utf-8"))
    cfg["mode"] = "real"
    detected = {}
    if sas_dir is not None:
        source = Path(sas_dir).expanduser().resolve()
        if not source.is_dir():
            raise DataError(f"SAS directory does not exist: {source}")
        cfg["inputs"]["sas_globs"] = [str(source / "**/*")]
        cfg["inputs"]["clinical_csv"] = ""
        cfg["inputs"]["clinical_formats_json"] = ""
    if sustain_dir is not None:
        detected = discover_sustain(sustain_dir)
        cfg["inputs"].update({k: detected[k] for k in ["derivatives_root", "participants_tsv"]})
        cfg["inputs"]["imaging_csv"] = ""
    for key, value in [("id_map_csv", id_map), ("qc_csv", qc_csv)]:
        if value is not None:
            path = Path(value).expanduser().resolve()
            if not path.is_file():
                raise DataError(f"{key} not found: {path}")
            cfg["inputs"][key] = str(path)
    if target.is_file():
        backup = target.with_name(target.name + "." + datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ") + ".bak")
        shutil.copy2(target, backup)
    target.parent.mkdir(exist_ok=True)
    target.write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return {"config": str(target), "mode": "real", "sas_globs": cfg["inputs"]["sas_globs"],
            "imaging": detected or cfg["inputs"].get("derivatives_root", ""),
            "next": "wmh-hcy audit --config " + str(target)}


def image_audit(cfg: dict) -> dict:
    from .imaging import read_imaging

    clinical_path = resolve(cfg, cfg["inputs"]["clinical_csv"]) if cfg["inputs"].get("clinical_csv") \
        else outdir(cfg) / "extracted/clinical_raw.csv"
    if not clinical_path.is_file():
        raise DataError("Run extract before image-audit to obtain clinical patient IDs")
    clinical = read_csv(clinical_path)
    unique_ids(clinical, "code_n", "clinical whitelist CSV")
    frame = read_imaging(cfg)
    matched = frame.patient_id.isin(clinical.code_n)
    ids = set(frame.patient_id)
    summary = {"mode": cfg["mode"], "clinical_unique_ids": len(clinical),
               "mapped_image_unique_ids": len(frame), "clinical_image_id_intersection": int(matched.sum()),
               "clinical_without_mapped_image": int((~clinical.code_n.isin(ids)).sum()),
               "images_without_clinical_id": int((~matched).sum()),
               "whole_wmh_available": int((np.isfinite(frame.wmh_ml) & frame.wmh_ml.ge(0)).sum()),
               "true_icv_available": int((np.isfinite(frame.icv_ml) & frame.icv_ml.gt(0)).sum()),
               "wmh_icv_qc_eligible_matched": int((matched & frame.image_valid).sum()),
               "t1_qc_pass_matched": int((matched & frame.t1_qc.eq("pass") & frame.gm119_ml.notna()).sum()),
               "acute_lesion_qc_pass_matched": int((matched & frame.lesion_qc.eq("pass") & frame.lesion_ml.notna()).sum()),
               "wmh_qc_states": frame.wmh_qc.value_counts().to_dict(),
               "icv_qc_states": frame.icv_qc.value_counts().to_dict()}
    unmatched_file = outdir(cfg) / "unmapped_images.csv"
    summary["unmapped_image_records"] = len(pd.read_csv(unmatched_file)) if unmatched_file.is_file() else 0
    summary["status"] = "READY_FOR_PREPARE" if summary["wmh_icv_qc_eligible_matched"] > 0 else "REVIEW_REQUIRED"
    dump_json(outdir(cfg) / "audit/imaging_summary.json", summary)
    rows = [{"measure": k, "count": v} for k, v in summary.items() if isinstance(v, int)]
    pd.DataFrame(rows).to_csv(outdir(cfg) / "audit/imaging_counts.csv", index=False)
    return summary


def run_pipeline(cfg: dict, through: str = "prepare", hypothesis: str | None = None) -> dict:
    from .analysis import analyse
    from .audit import clinical_audit
    from .cohorts import prepare
    from .reporting import report
    from .sas_extract import extract

    stop = {"audit": 0, "prepare": 3, "analyse": 4, "report": 5}[through]
    stages = [
        ("audit", lambda: clinical_audit(cfg)),
        ("extract", lambda: str(resolve(cfg, cfg["inputs"]["clinical_csv"])) if cfg["inputs"].get("clinical_csv") else str(extract(cfg))),
        ("image_audit", lambda: image_audit(cfg)),
        ("prepare", lambda: prepare(cfg)),
        ("analyse", lambda: analyse(cfg, only=hypothesis)),
        ("report", lambda: str(report(cfg))),
    ]
    state = {"started_utc": datetime.now(UTC).isoformat(), "config": cfg["_config"],
             "mode": cfg["mode"], "requested_through": through, "stages": {}, "status": "RUNNING"}
    target = outdir(cfg) / "pipeline_status.json"
    dump_json(target, state)
    try:
        for name, func in stages[:stop+1]:
            print("Stage:", name, flush=True)
            result = func()
            state["stages"][name] = result
            if isinstance(result, dict) and result.get("status") == "REVIEW_REQUIRED":
                state["status"] = "REVIEW_REQUIRED"
                dump_json(target, state)
                return state
            dump_json(target, state)
        state["status"] = "COMPLETED"
        if "analyse" in state["stages"]:
            estimates = state["stages"]["analyse"]["analyses"]
            primary = {"H1": "01_structure", "H2": "02_recurrence", "H3": "03_month3_update", "H4": "04_function"}
            if estimates.get(primary[hypothesis or "H2"], {}).get("status") != "ESTIMATED":
                state["status"] = "ANALYSIS_REVIEW_REQUIRED"
    except (DataError, OSError, ValueError) as exc:
        state.update(status="FAILED", error=str(exc))
        dump_json(target, state)
        raise
    state["finished_utc"] = datetime.now(UTC).isoformat()
    dump_json(target, state)
    return state


def default_config(command: str) -> str:
    local = Path("config/workstation.local.yml")
    return str(local) if command == "configure" or local.is_file() else "config/analysis.yml"


def load_workstation(config: str | None, command: str) -> dict:
    return load_config(config or default_config(command))
