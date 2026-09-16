"""Adapter for existing Substain products; never reruns image processing."""
from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from .common import DataError, dump_json, outdir, read_csv, read_json, resolve, unique_ids

VOLUMES = ["wmh_ml", "wmh_raw_ml", "icv_ml", "lesion_ml", "gm119_ml"]
IMAGE_ELIGIBILITY_RULE = "available_volumes_no_review_gate_20260916"


def available_volume(values: pd.Series, *, allow_zero: bool = False) -> pd.Series:
    """Basic numeric availability, independent of any image review label."""
    return np.isfinite(values) & (values.ge(0) if allow_zero else values.gt(0))


def mask_volume_ml(path: Path) -> float:
    img = nib.load(str(path))
    if len(img.shape) != 3:
        raise DataError(f"{path.name}: a 3D binary mask is required")
    unit = img.header.get_xyzt_units()[0]
    factor = {"mm": 1.0, "meter": 1000.0, "micron": 0.001}.get(unit)
    if factor is None:
        raise DataError(f"{path.name}: physical spatial unit is unknown")
    data = np.asanyarray(img.dataobj)
    if not np.isfinite(data).all() or not np.isin(data, [0, 1]).all():
        raise DataError(f"{path.name}: not a finite binary mask; no automatic thresholding")
    det = abs(np.linalg.det(img.affine[:3, :3]))
    if not np.isfinite(det) or det <= 0:
        raise DataError(f"{path.name}: invalid voxel-to-world geometry")
    return float(np.count_nonzero(data) * det * factor**3 / 1000)


def remap_path(cfg: dict, value: str, default: Path) -> Path:
    if not value:
        return default
    for old, new in cfg["imaging"].get("path_prefix_map", {}).items():
        if value == old or value.startswith(old.rstrip("/") + "/"):
            return Path(new + value[len(old):])
    return Path(value)


def read_subject(cfg: dict, subject: Path) -> dict:
    w = read_json(subject / "wmh/wmh_features.json")
    s = read_json(subject / "status/wmh.json")
    detail = s.get("details", {})
    corr = w.get("contralateral_correction", {}) or detail
    row = {"participant_id": str(w.get("participant_id", subject.name.removeprefix("sub-"))),
           "wmh_ml": corr.get("wmh_volume_after_correction_ml", np.nan),
           "wmh_raw_ml": corr.get("wmh_volume_before_correction_ml", np.nan),
           "wmh_source": str(subject / "wmh/wmh_features.json") if w else str(subject / "status/wmh.json")}
    corrected = remap_path(cfg, corr.get("final_wmh", detail.get("corrected_wmh", "")),
                          subject / "wmh/contralateral/wmh_corrected_mask.nii.gz")
    original = remap_path(cfg, corr.get("original_wmh", detail.get("segmentation", "")),
                         subject / "wmh/wmh_mask.nii.gz")
    for column, mask in [("wmh_ml", corrected), ("wmh_raw_ml", original)]:
        if pd.isna(row[column]) and mask.is_file():
            try:
                row[column] = mask_volume_ml(mask)
                row[f"{column}_source"] = str(mask)
            except (DataError, OSError, ValueError) as exc:
                if column == "wmh_ml":
                    raise
                row["wmh_raw_issue"] = str(exc)
    # No sum over the twenty atlas-covered features is used for whole-brain WMH.
    t = read_json(subject / "t1/t1_features.json")
    row["icv_ml"] = t.get("dlicv_icv", {}).get("icv_label702_ml", np.nan)
    gm = t.get("gm119_ml", {})
    row["gm119_ml"] = sum(gm.values()) if len(gm) == 119 else np.nan
    row["icv_source"] = str(subject / "t1/t1_features.json")
    raw_t1 = subject / "t1/nichart_tool_output/DLMUSE_Volumes.csv"
    if raw_t1.is_file():
        volumes = pd.read_csv(raw_t1)
        if len(volumes) != 1:
            raise DataError("Subject DLMUSE volume CSV must have exactly one row")
        if pd.isna(row["icv_ml"]) and "702" in volumes:
            row["icv_ml"] = float(volumes["702"].iloc[0]) / 1000
            row["icv_source"] = str(raw_t1)
        label_file = cfg["inputs"].get("gm119_labels_json")
        if pd.isna(row["gm119_ml"]) and label_file:
            labels = read_json(resolve(cfg, label_file))["labels"]
            columns = [str(x) for x in labels]
            if len(set(columns)) != 119:
                raise DataError("GM119 membership must contain exactly 119 distinct official labels")
            if set(columns) <= set(volumes.columns):
                row["gm119_ml"] = float(volumes[columns].iloc[0].sum()) / 1000
    lesion = subject / "lesion/lesion_space-FLAIR.nii.gz"
    row["lesion_ml"] = np.nan
    if lesion.is_file():
        try:
            row["lesion_ml"] = mask_volume_ml(lesion)
        except (DataError, OSError, ValueError) as exc:
            row["lesion_issue"] = str(exc)
    return row


def read_imaging(cfg: dict) -> pd.DataFrame:
    supplied = cfg["inputs"].get("imaging_csv")
    issues = []
    if supplied:
        frame = read_csv(resolve(cfg, supplied))
    else:
        root_string = cfg["inputs"].get("derivatives_root")
        if not root_string or not resolve(cfg, root_string).is_dir():
            raise DataError("Set inputs.derivatives_root or inputs.imaging_csv")
        root = resolve(cfg, root_string)
        records = []
        for subject in sorted(root.glob("sub-*")):
            if not subject.is_dir():
                continue
            try:
                records.append(read_subject(cfg, subject))
            except (DataError, OSError, ValueError, KeyError) as exc:
                issues.append({"participant_id": subject.name, "reason": str(exc)})
        frame = pd.DataFrame(records)
        pt = cfg["inputs"].get("participants_tsv")
        if pt and not frame.empty:
            participants = pd.read_csv(resolve(cfg, pt), sep="\t", dtype=str, keep_default_na=False)
            unique_ids(participants, "participant_id", "participants.tsv")
            masks = participants.set_index("participant_id").get("lesion_mask", pd.Series(dtype=str))
            for idx in frame.index[frame.lesion_ml.isna()]:
                value = masks.get(frame.at[idx, "participant_id"], "")
                if value:
                    path = remap_path(cfg, value, Path(value))
                    if path.is_file():
                        try:
                            frame.at[idx, "lesion_ml"] = mask_volume_ml(path)
                        except (DataError, OSError) as exc:
                            issues.append({"participant_id": frame.at[idx, "participant_id"],
                                           "reason": str(exc)})
    if frame.empty:
        raise DataError("No readable image records")
    for column in ["wmh_raw_issue", "lesion_issue"]:
        if column in frame:
            for _, row in frame.loc[frame[column].notna()].iterrows():
                issues.append({"participant_id": row.participant_id,
                               "optional_input": column, "reason": str(row[column])})
    unique_ids(frame, "participant_id", "imaging")
    for c in VOLUMES:
        frame[c] = pd.to_numeric(frame[c], errors="coerce") if c in frame else np.nan
    # Review columns in an existing CSV are retained as metadata only.
    # Legacy require_qc / qc_csv settings and external review tables are not used.
    mapping = cfg["inputs"].get("id_map_csv")
    if mapping:
        ids = read_csv(resolve(cfg, mapping))
        unique_ids(ids, "participant_id", "ID map")
        unique_ids(ids, "code_n", "ID map")
        frame = frame.merge(ids[["participant_id", "code_n"]], on="participant_id", how="left",
                            validate="one_to_one").rename(columns={"code_n": "patient_id"})
    else:
        frame["patient_id"] = frame.participant_id.astype(str)
    unmatched = frame.patient_id.isna() | frame.patient_id.eq("")
    frame.loc[unmatched].to_csv(outdir(cfg) / "unmapped_images.csv", index=False)
    frame = frame.loc[~unmatched].copy()
    unique_ids(frame, "patient_id", "mapped images")
    frame["image_valid"] = (available_volume(frame.wmh_ml, allow_zero=True)
                             & available_volume(frame.icv_ml)
                             & (frame.wmh_ml <= frame.icv_ml))
    dump_json(outdir(cfg) / "image_adapter_issues.json", issues)
    frame.to_csv(outdir(cfg) / "imaging.csv", index=False)
    return frame
