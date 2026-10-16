"""Revise downloaded aggregate publication assets without accessing patient data."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .baseline_revision import revise_baseline_exports
from .cli import dump, figure_legends, safe_output, write_gallery
from .figures import figure2, primary_probability_rows, save_figure
from .inputs import _probabilities
from .theme import FONT_SIZE, resolved_font

TRAJECTORY_LEGEND = (
    "Figure X. Associations of baseline brain volumes with longitudinal transitions in functional status. "
    "(A) White matter hyperintensity (WMH) volume. (B) Grey matter volume (GMV). Points indicate adjusted "
    "relative risk ratios (RRRs), and horizontal bars indicate 95% confidence intervals, per 1 SD higher "
    "log(1 + WMH volume) or untransformed GMV. RRRs compare the probability of each destination state "
    "with the probability of remaining in the starting state during an adjacent follow-up interval. "
    "Functional independence and dependence were defined as modified Rankin Scale scores of 0–2 and "
    "3–5, respectively; death was an absorbing state. Estimates were obtained from a discrete-time "
    "multinomial transition model across adjacent visits at 3, 12, 24, 36, 48, and 60 months; only "
    "intervals with known starting and ending states contributed to estimation. "
    "Both imaging measures were included jointly, with adjustment for age, sex, education, smoking, "
    "alcohol consumption, hypertension, diabetes, previous stroke, pre-stroke mRS category, "
    "admission National Institutes of Health Stroke Scale score, Trial of Org 10172 in Acute Stroke "
    "Treatment subtype, acute infarct volume, three-month mRS category, and intracranial volume. "
    "The model included follow-up interval and interval-by-starting-state interactions, allowing "
    "interval-specific intercepts to differ by starting state. Imaging associations were also allowed "
    "to vary by starting state. Estimates were pooled across 50 covariate-imputed datasets using "
    "Rubin’s rules, with patient-clustered covariance estimates. The horizontal axis is logarithmic; "
    "the dashed line indicates RRR = 1. This was an exploratory secondary analysis.\n"
)


def _hashes(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def _verify_export(v3: Path) -> dict:
    status = json.loads((v3 / "status.json").read_text(encoding="utf-8"))
    manifest = json.loads((v3 / "manifest.json").read_text(encoding="utf-8"))
    if status.get("status") != "COMPLETED" or status.get("mode") != "real":
        raise ValueError("Revision requires a completed real publication export")
    if manifest.get("mode") != "real" or not manifest.get("no_refitting"):
        raise ValueError("Source manifest does not identify a saved-results publication export")
    recorded = [item["artifact"] for item in manifest["outputs"]]
    required = {"source_data/primary_probability_grid.csv"}
    required.update(f"tables/{stem}.{suffix}" for stem in ("table1", "supplement_table1")
                    for suffix in ("csv", "docx", "html", "tex"))
    if "trajectory" in manifest["sources"]:
        required.add("source_data/trajectory_effects.csv")
    if not recorded or len(set(recorded)) != len(recorded) or not required.issubset(recorded):
        raise ValueError("Source manifest lacks unique checksum coverage of required aggregate inputs")
    local_word_changes = []
    for item in manifest["outputs"]:
        path = (v3 / item["artifact"]).resolve()
        if v3.resolve() not in path.parents:
            raise ValueError("Manifest artifact is outside the publication directory")
        if not path.is_file():
            raise ValueError(f"Source export checksum mismatch: {item['artifact']}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != item["sha256"]:
            if item["artifact"] not in {"tables/table1.docx", "tables/supplement_table1.docx"}:
                raise ValueError(f"Source export checksum mismatch: {item['artifact']}")
            # Word may have been edited locally. The baseline helper must independently
            # verify every Word table cell against the checksum-verified CSV before edits.
            local_word_changes.append({"artifact": item["artifact"],
                                       "workstation_sha256": item["sha256"], "local_sha256": actual})
    manifest["_local_word_changes"] = local_word_changes
    return manifest


def revise_export(source: Path, out: Path) -> Path:
    """Copy a complete export, then replace only requested table and Figure 2 assets."""
    source = source.resolve()
    out = safe_output(out, [source])
    if any(p.is_symlink() for p in source.rglob("*")):
        raise ValueError("Export revision requires ordinary files, not symbolic links")
    before = _hashes(source)
    manifest = _verify_export(source / "v3")
    if "heatmap" not in manifest["styles"]:
        raise ValueError("The source export does not contain the requested Figure 2 heatmap")
    frame = pd.read_csv(source / "v3/source_data/primary_probability_grid.csv")
    _probabilities(frame, ["wmh_ml", "gm119_ml"], "saved primary probability grid")
    if (frame.lower < 0).any() or (frame.upper > 1).any():
        raise ValueError("Saved probability intervals fall outside [0, 1]")
    frame = primary_probability_rows({"probabilities": frame})
    if _hashes(source) != before:
        raise ValueError("Source export changed while reading")
    shutil.copytree(source, out, dirs_exist_ok=True)
    dump(out / "revision_status.json", {"status": "BUILDING", "mode": "real"})
    try:
        v3 = out / "v3"
        table_audit = revise_baseline_exports(source / "v3/tables", v3 / "tables")
        save_figure(figure2({"probabilities": frame}, "heatmap"), v3 / "heatmap", "figure2",
                    preview_only=manifest.get("preview_only", False), demo=False)
        frame.loc[frame.state.isin(["dependent", "dead"])].to_csv(
            v3 / "source_data/figure2_heatmap_probabilities.csv", index=False)
        (v3 / "figure_legends.txt").write_text(
            figure_legends(manifest["styles"], manifest["effect_views"], "trajectory" in manifest["sources"]),
            encoding="utf-8")
        if "trajectory" in manifest["sources"]:
            mi_counts = pd.read_csv(source / "v3/source_data/trajectory_effects.csv", usecols=["m"])["m"]
            if mi_counts.isna().any() or mi_counts.nunique() != 1 or mi_counts.iloc[0] < 1:
                raise ValueError("Trajectory legend requires one saved imputation count")
            m = int(mi_counts.iloc[0])
            if m != mi_counts.iloc[0]:
                raise ValueError("Trajectory imputation count must be an integer")
            legend = TRAJECTORY_LEGEND.replace("50 covariate-imputed", f"{m} covariate-imputed")
            (out / "trajectory_rrr_legend.txt").write_text(legend, encoding="utf-8")
        write_gallery(v3, manifest["styles"], manifest["effect_views"], False,
                      "trajectory" in manifest["sources"])
        if _hashes(source) != before:
            raise ValueError("Source export changed during revision")
        after = _hashes(out)
        changed = sorted(name for name, digest in before.items() if after.get(name) != digest)
        allowed = {"v3/style_gallery.html", "v3/figure_legends.txt"}
        allowed.update(f"v3/heatmap/figure2.{suffix}" for suffix in ("pdf", "svg", "png", "tiff"))
        allowed.update(f"v3/tables/{stem}.{suffix}" for stem in ("table1", "supplement_table1")
                       for suffix in ("csv", "docx", "html", "tex"))
        if set(changed) - allowed:
            raise ValueError("Revision changed an unrequested existing artifact")
        # Retain the workstation manifest as provenance and refresh copied-output hashes.
        shutil.copy2(source / "v3/manifest.json", v3 / "manifest_before_revision.json")
        font = resolved_font()
        manifest["revision"] = {
            "source_export": str(source), "no_refitting": True,
            "baseline_source": "exported integer category frequencies",
            "figure2_states": ["dependent", "dead"],
            "source_docx_checksum_differences": manifest.pop("_local_word_changes"),
            "environment_applies_to": ["heatmap/figure2.*", "tables/table1.*", "tables/supplement_table1.*"],
            "python": platform.python_version(),
            "packages": {name: importlib.metadata.version(name) for name in
                         ("numpy", "pandas", "matplotlib", "seaborn", "SciencePlots", "python-docx")},
            "font": {"family": font["family"], "size_pt": FONT_SIZE,
                     "sha256": hashlib.sha256(Path(font["path"]).read_bytes()).hexdigest()},
        }
        manifest["outputs"] = [
            {"artifact": p.relative_to(v3).as_posix(), "bytes": p.stat().st_size,
             "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in sorted(v3.rglob("*")) if p.is_file() and p.name not in {"manifest.json", "status.json"}
        ]
        dump(v3 / "manifest.json", manifest)
        changed.append("v3/manifest.json")
        dump(out / "revision_status.json", {"status": "COMPLETED", "mode": "real", "no_refitting": True})
        dump(out / "revision_manifest.json", {
            "created_utc": datetime.now(UTC).isoformat(), "source_export": str(source),
            "source_file_sha256": before, "source_unchanged": True, "no_refitting": True,
            "patient_data_read": False, "baseline_revision": table_audit,
            "changed_existing_assets": sorted(changed), "unchanged_existing_assets": len(before) - len(changed),
            "code_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted(Path(__file__).parent.glob("*.py"))},
            "output_file_sha256": _hashes(out),
        })
    except Exception as exc:
        dump(out / "revision_status.json", {"status": "FAILED", "reason": str(exc)})
        raise
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Revise aggregate publication exports without patient data")
    parser.add_argument("--export-root", type=Path, required=True, help="Downloaded folder containing v3/")
    parser.add_argument("--out", type=Path, required=True, help="New empty revision folder")
    args = parser.parse_args(argv)
    try:
        out = revise_export(args.export_root, args.out)
    except (ValueError, OSError, KeyError, ImportError) as exc:
        print(f"Stopped: {exc}", file=sys.stderr)
        return 2
    print(f"Completed requested figure/table revision: {out}")
    return 0
