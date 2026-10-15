"""Export one descriptive Sankey from a frozen cohort; never rebuild v3 assets."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from .figures import save_figure
from .inputs import source_manifest
from .sankey_data import MONTHS, SOURCE_FILES, load_sankey
from .sankey_figure import sankey_figure
from .theme import FONT_SIZE, resolved_font


def _dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _output(path: Path, run: Path) -> Path:
    path, run = path.resolve(), run.resolve()
    if path == run or run in path.parents or path in run.parents:
        raise ValueError("Sankey output must not overlap the saved input run")
    if any(part in {"01_recovery", "01_recovery_path", "01_recovery_trajectory"} for part in path.parts):
        raise ValueError("Sankey output must be outside historical study directories")
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError("Output is not empty; select a new directory and preserve the approved v3 files")
    return path


def _preview(out: Path, *, n: int, demo: bool) -> None:
    from PIL import Image
    with Image.open(out / "mrs_sankey.png") as raster:
        width = raster.width / raster.info.get("dpi", (300, 300))[0] * 96
    notice = ("合成数据样式预览，不能作为研究结果。" if demo else
              f"读取同一保存队列的{n:,}人；缺失保留为独立节点，死亡为吸收状态。")
    html = f"""<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<title>Recovery mRS Sankey</title><style>
body{{font:16px Arial,sans-serif;line-height:1.6;margin:32px;color:#28313b;max-width:1000px}}
img{{display:block;max-width:100%;height:auto;margin:20px 0}}a{{color:#356a92}}
</style><h1>mRS变化桑基图</h1><p>{notice}</p>
<p>3、12、24、36、48、60个月；列间距表示相邻访视，不按月份间隔比例绘制。</p>
<a href="mrs_sankey.png"><img src="mrs_sankey.png" style="width:{width:.1f}px" alt="mRS Sankey"></a>
<p><a href="mrs_sankey.pdf">PDF</a> · <a href="mrs_sankey.svg">SVG</a> ·
<a href="figure_legend.txt">英文图注</a> · <a href="audit.json">人数及流量核对</a></p>
<p>按投稿物理宽度预览（96 CSS dpi）。已确认的v3图形不在此脚本的输出范围内。</p></html>"""
    (out / "preview.html").write_text(html, encoding="utf-8")


def build_sankey(run: Path, out: Path, *, preview_only: bool = False, demo: bool = False) -> Path:
    run = run.resolve()
    out = _output(out, run)
    out.mkdir(parents=True, exist_ok=True)
    status = {"study": "recovery_mrs_sankey", "scope": "descriptive_sankey_only",
              "mode": "synthetic" if demo else "real", "status": "READING"}
    _dump(out / "status.json", status)
    try:
        before = source_manifest(run, list(SOURCE_FILES))
        data = load_sankey(run, allow_demo=demo)
        if source_manifest(run, list(SOURCE_FILES)) != before:
            raise ValueError("Input artifacts changed while reading; retry with a stable saved run")
        audit = {"cohort_n": data["n"], "months": list(MONTHS), "source_run": run.name,
                 "denominator": "all baseline-eligible patients at every visit and interval",
                 "node_count": len(data["nodes"]), "link_count_including_zeros": len(data["links"]),
                 "positive_link_count": int(data["links"].n.gt(0).sum()),
                 **data["validation"], "no_refitting": True, "patient_data_exported": False}
        _dump(out / "audit.json", audit)
        data["nodes"].to_csv(out / "mrs_sankey_nodes.csv", index=False)
        data["links"].to_csv(out / "mrs_sankey_transitions.csv", index=False)
        save_figure(sankey_figure(data["nodes"], data["links"]), out, "mrs_sankey",
                    preview_only=preview_only, demo=demo)
        caption = (
            f"Descriptive mRS transitions among all {data['n']:,} patients independent at three months. "
            "Columns show 3, 12, 24, 36, 48 and 60 months. Node heights and ribbon widths are proportional "
            "to numbers of patients in the fixed starting cohort. Each ribbon represents the observed "
            "transition between two adjacent displayed visits; it does not identify a latent trajectory "
            "class or show an individually traceable path through all visits. Visit columns are equally "
            "spaced, irrespective of elapsed time. mRS 6 denotes death and is absorbing. When a saved "
            "functional state confirms death but the raw mRS is missing, mRS is represented as 6; no "
            "living mRS value is carried forward. Missing scores without confirmed prior death form a "
            "separate category, and patients can return from missing to an observed state. All counts "
            "are unweighted, before covariate imputation; five-year outcome availability does not "
            "determine inclusion. No model is fitted or re-estimated. mRS, modified Rankin scale.\n"
        )
        if demo:
            caption = "SYNTHETIC PREVIEW — NOT RESEARCH RESULTS\n\n" + caption
        (out / "figure_legend.txt").write_text(caption, encoding="utf-8")
        _preview(out, n=data["n"], demo=demo)
        if source_manifest(run, list(SOURCE_FILES)) != before:
            raise ValueError("Input artifacts changed during export; output is not certified")
        modules = ("sankey_cli.py", "sankey_data.py", "sankey_figure.py", "theme.py", "figures.py", "inputs.py")
        module_root = Path(__file__).parent
        font = resolved_font()
        _dump(out / "manifest.json", {
            "mode": status["mode"], "scope": status["scope"], "source_run": str(run), "sources": before,
            "python": platform.python_version(),
            "packages": {name: importlib.metadata.version(name)
                         for name in ("numpy", "pandas", "matplotlib", "seaborn", "SciencePlots")},
            "font": {"family": font["family"], "size_pt": FONT_SIZE,
                     "file_sha256": hashlib.sha256(Path(font["path"]).read_bytes()).hexdigest()},
            "code_sha256": {name: hashlib.sha256((module_root / name).read_bytes()).hexdigest()
                            for name in modules},
            "no_refitting": True, "patient_data_exported": False, "preview_only": preview_only,
            "outputs": [{"artifact": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                        for p in sorted(out.iterdir()) if p.is_file() and p.name != "status.json"],
        })
        _dump(out / "status.json", {**status, "status": "COMPLETED", "cohort_n": data["n"]})
    except (OSError, ValueError, KeyError, ImportError) as exc:
        failure = {**status, "status": "FAILED", "reason": str(exc), "error_type": type(exc).__name__}
        _dump(out / "status.json", failure)
        _dump(out / "failure.json", failure)
        raise
    return out


def create_synthetic_run(root: Path) -> Path:
    """Invent patient histories solely for a visibly labelled rendering example."""
    rng = np.random.default_rng(20261015)
    n = 360
    run = root / "SYNTHETIC_SANKEY_PREVIEW"
    run.mkdir(parents=True)
    frame = pd.DataFrame({"patient_id": [f"SYNTHETIC_SANKEY_{i:04d}" for i in range(n)]})
    latent = rng.choice([0, 1, 2], size=n, p=[.35, .49, .16])
    state_counts = {}
    for index, month in enumerate(MONTHS):
        if index:
            alive = latent != 6
            movement = rng.choice([-1, 0, 1, 2], size=n, p=[.19, .61, .16, .04])
            latent[alive] = np.clip(latent[alive] + movement[alive], 0, 5)
            latent[alive & (rng.random(n) < .012 + .009 * index)] = 6
        observed = latent.astype(float).copy()
        observed[(latent != 6) & (rng.random(n) < .018 * index)] = np.nan
        state = np.select([observed <= 2, (observed >= 3) & (observed <= 5), observed == 6],
                          [0., 1., 2.], default=np.nan)
        state_counts[str(month)] = {str(k): int((state == k).sum()) for k in range(3)}
        state_counts[str(month)]["unknown"] = int(np.isnan(state).sum())
        observed[(state == 2) & (rng.random(n) < .5)] = np.nan
        frame[f"mrs{month}"] = observed
        frame[f"state{month}"] = state
        frame[f"state{month}_conflict"] = False
    frame.to_csv(run / "base_eligible.csv", index=False)
    _dump(run / "audit.json", {"base_eligible_n": n, "state_counts": state_counts,
                               "future_recurrence_fields_requested": []})
    _dump(run / "status.json", {"study": "recovery_trajectory", "run": run.name, "path": str(run.resolve()),
                                "mode": "synthetic", "status": "COMPLETED", "analysis_scope": "trajectory_only",
                                "notice": "Invented rendering fixture; no model was fitted."})
    return run


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="One mRS Sankey from all saved recovery patients; no refits")
    parser.add_argument("--trajectory-run", type=Path, help="Exact saved recovery_trajectory run")
    parser.add_argument("--out", type=Path, help="New empty output directory, separate from approved v3 figures")
    parser.add_argument("--preview-only", action="store_true", help="Skip TIFF; retain PDF, SVG and PNG")
    parser.add_argument("--demo", action="store_true", help="Clearly labelled synthetic layout preview only")
    args = parser.parse_args(argv)
    if args.demo and args.trajectory_run:
        parser.error("--demo cannot read a real saved run")
    if not args.demo and args.trajectory_run is None:
        parser.error("--trajectory-run is required, or use --demo for a synthetic preview")
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    out = args.out or Path("outputs/publication/recovery_sankey") / (("demo_" if args.demo else "") + stamp)
    try:
        if args.demo:
            with tempfile.TemporaryDirectory(prefix="recovery_sankey_synthetic_") as temp:
                build_sankey(create_synthetic_run(Path(temp)), out, preview_only=args.preview_only, demo=True)
        else:
            build_sankey(args.trajectory_run, out, preview_only=args.preview_only)
    except (OSError, ValueError, KeyError, ImportError) as exc:
        print(f"Stopped: {exc}", file=sys.stderr)
        return 2
    print(f"Completed one descriptive mRS Sankey: {out.resolve()}")
    print("Open preview.html; approved v3 figures and model results were not modified.")
    return 0
