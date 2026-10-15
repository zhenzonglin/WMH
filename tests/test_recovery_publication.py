"""Publication smoke tests: saved evidence, not another statistical analysis."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("tableone")
pytest.importorskip("pubfig")
pytest.importorskip("docx")

from wmh_hcy.recovery_publication import cli
from wmh_hcy.recovery_publication.demo import create_demo
from wmh_hcy.recovery_publication.figures import STYLES, dependence_points, figure2, save_figure
from wmh_hcy.recovery_publication.inputs import load_recovery


def hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*") if p.is_file()}


def test_complete_export_preserves_inputs_and_never_exports_patient_rows(tmp_path):
    inputs = tmp_path / "synthetic_inputs"
    recovery, trajectory = create_demo(inputs)
    before = hashes(inputs)
    output = tmp_path / "publication"
    cli.build(recovery, trajectory, output, styles=list(STYLES), preview_only=True, demo=True)
    assert before == hashes(inputs)
    assert json.loads((output / "status.json").read_text())["status"] == "COMPLETED"
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["no_refitting"] is True
    assert manifest["patient_data_exported"] is False
    assert manifest["mode"] == "synthetic"
    assert "figure_legends.txt" in {p["artifact"] for p in manifest["outputs"]}
    assert len(list((output / "tables").glob("*.docx"))) == 7
    source = pd.read_csv(recovery / "eligible.csv")
    forbidden = source.patient_id.astype(str).tolist()
    for file in output.rglob("*"):
        if file.suffix in {".json", ".csv", ".html", ".txt", ".tex"}:
            content = file.read_text(encoding="utf-8-sig")
            assert not any(value in content for value in forbidden)
    for file in (output / "tables").glob("*.html"):
        assert "Synthetic style preview" in file.read_text(encoding="utf-8")
    for style in STYLES:
        assert (output / style / "figure2.pdf").exists()
        assert (output / style / "figure2.svg").exists()
        assert (output / style / "figure2.png").exists()
        assert not (output / style / "figure2.tiff").exists()
    assert not any(p.name in {"01_recovery", "01_recovery_path", "01_recovery_trajectory"}
                   for p in tmp_path.rglob("*"))


def test_all_styles_use_identical_saved_estimates_and_cis(tmp_path):
    recovery, _ = create_demo(tmp_path / "inputs")
    main = load_recovery(recovery, allow_demo=True)
    values = dependence_points(main)
    for style in STYLES:
        fig = figure2(main, style)
        error = fig.axes[0].containers[0]
        line = error.lines[0]
        actual = line.get_xdata() if style == "compact" else line.get_ydata()
        assert np.allclose(actual, values.probability * 100)
        segments = error.lines[2][0].get_segments()
        axis = 0 if style == "compact" else 1
        assert np.allclose([s[0, axis] for s in segments], values.lower * 100)
        assert np.allclose([s[1, axis] for s in segments], values.upper * 100)
        plt.close(fig)


def test_different_final_states_cannot_be_combined(tmp_path, monkeypatch):
    recovery, trajectory = create_demo(tmp_path / "inputs")
    original = cli.load_trajectory

    def bad(*args, **kwargs):
        data = original(*args, **kwargs)
        data["audit"]["state_counts"]["60"]["0"] -= 1
        data["audit"]["state_counts"]["60"]["1"] += 1
        return data

    monkeypatch.setattr(cli, "load_trajectory", bad)
    with pytest.raises(ValueError, match="five-year state counts"):
        cli.build(recovery, trajectory, tmp_path / "out", styles=["classic"], demo=True)
    assert not (tmp_path / "out").exists()


def test_tiff_is_single_frame_1200dpi_and_vectors_contain_text(tmp_path):
    from PIL import Image
    recovery, _ = create_demo(tmp_path / "inputs")
    data = load_recovery(recovery, allow_demo=True)
    save_figure(figure2(data, "classic"), tmp_path / "fig", "figure2", preview_only=False, demo=True)
    with Image.open(tmp_path / "fig/figure2.tiff") as image:
        assert image.n_frames == 1
        assert image.info["dpi"] == (1200., 1200.)
        assert abs(image.width / 1200 * 25.4 - 173) < .1
    assert "<text" in (tmp_path / "fig/figure2.svg").read_text(encoding="utf-8")


def test_output_refuses_overwrite_and_analysis_folders(tmp_path):
    protected = tmp_path / "outputs/01_recovery/runs/example"
    with pytest.raises(ValueError, match="historical study"):
        cli.safe_output(protected, [])
    existing = tmp_path / "existing"
    existing.mkdir()
    (existing / "keep.txt").write_text("keep")
    with pytest.raises(ValueError, match="not empty"):
        cli.safe_output(existing, [])
    assert (existing / "keep.txt").read_text() == "keep"


def test_publication_source_does_not_import_analysis_workflows():
    root = Path(__file__).parents[1] / "src/wmh_hcy/recovery_publication"
    for file in root.glob("*.py"):
        tree = ast.parse(file.read_text(encoding="utf-8"), feature_version=(3, 11))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert not (node.module or "").startswith("wmh_hcy.studies")
            elif isinstance(node, ast.Import):
                assert not any(a.name.startswith("wmh_hcy.studies") for a in node.names)


def test_no_real_inputs_accepted_with_demo():
    with pytest.raises(SystemExit):
        cli.main(["--demo", "--recovery-run", "not_allowed"])
