"""Standalone Sankey export cannot regenerate or overwrite approved v3 assets."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from wmh_hcy.recovery_publication import sankey_cli


def _hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*") if p.is_file()}


def test_one_sankey_preserves_sources_and_approved_figures_and_exports_only_aggregates(tmp_path):
    run = sankey_cli.create_synthetic_run(tmp_path / "inputs")
    approved = tmp_path / "approved_v3"
    approved.mkdir()
    (approved / "figure2.png").write_bytes(b"Approved illustration, do not modify")
    before, approved_before = _hashes(run), _hashes(approved)
    out = tmp_path / "new_sankey"
    sankey_cli.build_sankey(run, out, preview_only=True, demo=True)
    assert _hashes(run) == before and _hashes(approved) == approved_before
    status = json.loads((out / "status.json").read_text())
    assert status["status"] == "COMPLETED" and status["cohort_n"] == 360
    audit = json.loads((out / "audit.json").read_text())
    assert audit["months"] == [3, 12, 24, 36, 48, 60]
    assert audit["death_absorbing"] is True and audit["living_scores_carried_forward"] == 0
    assert [p.name for p in out.glob("*.png")] == ["mrs_sankey.png"]
    assert (out / "mrs_sankey.pdf").is_file() and (out / "mrs_sankey.svg").is_file()
    assert not list(out.glob("*.tiff"))
    nodes = pd.read_csv(out / "mrs_sankey_nodes.csv")
    links = pd.read_csv(out / "mrs_sankey_transitions.csv")
    assert len(nodes) == 48 and len(links) == 320
    assert nodes.groupby("month").n.sum().eq(360).all()
    assert links.groupby(["from_month", "to_month"]).n.sum().eq(360).all()
    for p in out.iterdir():
        if p.suffix in {".txt", ".csv", ".json", ".html", ".svg"}:
            assert "SYNTHETIC_SANKEY_000" not in p.read_text(encoding="utf-8")
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["no_refitting"] and not manifest["patient_data_exported"]
    assert manifest["font"]["size_pt"] == 9


def test_refuses_existing_output_and_does_not_modify_approved_content(tmp_path):
    run = sankey_cli.create_synthetic_run(tmp_path / "inputs")
    out = tmp_path / "v3"
    out.mkdir()
    (out / "keep.png").write_bytes(b"approved")
    before = _hashes(out)
    with pytest.raises(ValueError, match="not empty"):
        sankey_cli.build_sankey(run, out, demo=True)
    assert _hashes(out) == before


def test_failed_input_is_audited_without_a_figure(tmp_path):
    run = sankey_cli.create_synthetic_run(tmp_path / "inputs")
    frame = pd.read_csv(run / "base_eligible.csv")
    frame.loc[0, "mrs3"] = 99
    frame.to_csv(run / "base_eligible.csv", index=False)
    with pytest.raises(ValueError, match="Invalid mRS"):
        sankey_cli.build_sankey(run, tmp_path / "out", demo=True)
    failure = json.loads((tmp_path / "out/failure.json").read_text())
    assert failure["status"] == "FAILED" and "month 3" in failure["reason"]
    assert not list((tmp_path / "out").glob("*.png"))


def test_changed_input_is_rejected_before_rendering(tmp_path, monkeypatch):
    run = sankey_cli.create_synthetic_run(tmp_path / "inputs")
    original = sankey_cli.load_sankey

    def change(*args, **kwargs):
        data = original(*args, **kwargs)
        path = run / "audit.json"
        path.write_text(path.read_text() + "\n", encoding="utf-8")
        return data

    monkeypatch.setattr(sankey_cli, "load_sankey", change)
    with pytest.raises(ValueError, match="changed while reading"):
        sankey_cli.build_sankey(run, tmp_path / "out", demo=True)
    assert not list((tmp_path / "out").glob("*.png"))


def test_sankey_modules_do_not_import_model_workflows():
    root = Path(sankey_cli.__file__).parent
    for file in root.glob("sankey*.py"):
        tree = ast.parse(file.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert not (node.module or "").startswith("wmh_hcy.studies")
            if isinstance(node, ast.Import):
                assert not any(alias.name.startswith("wmh_hcy.studies") for alias in node.names)


def test_demo_cannot_accept_real_input():
    with pytest.raises(SystemExit):
        sankey_cli.main(["--demo", "--trajectory-run", "saved_run"])


def test_output_cannot_overlap_input_or_historical_studies(tmp_path):
    run = tmp_path / "input"
    for out in (run, run / "child", tmp_path, tmp_path / "01_recovery/runs/new"):
        with pytest.raises(ValueError):
            sankey_cli._output(out, run)
