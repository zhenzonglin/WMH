"""A downloaded-asset revision must preserve sources and unrequested results."""
from __future__ import annotations

import csv
import hashlib
import json

import pandas as pd
import pytest

pytest.importorskip("tableone")
pytest.importorskip("scienceplots")

from wmh_hcy.recovery_publication.cli import build, dump
from wmh_hcy.recovery_publication.demo import create_demo
from wmh_hcy.recovery_publication.export_revision import _hashes, _verify_export, revise_export


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    root = tmp_path_factory.mktemp("downloaded_test_export")
    recovery, trajectory = create_demo(root / "inputs")
    source = root / "export"
    build(recovery, trajectory, source / "v3", styles=["heatmap"], effect_views=["forest"],
          preview_only=True, demo=True)
    # Exercise the real-export identity gate using a test-only generated package.
    for stem in ("table1", "supplement_table1"):
        path = source / "v3/tables" / f"{stem}.csv"
        with path.open(encoding="utf-8-sig", newline="") as stream:
            rows = [row for row in csv.reader(stream) if not row[0].startswith("# ")]
        with path.open("w", encoding="utf-8-sig", newline="") as stream:
            csv.writer(stream).writerows(rows)
    for name in ("status", "manifest"):
        path = source / "v3" / f"{name}.json"
        value = json.loads(path.read_text())
        value["mode"] = "real"
        if name == "manifest":
            for row in value["outputs"]:
                row["sha256"] = hashlib.sha256((source / "v3" / row["artifact"]).read_bytes()).hexdigest()
        dump(path, value)
    (source / "sankey").mkdir()
    (source / "sankey/unchanged.txt").write_text("Saved unrelated asset")
    return source


def test_revision_preserves_sources_unrequested_assets_and_saved_inference(exported, tmp_path):
    before = _hashes(exported)
    out = revise_export(exported, tmp_path / "revised")
    assert _hashes(exported) == before
    changed = json.loads((out / "revision_manifest.json").read_text())["changed_existing_assets"]
    assert all(name in {"v3/manifest.json", "v3/figure_legends.txt", "v3/style_gallery.html"}
               or name.startswith(("v3/heatmap/figure2.", "v3/tables/table1.",
                                   "v3/tables/supplement_table1.")) for name in changed)
    for name in ("v3/source_data/analysis_statuses.csv", "v3/source_data/primary_probability_grid.csv",
                 "v3/effects/forest/trajectory_rrr.png", "sankey/unchanged.txt"):
        assert (out / name).read_bytes() == (exported / name).read_bytes()
    displayed = pd.read_csv(out / "v3/source_data/figure2_heatmap_probabilities.csv")
    assert len(displayed) == 18 and set(displayed.state) == {"dependent", "dead"}
    legend = (out / "trajectory_rrr_legend.txt").read_text(encoding="utf-8")
    assert "(A) White matter" in legend and "(B) Grey matter" in legend
    assert "interval-by-starting-state" in legend
    _verify_export(out / "v3")
    audit = json.loads((out / "revision_manifest.json").read_text())
    for name, digest in audit["output_file_sha256"].items():
        assert _hashes(out)[name] == digest
    assert json.loads((out / "revision_status.json").read_text())["status"] == "COMPLETED"


def test_modified_source_export_is_rejected(exported, tmp_path):
    import shutil
    source = tmp_path / "corrupt"
    shutil.copytree(exported, source)
    with (source / "v3/source_data/primary_probability_grid.csv").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="checksum"):
        revise_export(source, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_revision_refuses_original_or_existing_output(exported, tmp_path):
    before = _hashes(exported)
    with pytest.raises(ValueError, match="overlaps"):
        revise_export(exported, exported)
    out = tmp_path / "nonempty"
    out.mkdir()
    (out / "manual.txt").write_text("preserve")
    with pytest.raises(ValueError, match="not empty"):
        revise_export(exported, out)
    assert _hashes(exported) == before


def test_omitting_critical_checksum_cannot_bypass_integrity(exported, tmp_path):
    import shutil
    source = tmp_path / "omitted"
    shutil.copytree(exported, source)
    path = source / "v3/manifest.json"
    manifest = json.loads(path.read_text())
    manifest["outputs"] = [row for row in manifest["outputs"]
                           if row["artifact"] != "source_data/primary_probability_grid.csv"]
    dump(path, manifest)
    with pytest.raises(ValueError, match="checksum coverage"):
        revise_export(source, tmp_path / "out")
    assert not (tmp_path / "out").exists()
