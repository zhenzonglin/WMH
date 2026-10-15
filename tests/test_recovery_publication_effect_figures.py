"""Saved estimates, compact academic encodings, and fixed combined-figure scope."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.container import BarContainer

from wmh_hcy.recovery_publication.demo import create_demo
from wmh_hcy.recovery_publication.effect_figures import (
    COMBINED_ANALYSES,
    EFFECT_VIEWS,
    combined_effect_figure,
    combined_effect_rows,
    trajectory_effect_figure,
)
from wmh_hcy.recovery_publication.inputs import load_recovery, load_trajectory
from wmh_hcy.recovery_publication.theme import resolved_font


@pytest.fixture

def saved(tmp_path):
    main_run, trajectory_run = create_demo(tmp_path)
    return load_recovery(main_run, allow_demo=True), load_trajectory(trajectory_run, allow_demo=True)


def _artists(fig, *names):
    return [artist for artist in fig.findobj() if artist.get_gid() in names]


def _panel_expected(main):
    return combined_effect_rows(main).query("outcome == 'dependent'").set_index(["exposure", "analysis"]).loc[
        [(exposure, analysis) for exposure in ("wmh_ml", "gm119_ml") for analysis in COMBINED_ANALYSES]]


def test_combined_scope_is_exactly_four_analyses_and_sixteen_saved_rrrs(saved):
    main, _ = saved
    frame = combined_effect_rows(main)
    assert len(frame) == 16
    assert frame.analysis.drop_duplicates().tolist() == list(COMBINED_ANALYSES)
    assert set(frame.outcome) == {"dependent", "dead"}
    assert set(frame.exposure) == {"wmh_ml", "gm119_ml"}
    assert not frame.analysis.isin(["month12", "month24", "raw_wmh", "dependent_or_dead",
                                    "complete_case", "actual_improvement"]).any()
    for row in frame.itertuples(index=False):
        source = main["coefficients"][row.analysis].set_index("term").loc[f"{row.outcome}:{row.exposure}"]
        assert row.n == main["result_details"][row.analysis]["n"]
        np.testing.assert_array_equal([row.ratio, row.ratio_lower, row.ratio_upper, row.p],
                                      source[["ratio", "ratio_lower", "ratio_upper", "p"]].to_numpy(float))


def test_nonestimable_rows_never_reuse_stale_or_zero_estimates(saved):
    main, _ = saved
    main["result_details"]["month36"]["status"] = "NOT_ESTIMABLE"
    rows = combined_effect_rows(main).query("analysis == 'month36'")
    assert len(rows) == 4 and rows.status.eq("NOT_ESTIMABLE").all()
    assert rows[["ratio", "ratio_lower", "ratio_upper", "p"]].isna().all().all()
    del main["result_details"]["month48"]
    assert combined_effect_rows(main).query("analysis == 'month48'").status.eq("NOT_AVAILABLE").all()


@pytest.mark.parametrize("view", EFFECT_VIEWS)
def test_academic_combined_layout_typography_primary_hierarchy_and_no_footers(saved, view):
    main, _ = saved
    before = {name: frame.copy(deep=True) for name, frame in main["coefficients"].items()}
    fig = combined_effect_figure(main, view)
    assert fig.get_label() == "combined_rrr"
    assert fig.get_layout_engine() is not None
    np.testing.assert_allclose(fig.get_size_inches()[0] * 25.4, 178)
    fig.canvas.draw()
    visible = [artist for artist in fig.findobj() if hasattr(artist, "get_text")
               and artist.get_visible() and artist.get_text()]
    assert all(artist.get_fontsize() == 9 for artist in visible)
    assert all(resolved_font()["family"] in artist.get_fontfamily() for artist in visible)
    text = " ".join(artist.get_text() for artist in visible)
    assert "5 y (primary)" in text and "5 y (weighted)" in text
    assert "PRIMARY" not in text and "SENSITIVITY" not in text and "n=" not in text
    assert "12 months" not in text and "24 months" not in text and "OR (" not in text
    assert "Reference" not in text and "Imaging measures and" not in text
    primary_labels = _artists(fig, "primary-row-label")
    assert primary_labels and all(label.get_fontweight() == "bold" for label in primary_labels)
    if view == "heatmap":
        assert len(_artists(fig, "primary-annotation")) == 4
        assert len(_artists(fig, "primary-row-outline")) == 1
        assert all(artist.get_fontweight() == "bold" for artist in _artists(fig, "primary-annotation"))
    else:
        assert len(fig.axes) == 2
        assert len(_artists(fig, "primary-highlight")) == 4
        assert len(_artists(fig, "saved-ci", "saved-ci-transformed")) == 16
        assert not _artists(fig, "primary-annotation", "effect-annotation")
    for name, frame in before.items():
        pd.testing.assert_frame_equal(frame, main["coefficients"][name])
    plt.close(fig)


def test_forest_keeps_all_saved_ratios_and_log_axes_in_two_outcome_panels(saved):
    main, _ = saved
    fig = combined_effect_figure(main, "forest")
    np.testing.assert_allclose(fig.get_size_inches() * 25.4, [178, 105])
    assert len(fig.axes) == 2 and all(ax.get_xscale() == "log" for ax in fig.axes)
    np.testing.assert_allclose(fig.axes[0].get_xlim(), fig.axes[1].get_xlim())
    expected = _panel_expected(main)
    ax = fig.axes[0]
    markers = [line for line in ax.lines if line.get_gid() in {"primary-estimate", "effect-estimate"}]
    np.testing.assert_allclose([line.get_xdata()[0] for line in markers], expected.ratio)
    intervals = [collection.get_segments()[0][:, 0] for collection in ax.collections
                 if collection.get_gid() == "saved-ci"]
    np.testing.assert_allclose(intervals, expected[["ratio_lower", "ratio_upper"]])
    assert len(_artists(fig, "primary-estimate")) == 4 and len(_artists(fig, "effect-estimate")) == 12
    assert markers[0].get_marker() == "s" and markers[1].get_marker() == "o"
    assert markers[0].get_markersize() == 4.5 and markers[1].get_markersize() == 3.5
    plt.close(fig)


def test_diverging_uses_exact_transformed_estimates_and_intervals_without_numeric_columns(saved):
    main, _ = saved
    fig = combined_effect_figure(main, "diverging")
    ax = fig.axes[0]
    expected = _panel_expected(main)
    bars = [container.patches[0] for container in ax.containers if isinstance(container, BarContainer)]
    np.testing.assert_allclose([bar.get_width() for bar in bars], 100 * (expected.ratio - 1))
    np.testing.assert_allclose([bar.get_x() for bar in bars], 0.)
    intervals = [collection.get_segments()[0][:, 0] for collection in ax.collections
                 if collection.get_gid() == "saved-ci-transformed"]
    np.testing.assert_allclose(intervals, 100 * (expected[["ratio_lower", "ratio_upper"]] - 1))
    assert len(_artists(fig, "primary-bar")) == 4 and len(_artists(fig, "effect-bar")) == 12
    assert any(text.get_text() == "Change in outcome ratio (%)" for text in fig.texts)
    plt.close(fig)


def test_seaborn_heatmap_is_centered_log_rrr_with_point_only_annotations(saved):
    main, _ = saved
    fig = combined_effect_figure(main, "heatmap")
    ax = fig.axes[0]
    mesh = ax.collections[0]
    assert mesh.norm.vcenter == 0.
    expected = combined_effect_rows(main)
    np.testing.assert_allclose(mesh.get_array(), np.log(expected.ratio.to_numpy()).reshape(4, 4))
    annotations = _artists(fig, "primary-annotation", "effect-annotation")
    assert len(annotations) == 16
    assert [text.get_text() for text in annotations] == [f"{value:.2f}" for value in expected.ratio]
    assert all("(" not in text.get_text() for text in annotations)
    plt.close(fig)


@pytest.mark.parametrize("view", EFFECT_VIEWS)
def test_missing_secondary_remains_visible_without_a_zero_point(saved, view):
    main, _ = saved
    main["result_details"]["month36"]["status"] = "NOT_ESTIMABLE"
    fig = combined_effect_figure(main, view)
    assert len(_artists(fig, "effect-missing")) == 4
    if view != "heatmap":
        assert len(_artists(fig, "saved-ci", "saved-ci-transformed")) == 12
    else:
        assert np.ma.getmaskarray(fig.axes[0].collections[0].get_array()).sum() == 4
    plt.close(fig)


@pytest.mark.parametrize("view", EFFECT_VIEWS)
def test_trajectory_keeps_eight_saved_transitions_in_compact_layout(saved, view):
    _, trajectory = saved
    before = trajectory["effects"].copy(deep=True)
    fig = trajectory_effect_figure(trajectory, view)
    assert fig.get_label() == "trajectory_rrr"
    np.testing.assert_allclose(fig.get_size_inches() * 25.4, [178, 82])
    if view == "heatmap":
        assert len(_artists(fig, "effect-annotation")) == 8
        expected = before.set_index(["transition", "exposure"]).loc[
            [(transition, exposure) for transition in before.transition.drop_duplicates()
             for exposure in ("wmh_ml", "gm119_ml")], "ratio"].to_numpy().reshape(4, 2)
        np.testing.assert_allclose(fig.axes[0].collections[0].get_array(), np.log(expected))
    else:
        assert len(fig.axes) == 2
        assert len(_artists(fig, "saved-ci", "saved-ci-transformed")) == 8
    assert all("Reference" not in text.get_text() and "P(destination)" not in text.get_text() for text in fig.texts)
    pd.testing.assert_frame_equal(trajectory["effects"], before)
    plt.close(fig)


def test_unknown_views_and_missing_primary_reject(saved):
    main, trajectory = saved
    with pytest.raises(ValueError, match="Unknown effect figure type"):
        combined_effect_figure(main, "palette")
    with pytest.raises(ValueError, match="Unknown effect figure type"):
        trajectory_effect_figure(trajectory, "palette")
    main["result_details"]["primary"]["status"] = "NOT_ESTIMABLE"
    with pytest.raises(ValueError, match="primary analysis"):
        combined_effect_rows(main)
