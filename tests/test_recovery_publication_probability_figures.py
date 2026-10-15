"""Figure semantics against saved synthetic artifacts; no statistical refits."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PathCollection, PolyCollection, QuadMesh
from matplotlib.container import BarContainer, ErrorbarContainer
from matplotlib.text import Text

from wmh_hcy.recovery_publication import figures
from wmh_hcy.recovery_publication.demo import create_demo
from wmh_hcy.recovery_publication.figures import (
    QUANTILES,
    STATES,
    STYLES,
    VISITS,
    dependence_points,
    figure2,
    primary_probability_rows,
    trajectory_figure,
)
from wmh_hcy.recovery_publication.inputs import load_recovery, load_trajectory
from wmh_hcy.recovery_publication.theme import resolved_font


@pytest.fixture
def saved(tmp_path):
    recovery, trajectory = create_demo(tmp_path / "synthetic_inputs")
    return load_recovery(recovery, allow_demo=True), load_trajectory(trajectory, allow_demo=True)


def _assert_stack(ax, source, key, order):
    assert len(ax.containers) == 3
    assert all(isinstance(item, BarContainer) for item in ax.containers)
    matrix = source.pivot(index=key, columns="state", values="probability").loc[list(order), list(STATES)]
    bottom = np.zeros(3)
    for container, state in zip(ax.containers, STATES, strict=True):
        heights = np.array([patch.get_height() for patch in container.patches])
        np.testing.assert_allclose(heights, matrix[state] * 100)
        np.testing.assert_allclose([patch.get_y() for patch in container.patches], bottom)
        assert all(not patch.get_hatch() for patch in container.patches)
        bottom += heights
        assert container.errorbar is None
    np.testing.assert_allclose(bottom, 100, atol=.01)
    assert ax.get_ylim() == (0, 100) and not ax.texts


def _mesh(ax):
    meshes = [item for item in ax.collections if isinstance(item, QuadMesh)]
    assert len(meshes) == 1
    return meshes[0]


def _bands(ax):
    return [item for item in ax.collections if isinstance(item, PolyCollection)]


def _assert_band(band, source):
    vertices = band.get_paths()[0].vertices
    for row in source.itertuples(index=False):
        ys = vertices[np.isclose(vertices[:, 0], row.month), 1]
        np.testing.assert_allclose([ys.min(), ys.max()], [row.lower * 100, row.upper * 100])


def test_three_layout_names_are_chart_types_not_palette_names():
    assert set(STYLES) == {"estimates", "heatmap", "composition"}


def test_full_primary_rows_preserve_all_saved_values_and_add_only_display_labels(saved):
    main, _ = saved
    before = main["probabilities"].copy(deep=True)
    full = primary_probability_rows(main)
    assert len(full) == 27 and full.scenario.nunique() == 9
    assert set(full.wmh_quantile) == set(full.gm_quantile) == {25, 50, 75}
    keys = ["wmh_ml", "gm119_ml", "state"]
    pd.testing.assert_frame_equal(full[before.columns].sort_values(keys).reset_index(drop=True),
                                  before.sort_values(keys).reset_index(drop=True))
    pd.testing.assert_frame_equal(main["probabilities"], before)


def test_primary_estimates_are_three_unconnected_original_volume_points_with_saved_cis(saved):
    main, _ = saved
    fig = figure2(main, "estimates")
    try:
        ax = fig.axes[0]
        assert len(fig.axes) == 1 and len(ax.containers) == 1
        source = dependence_points(main)
        errorbar = ax.containers[0]
        assert isinstance(errorbar, ErrorbarContainer) and errorbar.lines[0] is None
        segments = errorbar.lines[2][0].get_segments()
        np.testing.assert_allclose([s[0, 1] for s in segments], source.lower * 100)
        np.testing.assert_allclose([s[1, 1] for s in segments], source.upper * 100)
        points = [item for item in ax.collections if isinstance(item, PathCollection)]
        assert len(points) == 1
        np.testing.assert_allclose(points[0].get_offsets(), np.column_stack((source.wmh_ml, source.probability * 100)))
        assert not ax.texts and not fig.texts
    finally:
        plt.close(fig)


def test_primary_seaborn_heatmap_shows_nine_saved_probabilities_without_interval_text(saved):
    main, _ = saved
    fig = figure2(main, "heatmap")
    try:
        ax = fig.axes[0]
        source = primary_probability_rows(main)
        source = source.loc[source.state.eq("dependent")]
        matrix = source.pivot(index="gm_quantile", columns="wmh_quantile", values="probability")
        expected = matrix.loc[list(QUANTILES[::-1]), list(QUANTILES)] * 100
        np.testing.assert_allclose(_mesh(ax).get_array().reshape(3, 3), expected)
        assert len(ax.texts) == 9
        for text in ax.texts:
            col, row = [int(v - .5) for v in text.get_position()]
            assert text.get_text() == f"{expected.iloc[row, col]:.1f}"
    finally:
        plt.close(fig)


def test_primary_composition_uses_all_states_without_patterns_or_invented_cis(saved):
    main, _ = saved
    fig = figure2(main, "composition")
    try:
        assert len(fig.axes) == 3
        source = primary_probability_rows(main)
        for ax, q in zip(fig.axes, QUANTILES, strict=True):
            _assert_stack(ax, source.loc[source.gm_quantile.eq(q)], "wmh_quantile", QUANTILES)
        assert [t.get_text() for t in fig.legends[0].get_texts()] == [s.capitalize() for s in STATES]
    finally:
        plt.close(fig)


def test_trajectory_seaborn_lines_preserve_saved_points_and_pointwise_bands(saved, monkeypatch):
    _, data = saved
    calls = []
    original = figures.sns.lineplot

    def spy(*args, **kwargs):
        calls.append(kwargs.copy())
        return original(*args, **kwargs)

    monkeypatch.setattr(figures.sns, "lineplot", spy)
    fig = trajectory_figure(data, "estimates")
    try:
        assert len(fig.axes) == 6 and len(calls) == 12
        assert all(call["estimator"] is None and call["errorbar"] is None for call in calls)
        for index, ax in enumerate(fig.axes):
            row, col = divmod(index, 3)
            exposure = "wmh" if row == 0 else "gm"
            assert len(ax.lines) == len(_bands(ax)) == 2
            for line, band, q in zip(ax.lines, _bands(ax), (25, 75), strict=True):
                source = data["probabilities"].loc[
                    data["probabilities"].scenario.eq(f"{exposure}_q{q}")
                    & data["probabilities"].state.eq(STATES[col])].sort_values("month")
                np.testing.assert_allclose(line.get_xdata(), VISITS)
                np.testing.assert_allclose(line.get_ydata(), source.probability * 100)
                _assert_band(band, source)
                assert line.get_linestyle() in {"-", "--"}
            assert ax.lines[0].get_marker() != ax.lines[1].get_marker()
    finally:
        plt.close(fig)


def test_trajectory_heatmaps_cover_all_54_saved_point_probabilities(saved):
    _, data = saved
    fig = trajectory_figure(data, "heatmap")
    try:
        axes = [ax for ax in fig.axes if len(ax.texts) == 9]
        assert len(axes) == 6
        for index, ax in enumerate(axes):
            row, col = divmod(index, 2)
            exposure = "wmh" if col == 0 else "gm"
            source = data["probabilities"].loc[
                data["probabilities"].scenario.str.startswith(exposure + "_")
                & data["probabilities"].state.eq(STATES[row])].copy()
            source["q"] = source.scenario.str.extract(r"q(\d+)$", expand=False).astype(int)
            expected = source.pivot(index="q", columns="month", values="probability").loc[
                list(QUANTILES[::-1]), list(VISITS)] * 100
            np.testing.assert_allclose(_mesh(ax).get_array().reshape(3, 3), expected)
            for text in ax.texts:
                col_index, row_index = [int(v - .5) for v in text.get_position()]
                assert text.get_text() == f"{expected.iloc[row_index, col_index]:.1f}"
        for row in range(3):
            assert _mesh(axes[row * 2]).norm is _mesh(axes[row * 2 + 1]).norm
    finally:
        plt.close(fig)


def test_trajectory_composition_is_four_panels_with_three_unit_sum_visit_stacks(saved):
    _, data = saved
    fig = trajectory_figure(data, "composition")
    try:
        assert len(fig.axes) == 4
        for ax, scenario in zip(fig.axes, ("wmh_q25", "wmh_q75", "gm_q25", "gm_q75"), strict=True):
            source = data["probabilities"].loc[data["probabilities"].scenario.eq(scenario)]
            _assert_stack(ax, source, "month", VISITS)
    finally:
        plt.close(fig)


def test_trajectory_bands_do_not_clip_saved_normal_intervals_outside_probability_domain(saved):
    _, source = saved
    data = {**source, "probabilities": source["probabilities"].copy()}
    frame = data["probabilities"]
    frame.loc[frame.scenario.eq("wmh_q25") & frame.state.eq("dependent") & frame.month.eq(36), "lower"] = -.025
    frame.loc[frame.scenario.eq("gm_q75") & frame.state.eq("independent") & frame.month.eq(60), "upper"] = 1.04
    fig = trajectory_figure(data, "estimates")
    try:
        assert fig.axes[1].get_ylim()[0] < -2.5
        assert fig.axes[3].get_ylim()[1] > 104
        assert _bands(fig.axes[1])[0].get_paths()[0].vertices[:, 1].min() == -2.5
    finally:
        plt.close(fig)


@pytest.mark.parametrize("style", list(STYLES))
def test_probability_layouts_do_not_modify_inputs_and_render(saved, style):
    main, data = saved
    before = [main["probabilities"].copy(deep=True), data["probabilities"].copy(deep=True)]
    for function, source in ((figure2, main), (trajectory_figure, data)):
        fig = function(source, style)
        try:
            fig.canvas.draw()
        finally:
            plt.close(fig)
    pd.testing.assert_frame_equal(main["probabilities"], before[0])
    pd.testing.assert_frame_equal(data["probabilities"], before[1])


@pytest.mark.parametrize("style", list(STYLES))
@pytest.mark.parametrize("kind", ["primary", "trajectory"])
def test_journal_scale_uses_resolved_font_uniform_size_and_no_explanatory_titles(saved, style, kind):
    main, data = saved
    fig = figure2(main, style) if kind == "primary" else trajectory_figure(data, style)
    sizes = {"primary": {"estimates": (85, 68), "heatmap": (110, 90), "composition": (178, 90)},
             "trajectory": {"estimates": (178, 110), "heatmap": (178, 160), "composition": (178, 120)}}
    try:
        fig.canvas.draw()
        np.testing.assert_allclose(fig.get_size_inches() * 25.4, sizes[kind][style])
        for text in fig.findobj(Text):
            if not text.get_visible() or not text.get_text():
                continue
            assert text.get_fontfamily() == [resolved_font()["family"]]
            assert text.get_fontsize() == 9
        assert not fig.texts and fig._suptitle is None
    finally:
        plt.close(fig)


@pytest.mark.parametrize("kind,style", [("primary", "composition"), ("trajectory", "estimates")])
def test_direct_export_preserves_dimensions_and_demo_marker_does_not_overlap_legend(saved, tmp_path, kind, style):
    main, data = saved
    fig = figure2(main, style) if kind == "primary" else trajectory_figure(data, style)
    before = fig.get_size_inches().copy()
    output = figures.save_figure(fig, tmp_path / kind, style, preview_only=True, demo=True)
    assert {path.suffix for path in output} == {".pdf", ".svg", ".png"}
    np.testing.assert_array_equal(fig.get_size_inches(), before)
    FigureCanvasAgg(fig)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    assert fig._suptitle.get_text() == "Synthetic preview"
    assert fig._suptitle.get_fontsize() == 9 and fig._suptitle.get_fontweight() == "normal"
    assert not fig._suptitle.get_window_extent(renderer).overlaps(fig.legends[0].get_window_extent(renderer))
    assert "<text" in next(path for path in output if path.suffix == ".svg").read_text(encoding="utf-8")
