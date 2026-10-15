"""Aggregate alluvial conservation, missingness, absorbing death, and geometry."""
from __future__ import annotations

from itertools import pairwise, product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.patches import PathPatch, Rectangle

from wmh_hcy.recovery_publication.sankey_figure import (
    FLOW_HEIGHT,
    MONTHS,
    NODE_WIDTH,
    STATE_ORDER,
    sankey_figure,
)
from wmh_hcy.recovery_publication.theme import resolved_font


def aggregate_visits(states: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    n = states.shape[1]
    nodes = [{"month": month, "state": state, "n": int((states[index] == state).sum()),
              "proportion": float((states[index] == state).sum() / n)}
             for index, month in enumerate(MONTHS) for state in STATE_ORDER]
    links = []
    for index, (start, stop) in enumerate(pairwise(MONTHS)):
        for source, target in product(STATE_ORDER, repeat=2):
            count = int(((states[index] == source) & (states[index+1] == target)).sum())
            links.append({"from_month": start, "to_month": stop, "source": source, "target": target,
                          "n": count, "proportion": count/n})
    return pd.DataFrame(nodes), pd.DataFrame(links)


@pytest.fixture
def data():
    n = 60
    states = np.tile(np.array([str(i % 3) for i in range(n)], dtype=object), (len(MONTHS), 1))
    for visit in range(1, len(MONTHS)):
        states[visit, :visit] = "6"
        states[visit, 10+visit:13+visit] = "3"
        states[visit, 20:20+visit] = "4"
        states[visit, 30] = "5"  # Positive rare state must not disappear.
        states[visit, n-visit:] = "missing"
    return aggregate_visits(states)


def _patches(fig):
    ax = fig.axes[0]
    nodes = {patch.get_gid(): patch for patch in ax.patches if isinstance(patch, Rectangle)}
    links = {patch.get_gid(): patch for patch in ax.patches if isinstance(patch, PathPatch)}
    return nodes, links


def test_static_six_visit_diagram_keeps_full_cohort_and_omits_zero_nodes(data):
    nodes, links = data
    before_nodes, before_links = nodes.copy(deep=True), links.copy(deep=True)
    fig = sankey_figure(nodes, links)
    drawn_nodes, drawn_links = _patches(fig)
    assert fig.get_label() == "mrs_sankey"
    np.testing.assert_allclose(fig.get_size_inches()*25.4, [178, 115])
    assert len(drawn_nodes) == nodes.n.gt(0).sum()
    assert len(drawn_links) == links.n.gt(0).sum()
    assert "node-3-6" not in drawn_nodes and "node-12-5" in drawn_nodes
    assert drawn_nodes["node-12-5"].get_height() > 0
    assert drawn_nodes["node-12-missing"].get_height() > 0
    assert [tick.get_text() for tick in fig.axes[0].get_xticklabels()] == list(map(str, MONTHS))
    np.testing.assert_allclose(fig.axes[0].get_xticks(), np.arange(6))
    assert fig.axes[0].get_xlabel() == "Visit (months)"
    pd.testing.assert_frame_equal(nodes, before_nodes)
    pd.testing.assert_frame_equal(links, before_links)
    plt.close(fig)


def test_ribbon_geometry_preserves_both_endpoint_counts_and_only_adjacent_links(data):
    nodes, links = data
    fig = sankey_figure(nodes, links)
    drawn_nodes, drawn_links = _patches(fig)
    cohort_n = int(nodes.loc[nodes.month.eq(3), "n"].sum())
    for row in nodes.loc[nodes.n.gt(0)].itertuples(index=False):
        patch = drawn_nodes[f"node-{row.month}-{row.state}"]
        np.testing.assert_allclose(patch.get_height(), FLOW_HEIGHT * row.n / cohort_n, atol=1e-14)
        assert patch.get_y() >= 0 and patch.get_height() > 0
    for row in links.loc[links.n.gt(0)].itertuples(index=False):
        patch = drawn_links[f"link-{row.from_month}-{row.to_month}-{row.source}-{row.target}"]
        vertices = patch.get_path().vertices
        assert np.isfinite(vertices).all() and (vertices[:, 1] >= -1e-14).all()
        np.testing.assert_allclose(vertices[3, 0]-vertices[0, 0], 1-NODE_WIDTH)
        np.testing.assert_allclose(vertices[7, 1]-vertices[0, 1], FLOW_HEIGHT * row.n/cohort_n, atol=1e-14)
        np.testing.assert_allclose(vertices[4, 1]-vertices[3, 1], FLOW_HEIGHT * row.n/cohort_n, atol=1e-14)
        assert np.all(np.diff(vertices[:4, 0]) >= 0)
    # Each interval has N flows, not merely an observed-state denominator.
    for start, stop in pairwise(MONTHS):
        total_height = sum(patch.get_path().vertices[7, 1]-patch.get_path().vertices[0, 1]
                           for key, patch in drawn_links.items() if key.startswith(f"link-{start}-{stop}-"))
        np.testing.assert_allclose(total_height, FLOW_HEIGHT, atol=1e-14)
    plt.close(fig)


def test_state_order_and_uniform_typography_have_no_node_value_clutter(data):
    nodes, links = data
    fig = sankey_figure(nodes, links)
    drawn_nodes, _ = _patches(fig)
    for month in MONTHS:
        top_down = [drawn_nodes[f"node-{month}-{state}"].get_y()
                    for state in STATE_ORDER if f"node-{month}-{state}" in drawn_nodes]
        assert all(a > b for a, b in pairwise(top_down))
    legend = fig.axes[0].get_legend()
    assert [text.get_text() for text in legend.get_texts()] == [*map(str, range(7)), "Missing"]
    assert legend.get_title().get_text() == "mRS"
    fig.canvas.draw()
    visible = [artist for artist in fig.findobj() if hasattr(artist, "get_fontsize")
               and artist.get_visible() and artist.get_text()]
    assert all(artist.get_fontsize() == 9 for artist in visible)
    assert all(resolved_font()["family"] in artist.get_fontfamily() for artist in visible)
    assert not fig.texts and not fig.axes[0].texts
    plt.close(fig)


@pytest.mark.parametrize("target", ["0", "5", "missing"])
def test_confirmed_death_cannot_flow_to_alive_or_missing(target):
    states = np.full((6, 10), "0", dtype=object)
    states[1:, 0] = "6"
    states[2, 0] = target
    nodes, links = aggregate_visits(states)
    with pytest.raises(ValueError, match="death must remain absorbing"):
        sankey_figure(nodes, links)


def test_crossing_adjacent_flows_are_retained_without_inventing_long_paths():
    states = np.array([["0", "1", "0", "1"], ["1", "0", "0", "1"],
                       ["0", "1", "1", "0"], ["1", "0", "1", "0"],
                       ["0", "1", "0", "1"], ["1", "0", "0", "1"]], dtype=object)
    nodes, links = aggregate_visits(states)
    fig = sankey_figure(nodes, links)
    _, ribbons = _patches(fig)
    assert "link-3-12-0-1" in ribbons and "link-3-12-1-0" in ribbons
    assert len(ribbons) == links.n.gt(0).sum()
    plt.close(fig)


def test_single_count_state_in_large_cohort_retains_exact_nonzero_vector_geometry():
    states = np.full((6, 7529), "0", dtype=object)
    states[1:, 0] = "6"
    states[1:, 1] = "5"
    states[1:, 2] = "missing"
    nodes, links = aggregate_visits(states)
    fig = sankey_figure(nodes, links)
    drawn_nodes, drawn_links = _patches(fig)
    rare = drawn_nodes["node-12-5"]
    assert rare.get_height() > 0 and rare.get_linewidth() == .30
    np.testing.assert_allclose(rare.get_height(), FLOW_HEIGHT / 7529, atol=1e-14)
    vertices = drawn_links["link-3-12-0-5"].get_path().vertices
    np.testing.assert_allclose(vertices[4, 1]-vertices[3, 1], FLOW_HEIGHT / 7529, atol=1e-14)
    plt.close(fig)


@pytest.mark.parametrize("mutation", ["node_total", "node_duplicate", "link_duplicate", "missing_grid",
                                     "nonadjacent", "negative", "fractional", "flow", "proportion"])
def test_inconsistent_aggregate_evidence_is_rejected(data, mutation):
    nodes, links = data
    if mutation == "node_total":
        nodes.loc[0, "n"] += 1
    elif mutation == "node_duplicate":
        nodes = pd.concat([nodes, nodes.iloc[:1]])
    elif mutation == "link_duplicate":
        links = pd.concat([links, links.iloc[:1]])
    elif mutation == "missing_grid":
        links = links.iloc[1:]
    elif mutation == "nonadjacent":
        links.loc[0, "to_month"] = 24
    elif mutation == "negative":
        links.loc[0, "n"] = -1
    elif mutation == "fractional":
        links["n"] = links.n.astype(float)
        links.loc[0, "n"] = .5
    elif mutation == "flow":
        # Preserve total and proportions while breaking source-node conservation.
        a = links.index[(links.from_month == 3) & (links.source == "0") & links.n.gt(0)][0]
        b = links.index[(links.from_month == 3) & (links.source == "1") & (links.target == "0")][0]
        links.loc[a, "n"] -= 1
        links.loc[b, "n"] += 1
        links["proportion"] = links.n / 60
    else:
        links.loc[0, "proportion"] = .9
    with pytest.raises(ValueError):
        sankey_figure(nodes, links)
