"""Adjacent-visit alluvial display from aggregate mRS counts only.

Ribbons represent saved two-visit transitions. Joining their ports does not
reconstruct individual trajectories across three or more visits. Every column
and interval uses the same full-cohort denominator, including missing states.
"""
from __future__ import annotations

from itertools import pairwise, product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.path import Path

from .theme import publication_style

MONTHS = (3, 12, 24, 36, 48, 60)
STATE_ORDER = (*map(str, range(7)), "missing")
NODE_GAP = .010
NODE_WIDTH = .075
FLOW_HEIGHT = 1. - NODE_GAP * (len(STATE_ORDER) - 1)


def _counts(frame: pd.DataFrame, columns: tuple[str, ...], artifact: str) -> None:
    for column in columns:
        try:
            values = pd.to_numeric(frame[column], errors="raise")
        except (ValueError, TypeError) as exc:
            raise ValueError(f"Sankey {artifact} contains a nonnumeric count") from exc
        if (not np.isfinite(values).all() or values.lt(0).any()
                or values.ne(np.floor(values)).any()):
            raise ValueError(f"Sankey {artifact} contains a missing, negative, or fractional count")
        frame[column] = values.astype(int)


def _validate(nodes: pd.DataFrame, links: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    node_columns = ["month", "state", "n", "proportion"]
    link_columns = ["from_month", "to_month", "source", "target", "n", "proportion"]
    if not set(node_columns).issubset(nodes) or not set(link_columns).issubset(links):
        raise ValueError("Sankey source lacks required aggregate columns")
    # Whitelist aggregate fields; no patient-level identifiers or paths are read.
    nodes, links = nodes[node_columns].copy(), links[link_columns].copy()
    _counts(nodes, ("month", "n"), "nodes")
    _counts(links, ("from_month", "to_month", "n"), "links")
    nodes["state"] = nodes.state.astype(str)
    links["source"], links["target"] = links.source.astype(str), links.target.astype(str)
    if (nodes.duplicated(["month", "state"]).any()
            or set(zip(nodes.month, nodes.state, strict=True)) != set(product(MONTHS, STATE_ORDER))):
        raise ValueError("Sankey nodes must contain the complete six-visit, eight-state grid")
    expected_links = {(start, stop, source, target) for start, stop in pairwise(MONTHS)
                      for source, target in product(STATE_ORDER, repeat=2)}
    if (links.duplicated(["from_month", "to_month", "source", "target"]).any()
            or set(zip(links.from_month, links.to_month, links.source, links.target, strict=True))
            != expected_links):
        raise ValueError("Sankey links must contain all state pairs for adjacent visits only")
    totals = nodes.groupby("month").n.sum()
    if len(totals.unique()) != 1 or int(totals.iloc[0]) <= 0:
        raise ValueError("Sankey node totals must share a positive full-cohort denominator")
    cohort_n = int(totals.iloc[0])
    for frame, name in ((nodes, "nodes"), (links, "links")):
        try:
            proportions = pd.to_numeric(frame.proportion, errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Sankey {name} has an invalid full-cohort proportion") from exc
        if (not np.isfinite(proportions).all() or not proportions.between(0, 1).all()
                or not np.allclose(proportions, frame.n / cohort_n, atol=1e-8, rtol=0)):
            raise ValueError(f"Sankey {name} proportions do not use the full-cohort denominator")
        frame["proportion"] = proportions
    indexed = nodes.set_index(["month", "state"]).n
    for start, stop in pairwise(MONTHS):
        interval = links.loc[links.from_month.eq(start) & links.to_month.eq(stop)]
        outgoing = interval.groupby("source").n.sum().reindex(STATE_ORDER)
        incoming = interval.groupby("target").n.sum().reindex(STATE_ORDER)
        if (int(interval.n.sum()) != cohort_n
                or not np.array_equal(outgoing, indexed.loc[start].reindex(STATE_ORDER))
                or not np.array_equal(incoming, indexed.loc[stop].reindex(STATE_ORDER))):
            raise ValueError("Sankey adjacent flows do not reconcile with both node columns")
    if links.loc[links.source.eq("6") & links.target.ne("6"), "n"].gt(0).any():
        raise ValueError("Sankey confirmed death must remain absorbing, including at missing visits")
    return nodes, links, cohort_n


def _colours() -> dict[str, object]:
    living = sns.blend_palette(["#62869E", "#B5C7CA", "#DCC39C", "#BF806A"], n_colors=6)
    return {**dict(zip(map(str, range(6)), living, strict=True)), "6": "#626B75", "missing": "#D9DDE1"}


def _node_positions(nodes: pd.DataFrame, cohort_n: int) -> dict[tuple[int, str], tuple[float, float]]:
    counts = nodes.set_index(["month", "state"]).n
    positions = {}
    for month in MONTHS:
        cursor = 1.
        for state in STATE_ORDER:
            height = FLOW_HEIGHT * int(counts.loc[month, state]) / cohort_n
            lower = max(0., cursor - height)
            positions[month, state] = (lower, cursor)
            cursor = lower - NODE_GAP
    return positions


def _ribbon(x_start: float, x_stop: float, source: tuple[float, float], target: tuple[float, float]) -> Path:
    source_low, source_high = source
    target_low, target_high = target
    bend = .42 * (x_stop - x_start)
    vertices = [(x_start, source_low), (x_start+bend, source_low), (x_stop-bend, target_low),
                (x_stop, target_low), (x_stop, target_high), (x_stop-bend, target_high),
                (x_start+bend, source_high), (x_start, source_high), (x_start, source_low)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    return Path(vertices, codes)


def sankey_figure(nodes: pd.DataFrame, links: pd.DataFrame):
    """Draw aggregate adjacent-visit flows without fitting or completing missing visits."""
    nodes, links, cohort_n = _validate(nodes, links)
    positions = _node_positions(nodes, cohort_n)
    colours = _colours()
    node_counts = nodes.set_index(["month", "state"]).n
    link_counts = links.set_index(["from_month", "to_month", "source", "target"]).n
    unit = FLOW_HEIGHT / cohort_n
    with publication_style("estimates"):
        fig, ax = plt.subplots(figsize=(178/25.4, 115/25.4), layout="constrained")
        fig.get_layout_engine().set(w_pad=.04, h_pad=.035)
        ribbons = []
        for column, (start, stop) in enumerate(pairwise(MONTHS)):
            used_source = dict.fromkeys(STATE_ORDER, 0)
            used_target = dict.fromkeys(STATE_ORDER, 0)
            for source, target in product(STATE_ORDER, repeat=2):
                n = int(link_counts.loc[start, stop, source, target])
                if n == 0:
                    continue
                source_low, source_top = positions[start, source]
                target_low, target_top = positions[stop, target]
                source_port = (max(source_low, source_top-(used_source[source]+n)*unit),
                               source_top-used_source[source]*unit)
                target_port = (max(target_low, target_top-(used_target[target]+n)*unit),
                               target_top-used_target[target]*unit)
                path = _ribbon(column+NODE_WIDTH/2, column+1-NODE_WIDTH/2, source_port, target_port)
                patch = PathPatch(path, facecolor=colours[source], edgecolor="none",
                                  alpha=.60 if source == "missing" else .28, zorder=1)
                patch.set_gid(f"link-{start}-{stop}-{source}-{target}")
                ribbons.append((n, patch))
                used_source[source] += n
                used_target[target] += n
        # Thin flows are painted last so larger ribbons cannot conceal them.
        for _, patch in sorted(ribbons, key=lambda item: -item[0]):
            ax.add_patch(patch)
        for column, month in enumerate(MONTHS):
            for state in STATE_ORDER:
                if int(node_counts.loc[month, state]) == 0:
                    continue
                low, high = positions[month, state]
                node = Rectangle((column-NODE_WIDTH/2, low), NODE_WIDTH, high-low,
                                 facecolor=colours[state], edgecolor=colours[state], linewidth=.30, zorder=3)
                node.set_gid(f"node-{month}-{state}")
                ax.add_patch(node)
        handles = [Patch(facecolor=colours[state], edgecolor="none", label="Missing" if state == "missing" else state)
                   for state in STATE_ORDER]
        ax.legend(handles=handles, title="mRS", ncol=8, loc="lower center", bbox_to_anchor=(.5, 1.015),
                  fontsize=9, title_fontsize=9, frameon=False, columnspacing=.8, handlelength=1.,
                  handletextpad=.4, borderpad=0, labelspacing=.35)
        ax.set(xlim=(-.12, len(MONTHS)-1+.12), ylim=(-.01, 1.01), yticks=[],
               xticks=np.arange(len(MONTHS)), xticklabels=[str(month) for month in MONTHS])
        ax.set_xlabel("Visit (months)", fontsize=9)
        ax.tick_params(axis="x", length=0, labelsize=9, pad=4)
        ax.grid(False)
        sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
        fig.set_label("mrs_sankey")
        return fig
