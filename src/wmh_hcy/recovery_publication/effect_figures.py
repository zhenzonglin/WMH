"""Compact Seaborn/SciencePlots displays of saved RRRs; no fitting or tests.

The combined figure is restricted to primary five-year, three-year,
four-year, and observation-weighted five-year analyses. Confidence intervals
are plotted unchanged in forests, monotonically transformed in diverging
charts, and retained in companion tables for RRR-only heatmap annotations.
"""
from __future__ import annotations

from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator, ScalarFormatter

from .theme import publication_style

EFFECT_VIEWS = {
    "forest": {"description": "Grouped forest plot of saved RRRs and 95% CIs on a logarithmic axis"},
    "heatmap": {"description": "Seaborn heatmap of log RRR; RRR annotations and companion-table 95% CIs"},
    "diverging": {"description": "Grouped bars of 100 × (RRR − 1), with transformed saved 95% CIs"},
}
COMBINED_ANALYSES = ("primary", "month36", "month48", "observation_weighted")
EXPOSURES = ("wmh_ml", "gm119_ml")
OUTCOMES = ("dependent", "dead")
TRANSITIONS = ("independent_to_dependent", "dependent_to_independent",
               "independent_to_dead", "dependent_to_dead")
ANALYSIS_LABELS = {"primary": "5 y (primary)", "month36": "3 y", "month48": "4 y",
                   "observation_weighted": "5 y (weighted)"}
EXPOSURE_LABELS = {"wmh_ml": "WMH", "gm119_ml": "GMV"}
OUTCOME_LABELS = {"dependent": "Dependence", "dead": "Death"}
EFFECT_COLUMNS = ("ratio", "ratio_lower", "ratio_upper", "p")
COMBINED_Y = np.array([8., 7., 6., 5., 3., 2., 1., 0.])


def combined_effect_rows(main: dict) -> pd.DataFrame:
    """Return 16 selected rows without filtering on P values or reusing failures."""
    rows = []
    for analysis in COMBINED_ANALYSES:
        detail = main["result_details"].get(analysis)
        if analysis == "primary" and (detail is None or detail.get("status") != "ESTIMATED"):
            raise ValueError("A saved estimated primary analysis is required for the combined RRR figure")
        status = detail.get("status", "NOT_AVAILABLE") if detail else "NOT_AVAILABLE"
        n = detail.get("n", np.nan) if detail else np.nan
        coefficients = None
        if status == "ESTIMATED":
            if analysis not in main["coefficients"]:
                raise ValueError(f"Saved imaging coefficients are missing for {analysis}")
            coefficients = main["coefficients"][analysis].set_index("term")
        for outcome, exposure in product(OUTCOMES, EXPOSURES):
            row = {"analysis": analysis, "n": n, "outcome": outcome, "exposure": exposure,
                   "status": status, **{name: np.nan for name in EFFECT_COLUMNS}}
            if coefficients is not None:
                term = f"{outcome}:{exposure}"
                if term not in coefficients.index or coefficients.index.duplicated().any():
                    raise ValueError(f"Saved RRR coefficient set is incomplete or duplicated for {analysis}")
                row.update({name: float(coefficients.loc[term, name]) for name in EFFECT_COLUMNS})
            rows.append(row)
    frame = pd.DataFrame(rows)
    _validate_estimated(frame)
    return frame


def _validate_estimated(frame: pd.DataFrame) -> None:
    estimated = frame.loc[frame.status.eq("ESTIMATED")]
    values = estimated[list(EFFECT_COLUMNS)].to_numpy(float)
    if (not np.isfinite(values).all() or (values[:, :3] <= 0).any()
            or (estimated.ratio_lower > estimated.ratio).any()
            or (estimated.ratio > estimated.ratio_upper).any()
            or not estimated.p.between(0, 1).all()):
        raise ValueError("Malformed saved RRR, confidence interval, or P value")


def _trajectory_rows(data: dict) -> pd.DataFrame:
    frame = data["effects"].copy()
    expected = set(product(EXPOSURES, TRANSITIONS))
    if (frame.duplicated(["exposure", "transition"]).any()
            or set(zip(frame.exposure, frame.transition, strict=True)) != expected):
        raise ValueError("Trajectory display requires eight saved exposure/transition contrasts")
    if "status" not in frame:
        frame["status"] = "ESTIMATED"
    _validate_estimated(frame)
    for column in EFFECT_COLUMNS:
        frame.loc[frame.status.ne("ESTIMATED"), column] = np.nan
    return frame


def _limits(frame: pd.DataFrame, view: str) -> tuple[float, float]:
    frame = frame.loc[frame.status.eq("ESTIMATED")]
    if frame.empty:
        return (.5, 2.) if view == "forest" else (-50., 100.)
    lower, upper = min(float(frame.ratio_lower.min()), 1.), max(float(frame.ratio_upper.max()), 1.)
    if view == "forest":
        return lower / 1.12, upper * 1.12
    lower, upper = 100 * (lower - 1), 100 * (upper - 1)
    pad = max((upper - lower) * .10, 5.)
    return max(-100., lower - pad), upper + pad


def _axes_finish(ax, view: str, limits: tuple[float, float]) -> None:
    reference = ax.axvline(1. if view == "forest" else 0., color="#777777", lw=.65, ls="--", zorder=0)
    reference.set_gid("rrr-null-line")
    ax.set_xlim(limits)
    ax.grid(False)
    sns.despine(ax=ax, top=True, right=True, left=True)
    ax.tick_params(axis="y", length=0, labelsize=9)
    ax.tick_params(axis="x", labelsize=9, width=.6, length=3)
    if view == "forest":
        ax.set_xscale("log")
        ticks = [v for v in (.125, .25, .5, .75, 1., 1.5, 2., 3., 4., 6., 8.) if limits[0] < v < limits[1]]
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.minorticks_off()
    else:
        ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.get_offset_text().set_fontsize(9)


def _draw_interval(ax, row, y: float, colour: str, view: str, *, primary: bool = False) -> None:
    if row.status != "ESTIMATED":
        label = "Not estimable" if row.status == "NOT_ESTIMABLE" else "Not available"
        text = ax.text(.98, y, label, transform=ax.get_yaxis_transform(), va="center", ha="right",
                       fontsize=9, color="#777777")
        text.set_gid("effect-missing")
        return
    if primary:
        background = ax.axhspan(y-.38, y+.38, facecolor="#E9EDF2", edgecolor="none", alpha=.55, zorder=-1)
        background.set_gid("primary-highlight")
    values = np.array([row.ratio, row.ratio_lower, row.ratio_upper])
    if view == "diverging":
        values = 100 * (values - 1)
    estimate, lower, upper = values
    errors = np.array([[estimate-lower], [upper-estimate]])
    if view == "forest":
        interval = ax.errorbar(estimate, y, xerr=errors, fmt="s" if primary else "o", color=colour,
                               markersize=4.5 if primary else 3.5, elinewidth=.85, capsize=2,
                               capthick=.75, markeredgewidth=.4, zorder=3)
        interval.lines[0].set_gid("primary-estimate" if primary else "effect-estimate")
    else:
        bar = ax.barh(y, estimate, left=0., height=.48 if primary else .40, color=colour,
                       alpha=.95 if primary else .72, edgecolor="none", zorder=2)
        bar.patches[0].set_gid("primary-bar" if primary else "effect-bar")
        interval = ax.errorbar(estimate, y, xerr=errors, fmt="none", ecolor="#414141",
                               elinewidth=.85, capsize=2, capthick=.75, zorder=3)
    interval.lines[2][0].set_gid("saved-ci-transformed" if view == "diverging" else "saved-ci")


def _combined_panel_rows(frame: pd.DataFrame, outcome: str) -> pd.DataFrame:
    return frame.loc[frame.outcome.eq(outcome)].set_index(["exposure", "analysis"]).loc[
        list(product(EXPOSURES, COMBINED_ANALYSES))].reset_index()


def _combined_intervals(frame: pd.DataFrame, view: str, theme: dict):
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(178/25.4, 105/25.4), layout="constrained")
    fig.get_layout_engine().set(w_pad=.045, h_pad=.035, wspace=.07)
    limits = _limits(frame, view)
    tick_positions = [8.85, *COMBINED_Y[:4], 3.85, *COMBINED_Y[4:]]
    tick_labels = ["WMH", *(ANALYSIS_LABELS[k] for k in COMBINED_ANALYSES),
                   "GMV", *(ANALYSIS_LABELS[k] for k in COMBINED_ANALYSES)]
    for index, (ax, outcome) in enumerate(zip(axes, OUTCOMES, strict=True)):
        rows = _combined_panel_rows(frame, outcome)
        for row, y in zip(rows.itertuples(index=False), COMBINED_Y, strict=True):
            _draw_interval(ax, row, y, theme["colour"] if row.exposure == "wmh_ml" else theme["second"],
                           view, primary=row.analysis == "primary")
        _axes_finish(ax, view, limits)
        ax.set(ylim=(-.65, 9.30), yticks=tick_positions, yticklabels=tick_labels)
        ax.set_title(f"{chr(65+index)}  {OUTCOME_LABELS[outcome]}", loc="left", fontsize=9, weight="bold")
    for tick in axes[0].get_yticklabels():
        if tick.get_text() in {"WMH", "GMV", "5 y (primary)"}:
            tick.set_weight("bold")
        tick.set_gid("primary-row-label" if tick.get_text() == "5 y (primary)" else "effect-row-label")
    fig.supxlabel("RRR" if view == "forest" else "Change in outcome ratio (%)", fontsize=9)
    return fig


def _heatmap(frame: pd.DataFrame, *, trajectory: bool, theme: dict):
    row_keys = list(TRANSITIONS) if trajectory else list(COMBINED_ANALYSES)
    column_keys = list(EXPOSURES) if trajectory else list(product(OUTCOMES, EXPOSURES))
    values = np.full((len(row_keys), len(column_keys)), np.nan)
    labels = np.full(values.shape, "", dtype=object)
    unavailable = {}
    for i, key in enumerate(row_keys):
        for j, column in enumerate(column_keys):
            if trajectory:
                row = frame.loc[frame.transition.eq(key) & frame.exposure.eq(column)].iloc[0]
            else:
                outcome, exposure = column
                row = frame.loc[frame.analysis.eq(key) & frame.outcome.eq(outcome)
                                & frame.exposure.eq(exposure)].iloc[0]
            if row.status == "ESTIMATED":
                values[i, j] = np.log(row.ratio)
                labels[i, j] = f"{row.ratio:.2f}"
            else:
                unavailable[i, j] = "Not estimable" if row.status == "NOT_ESTIMABLE" else "Not available"
    finite = values[np.isfinite(values)]
    bound = max(float(np.abs(finite).max()) if len(finite) else .1, .1)
    norm = TwoSlopeNorm(vmin=-bound, vcenter=0., vmax=bound)
    cmap = sns.blend_palette([theme["colour"], "#F8F8F8", theme["second"]], as_cmap=True)
    fig, ax = plt.subplots(figsize=(178/25.4, (82 if trajectory else 88)/25.4), layout="constrained")
    fig.get_layout_engine().set(w_pad=.045, h_pad=.035)
    ylabels = ([key.replace("_to_", " → ").replace("_", " ").capitalize() for key in row_keys]
               if trajectory else [ANALYSIS_LABELS[key] for key in row_keys])
    xlabels = ([EXPOSURE_LABELS[key] for key in column_keys] if trajectory else
               [f"{EXPOSURE_LABELS[exposure]}\n{OUTCOME_LABELS[outcome]}" for outcome, exposure in column_keys])
    ax.set_facecolor("#EEEEEE")
    sns.heatmap(values, mask=~np.isfinite(values), annot=labels, fmt="", ax=ax, cmap=cmap, norm=norm,
                linewidths=1.3, linecolor="white", xticklabels=xlabels, yticklabels=ylabels,
                annot_kws={"size": 9}, cbar_kws={"orientation": "horizontal", "fraction": .055,
                                                "pad": .10, "shrink": .72, "aspect": 35})
    ax.xaxis.tick_top()
    ax.tick_params(axis="both", length=0, labelsize=9, pad=6)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    for text in ax.texts:
        primary = not trajectory and text.get_position()[1] < 1
        text.set_weight("bold" if primary else "normal")
        text.set_gid("primary-annotation" if primary else "effect-annotation")
    for (i, j), label in unavailable.items():
        text = ax.text(j+.5, i+.5, label.replace(" ", "\n", 1), ha="center", va="center", fontsize=9)
        text.set_gid("effect-missing")
    if not trajectory:
        outline = Rectangle((0, 0), len(column_keys), 1., fill=False, edgecolor="#626A73", lw=.75)
        outline.set_gid("primary-row-outline")
        ax.add_patch(outline)
        ax.get_yticklabels()[0].set_weight("bold")
        ax.get_yticklabels()[0].set_gid("primary-row-label")
    colourbar = ax.collections[0].colorbar
    colourbar.set_ticks([-bound, 0., bound], labels=[f"{np.exp(-bound):.2f}", "1.00", f"{np.exp(bound):.2f}"])
    colourbar.set_label("RRR", fontsize=9)
    colourbar.ax.tick_params(labelsize=9, width=.6, length=2)
    return fig


def _check_view(view: str) -> None:
    if view not in EFFECT_VIEWS:
        raise ValueError(f"Unknown effect figure type: {view}")


def combined_effect_figure(main: dict, view: str):
    """One combined figure; all four saved analyses retained without refitting."""
    _check_view(view)
    frame = combined_effect_rows(main)
    with publication_style("estimates") as theme:
        fig = _heatmap(frame, trajectory=False, theme=theme) if view == "heatmap" else _combined_intervals(
            frame, view, theme)
        fig.set_label("combined_rrr")
        return fig


def trajectory_effect_figure(data: dict, view: str):
    """Four transitions per exposure; staying in origin is defined in the caption."""
    _check_view(view)
    frame = _trajectory_rows(data)
    with publication_style("estimates") as theme:
        if view == "heatmap":
            fig = _heatmap(frame, trajectory=True, theme=theme)
        else:
            fig, axes = plt.subplots(1, 2, sharey=True, figsize=(178/25.4, 82/25.4), layout="constrained")
            fig.get_layout_engine().set(w_pad=.045, h_pad=.035, wspace=.07)
            limits = _limits(frame, view)
            y = np.arange(3, -1, -1, dtype=float)
            labels = [key.replace("_to_", " → ").replace("_", " ").capitalize() for key in TRANSITIONS]
            for index, (ax, exposure) in enumerate(zip(axes, EXPOSURES, strict=True)):
                rows = frame.loc[frame.exposure.eq(exposure)].set_index("transition").loc[list(TRANSITIONS)]
                for row, position in zip(rows.itertuples(index=False), y, strict=True):
                    _draw_interval(ax, row, position, theme["colour"] if index == 0 else theme["second"], view)
                _axes_finish(ax, view, limits)
                ax.set(ylim=(-.55, 3.55), yticks=y, yticklabels=labels)
                ax.set_title(f"{chr(65+index)}  {EXPOSURE_LABELS[exposure]}", loc="left", fontsize=9, weight="bold")
            fig.supxlabel("RRR" if view == "forest" else "Change in outcome ratio (%)", fontsize=9)
        fig.set_label("trajectory_rrr")
        return fig
