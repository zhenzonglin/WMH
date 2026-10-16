"""Seaborn/SciencePlots layouts from saved probabilities; no model refitting."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, PercentFormatter

from .theme import PALETTE, STATE_COLOURS, publication_style

STYLES = {
    "estimates": {**PALETTE,
                  "description": "Saved estimates with confidence intervals or pointwise confidence bands"},
    "heatmap": {**PALETTE,
                "description": "Probability heatmaps; exact 95% intervals in companion tables"},
    "composition": {**PALETTE,
                    "description": "Stacked state probabilities; exact 95% intervals in companion tables"},
}
STATES = ("independent", "dependent", "dead")
QUANTILES = (25, 50, 75)
VISITS = (36, 48, 60)
MM = 1 / 25.4


def save_figure(fig, outdir: Path, stem: str, *, preview_only: bool, demo: bool) -> list[Path]:
    """Preserve canvas dimensions in vector PDF/SVG, PNG 300 dpi and TIFF 1200 dpi.

    GridSpec/constrained layout belongs to the figure, not the exporter. No
    post-export restyling, bounding-box crop, or automatic physical resizing.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    with publication_style("estimates"):
        if demo:
            fig.suptitle("Synthetic preview", x=.995, ha="right", fontsize=9, color=".45", fontweight="normal")
        files = []
        for suffix in (".pdf", ".svg", ".png"):
            path = outdir / (stem + suffix)
            fig.savefig(path, dpi=300, facecolor="white", transparent=False, bbox_inches=None)
            files.append(path)
        if not preview_only:
            path = outdir / (stem + ".tiff")
            fig.savefig(path, dpi=1200, facecolor="white", bbox_inches=None,
                        pil_kwargs={"compression": "tiff_lzw"})
            files.append(path)
    plt.close(fig)
    return files


def primary_probability_rows(data: dict) -> pd.DataFrame:
    """Copy all 27 saved rows, adding wmh_quantile, gm_quantile and scenario."""
    frame = data["probabilities"].copy()
    for exposure, label in (("wmh_ml", "wmh_quantile"), ("gm119_ml", "gm_quantile")):
        values = sorted(frame[exposure].unique())
        if len(values) != 3:
            raise ValueError(f"Expected exactly three saved {exposure} display quantiles")
        frame[label] = frame[exposure].map(dict(zip(values, QUANTILES, strict=True)))
    if len(frame) != 27 or frame.duplicated(["wmh_quantile", "gm_quantile", "state"]).any():
        raise ValueError("Expected a complete 27-row primary probability grid")
    frame["scenario"] = [f"wmh_q{w}_gm_q{g}" for w, g in zip(
        frame.wmh_quantile, frame.gm_quantile, strict=True)]
    return frame.sort_values(["gm_quantile", "wmh_quantile", "state"]).reset_index(drop=True)


def dependence_points(data: dict) -> pd.DataFrame:
    """Three saved dependence estimates at median GMV, ordered by original WMH mL."""
    frame = primary_probability_rows(data)
    return frame.loc[frame.state.eq("dependent") & frame.gm_quantile.eq(50)].sort_values(
        "wmh_ml").reset_index(drop=True)


def _probability_axis(ax, *, upper=100., lower=0.):
    ax.set_ylim(lower, upper)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.yaxis.set_major_locator(MaxNLocator(4, integer=True))
    ax.minorticks_off()
    sns.despine(ax=ax)


def _state_legend(fig):
    handles = [Patch(facecolor=STATE_COLOURS[s], edgecolor="none", label=s.capitalize()) for s in STATES]
    fig.legend(handles=handles, frameon=False, loc="outside upper center", ncol=3)


def _composition_bars(ax, frame, key, order):
    """Saved occupancy only; no hatching, numeric clutter, or invented intervals."""
    values = frame.pivot(index=key, columns="state", values="probability").loc[list(order), list(STATES)]
    if not np.allclose(values.sum(axis=1), 1, atol=1e-4, rtol=0):
        raise ValueError("Composition probabilities must sum to one")
    bottom = np.zeros(len(order))
    for state in STATES:
        height = values[state].to_numpy(float) * 100
        ax.bar(np.arange(len(order)), height, width=.68, bottom=bottom, color=STATE_COLOURS[state],
               edgecolor="white", linewidth=.5, label=state.capitalize())
        bottom += height
    ax.set_xlim(-.6, len(order) - .4)
    _probability_axis(ax)


def _heatmap(ax, frame, *, row_key, col_key, rows, cols, norm):
    """Seaborn displays saved point values; interval values stay in companion tables."""
    values = frame.pivot(index=row_key, columns=col_key, values="probability").loc[list(rows), list(cols)]
    sns.heatmap(values * 100, ax=ax, annot=True, fmt=".1f", annot_kws={"fontsize": 9},
                cmap=sns.light_palette("#456F89", as_cmap=True), norm=norm, cbar=False,
                linewidths=.4, linecolor="white", xticklabels=list(cols), yticklabels=list(rows))
    ax.tick_params(axis="both", length=0, labelrotation=0)
    ax.minorticks_off()
    return ax.collections[0]


def _colourbar(fig, image, ax):
    bar = fig.colorbar(image, cax=ax)
    bar.ax.yaxis.set_major_locator(MaxNLocator(3))
    bar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    bar.ax.minorticks_off()
    bar.outline.set_visible(False)
    return bar


def figure2(data: dict, style: str):
    """Three distinct primary probability presentations, all using saved estimates."""
    if style not in STYLES:
        raise ValueError(f"Unknown style: {style}")
    with publication_style(style) as theme:
        frame = primary_probability_rows(data)
        if style == "estimates":
            points = dependence_points(data)
            fig, ax = plt.subplots(figsize=(85 * MM, 68 * MM), layout="constrained")
            x, p = points.wmh_ml.to_numpy(float), points.probability.to_numpy(float) * 100
            low, high = points.lower.to_numpy(float) * 100, points.upper.to_numpy(float) * 100
            ax.errorbar(x, p, yerr=np.vstack((p - low, high - p)), fmt="none",
                        ecolor=theme["colour"], elinewidth=1.0, capsize=2.5, zorder=2)
            sns.scatterplot(x=x, y=p, ax=ax, s=25, color=theme["colour"], edgecolor="none", zorder=3)
            pad = max((x[-1] - x[0]) * .2, .1)
            ax.set(xlim=(x[0] - pad, x[-1] + pad), xticks=x,
                   xticklabels=[f"{v:.1f}" for v in x], xlabel="WMH volume (mL)",
                   ylabel="Dependence probability (%)")
            _probability_axis(ax, upper=min(100., high.max() * 1.15))
            return fig
        if style == "heatmap":
            fig = plt.figure(figsize=(178 * MM, 90 * MM), layout="constrained")
            grid = fig.add_gridspec(1, 3, width_ratios=(1, 1, .045))
            selected = frame.loc[frame.state.isin(("dependent", "dead"))]
            # One absolute-probability scale makes colours comparable across outcomes.
            norm = Normalize(vmin=0, vmax=selected.probability.max() * 100)
            for col, (state, title) in enumerate((("dependent", "Dependence"), ("dead", "Death"))):
                ax = fig.add_subplot(grid[0, col])
                values = selected.loc[selected.state.eq(state)]
                image = _heatmap(ax, values, row_key="gm_quantile", col_key="wmh_quantile",
                                 rows=QUANTILES[::-1], cols=QUANTILES, norm=norm)
                wmh = values.groupby("wmh_quantile").wmh_ml.first()
                gm = values.groupby("gm_quantile").gm119_ml.first()
                ax.set_xticklabels([f"Q{q}\n{wmh.loc[q]:.1f}" for q in QUANTILES])
                ax.set_yticklabels([f"Q{q}\n{gm.loc[q]:.1f}" for q in QUANTILES[::-1]], rotation=0)
                ax.set(xlabel="WMH volume (mL)", ylabel="GMV (mL)")
                ax.set_title(f"({chr(65 + col)}) {title}", loc="left")
            _colourbar(fig, image, fig.add_subplot(grid[0, 2]))
            return fig
        fig, axes = plt.subplots(1, 3, figsize=(178 * MM, 90 * MM), layout="constrained", sharey=True)
        for col, (ax, gm_q) in enumerate(zip(axes, QUANTILES, strict=True)):
            selected = frame.loc[frame.gm_quantile.eq(gm_q)]
            _composition_bars(ax, selected, "wmh_quantile", QUANTILES)
            wmh = selected.groupby("wmh_quantile").wmh_ml.first()
            ax.set_xticks(range(3), [f"{wmh.loc[q]:.1f}" for q in QUANTILES])
            ax.set_xlabel("WMH volume (mL)")
            ax.set_title(f"{chr(65 + col)}  GMV Q{gm_q}", loc="left")
            ax.set_ylabel("State probability (%)" if col == 0 else "")
        _state_legend(fig)
        return fig


def _trajectory_rows(data):
    frame = data["probabilities"].copy()
    frame["quantile"] = frame.scenario.str.extract(r"_q(25|50|75)$", expand=False).astype(int)
    frame["exposure"] = frame.scenario.str.split("_", n=1).str[0]
    return frame


def _trajectory_heatmap(frame):
    fig = plt.figure(figsize=(178 * MM, 160 * MM), layout="constrained")
    grid = fig.add_gridspec(3, 3, width_ratios=(1, 1, .045))
    for row, state in enumerate(STATES):
        all_states = frame.loc[frame.state.eq(state)]
        norm = Normalize(vmin=max(0., all_states.probability.min() * 95),
                         vmax=min(100., all_states.probability.max() * 105))
        for col, exposure in enumerate(("wmh", "gm")):
            ax = fig.add_subplot(grid[row, col])
            selected = all_states.loc[all_states.exposure.eq(exposure)]
            image = _heatmap(ax, selected, row_key="quantile", col_key="month",
                             rows=QUANTILES[::-1], cols=VISITS, norm=norm)
            name = "WMH" if exposure == "wmh" else "GMV"
            ax.set_title(f"{chr(65 + row * 2 + col)}  {name} · {state.capitalize()}", loc="left")
            ax.set_yticklabels([f"Q{q}" for q in QUANTILES[::-1]], rotation=0)
            ax.set_ylabel("Display quantile" if col == 0 else "")
            ax.set_xlabel("Months after stroke" if row == 2 else "")
        _colourbar(fig, image, fig.add_subplot(grid[row, 2]))
    return fig


def trajectory_figure(data: dict, style: str):
    """Saved visit occupancies, never a second fit, reaggregation, or CI procedure."""
    if style not in STYLES:
        raise ValueError(f"Unknown style: {style}")
    with publication_style(style) as theme:
        frame = _trajectory_rows(data)
        if style == "heatmap":
            return _trajectory_heatmap(frame)
        if style == "composition":
            fig, axes = plt.subplots(2, 2, figsize=(178 * MM, 120 * MM), layout="constrained", sharey=True)
            for index, (exposure, quantile) in enumerate((("wmh", 25), ("wmh", 75), ("gm", 25), ("gm", 75))):
                row, col = divmod(index, 2)
                ax = axes[row, col]
                selected = frame.loc[frame.scenario.eq(f"{exposure}_q{quantile}")]
                _composition_bars(ax, selected, "month", VISITS)
                name = "WMH" if exposure == "wmh" else "GMV"
                ax.set_title(f"{chr(65 + index)}  {name} Q{quantile}", loc="left")
                ax.set_xticks(range(3), VISITS)
                ax.set_ylabel("State probability (%)" if col == 0 else "")
                ax.set_xlabel("Months after stroke" if row == 1 else "")
            _state_legend(fig)
            return fig
        fig, axes = plt.subplots(2, 3, figsize=(178 * MM, 110 * MM), layout="constrained", sharex=True)
        for col, state in enumerate(STATES):
            all_states = frame.loc[frame.state.eq(state)]
            for row, exposure in enumerate(("wmh", "gm")):
                ax = axes[row, col]
                selected = all_states.loc[all_states.exposure.eq(exposure)]
                for q, colour, marker, line in ((25, theme["colour"], "o", "-"),
                                                (75, theme["second"], "s", "--")):
                    values = selected.loc[selected["quantile"].eq(q)].sort_values("month").copy()
                    values["probability_percent"] = values.probability * 100
                    sns.lineplot(data=values, x="month", y="probability_percent", ax=ax,
                                 estimator=None, errorbar=None, sort=False, legend=False,
                                 marker=marker, linestyle=line, color=colour, linewidth=1.1, markersize=4)
                    ax.fill_between(values.month.to_numpy(), values.lower.to_numpy() * 100,
                                    values.upper.to_numpy() * 100, color=colour, alpha=.14, linewidth=0)
                low, high = all_states.lower.min() * 100 - 2, all_states.upper.max() * 100 + 2
                _probability_axis(ax, lower=low if state == "independent" else min(0., low), upper=high)
                ax.set(xticks=VISITS, xlim=(34, 62), ylabel="", xlabel="")
                if col == 0:
                    ax.set_ylabel(("WMH" if row == 0 else "GMV") + "\nProbability (%)")
                if row == 1:
                    ax.set_xlabel("Months after stroke")
                ax.set_title(f"{chr(65 + row * 3 + col)}  {state.capitalize()}", loc="left")
        handles = [Line2D([], [], color=theme["colour"], marker="o", linestyle="-", label="Q25"),
                   Line2D([], [], color=theme["second"], marker="s", linestyle="--", label="Q75")]
        fig.legend(handles=handles, frameon=False, loc="outside upper center", ncol=2)
        return fig
