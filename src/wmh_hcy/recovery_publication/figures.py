"""Publication layouts; never compute model estimates or invent unsaved CIs."""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.ticker import MaxNLocator, PercentFormatter, ScalarFormatter

STYLES = {
    "classic": {"colour": "#222222", "second": "#777777", "marker": "o",
                "description": "Monochrome, vertical point estimates and 95% confidence intervals"},
    "colour": {"colour": "#0072B2", "second": "#D55E00", "marker": "o",
               "description": "Blue and vermilion, vertical estimates; shape and line cues retained"},
    "compact": {"colour": "#00796B", "second": "#6A3D9A", "marker": "s",
                "description": "Horizontal interval layout with exact values in a separate column"},
}
FLOW_LABELS = {
    "start": "Registered patients", "clinical": "Registered patients", "base": "Eligible patients",
    "adult_ischemic_stroke": "Adults with ischaemic stroke",
    "available_true_icv": "Intracranial volume available",
    "available_wmh": "WMH volume available", "available_gm119": "Grey matter volume available",
    "available_lesion": "Acute lesion volume available",
    "alive_at_month3": "Alive at three months",
    "independent_at_month3": "Independent at three months\n(mRS 0–2)",
    "known_five_year_functional_state": "Known five-year outcome\nPrimary analysis",
}
EXCLUSION_LABELS = {
    "adult_ischemic_stroke": "Did not meet adult ischaemic\nstroke criteria",
    "available_true_icv": "Intracranial volume unavailable",
    "alive_at_month3": "Did not meet three-month\nsurvival criterion",
    "independent_at_month3": "Did not meet three-month\nindependence criterion",
    "known_five_year_functional_state": "Five-year state unknown",
}


@contextmanager
def publication_style(name: str):
    if name not in STYLES:
        raise ValueError(f"Unknown style: {name}")
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9, "axes.labelsize": 9,
                         "axes.titlesize": 10, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "legend.fontsize": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": .7, "lines.linewidth": 1.1, "pdf.fonttype": 42,
                         "ps.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white",
                         "axes.facecolor": "white", "savefig.facecolor": "white", "axes.grid": False}):
        yield STYLES[name]


def save_figure(fig, outdir: Path, stem: str, *, preview_only: bool, demo: bool) -> list[Path]:
    """Vector PDF/SVG, 300 dpi PNG; optionally single-frame 1200 dpi TIFF.

    BMJ-compatible export choices use an explicit custom pubfig spec; no
    unrelated journal preset is silently applied. Widths are design choices.
    """
    import pubfig as pf

    outdir.mkdir(parents=True, exist_ok=True)
    if demo:
        fig.text(.5, .995, "SYNTHETIC STYLE PREVIEW — NOT RESEARCH RESULTS", ha="center", va="top",
                 fontsize=8, color="#9A3412", weight="bold")
    width_mm, height_mm = fig.get_size_inches() * 25.4
    spec = pf.FigureSpec("bmj_project", font_family="DejaVu Sans", double_column_mm=173,
                         default_raster_dpi=300)
    files = []
    for suffix in (".pdf", ".svg", ".png"):
        path = outdir / (stem + suffix)
        pf.save_figure(fig, path, spec=spec, width=width_mm, height_mm=height_mm,
                       raster_dpi=300, transparent=False, trim=False)
        files.append(path)
    if not preview_only:
        path = outdir / (stem + ".tiff")
        fig.savefig(path, dpi=1200, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
        files.append(path)
    plt.close(fig)
    return files


def flow_rows(data: dict) -> list[dict]:
    """Retain saved exclusion steps; collapse only zero-exclusion intermediate steps."""
    frame = data["flow"]
    rows = []
    for i, row in enumerate(frame.itertuples(index=False)):
        if i == 0 or int(row.excluded_here) or row.step in {
                "independent_at_month3", "known_five_year_functional_state"}:
            label = FLOW_LABELS.get(row.step, row.step.replace("_", " ").capitalize())
            if i == 0:
                label = "CNSR-III registered participants"
            rows.append({"step": row.step, "label": label, "remaining": int(row.remaining),
                         "excluded_here": int(row.excluded_here)})
    if len(rows) > 9:
        raise ValueError("Flow has more than nine nonzero steps; review layout rather than omit exclusions")
    return rows


def figure1(data: dict, style: str):
    with publication_style(style) as theme:
        rows = flow_rows(data)
        fig = plt.figure(figsize=(173 / 25.4, max(180, 27 * len(rows)) / 25.4))
        # Explicit axes keep physical margins stable under pubfig export.
        ax = fig.add_axes((.03, .09, .95, .84))
        ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
        positions = np.linspace(.95, .08, len(rows))
        box_height = min(.095, .65 / len(rows))
        for i, (row, y) in enumerate(zip(rows, positions)):
            ax.add_patch(Rectangle((.035, y - box_height / 2), .54, box_height, facecolor="white",
                                   edgecolor=theme["colour"], linewidth=.9))
            ax.text(.305, y, f'{row["label"]}\nn = {row["remaining"]:,}',
                    ha="center", va="center", fontsize=8.5, linespacing=1.3)
            if i:
                previous = positions[i - 1]
                ax.add_patch(FancyArrowPatch((.305, previous - box_height / 2),
                                            (.305, y + box_height / 2), arrowstyle="-|>",
                                            mutation_scale=10, color="#444444", linewidth=.8))
                if row["excluded_here"]:
                    midpoint = (previous + y) / 2
                    ax.plot([.305, .63], [midpoint, midpoint], color="#555555", lw=.7)
                    reason = EXCLUSION_LABELS.get(row["step"], "Did not meet this criterion")
                    ax.text(.65, midpoint, f'{reason}\nn = {row["excluded_here"]:,}',
                            ha="left", va="center", fontsize=8, linespacing=1.3)
        counts = data["audit"]
        fig.text(.5, .035, f'Five-year states: independent {counts["independent"]:,}; '
                 f'dependent {counts["dependent"]:,}; dead {counts["dead"]:,}', ha="center", fontsize=8)
        return fig


def dependence_points(data: dict) -> pd.DataFrame:
    probs = data["probabilities"]
    gm_values = sorted(probs.gm119_ml.unique())
    if len(gm_values) != 3:
        raise ValueError("Expected exactly three saved GM119 display quantiles")
    points = probs.loc[(probs.state == "dependent") & (probs.gm119_ml == gm_values[1])].copy()
    return points.sort_values("wmh_ml").reset_index(drop=True)


def figure2(data: dict, style: str):
    points = dependence_points(data)
    p = points.probability.to_numpy() * 100
    lo = points.lower.to_numpy() * 100
    hi = points.upper.to_numpy() * 100
    error = np.array([p - lo, hi - p])
    labels = [f'{q}th percentile\n{v:.2f} mL' for q, v in zip((25, 50, 75), points.wmh_ml)]
    with publication_style(style) as theme:
        fig = plt.figure(figsize=(173 / 25.4, 112 / 25.4))
        if style == "compact":
            ax = fig.add_axes((.22, .27, .42, .54))
            y = np.arange(3)
            ax.errorbar(p, y, xerr=error, fmt=theme["marker"], capsize=4, color=theme["colour"], ms=6)
            ax.set(yticks=y, yticklabels=labels,
                   xlabel="Standardised five-year\ndependence probability (%)")
            ax.invert_yaxis()
            ax.set_xlim(0, max(hi) * 1.2)
            ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
            ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
            for yy, prob, lower, upper in zip(y, p, lo, hi):
                ax.text(1.06, yy, f"{prob:.2f}% ({lower:.2f}–{upper:.2f})",
                        transform=ax.get_yaxis_transform(), va="center", fontsize=9)
            ax.text(1.06, 1.08, "Probability (95% CI)", transform=ax.transAxes, fontsize=9)
        else:
            ax = fig.add_axes((.16, .27, .79, .56))
            x = points.wmh_ml.to_numpy()
            ax.errorbar(x, p, yerr=error, fmt=theme["marker"], capsize=4, color=theme["colour"], ms=6)
            pad = max((x[-1] - x[0]) * .17, .1)
            ax.set(xlim=(x[0] - pad, x[-1] + pad), ylim=(0, max(hi) * 1.25),
                   xticks=x, xticklabels=labels, xlabel="Baseline WMH volume",
                   ylabel="Standardised five-year\ndependence probability (%)")
            ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
            for xx, prob, upper in zip(x, p, hi):
                ax.annotate(f"{prob:.2f}%", (xx, upper), xytext=(0, 7), textcoords="offset points",
                            ha="center", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        fig.text(.15, .92, f'GMV fixed at the median ({points.gm119_ml.iloc[0]:.2f} mL)', fontsize=9)
        fig.text(.15, .04, "Points are fixed display values; bars show pointwise 95% CIs.\nNo cut-off is implied.",
                 fontsize=8)
        return fig


def trajectory_figure(data: dict, style: str):
    with publication_style(style) as theme:
        fig = plt.figure(figsize=(173 / 25.4, 145 / 25.4))
        axes = np.array([[fig.add_axes((.12 + c*.30, .59 - r*.41, .24, .27))
                          for c in range(3)] for r in range(2)])
        for r, exposure in enumerate(("wmh", "gm")):
            for c, state in enumerate(("independent", "dependent", "dead")):
                ax = axes[r, c]
                for q, colour, marker, linestyle in ((25, theme["colour"], "o", "-"),
                                                     (75, theme["second"], "s", "--")):
                    values = data["probabilities"].loc[
                        (data["probabilities"].scenario == f"{exposure}_q{q}")
                        & (data["probabilities"].state == state)].sort_values("month")
                    y = values.probability.to_numpy() * 100
                    ax.errorbar(values.month, y,
                                yerr=np.vstack((y-values.lower.to_numpy()*100,
                                                values.upper.to_numpy()*100-y)),
                                fmt=marker, linestyle=linestyle, color=colour, capsize=2, ms=4,
                                label=f"{q}th percentile")
                ax.set(xticks=(36, 48, 60), xlim=(33, 63))
                ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
                ax.yaxis.set_major_locator(MaxNLocator(4, integer=True))
                ax.set_title(f'{chr(65 + r*3+c)}  {state.capitalize()}', loc="left")
                if c == 0:
                    ax.set_ylabel(("WMH" if r == 0 else "GMV") + " scenarios\nState probability (%)")
                if r == 1:
                    ax.set_xlabel("Months after stroke")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, frameon=False, loc="upper center", bbox_to_anchor=(.55, .965), ncol=2)
        fig.text(.11, .045, "Other imaging measure fixed at its median; pointwise 95% CIs.\n"
                 "Lines connect saved visits only; panel y-axis ranges differ.", fontsize=8)
        return fig


def transition_figure(data: dict, style: str):
    transitions = ["independent_to_dependent", "dependent_to_independent",
                   "independent_to_dead", "dependent_to_dead"]
    labels = ["Independent → dependent", "Dependent → independent", "Independent → dead", "Dependent → dead"]
    with publication_style(style) as theme:
        fig = plt.figure(figsize=(173 / 25.4, 118 / 25.4))
        axes = [fig.add_axes((.31 + i*.35, .28, .28, .55)) for i in range(2)]
        vals = data["effects"]
        low = min(vals.ratio_lower.min(), .5) * .9
        high = max(vals.ratio_upper.max(), 2.) * 1.1
        for i, (ax, exposure) in enumerate(zip(axes, ("wmh_ml", "gm119_ml"))):
            frame = vals.loc[vals.exposure == exposure].set_index("transition").loc[transitions]
            y = np.arange(4)
            ratio = frame.ratio.to_numpy()
            ax.axvline(1, color="#888888", linestyle="--", lw=.8)
            ax.errorbar(ratio, y,
                        xerr=np.vstack((ratio-frame.ratio_lower.to_numpy(), frame.ratio_upper.to_numpy()-ratio)),
                        fmt=theme["marker"], color=theme["colour"] if i == 0 else theme["second"], capsize=3)
            ax.set(xscale="log", xlim=(low, high), yticks=y, yticklabels=labels,
                   title="WMH" if i == 0 else "GMV", xlabel="RRR (95% CI) per SD")
            ax.set_xticks([v for v in (.5, .75, 1., 1.5, 2., 3.) if low < v < high])
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.minorticks_off()
            ax.invert_yaxis()
            if i:
                ax.tick_params(axis="y", labelleft=False)
        fig.text(.06, .07, "Reference: remaining in the starting state. RRR compares P(destination)/P(stay).\n"
                 "WMH: log(1 + volume) per SD; GMV: raw volume per SD. Not a hazard ratio.", fontsize=8)
        return fig
