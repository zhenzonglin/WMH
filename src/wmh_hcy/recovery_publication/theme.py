"""Journal-scale typography and Seaborn/SciencePlots style, with no estimators."""
from __future__ import annotations

from contextlib import contextmanager
from functools import lru_cache

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401 -- registers the science and no-latex styles
import seaborn as sns
from matplotlib import font_manager

FONT_SIZE = 9
PALETTE = {"colour": "#4C72B0", "second": "#C56E5A", "marker": "o"}
STATE_COLOURS = {"independent": "#6B8BA4", "dependent": "#D5A067", "dead": "#89919B"}


@lru_cache(maxsize=1)
def resolved_font() -> dict[str, str]:
    """Use a metrically similar sans-serif fallback on Linux; record the actual font."""
    for family in ("Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"):
        try:
            path = font_manager.findfont(family, fallback_to_default=False)
        except ValueError:
            continue
        return {"family": font_manager.FontProperties(fname=path).get_name(), "path": path}
    raise ValueError("No supported sans-serif font is available")


@contextmanager
def publication_style(name: str = "estimates"):
    if name not in {"estimates", "heatmap", "composition"}:
        raise ValueError(f"Unknown probability representation: {name}")
    # SciencePlots supplies the baseline line/tick styling. Seaborn controls
    # aesthetics and palette; explicit rc overrides preserve journal typography.
    with (
        plt.style.context(["science", "no-latex"]),
        sns.axes_style("ticks"),
        sns.plotting_context("paper"),
        plt.rc_context({
            "font.family": resolved_font()["family"], "font.size": FONT_SIZE,
            "axes.labelsize": FONT_SIZE, "axes.titlesize": FONT_SIZE,
            "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
            "legend.fontsize": FONT_SIZE, "legend.title_fontsize": FONT_SIZE,
            "figure.titlesize": FONT_SIZE, "axes.titleweight": "normal",
            "axes.labelweight": "normal", "font.weight": "normal",
            "axes.spines.top": False, "axes.spines.right": False,
            "axes.linewidth": .6, "axes.edgecolor": "#555555", "axes.labelcolor": "#222222",
            "text.color": "#222222", "xtick.color": "#444444", "ytick.color": "#444444",
            "xtick.direction": "out", "ytick.direction": "out",
            "xtick.top": False, "ytick.right": False,
            "xtick.minor.visible": False, "ytick.minor.visible": False,
            "xtick.major.size": 2.5, "ytick.major.size": 2.5,
            "xtick.major.width": .6, "ytick.major.width": .6,
            "lines.linewidth": 1., "lines.markersize": 3.5,
            "axes.grid": False, "legend.frameon": False,
            "grid.color": "#E7E9EC", "grid.linewidth": .5,
            "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
            "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "text.usetex": False,
            "savefig.bbox": None,
        }),
    ):
        yield dict(PALETTE)
