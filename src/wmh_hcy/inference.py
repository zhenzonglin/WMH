"""Full-covariance contrasts and fixed-family directional hypothesis decisions."""
from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .common import DataError
from .imputation import pool_scalar


def contrast_vector(data, spec, columns, low=10., high=15.):
    a, b = data.iloc[:1].copy(), data.iloc[:1].copy()
    a[spec.hcy], b[spec.hcy] = low, high
    delta = (spec.transform(b)-spec.transform(a)).iloc[0]
    return np.array([0. if name == "const" else delta[name] for name in columns])


def pooled_contrast(data, spec, fits, kind="cox", low=10., high=15.):
    support = data[spec.hcy].quantile([.05, .95]).to_numpy()
    if not support[0] <= low < high <= support[1]:
        raise DataError(f"Prespecified 10/15 contrast outside central observed {spec.hcy} support: {support.tolist()}")
    estimates, variances = [], []
    for fit in fits:
        if kind == "cox":
            columns, beta, cov = fit.columns, fit.params, fit.covariance
        else:
            columns, beta, cov = list(fit.term), fit.estimate.to_numpy(), fit.attrs["covariance"]
        c = contrast_vector(data, spec, columns, low, high)
        estimates.append(float(c @ beta))
        variances.append(float(c @ cov @ c))
    result = pool_scalar(estimates, variances)
    result.update(hcy_low=low, hcy_high=high,
                  ratio=float(np.exp(result["estimate"])),
                  ratio_lower=float(np.exp(result["lower"])),
                  ratio_upper=float(np.exp(result["upper"])))
    return result


def structural_curve(frames, spec, tables, points=25):
    """Standardized mean log1p-WMH with full HC3 covariance, then Rubin pooling."""
    grid = np.linspace(*frames[0].hcy.quantile([.05, .95]), points)
    rows = []
    for hcy in grid:
        estimates, variances = [], []
        for data, fit in zip(frames, tables, strict=True):
            new = data.copy()
            new["hcy"] = hcy
            matrix = spec.transform(new)
            vector = np.array([1. if name == "const" else matrix[name].mean() for name in fit.term])
            estimates.append(float(vector @ fit.estimate.to_numpy()))
            variances.append(float(vector @ fit.attrs["covariance"] @ vector))
        rows.append({"hcy": hcy, **pool_scalar(estimates, variances)})
    return pd.DataFrame(rows)


def decide_hypotheses(rows):
    """H2 is primary; H1/H3/H4 always form the same three-test Holm family."""
    frame = pd.DataFrame(rows).set_index("hypothesis").reindex(["H1", "H2", "H3", "H4"])
    secondary = ["H1", "H3", "H4"]
    adjusted = multipletests(frame.loc[secondary, "p"].fillna(1), method="holm")[1]
    frame["p_decision"] = frame.p
    frame.loc[secondary, "p_decision"] = adjusted
    frame["decision"] = "NOT_SUPPORTED"
    significant = frame.p_decision.lt(.05)
    frame.loc[significant & frame.estimate.gt(0), "decision"] = "SUPPORTED_POSITIVE"
    frame.loc[significant & frame.estimate.lt(0), "decision"] = "SUPPORTED_OPPOSITE"
    frame.loc[frame.p.isna(), ["p_decision", "decision"]] = [np.nan, "NOT_ESTIMABLE"]
    frame["family"] = ["secondary_Holm_3", "primary_alpha_0.05", "secondary_Holm_3", "secondary_Holm_3"]
    return frame.reset_index()
