"""Continuous displays and covariance-aware contrasts; no additional model fits."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import expit, logit

from ..imputation import pool_scalar
from .models import state_probability_gradients
from .pooling import contrast


def bp_effects(datasets, design, fits, continuous=False):
    observed = datasets[0]
    limits = np.quantile(observed.wmh_ml, [.125, .375, .625, .875])
    grid = np.unique(np.r_[np.linspace(*observed.sbp3.quantile([.025, .975]), 101), 140])
    rows = []
    for index, percentile in enumerate((25, 50, 75)):
        w = float(observed.wmh_ml.quantile(percentile/100))
        local = observed.loc[observed.wmh_ml.between(limits[index], limits[index+1]), "sbp3"]
        low, high = local.quantile([.025, .975]) if continuous else (local.min(), local.max())
        for sbp in grid if continuous else (120, 130, 150):
            supported = len(local) >= 20 and low <= min(sbp, 140) and high >= max(sbp, 140)
            row = {"wmh_percentile": percentile, "wmh_ml": w, "sbp": float(sbp), "reference_sbp": 140,
                   "local_n": len(local), "support_min": float(low), "support_max": float(high),
                   "display": "continuous_central95" if continuous else "fixed_contrast_observed_range"}
            if not supported:
                rows.append({**row, "status": "OUTSIDE_OBSERVED_SUPPORT"})
                continue
            vectors = [design.contrast(d, {"sbp3": sbp, "wmh_ml": w}, {"sbp3": 140, "wmh_ml": w}) for d in datasets]
            r = contrast(fits, vectors)
            rows.append({**row, **r, "HR": np.exp(r["estimate"]), "HR_lower": np.exp(r["lower"]),
                         "HR_upper": np.exp(r["upper"]), "status": "ESTIMATED"})
    return pd.DataFrame(rows)


def bp_distribution(data):
    edges = np.linspace(data.sbp3.min(), data.sbp3.max(), 26)
    limits = np.quantile(data.wmh_ml, [.125, .375, .625, .875])
    rows = []
    for index, percentile in enumerate((25, 50, 75)):
        values = data.loc[data.wmh_ml.between(limits[index], limits[index+1]), "sbp3"]
        counts, _ = np.histogram(values, bins=edges)
        rows.extend({"wmh_percentile": percentile, "sbp_left": edges[j], "sbp_right": edges[j+1], "n": int(n)}
                    for j, n in enumerate(counts))
    return pd.DataFrame(rows)


def pool_probability(values, variances):
    r = pool_scalar(values, variances)
    q = np.clip(r["estimate"], 1e-8, 1-1e-8)
    margin = stats.t.ppf(.975, r["df"])*r["se"]/(q*(1-q))
    return {"probability": r["estimate"], "lower": float(expit(logit(q)-margin)),
            "upper": float(expit(logit(q)+margin))}


def standardized_difference(fitted, x_exposed, x_reference, state=1):
    """Difference of two standardized probabilities with their shared-fit covariance."""
    a, ga = state_probability_gradients(fitted, x_exposed)[state]
    b, gb = state_probability_gradients(fitted, x_reference)[state]
    gradient = ga-gb
    return a-b, float(gradient @ fitted.covariance @ gradient)


def kidney_effects(datasets, design, fits):
    """Persistent versus both-low albuminuria across continuously varying WMH."""
    data = datasets[0]
    groups = [data.loc[data.albuminuria.eq(g), "wmh_ml"] for g in (0, 3)]
    lo = max(v.quantile(.025) for v in groups)
    hi = min(v.quantile(.975) for v in groups)
    if min(map(len, groups)) < 20 or not np.isfinite([lo, hi]).all() or lo >= hi:
        return pd.DataFrame([{"status": "INSUFFICIENT_COMMON_WMH_SUPPORT", "both_low_n": len(groups[0]),
                              "persistent_n": len(groups[1])}])
    grid = np.unique(np.r_[np.expm1(np.linspace(np.log1p(lo), np.log1p(hi), 41)),
                           data.wmh_ml.quantile([.25, .5, .75])])
    rows = []
    for w in grid:
        row = {"wmh_ml": float(w), "support_min": float(lo), "support_max": float(hi),
               "both_low_n": len(groups[0]), "persistent_n": len(groups[1])}
        if not lo <= w <= hi:
            rows.append({**row, "status": "OUTSIDE_OBSERVED_SUPPORT"})
            continue
        values, variances, differences, diffvars, vectors = [], [], [], [], []
        for d, fitted in zip(datasets, fits, strict=True):
            x0 = design.transform(d.assign(albuminuria=0, wmh_ml=w))
            x3 = design.transform(d.assign(albuminuria=3, wmh_ml=w))
            p0, g0 = state_probability_gradients(fitted, x0)[1]
            p3, g3 = state_probability_gradients(fitted, x3)[1]
            values.append([p0, p3])
            variances.append([g0 @ fitted.covariance @ g0, g3 @ fitted.covariance @ g3])
            difference, variance = standardized_difference(fitted, x3, x0)
            differences.append(difference)
            diffvars.append(variance)
            delta = design.contrast(d, {"albuminuria": 3, "wmh_ml": w}, {"albuminuria": 0, "wmh_ml": w})
            # Dependence coefficients are the first multinomial equation; death stays in the full covariance.
            vectors.append(np.r_[delta, np.zeros(len(delta))])
        values, variances = np.asarray(values), np.asarray(variances)
        for j, label in enumerate(("both_low", "persistent")):
            row.update({f"{label}_{k}": v for k, v in pool_probability(values[:, j], variances[:, j]).items()})
        rd = pool_scalar(differences, diffvars)
        row.update({f"difference_{k}": rd[k] for k in ("estimate", "lower", "upper")})
        ratio = contrast(fits, vectors)
        row.update(relative_probability_ratio=float(np.exp(ratio["estimate"])),
                   ratio_lower=float(np.exp(ratio["lower"])), ratio_upper=float(np.exp(ratio["upper"])),
                   status="ESTIMATED", ci="pointwise delta + Rubin; fixed covariate distribution")
        rows.append(row)
    return pd.DataFrame(rows)
