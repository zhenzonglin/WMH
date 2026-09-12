"""Delayed-entry cause-specific cumulative incidence and patient bootstrap + MI uncertainty."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .common import DataError
from .design import Design
from .imputation import pool_scalar
from .models import CoxFit, fit_cause


def cumulative_incidence(data: pd.DataFrame, spec: Design, stroke: CoxFit, death: CoxFit,
                         horizon: float = 365) -> np.ndarray:
    x = spec.transform(data).to_numpy()
    scores = [np.exp(x @ f.params) for f in [stroke, death]]
    times = np.union1d(stroke.times, death.times)
    times = times[times <= horizon]
    maps = [dict(zip(f.times, f.increments, strict=True)) for f in [stroke, death]]
    survival = np.ones(len(data))
    cif = np.zeros(len(data))
    for t in times:
        active = data.entry.to_numpy() < t
        a = maps[0].get(t, 0.0) * scores[0] * active
        b = maps[1].get(t, 0.0) * scores[1] * active
        total = a+b
        # Exponential integration of each pooled baseline-hazard increment.
        # Stable under tied daily events and nonnegative for all positive hazards.
        fraction = np.divide(a, total, out=np.zeros_like(a), where=total > 0)
        cif += survival * (-np.expm1(-total)) * fraction
        survival *= np.exp(-total)
    if not np.isfinite(cif).all() or (cif < 0).any() or (cif > 1+1e-10).any():
        raise DataError("Invalid cumulative incidence")
    return cif


def make_grid(data: pd.DataFrame, cfg: dict, hcol: str = "hcy") -> pd.DataFrame:
    ranks = data.wmh_ml.rank(pct=True)
    hgrid = np.linspace(*data[hcol].quantile([0.05, 0.95]), cfg["analysis"]["risk_hcy_grid_points"])
    refs = cfg["analysis"]["risk_hcy_values"]
    rows = []
    for wq in cfg["analysis"]["wmh_quantiles"]:
        near = data.loc[(ranks-wq).abs() <= 0.10, hcol]
        lower, upper = near.quantile([0.05, 0.95])
        w = float(data.wmh_ml.quantile(wq))
        for h in sorted(set(hgrid.tolist()+refs)):
            rows.append({"wmh_quantile": wq, "wmh_ml": w, "hcy": float(h),
                         "supported": bool(lower <= h <= upper), "reference": h in refs})
    return pd.DataFrame(rows)


def standardized_risks(data: pd.DataFrame, spec: Design, stroke: CoxFit, death: CoxFit,
                       grid: pd.DataFrame, horizon: float) -> np.ndarray:
    estimates = []
    for row in grid.itertuples():
        if not row.supported:
            estimates.append(np.nan)
            continue
        new = data.copy()
        new[spec.hcy] = row.hcy
        new[spec.wmh_column] = row.wmh_ml
        estimates.append(float(cumulative_incidence(new, spec, stroke, death, horizon).mean()))
    return np.asarray(estimates)


def risk_analysis(completed: list[pd.DataFrame], spec: Design, fits: list[tuple[CoxFit, CoxFit]],
                  cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    grid = make_grid(completed[0], cfg, spec.hcy)
    predictions = [standardized_risks(d, spec, s, death, grid, cfg["analysis"]["horizon"])
                   for d, (s, death) in zip(completed, fits, strict=True)]
    curve = grid.copy()
    with np.errstate(invalid="ignore"):
        curve["risk"] = np.mean(predictions, axis=0)
    # Confidence intervals at prespecified 10/15 µmol/L contrasts, not thousands of searched cutpoints.
    landmarks = grid.loc[grid.reference & grid.supported].reset_index(drop=True)
    b = cfg["analysis"]["bootstrap_per_imputation"]
    point = []
    variances = []
    rdpoint = []
    rdvar = []
    failures = []
    rng = np.random.default_rng(cfg["analysis"]["seed"] + 471)
    pair_keys = []
    for w in cfg["analysis"]["wmh_quantiles"]:
        pair = landmarks.index[landmarks.wmh_quantile.eq(w)].tolist()
        if len(pair) == 2:
            pair_keys.append((w, pair[0], pair[1]))
    if b < 2 or landmarks.empty:
        return curve, pd.DataFrame(), pd.DataFrame(), {"status": "NO_BOOTSTRAP_OR_SUPPORTED_LANDMARKS"}
    for m, (d, (s, death)) in enumerate(zip(completed, fits, strict=True)):
        estimates = standardized_risks(d, spec, s, death, landmarks, cfg["analysis"]["horizon"])
        samples = []
        for rep in range(b):
            sample = d.iloc[rng.integers(0, len(d), len(d))].reset_index(drop=True)
            try:
                sf = fit_cause(sample, spec, 1)
                df = fit_cause(sample, spec, 2)
                samples.append(standardized_risks(sample, spec, sf, df, landmarks,
                                                   cfg["analysis"]["horizon"]))
            except (DataError, ValueError, np.linalg.LinAlgError) as exc:
                failures.append({"imputation": m, "replicate": rep, "reason": str(exc)})
        if len(samples) < max(2, int(np.ceil(0.9*b))):
            return curve, pd.DataFrame(), pd.DataFrame(), {
                "status": "BOOTSTRAP_UNSTABLE", "failures": failures, "valid": len(samples), "requested": b}
        samples = np.asarray(samples)
        eps = 1e-10
        def logit(z, clip=eps):
            return np.log(np.clip(z, clip, 1-clip)/(1-np.clip(z, clip, 1-clip)))
        point.append(logit(estimates))
        variances.append(np.var(logit(samples), axis=0, ddof=1))
        rdpoint.append([estimates[j]-estimates[i] for _, i, j in pair_keys])
        rdvar.append([np.var(samples[:, j]-samples[:, i], ddof=1) for _, i, j in pair_keys])
    ci = []
    for j, row in enumerate(landmarks.to_dict("records")):
        r = pool_scalar([v[j] for v in point], [v[j] for v in variances])
        expit = lambda z: float(1/(1+np.exp(-z)))
        ci.append({**row, "risk": expit(r["estimate"]), "lower": expit(r["lower"]),
                   "upper": expit(r["upper"]), "m": r["m"]})
    differences = []
    for j, (w, a, z) in enumerate(pair_keys):
        r = pool_scalar([v[j] for v in rdpoint], [v[j] for v in rdvar])
        differences.append({"wmh_quantile": w, "hcy_low": landmarks.loc[a, "hcy"],
                            "hcy_high": landmarks.loc[z, "hcy"], **r})
    return curve, pd.DataFrame(ci), pd.DataFrame(differences), {
        "status": "ESTIMATED", "bootstrap_per_imputation": b, "failures": failures,
        "method": "MI then within-completed-data patient bootstrap; Rubin pooling. Approximate congeniality."}
