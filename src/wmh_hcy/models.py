"""Frozen Python regressions and explicit convergence/estimability diagnostics."""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from .common import DataError
from .design import Design
from .imputation import pool_scalar


@dataclass
class CoxFit:
    columns: list[str]
    params: np.ndarray
    covariance: np.ndarray
    times: np.ndarray
    increments: np.ndarray
    diagnostics: dict


def fit_cause(data: pd.DataFrame, spec: Design, cause: int) -> CoxFit:
    x = spec.transform(data)
    status = data.event_type.eq(cause).astype(int)
    events = int(status.sum())
    if events == 0:
        raise DataError(f"No events for cause {cause}; Cox model not estimable")
    rank = int(np.linalg.matrix_rank(x.to_numpy()))
    if rank < x.shape[1]:
        raise DataError("Rank-deficient model matrix; no automatic covariate deletion")
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        model = PHReg(data.exit.to_numpy(), x.to_numpy(), status=status.to_numpy(),
                      entry=np.nextafter(data.entry.to_numpy(float), np.inf), ties="efron", missing="raise")
        result = model.fit(maxiter=150, disp=False)
    if any(issubclass(w.category, ConvergenceWarning) for w in captured):
        raise DataError("Cox convergence failure: " + "; ".join(str(w.message) for w in captured))
    params = np.asarray(result.params)
    cov = np.asarray(result.cov_params())
    if not np.isfinite(params).all() or not np.isfinite(cov).all() or (np.diag(cov) < 0).any():
        raise DataError("Invalid Cox coefficients or covariance")
    eta = x.to_numpy() @ params
    risk_score = np.exp(eta)
    times, counts = np.unique(data.loc[status.eq(1), "exit"], return_counts=True)
    jumps = []
    for t, n in zip(times, counts, strict=True):
        at_risk = data.entry.lt(t) & data.exit.ge(t)
        jumps.append(float(n / risk_score[at_risk].sum()))
    residuals = np.asarray(result.schoenfeld_residuals)[status.eq(1)]
    log_t = np.log1p(data.loc[status.eq(1), "exit"])
    ph = {}
    for c in ["H", "W", "H_x_W"]:
        if c not in x:
            continue
        j = list(x).index(c)
        rho, p = stats.spearmanr(log_t, residuals[:, j])
        ph[c] = {"rho": float(rho), "p_descriptive": float(p)}
    scores = np.asarray(result.score_residuals)
    scores = np.nan_to_num(scores)
    influence = scores @ cov
    max_influence = np.max(np.abs(influence) / np.sqrt(np.diag(cov)), axis=0)
    standardized = (x - x.mean()) / x.std(ddof=0)
    condition = float(np.linalg.cond(standardized.to_numpy()))
    diag = {
        "n": len(data), "events": events, "parameters": x.shape[1], "rank": rank,
        "events_per_parameter_descriptive": events/x.shape[1], "scaled_condition_number": condition,
        "max_abs_gradient": float(np.max(np.abs(model.score(params)))),
        "max_standardized_case_influence": dict(zip(x.columns, max_influence, strict=True)),
        "schoenfeld_time_correlations": ph,
        "ph_note": "Descriptive diagnostic; not the Grambsch-Therneau omnibus test. See early/late fits.",
        "warnings": [str(w.message) for w in captured],
        "ties": "efron", "entry_convention": "(entry,exit] implemented with nextafter entry",
    }
    return CoxFit(list(x), params, cov, times, np.asarray(jumps), diag)


def pool_coefficients(fits: list[CoxFit], name: str) -> pd.DataFrame:
    rows = []
    for j, term in enumerate(fits[0].columns):
        r = pool_scalar([f.params[j] for f in fits], [f.covariance[j, j] for f in fits])
        r.update({"analysis": name, "term": term, "HR": np.exp(r["estimate"]),
                  "HR_lower": np.exp(r["lower"]), "HR_upper": np.exp(r["upper"])})
        rows.append(r)
    return pd.DataFrame(rows)


def cross_sectional(data: pd.DataFrame, spec: Design) -> pd.DataFrame:
    x = spec.transform(data).drop(columns=["W", "W_rcs", "H_x_W"])
    x = sm.add_constant(x, has_constant="add")
    fit = sm.OLS(np.log1p(data.wmh_ml), x, missing="raise").fit(cov_type="HC3")
    ci = fit.conf_int()
    table = pd.DataFrame({"term": x.columns, "estimate": fit.params.to_numpy(),
                         "variance": np.diag(fit.cov_params()), "lower": ci[0].to_numpy(),
                         "upper": ci[1].to_numpy(), "p": fit.pvalues.to_numpy()})
    table.attrs["covariance"] = np.asarray(fit.cov_params())
    return table


def ordinal(data: pd.DataFrame, spec: Design, outcome: str = "mrs12") -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    x = spec.transform(data)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        fit = OrderedModel(data[outcome].astype(int), x, distr="logit").fit(method="bfgs", maxiter=500, disp=False)
    if not fit.mle_retvals.get("converged", False):
        raise DataError("Ordinal logistic model did not converge")
    ncoef = x.shape[1]
    estimates = np.asarray(fit.params)[:ncoef]
    covariance = np.asarray(fit.cov_params())[:ncoef, :ncoef]
    table = pd.DataFrame({"term": list(x), "estimate": estimates, "variance": np.diag(covariance)})
    thresholds = []
    for cut in range(6):
        y = data[outcome].gt(cut).astype(int)
        if y.nunique() < 2:
            continue
        try:
            binary = sm.GLM(y, sm.add_constant(x), family=sm.families.Binomial()).fit()
            for term in ["H", "W", "H_x_W"]:
                thresholds.append({"threshold": cut, "term": term,
                                   "estimate": float(binary.params[term]), "se": float(binary.bse[term])})
        except (ValueError, np.linalg.LinAlgError) as exc:
            thresholds.append({"threshold": cut, "failure": str(exc)})
    diag = {"n": len(data), "converged": True, "outcome": outcome,
            "mrs_distribution": data[outcome].value_counts().to_dict(),
            "warnings": [str(w.message) for w in captured],
            "proportional_odds_check": "Inspect threshold-specific coefficients; no automated claim of PO validity"}
    return table, diag, pd.DataFrame(thresholds)


def pool_regression(tables: list[pd.DataFrame], name: str, exponentiate: bool = False) -> pd.DataFrame:
    rows = []
    for term in tables[0].term:
        sub = [t.set_index("term").loc[term] for t in tables]
        r = pool_scalar([v.estimate for v in sub], [v.variance for v in sub])
        r.update({"analysis": name, "term": term})
        if exponentiate:
            r.update({"OR": np.exp(r["estimate"]), "OR_lower": np.exp(r["lower"]),
                      "OR_upper": np.exp(r["upper"])})
        rows.append(r)
    return pd.DataFrame(rows)
