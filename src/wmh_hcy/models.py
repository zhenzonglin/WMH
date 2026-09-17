"""Frozen Python regressions and explicit convergence/estimability diagnostics."""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg, PHRegResults
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


def fit_phreg_checked(model: PHReg, x: np.ndarray, events: int) -> tuple[PHRegResults, dict]:
    """Same unpenalized likelihood; scaled BFGS initialization only if Newton fails."""
    scale = x.std(axis=0)
    if (scale <= 0).any():
        raise DataError("Rank-deficient centered Cox matrix; no automatic covariate deletion")

    def validate(result):
        beta = np.asarray(result.params)
        cov = np.asarray(result.cov_params())
        if not np.isfinite(beta).all() or not np.isfinite(cov).all() or (np.diag(cov) <= 0).any():
            raise DataError("Invalid Cox coefficients or covariance")
        score = model.score(beta) / scale
        information = -model.hessian(beta) / np.outer(scale, scale)
        if not np.isfinite(score).all() or not np.isfinite(information).all():
            raise DataError("Non-finite Cox score or information matrix")
        eigenvalues = np.linalg.eigvalsh((information + information.T) / 2)
        if eigenvalues[0] <= 0 or eigenvalues[-1] / eigenvalues[0] > 1e12:
            raise DataError("Singular or ill-conditioned Cox information matrix")
        scaled_score = float(np.max(np.abs(score)) / events)
        if scaled_score > 1e-6:
            raise DataError("Cox convergence failure: standardized score per event exceeds 1e-6")
        return {"max_scaled_score_per_event": scaled_score,
                "scaled_information_condition": float(eigenvalues[-1] / eigenvalues[0])}

    attempts = []
    for method in ("newton", "scaled_bfgs_then_newton"):
        try:
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter("always")
                if method == "newton":
                    result = model.fit(maxiter=150, disp=False)
                else:
                    # Centering/scaling are invertible reparameterizations of the Cox likelihood.
                    z = (x - x.mean(axis=0)) / scale
                    scaled_model = PHReg(model.endog, z, status=model.status, entry=model.entry,
                                         ties=model.ties, missing="raise")
                    initial = scaled_model.fit(method="bfgs", maxiter=1000, gtol=1e-8, disp=False)
                    if not np.isfinite(initial.params).all():
                        raise DataError("Invalid Cox initialization")
                    # A warning from the initializer is allowed only if polishing passes every check.
                    captured.clear()
                    polished = scaled_model.fit(start_params=initial.params, method="newton",
                                                maxiter=150, disp=False)
                    result = PHRegResults(model, np.asarray(polished.params) / scale,
                                          np.asarray(polished.cov_params()) / np.outer(scale, scale))
                if any(issubclass(w.category, ConvergenceWarning) for w in captured):
                    raise DataError("Cox convergence failure")
                diagnostics = validate(result)
            diagnostics.update(optimizer=method, failed_attempts=attempts,
                               warnings=[str(w.message) for w in captured])
            return result, diagnostics
        except (DataError, ValueError, np.linalg.LinAlgError) as exc:
            attempts.append({"optimizer": method, "reason": str(exc)})
    raise DataError("Cox fit failed with both solvers: " + "; ".join(r["reason"] for r in attempts))


def fit_cause(data: pd.DataFrame, spec: Design, cause: int) -> CoxFit:
    x = spec.transform(data)
    status = data.event_type.eq(cause).astype(int)
    events = int(status.sum())
    if events == 0:
        raise DataError(f"No events for cause {cause}; Cox model not estimable")
    rank = int(np.linalg.matrix_rank(x.to_numpy()))
    if rank < x.shape[1]:
        raise DataError("Rank-deficient model matrix; no automatic covariate deletion")
    model = PHReg(data.exit.to_numpy(), x.to_numpy(), status=status.to_numpy(),
                  entry=np.nextafter(data.entry.to_numpy(float), np.inf), ties="efron", missing="raise")
    result, numerical = fit_phreg_checked(model, x.to_numpy(), events)
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
        "warnings": numerical["warnings"], "numerical_fit": numerical,
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
