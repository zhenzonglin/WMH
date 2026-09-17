"""No variable selection, penalization, death regression or bootstrap."""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.special import softmax
from statsmodels.duration.hazard_regression import PHReg

from ..common import DataError
from ..models import fit_phreg_checked


@dataclass
class Fit:
    terms: list[str]
    params: np.ndarray
    covariance: np.ndarray
    diagnostics: dict


def fit(data, design, weights=None, observation=None):
    spec = design.spec
    frame = design.transform(data)
    x, y = frame.to_numpy(), data[spec.outcome].to_numpy(float)
    if not np.isfinite(y).all():
        raise DataError("Missing/nonfinite outcome cannot be imputed or silently dropped")
    if np.linalg.matrix_rank(x) < x.shape[1]:
        raise DataError("Rank deficient design; no automatic covariate deletion")
    variable = x.std(axis=0) > 0
    condition = np.linalg.cond((x[:, variable]-x[:, variable].mean(axis=0))/x[:, variable].std(axis=0))
    diagnostics = {"n": len(data), "parameters": x.shape[1], "design_condition": float(condition)}
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        if spec.family == "cox":
            events = int(y.sum())
            if events <= 0 or not data.exit.gt(data.entry).all():
                raise DataError("No events or invalid follow-up interval")
            model = PHReg(data.exit.to_numpy(), x, status=y, entry=np.nextafter(data.entry.to_numpy(float), np.inf),
                          ties="efron", missing="raise")
            result, numerical = fit_phreg_checked(model, x, events)
            diagnostics.update(events=events, numerical=numerical)
            residual = np.asarray(result.schoenfeld_residuals)[y == 1]
            diagnostics["PH_descriptive"] = {
                t: float(stats.spearmanr(np.log1p(data.loc[y == 1, "exit"]), residual[:, j]).pvalue)
                for j, t in enumerate(frame) if t.startswith(("sbp3", "wmh_ml"))}
            scores = np.nan_to_num(result.score_residuals)
        elif spec.family == "ols":
            result = sm.OLS(y, x, missing="raise").fit(cov_type="HC3")
            diagnostics.update(max_leverage=float(result.get_influence().hat_matrix_diag.max()),
                               max_cooks_distance=float(np.nanmax(result.get_influence().cooks_distance[0])))
            scores = None
        elif spec.family == "binary":
            if set(np.unique(y)) != {0, 1}:
                raise DataError("Binary outcome requires both observed states")
            if weights is not None:
                raise DataError("Weighted binary likelihood is not a specified study analysis")
            result = sm.Logit(y, x, missing="raise").fit(method="newton", maxiter=150, disp=False, cov_type="HC0")
            if not result.mle_retvals.get("converged", False):
                raise DataError("Binomial model failed to converge")
            scores = None
        else:
            if set(np.unique(y)) != {0, 1, 2}:
                raise DataError("Multinomial model requires independent, dependent and death states")
            if weights is not None:
                return weighted_multinomial(x, y, list(frame), weights, diagnostics, observation)
            model = sm.MNLogit(y.astype(int), x, missing="raise")
            result = model.fit(method="newton", maxiter=150, disp=False)
            if not result.mle_retvals.get("converged", False):
                raise DataError("Multinomial model failed to converge")
            scores = result.model.score_obs(result.params)
    params = np.asarray(result.params).reshape(-1, order="F")
    covariance = np.asarray(result.cov_params())
    if not np.isfinite(params).all() or not np.isfinite(covariance).all():
        raise DataError("Nonfinite model estimates")
    if np.linalg.eigvalsh(covariance).min() <= 0:
        raise DataError("Non-positive-definite model covariance")
    terms = list(frame) if spec.family != "multinomial" else [
        f"{state}:{term}" for state in ("dependent", "dead") for term in frame]
    if scores is not None:
        diagnostics["max_standardized_case_influence"] = float(
            np.max(np.abs(scores @ covariance) / np.sqrt(np.diag(covariance))))
    diagnostics.update(converged=True, parameters=len(params),
                       warnings=[str(w.message) for w in captured])
    if spec.family == "multinomial":
        score = result.model.score(params)
        diagnostics["max_score_per_patient"] = float(np.max(np.abs(score))/len(y))
        if diagnostics["max_score_per_patient"] > 1e-5:
            raise DataError("Multinomial score check failed")
    return Fit(terms, params, covariance, diagnostics)


def weighted_multinomial(x, y, terms, weights, diagnostics, observation=None):
    """IP-observation weighted likelihood; stacked sandwich includes weight estimation."""
    from scipy.optimize import minimize

    n, k = x.shape
    weights = np.asarray(weights, float)
    onehot = np.eye(3)[y.astype(int)]

    def quantities(beta):
        matrix = beta.reshape(k, 2, order="F")
        probabilities = softmax(np.column_stack([np.zeros(n), x @ matrix]), axis=1)
        scores = np.concatenate([x*(onehot[:, j]-probabilities[:, j])[:, None] for j in (1, 2)], axis=1)
        return probabilities, scores

    def objective(beta):
        p, scores = quantities(beta)
        return -np.sum(weights*np.log(np.maximum(p[np.arange(n), y.astype(int)], 1e-300))), -(weights[:, None]*scores).sum(axis=0)

    result = minimize(objective, np.zeros(2*k), jac=True, method="BFGS", options={"maxiter": 1000, "gtol": 1e-7})
    p, scores = quantities(result.x)
    if np.max(np.abs(objective(result.x)[1]))/n > 1e-6:
        raise DataError("Weighted multinomial score check failed")
    information = np.block([[x.T @ (x*(weights*p[:, j]*((j == l)-p[:, l]))[:, None])
                             for l in (1, 2)] for j in (1, 2)])
    bread = np.linalg.inv(information)
    weighted_scores = weights[:, None]*scores
    cov = bread @ weighted_scores.T @ weighted_scores @ bread
    if observation is not None:
        z, observed, probability = observation
        observed = np.asarray(observed, bool)
        a_alpha = z.T @ (z*(probability*(1-probability))[:, None])
        cross = (scores*((1-probability[observed])/probability[observed])[:, None]).T @ z[observed]
        stacked_bread = np.block([[information, cross], [np.zeros((z.shape[1], 2*k)), a_alpha]])
        beta_scores = np.zeros((len(z), 2*k))
        beta_scores[observed] = weighted_scores
        stacked_score = np.column_stack([beta_scores, z*(observed-probability)[:, None]])
        inverse = np.linalg.inv(stacked_bread)
        all_cov = inverse @ stacked_score.T @ stacked_score @ inverse.T
        cov = all_cov[:2*k, :2*k]
    if not np.isfinite(cov).all() or np.linalg.eigvalsh(cov).min() <= 0:
        raise DataError("Invalid weighted multinomial covariance")
    diagnostics.update(converged=True, parameters=2*k, observation_weights=(
                           "stacked sandwich includes observation-weight estimation" if observation is not None else
                           "fixed-weight sandwich"),
                       weight_min=float(weights.min()), weight_max=float(weights.max()),
                       weight_ess=float(weights.sum()**2/np.sum(weights**2)))
    return Fit([f"{s}:{t}" for s in ("dependent", "dead") for t in terms], result.x, cov, diagnostics)


def state_probabilities(fitted, x):
    """Mean standardized probabilities and delta-method variance, all covariance blocks."""
    x = np.asarray(x)
    k = x.shape[1]
    p = softmax(np.column_stack([np.zeros(len(x)), x @ fitted.params.reshape(k, 2, order="F")]), axis=1)
    result = []
    for s in range(3):
        gradient = np.concatenate([(x*(p[:, s]*((s == j)-p[:, j]))[:, None]).mean(axis=0) for j in (1, 2)])
        result.append((float(p[:, s].mean()), float(gradient @ fitted.covariance @ gradient)))
    return result
