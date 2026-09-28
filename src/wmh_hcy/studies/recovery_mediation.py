"""G-computation of exploratory recurrence-process effects, with a strict data gate.

The gate deliberately declines mediation when complete interval histories or
the first-stroke observation window are unavailable. Association fits remain
valid outputs in that case. No exact long-term death dates are constructed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit

from ..common import DataError, dump_json
from .design import StudyDesign
from .imputation import impute
from .registry import primary_spec

UPPER = np.array([365, 730, 1095, 1460, 1825], float)


def vital_state(data: pd.DataFrame, month: int) -> np.ndarray:
    state = pd.to_numeric(data[f"state{month}"], errors="coerce")
    flag = pd.to_numeric(data.get(f"death{month}", pd.Series(np.nan, index=data.index)), errors="coerce")
    conflict = (flag.eq(1) & state.isin((0, 1))) | (flag.eq(2) & state.eq(2))
    if conflict.any():
        raise DataError(f"Contradictory death/functional status at month {month}")
    return np.where(flag.eq(1) | state.eq(2), 2,
                    np.where(flag.eq(2) | state.isin((0, 1)), 0, np.nan))


def feasibility(data: pd.DataFrame, chronology: dict) -> dict:
    try:
        missing = {f"vital{m}": int(np.isnan(vital_state(data, m)).sum()) for m in (12, 24, 36, 48, 60)}
        vital_error = None
    except DataError as exc:
        missing, vital_error = {}, str(exc)
    findings = {
        "interval_state_missing": missing,
        "early_censor": int((data.event_type.eq(0) & data.exit.lt(1825)).sum()),
        "unresolved_loss": chronology.get("unresolved_loss_n", int((data.event_type.eq(0) & data.exit.lt(1825)).sum())),
        "entry_at_or_after_first_interval": int(data.entry.ge(365).sum()),
        "ambiguous_event_death_order": chronology["death_interval_ambiguous_event_n"],
        "event_after_known_death": chronology["event_after_known_death_n"],
        "events": int(data.event_type.sum()),
        "dependent": int(data.state60.eq(1).sum()),
        "deaths": int(data.state60.eq(2).sum()),
    }
    reasons = []
    if vital_error:
        reasons.append(vital_error)
    if any(missing.values()):
        reasons.append("Unknown annual vital state/death interval; no vital-state imputation authorized")
    if data.state60.notna().sum() < 100:
        reasons.append("Too few known five-year functional outcomes")
    if findings["unresolved_loss"]:
        reasons.append("First-stroke observation ends before five years without a confirmed death interval")
    if findings["entry_at_or_after_first_interval"]:
        reasons.append("Some actual three-month visits occur after the first interval endpoint")
    if findings["ambiguous_event_death_order"] or findings["event_after_known_death"]:
        reasons.append("Recurrence and death order cannot be resolved from interval data")
    if min(findings["events"], findings["dependent"], findings["deaths"]) < 20:
        reasons.append("Too few events or final states for transition models")
    return {"status": "READY" if not reasons else "NOT_ESTIMABLE", "reasons": reasons,
            "checks": findings,
            "scope": "Exploratory interventional recurrence-process decomposition; not a proven causal mechanism"}


def interval_rows(data: pd.DataFrame, design: StudyDesign) -> tuple[np.ndarray, ...]:
    x = design.transform(data).to_numpy()
    event_day = data.y5_stroke_day.to_numpy(float)
    event = data.event_type.to_numpy(bool)
    states = np.column_stack([vital_state(data, m).astype(int) for m in (12, 24, 36, 48, 60)])
    death_x, death_y, recurrence_x, recurrence_y = [], [], [], []
    alive = np.ones(len(data), bool)
    recurred = np.zeros(len(data), bool)
    for j, upper in enumerate(UPPER):
        period = np.eye(len(UPPER)-1)[np.full(len(data), j)] if j < len(UPPER)-1 else np.zeros((len(data), len(UPPER)-1))
        # Period 5 is the reference; fixed coding is also used for standardization.
        death_now = states[:, j] == 2
        stroke_now = event & (event_day <= upper) & (event_day > (data.entry.to_numpy() if j == 0 else UPPER[j-1]))
        if np.any(alive & stroke_now & death_now):
            raise DataError("Event and death share an unresolved visit interval")
        death_x.append(np.column_stack([x[alive], recurred[alive].astype(float), period[alive]]))
        death_y.append(death_now[alive].astype(int))
        at_risk = alive & ~recurred & ~death_now
        recurrence_x.append(np.column_stack([x[at_risk], period[at_risk]]))
        recurrence_y.append(stroke_now[at_risk].astype(int))
        recurred |= stroke_now
        alive &= ~death_now
    observed = data.state60.isin((0, 1)).to_numpy()[alive]
    final_x = np.column_stack([x[alive], recurred[alive].astype(float)])
    return (np.concatenate(recurrence_x), np.concatenate(recurrence_y),
            np.concatenate(death_x), np.concatenate(death_y), final_x,
            observed, data.state60.to_numpy()[alive][observed].astype(int))


def _glm(x: np.ndarray, y: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    if np.linalg.matrix_rank(x) != x.shape[1] or np.unique(y).size != 2:
        raise DataError("Transition design rank/class failure; no terms removed")
    model = sm.GLM(y, x, family=sm.families.Binomial(), freq_weights=weights).fit(maxiter=150)
    if not model.converged or not np.isfinite(model.params).all():
        raise DataError("Transition GLM did not converge")
    if np.max(np.abs(model.params)) > 20:
        raise DataError("Transition GLM has unstable near-separation coefficients")
    return np.asarray(model.params)


def transition_models(data: pd.DataFrame, design: StudyDesign) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    recurrence_x, recurrence_y, death_x, death_y, final_x, observed, dependence_y = interval_rows(data, design)
    recurrence = _glm(recurrence_x, recurrence_y)
    death = _glm(death_x, death_y)
    if observed.all():
        weights = np.ones(observed.sum())
    else:
        observation = _glm(final_x, observed.astype(int))
        probabilities = expit(final_x @ observation)
        if probabilities.min() < .05:
            raise DataError("Five-year outcome observation positivity failure")
        weights = 1/probabilities[observed]
    dependence = _glm(final_x[observed], dependence_y, weights=weights)
    return recurrence, death, dependence


def gformula(x_outcome: np.ndarray, x_recurrence: np.ndarray,
             coefficients: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    """Exact probability recursion over no event, post-event and death states."""
    rec_beta, death_beta, dep_beta = coefficients
    n = len(x_outcome)
    never = np.ones(n)
    recurrent = np.zeros(n)
    dead = np.zeros(n)
    for j in range(len(UPPER)):
        period = np.zeros((n, len(UPPER)-1))
        if j < len(UPPER)-1:
            period[:, j] = 1
        death0 = expit(np.column_stack([x_outcome, np.zeros(n), period]) @ death_beta)
        death1 = expit(np.column_stack([x_outcome, np.ones(n), period]) @ death_beta)
        hazard = expit(np.column_stack([x_recurrence, period]) @ rec_beta)
        dead += never*death0 + recurrent*death1
        recurrent = recurrent*(1-death1) + never*(1-death0)*hazard
        never *= (1-death0)*(1-hazard)
    dep0 = expit(np.column_stack([x_outcome, np.zeros(n)]) @ dep_beta)
    dep1 = expit(np.column_stack([x_outcome, np.ones(n)]) @ dep_beta)
    dependent = never*dep0 + recurrent*dep1
    independent = never*(1-dep0) + recurrent*(1-dep1)
    result = np.array([independent.mean(), dependent.mean(), dead.mean()])
    if not np.isclose(result.sum(), 1, atol=1e-8):
        raise DataError("G-formula probabilities do not sum to one")
    return result


def one_completed(data: pd.DataFrame, design: StudyDesign, contrasts: dict) -> dict:
    coefficients = transition_models(data, design)
    estimates = {}
    for exposure, (high_risk, low_risk) in contrasts.items():
        x_high = design.transform(data.assign(**{exposure: high_risk})).to_numpy()
        x_low = design.transform(data.assign(**{exposure: low_risk})).to_numpy()
        # Outcome and death processes stay at high-risk exposure; only the
        # stochastic first-recurrence process is shifted to the low-risk level.
        high = gformula(x_high, x_high, coefficients)
        shifted = gformula(x_high, x_low, coefficients)
        estimates[exposure] = {"high_risk_regime": high, "recurrence_shifted_regime": shifted,
                               "indirect_probability_difference": high-shifted}
    return estimates


def estimate(data: pd.DataFrame, directory, settings: dict, gate: dict) -> dict:
    if gate["status"] != "READY":
        dump_json(directory / "mediation.json", gate)
        return gate
    design = StudyDesign.freeze(data, primary_spec("recovery").variant("mediation_design"))
    contrasts = {"wmh_ml": tuple(np.quantile(data.wmh_ml, [.75, .25])),
                 "gm119_ml": tuple(np.quantile(data.gm119_ml, [.25, .75]))}
    completed, imputation = impute(data, design, settings)
    estimates = [one_completed(d, design, contrasts) for d in completed]
    point = {name: np.mean([e[name]["indirect_probability_difference"] for e in estimates], axis=0)
             for name in contrasts}
    rng = np.random.default_rng(int(settings.get("seed", 20260917))+1907)
    bootstrap = int(settings.get("mediation_bootstrap", 100))
    if bootstrap < 50:
        raise DataError("Mediation interval requires at least 50 patient bootstrap replicates")
    boot_values = {name: [] for name in contrasts}
    failures = []
    for rep in range(bootstrap):
        sampled = data.iloc[rng.integers(0, len(data), len(data))].reset_index(drop=True)
        try:
            # Re-freeze and re-impute inside each patient resample.
            sampled_design = StudyDesign.freeze(sampled, design.spec)
            samples, _ = impute(sampled, sampled_design, {**settings, "seed": int(settings.get("seed", 20260917))+rep+1},
                                progress=lambda *_: None)
            values = [one_completed(d, sampled_design, contrasts) for d in samples]
            for name in contrasts:
                boot_values[name].append(np.mean([v[name]["indirect_probability_difference"] for v in values], axis=0))
        except (DataError, ValueError, np.linalg.LinAlgError) as exc:
            failures.append({"replicate": rep, "reason": str(exc)})
    if any(len(v) < max(50, int(.8*bootstrap)) for v in boot_values.values()):
        result = {"status": "NOT_ESTIMABLE", "reasons": ["Patient bootstrap had insufficient successful refits"],
                  "bootstrap_failures": failures, "gate": gate}
        dump_json(directory / "mediation.json", result)
        return result
    rows = []
    for name in contrasts:
        samples = np.asarray(boot_values[name])
        for index, label in enumerate(("independent", "dependent", "dead")):
            rows.append({"exposure": name, "outcome": label, "estimate": point[name][index],
                         "lower": np.quantile(samples[:, index], .025),
                         "upper": np.quantile(samples[:, index], .975),
                         "scale": "probability_difference"})
        adverse = samples[:, 1]+samples[:, 2]
        rows.append({"exposure": name, "outcome": "dependent_or_dead",
                     "estimate": point[name][1]+point[name][2],
                     "lower": np.quantile(adverse, .025), "upper": np.quantile(adverse, .975),
                     "scale": "probability_difference"})
    pd.DataFrame(rows).to_csv(directory / "mediation_effects.csv", index=False)
    result = {"status": "ESTIMATED_EXPLORATORY", "contrasts": contrasts, "imputation": imputation,
              "bootstrap_attempts": bootstrap, "bootstrap_success": len(boot_values["wmh_ml"]),
              "bootstrap_failures": failures,
              "interpretation": "Interventional recurrence-process contrast under full temporal/confounding/model assumptions; not causal proof or proportion mediated"}
    dump_json(directory / "mediation.json", result)
    return result
