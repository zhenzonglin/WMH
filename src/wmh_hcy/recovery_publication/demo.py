"""Deterministic synthetic saved artifacts for visual style previews only.

No model, imputation, or inferential procedure is run. Coefficients, intervals,
and P values are explicitly invented display fixtures, unrelated to the study.
``ESTIMATED`` and ``converged`` are schema tokens required by the saved-artifact
reader; they do not claim that these illustrative data were fitted.
"""
from __future__ import annotations

import json
from itertools import pairwise, product
from pathlib import Path

import numpy as np
import pandas as pd

from .inputs import COUNT_TRANSITIONS, MONTHS, PRIMARY_TERMS, STATES, TRANSITIONS

NOTICE = (
    "SYNTHETIC STYLE PREVIEW ONLY. Invented demonstration data, estimates, confidence "
    "intervals, P values, and fit-status tokens; no research results and no models fitted. "
    "The illustrative coefficient/probability fixtures were not estimated from the "
    "synthetic descriptive cohort. Do not use these files for scientific inference."
)
N = 120
M = 2


def _json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")


def _run(root: Path, name: str) -> Path:
    run = root / name / "runs" / "SYNTHETIC_STYLE_PREVIEW"
    marker = run / "SYNTHETIC_NOTICE.txt"
    if run.exists() and (not marker.is_file() or marker.read_text(encoding="utf-8") != NOTICE):
        raise ValueError("Refusing to overwrite an unrecognized synthetic demo directory")
    run.mkdir(parents=True, exist_ok=True)
    marker.write_text(NOTICE, encoding="utf-8")
    return run


def _state_history() -> np.ndarray:
    """Construct illustrative visits with recovery, missingness, and absorbing death."""
    states = np.zeros((N, len(MONTHS)))
    assignments = (
        (1, slice(0, 8), slice(112, 116), (slice(116, 118),)),
        (2, slice(3, 12), slice(110, 116), (slice(116, 120),)),
        (3, slice(6, 16), slice(108, 116), (slice(106, 108), slice(116, 120))),
        (4, slice(8, 20), slice(104, 116), (slice(100, 104), slice(116, 120))),
        (5, slice(12, 26), slice(100, 116), (slice(88, 100),)),
    )
    for visit, dependent, dead, unknown in assignments:
        states[dependent, visit] = 1
        states[dead, visit] = 2
        for selection in unknown:
            states[selection, visit] = np.nan
    return states


def _cohort(states: np.ndarray) -> pd.DataFrame:
    rng = np.random.default_rng(70133)
    data = pd.DataFrame({
        "patient_id": [f"SYNTHETIC_STYLE_{i:04d}" for i in range(N)],
        "state60": states[:, -1],
        "age": np.rint(np.clip(rng.normal(61, 10, N), 30, 86)),
        "wmh_ml": np.round(np.clip(rng.lognormal(2.25, .65, N), .5, 55), 3),
        "gm119_ml": np.round(np.clip(rng.normal(548, 58, N), 370, 715), 3),
        "icv_ml": np.round(np.clip(rng.normal(1375, 112, N), 1080, 1710), 3),
        "lesion_ml": np.round(np.clip(rng.lognormal(1.5, 1.1, N), .05, 95), 3),
        "nihss": np.clip(rng.poisson(4, N), 0, 18),
        "sex": rng.choice([1, 2], N, p=[.62, .38]),
        "smoking": rng.choice([1, 2, 3, 4], N, p=[.48, .06, .32, .14]),
        "drinking": rng.choice([1, 2, 3, 4], N, p=[.55, .22, .17, .06]),
        "hypertension": rng.choice([1, 2], N, p=[.34, .66]).astype(float),
        "diabetes": rng.choice([1, 2], N, p=[.74, .26]).astype(float),
        "prior_stroke": rng.choice([1, 2], N, p=[.79, .21]).astype(float),
        "education": rng.choice([1, 2, 3, 4, 5], N, p=[.12, .24, .30, .26, .08]).astype(float),
        "pre_mrs": rng.choice(range(6), N, p=[.70, .18, .06, .03, .02, .01]),
        "toast": rng.choice(range(1, 6), N, p=[.28, .07, .27, .03, .35]),
        "mrs3": np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2])[np.arange(N) % 10],
    })
    for column, count in (("hypertension", 3), ("diabetes", 2), ("prior_stroke", 2), ("education", 15)):
        data.loc[rng.choice(N, count, replace=False), column] = np.nan
    return data


def _effect(term: str, ratio: float, low: float, high: float, p: float) -> dict:
    return {"term": term, "estimate": float(np.log(ratio)), "lower": float(np.log(low)),
            "upper": float(np.log(high)), "se": float((np.log(high) - np.log(low)) / 3.92),
            "ratio": ratio, "ratio_lower": low, "ratio_upper": high, "p": p, "m": M,
            "scale": "SYNTHETIC DISPLAY FIXTURE ONLY"}


def _probability_rows(wmh: float, gm: float, probabilities: tuple[float, float, float],
                      **keys: object) -> list[dict]:
    rows = []
    for state, probability in zip(STATES, probabilities, strict=True):
        width = .025 if state == "independent" else .015
        rows.append({**keys, "wmh_ml": wmh, "gm119_ml": gm, "state": state,
                     "probability": probability, "lower": max(0., probability - width),
                     "upper": min(1., probability + width), "ci": NOTICE, "ci_method": NOTICE})
    return rows


def _recovery(run: Path, data: pd.DataFrame) -> None:
    data.to_csv(run / "eligible.csv", index=False)
    observed = data.loc[data.state60.notna()]
    observed_n = len(observed)
    _json(run / "status.json", {"study": "recovery", "run": run.name, "path": str(run),
          "mode": "synthetic", "status": "COMPLETED", "synthetic_notice": NOTICE})
    _json(run / "audit.json", {"eligible_n": N, "outcome_observed_n": observed_n,
          "independent": int(data.state60.eq(0).sum()), "dependent": int(data.state60.eq(1).sum()),
          "dead": int(data.state60.eq(2).sum()), "unknown": int(data.state60.isna().sum()),
          "clinical_n": 240, "synthetic_notice": NOTICE})
    flow = [("clinical_start", 240, 0), ("adult_ischemic_stroke", 220, 20),
            ("available_true_icv", 160, 60), ("alive_at_month3", 145, 15),
            ("independent_at_month3", N, 25),
            ("known_five_year_functional_state", observed_n, N - observed_n)]
    pd.DataFrame(flow, columns=["step", "remaining", "excluded_here"]).to_csv(
        run / "cohort_flow.csv", index=False)
    examples = [
        ("primary", observed_n, .024, None),
        ("complete_case", 87, None, "Synthetic example: nonfinite saved model estimates"),
        ("month12", 118, None, "Synthetic example: saved model did not converge"),
        ("month24", 116, None, "Synthetic example: saved model did not converge"),
        ("month36", 114, .043, None), ("month48", 112, .031, None),
        ("actual_improvement", 61, None, "Synthetic example: saved model did not converge"),
        ("raw_wmh", observed_n, .028, None),
        ("dependent_or_dead", observed_n, .018, None),
        ("observation_weighted", observed_n, .027, None),
    ]
    base = ((1.30, .98, 1.72, .07), (.88, .70, 1.11, .28),
            (1.18, .94, 1.48, .15), (.91, .73, 1.14, .42))
    results = []
    for index, (name, n, p, reason) in enumerate(examples):
        status = "ESTIMATED" if reason is None else "NOT_ESTIMABLE"
        result = {"analysis": name, "status": status, "n": n,
                  "family": "synthetic_style_preview", "tier": "primary" if name == "primary" else "secondary",
                  "reason": reason, "synthetic_notice": NOTICE}
        if reason is None:
            result.update({"p": p, "imputations": M})
        _json(run / name / "result.json", result)
        results.append(result)
        if reason is not None:
            continue
        if name == "dependent_or_dead":
            effects = [_effect("wmh_ml", 1.27, 1.02, 1.58, .032),
                       _effect("gm119_ml", .86, .69, 1.07, .18)]
        else:
            factor = 1. + index * .006
            effects = [_effect(term, ratio * factor, low * factor, high * factor, value)
                       for term, (ratio, low, high, value) in zip(PRIMARY_TERMS, base, strict=True)]
        pd.DataFrame(effects).to_csv(run / name / "coefficients.csv", index=False)
    pd.DataFrame(results).to_csv(run / "results.csv", index=False)
    _json(run / "primary" / "frozen_design.json", {
        "wmh_ml": {"scale": float(np.log1p(observed.wmh_ml).std(ddof=1))},
        "gm119_ml": {"scale": float(observed.gm119_ml.std(ddof=1))}, "synthetic_notice": NOTICE})
    rows = []
    for (wi, wmh), (gi, gm) in product(enumerate(observed.wmh_ml.quantile([.25, .5, .75])),
                                     enumerate(observed.gm119_ml.quantile([.25, .5, .75]))):
        dependent = .065 + .022 * wi - .012 * gi
        dead = .095 + .018 * wi - .014 * gi
        rows.extend(_probability_rows(float(wmh), float(gm), (1. - dependent - dead, dependent, dead)))
    pd.DataFrame(rows).to_csv(run / "primary" / "standardized_states.csv", index=False)


def _trajectory(run: Path, data: pd.DataFrame, states: np.ndarray) -> None:
    counts = []
    for index, (start, stop) in enumerate(pairwise(MONTHS)):
        origin, destination = states[:, index], states[:, index + 1]
        known = np.isfinite(origin) & np.isfinite(destination)
        at_risk = known & np.isin(origin, [0, 1])
        row = {"interval_index": index, "from_month": start, "to_month": stop,
               "eligible_start_n": int(np.isin(origin, [0, 1]).sum()),
               "known_adjacent_state_n": int(known.sum()), "transition_rows_n": int(at_risk.sum()),
               "unknown_start_n": int(np.isnan(origin).sum()), "unknown_end_n": int(np.isnan(destination).sum()),
               "conflict_n": 0, "death_start_n": int((origin == 2).sum()),
               "death_to_alive_n": int(((origin == 2) & np.isin(destination, [0, 1])).sum()),
               "death_absorbing_check": True}
        for transition in COUNT_TRANSITIONS:
            before, after = transition.split("_to_")
            row[transition + "_n"] = int(((origin == STATES.index(before)) &
                                           (destination == STATES.index(after))).sum())
        counts.append(row)
    count_frame = pd.DataFrame(counts)
    total = int(count_frame.transition_rows_n.sum())
    count_frame.to_csv(run / "transition_counts.csv", index=False)
    _json(run / "status.json", {"study": "recovery_trajectory", "run": run.name, "path": str(run),
          "mode": "synthetic", "status": "COMPLETED", "analysis_scope": "trajectory_only",
          "analysis": {"status": "ESTIMATED", "imputations": M, "transition_rows": total,
                       "trajectory_probability_draws": 16}, "synthetic_notice": NOTICE})
    _json(run / "audit.json", {"base_eligible_n": N, "transition_rows_n": total,
          "state_counts": {str(month): {**{str(state): int((states[:, index] == state).sum()) for state in range(3)},
                                         "unknown": int(np.isnan(states[:, index]).sum())}
                           for index, month in enumerate(MONTHS)},
          "future_recurrence_fields_requested": [], "synthetic_notice": NOTICE})
    base = ((1.26, 1.04, 1.53, .020), (.84, .67, 1.05, .13),
            (1.16, .94, 1.43, .17), (1.08, .82, 1.42, .58),
            (.89, .73, 1.09, .26), (1.32, 1.03, 1.69, .029),
            (.90, .73, 1.11, .32), (.95, .73, 1.24, .70))
    effects = []
    for (exposure, transition), values in zip(product(("wmh_ml", "gm119_ml"), TRANSITIONS), base, strict=True):
        effects.append({**_effect(transition, *values), "exposure": exposure, "transition": transition})
    pd.DataFrame(effects).to_csv(run / "transition_effects.csv", index=False)
    wmh_values, gm_values = data.wmh_ml.quantile([.25, .5, .75]), data.gm119_ml.quantile([.25, .5, .75])
    rows = []
    for exposure, quantile_index, month in product(("wmh", "gm"), range(3), (36, 48, 60)):
        time = (month - 36) / 12
        contrast = quantile_index - 1
        dependent = .058 + .005 * time + (.018 if exposure == "wmh" else -.011) * contrast
        dead = .060 + .023 * time + (.011 if exposure == "wmh" else -.012) * contrast
        wmh = float(wmh_values.iloc[quantile_index if exposure == "wmh" else 1])
        gm = float(gm_values.iloc[quantile_index if exposure == "gm" else 1])
        rows.extend(_probability_rows(wmh, gm, (1. - dependent - dead, dependent, dead),
                                     scenario=f"{exposure}_q{(25, 50, 75)[quantile_index]}", month=month))
    pd.DataFrame(rows).to_csv(run / "state_probabilities_36_48_60.csv", index=False)
    rows = []
    for index, month in enumerate(MONTHS):
        state = states[:, index]
        known_n = int(np.isfinite(state).sum())
        scores = np.array(data.mrs3, dtype=float)
        dependent = state == 1
        scores[dependent] = np.array([3, 3, 3, 4, 4, 5])[np.flatnonzero(dependent) % 6]
        scores[state == 2] = 6
        scores[np.isnan(state)] = np.nan
        for mrs in (*map(str, range(7)), "unknown"):
            n = int(np.isnan(scores).sum()) if mrs == "unknown" else int((scores == int(mrs)).sum())
            rows.append({"month": month, "mrs": mrs, "n": n, "denominator_known": known_n,
                         "proportion": np.nan if mrs == "unknown" else n / known_n})
    pd.DataFrame(rows).to_csv(run / "mrs_trajectory_distribution.csv", index=False)
    _json(run / "fit_diagnostics.json", [{"converged": True, "n": total, "synthetic_notice": NOTICE}
                                          for _ in range(M)])


def create_demo(root: Path) -> tuple[Path, Path]:
    """Write two unmistakably synthetic saved runs; return recovery and trajectory paths."""
    root = Path(root).resolve()
    recovery = _run(root, "synthetic_recovery")
    trajectory = _run(root, "synthetic_trajectory")
    states = _state_history()
    data = _cohort(states)
    _recovery(recovery, data)
    _trajectory(trajectory, data, states)
    return recovery, trajectory
