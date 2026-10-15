"""Synthetic, saved-artifact contracts; no patient data or analysis fitting."""
from __future__ import annotations

import json
from itertools import pairwise, product
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from wmh_hcy.recovery_publication.inputs import (
    CATEGORICAL,
    CONTINUOUS,
    COUNT_TRANSITIONS,
    MONTHS,
    PRIMARY_TERMS,
    STATES,
    TRANSITIONS,
    load_recovery,
    load_trajectory,
    source_manifest,
)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def make_recovery(root: Path, n: int = 90) -> Path:
    """Complete explicitly synthetic recovery fixture, also usable for style demos."""
    run = root / "synthetic_recovery" / "runs" / "20261015T120000000000Z"
    run.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    data = pd.DataFrame({"patient_id": [f"SYNTHETIC_{i:04d}" for i in range(n)],
                         "state60": np.tile([0., 0., 0., 0., 0., 0., 1., 1., 2., np.nan], (n + 9) // 10)[:n]})
    for name, center in zip(CONTINUOUS, (60, 15, 560, 1350, 8, 4), strict=True):
        data[name] = center + rng.uniform(0, 3, n)
    for name in CATEGORICAL:
        data[name] = np.arange(n) % (3 if name in {"pre_mrs", "mrs3"} else 2)
        if name not in {"pre_mrs", "mrs3"}:
            data[name] += 1
    data.to_csv(run / "eligible.csv", index=False)
    observed_n = int(data.state60.notna().sum())
    write_json(run / "status.json", {"study": "recovery", "run": run.name, "path": str(run),
               "mode": "synthetic", "status": "COMPLETED"})
    write_json(run / "audit.json", {"eligible_n": n, "outcome_observed_n": observed_n,
               "independent": int(data.state60.eq(0).sum()), "dependent": int(data.state60.eq(1).sum()),
               "dead": int(data.state60.eq(2).sum()), "unknown": int(data.state60.isna().sum())})
    pd.DataFrame([{"step": "base", "remaining": n, "excluded_here": 0},
                  {"step": "known_five_year_functional_state", "remaining": observed_n,
                   "excluded_here": n - observed_n}]).to_csv(run / "cohort_flow.csv", index=False)
    pd.DataFrame([{"analysis": "primary", "status": "ESTIMATED", "n": observed_n, "p": .04},
                  {"analysis": "complete_case", "status": "NOT_ESTIMABLE", "n": observed_n,
                   "reason": "Synthetic example of a saved nonestimable analysis", "p": np.nan}]).to_csv(
                       run / "results.csv", index=False)
    write_json(run / "primary" / "result.json", {"analysis": "primary", "status": "ESTIMATED",
               "n": observed_n, "p": .04, "imputations": 2})
    pd.DataFrame([{"term": term, "ratio": 1.2, "ratio_lower": 1.01, "ratio_upper": 1.4, "p": .04}
                  for term in PRIMARY_TERMS]).to_csv(run / "primary" / "coefficients.csv", index=False)
    write_json(run / "primary" / "frozen_design.json", {"wmh_ml": {"scale": .5},
               "gm119_ml": {"scale": 60}})
    rows = []
    known = data.loc[data.state60.notna()]
    for wmh, gm, (state, probability) in product(known.wmh_ml.quantile([.25, .5, .75]),
            known.gm119_ml.quantile([.25, .5, .75]), zip(STATES, (.8, .1, .1), strict=True)):
        rows.append({"wmh_ml": wmh, "gm119_ml": gm, "state": state, "probability": probability,
                     "lower": probability - .03, "upper": probability + .03,
                     "ci": "SYNTHETIC DISPLAY FIXTURE"})
    pd.DataFrame(rows).to_csv(run / "primary" / "standardized_states.csv", index=False)
    return run


def make_trajectory(root: Path, n: int = 90) -> Path:
    """Complete synthetic saved trajectory fixture; no model or imputation run."""
    run = root / "synthetic_trajectory" / "runs" / "20261015T120001000000Z"
    run.mkdir(parents=True, exist_ok=True)
    analysis = {"status": "ESTIMATED", "imputations": 2, "transition_rows": 5 * n,
                "trajectory_probability_draws": 16}
    write_json(run / "status.json", {"study": "recovery_trajectory", "run": run.name, "path": str(run),
               "mode": "synthetic", "status": "COMPLETED", "analysis_scope": "trajectory_only",
               "analysis": analysis})
    write_json(run / "audit.json", {"base_eligible_n": n, "transition_rows_n": 5 * n,
               "state_counts": {str(month): {"0": n, "1": 0, "2": 0, "unknown": 0} for month in MONTHS},
               "future_recurrence_fields_requested": []})
    effects = [{"exposure": exposure, "transition": transition, "ratio": 1.2,
                "ratio_lower": 1.01, "ratio_upper": 1.4, "p": .04, "m": 2}
               for exposure, transition in product(("wmh_ml", "gm119_ml"), TRANSITIONS)]
    pd.DataFrame(effects).to_csv(run / "transition_effects.csv", index=False)
    rows = []
    scenarios = {"wmh_q25": (8, 560), "wmh_q50": (13, 560), "wmh_q75": (19, 560),
                 "gm_q25": (13, 520), "gm_q50": (13, 560), "gm_q75": (13, 600)}
    for (scenario, (wmh, gm)), month, (state, probability) in product(
            scenarios.items(), (36, 48, 60), zip(STATES, (.8, .1, .1), strict=True)):
        rows.append({"scenario": scenario, "month": month, "state": state, "wmh_ml": wmh, "gm119_ml": gm,
                     "probability": probability, "lower": probability - .03, "upper": probability + .03,
                     "ci_method": "SYNTHETIC DISPLAY FIXTURE"})
    pd.DataFrame(rows).to_csv(run / "state_probabilities_36_48_60.csv", index=False)
    counts = []
    for index, (start, stop) in enumerate(pairwise(MONTHS)):
        counts.append({"interval_index": index, "from_month": start, "to_month": stop,
                       "eligible_start_n": n, "known_adjacent_state_n": n, "transition_rows_n": n,
                       "unknown_start_n": 0, "unknown_end_n": 0, "conflict_n": 0, "death_start_n": 0,
                       "death_to_alive_n": 0, "death_absorbing_check": True,
                       **{name + "_n": n if name == "independent_to_independent" else 0
                          for name in COUNT_TRANSITIONS}})
    pd.DataFrame(counts).to_csv(run / "transition_counts.csv", index=False)
    rows = []
    for month, mrs in product(MONTHS, (*map(str, range(7)), "unknown")):
        rows.append({"month": month, "mrs": mrs, "n": n if mrs == "0" else 0,
                     "denominator_known": n,
                     "proportion": np.nan if mrs == "unknown" else 1. if mrs == "0" else 0.})
    pd.DataFrame(rows).to_csv(run / "mrs_trajectory_distribution.csv", index=False)
    write_json(run / "fit_diagnostics.json", [{"converged": True, "n": 5 * n} for _ in range(2)])
    return run


def modify_json(path: Path, key: str, value: object) -> None:
    data = json.loads(path.read_text(encoding="utf-8"))
    data[key] = value
    write_json(path, data)


def test_valid_readers_are_read_only_and_aggregate_contracts_are_complete(tmp_path):
    recovery = make_recovery(tmp_path)
    trajectory = make_trajectory(tmp_path)
    before = source_manifest(recovery, ["eligible.csv", "status.json"])
    main = load_recovery(recovery, allow_demo=True)
    longitudinal = load_trajectory(trajectory, allow_demo=True)
    assert "patient_id" not in main["cohort"]
    assert len(main["cohort"]) == 90 and main["cohort"].state60.notna().sum() == 81
    assert len(main["probabilities"]) == 27 and len(longitudinal["probabilities"]) == 54
    assert "complete_case" not in main["coefficients"]
    assert "nonestimable" in main["result_details"]["complete_case"]["reason"]
    assert "16 draws" in longitudinal["warnings"][0]
    assert source_manifest(recovery, ["eligible.csv", "status.json"]) == before


@pytest.mark.parametrize("key,value", [("study", "recovery_path"), ("run", "wrong"),
                                     ("path", "unrelated"), ("status", "FAILED"), ("status", "PREPARED")])
def test_recovery_requires_correct_run_identity_and_result_state(tmp_path, key, value):
    run = make_recovery(tmp_path)
    modify_json(run / "status.json", key, value)
    with pytest.raises(ValueError):
        load_recovery(run, allow_demo=True)


def test_demo_requires_opt_in_and_cannot_read_real_run(tmp_path):
    run = make_recovery(tmp_path)
    with pytest.raises(ValueError, match="mode=real"):
        load_recovery(run)
    modify_json(run / "status.json", "mode", "real")
    assert load_recovery(run)["status"]["mode"] == "real"
    with pytest.raises(ValueError, match="cannot consume a real"):
        load_recovery(run, allow_demo=True)


def test_duplicate_patients_error_does_not_disclose_identifiers(tmp_path):
    run = make_recovery(tmp_path)
    data = pd.read_csv(run / "eligible.csv")
    data.loc[1, "patient_id"] = data.loc[0, "patient_id"]
    data.to_csv(run / "eligible.csv", index=False)
    with pytest.raises(ValueError, match="duplicate patient") as exc:
        load_recovery(run, allow_demo=True)
    assert "SYNTHETIC_0000" not in str(exc.value)


@pytest.mark.parametrize("mutation", ["state", "audit", "model_n", "missing_term", "duplicate_term",
                                     "ratio", "ci", "grid_missing", "grid_duplicate", "sum"])
def test_recovery_rejects_inconsistent_saved_evidence(tmp_path, mutation):
    run = make_recovery(tmp_path)
    if mutation == "state":
        path = run / "eligible.csv"
        frame = pd.read_csv(path)
        frame.loc[0, "state60"] = 3
    elif mutation == "audit":
        modify_json(run / "audit.json", "eligible_n", 89)
    elif mutation == "model_n":
        modify_json(run / "primary" / "result.json", "n", 80)
    elif mutation in {"missing_term", "duplicate_term", "ratio", "ci"}:
        path = run / "primary" / "coefficients.csv"
        frame = pd.read_csv(path)
        if mutation == "missing_term":
            frame = frame.iloc[1:]
        elif mutation == "duplicate_term":
            frame = pd.concat([frame, frame.iloc[:1]])
        elif mutation == "ratio":
            frame["ratio"] = frame.ratio.astype(object)
            frame.loc[0, "ratio"] = "invalid"
        else:
            frame.loc[0, "ratio_lower"] = 3
    else:
        path = run / "primary" / "standardized_states.csv"
        frame = pd.read_csv(path)
        if mutation == "grid_missing":
            frame = frame.iloc[3:]
        elif mutation == "grid_duplicate":
            frame = pd.concat([frame, frame.iloc[:1]])
        else:
            frame.loc[0, "probability"] = .79
    if mutation not in {"audit", "model_n"}:
        frame.to_csv(path, index=False)
    with pytest.raises(ValueError):
        load_recovery(run, allow_demo=True)


def test_failed_sensitivity_never_reuses_stale_coefficients(tmp_path):
    run = make_recovery(tmp_path)
    path = run / "complete_case"
    path.mkdir()
    (path / "coefficients.csv").write_text("invalid stale file", encoding="utf-8")
    frame = pd.read_csv(run / "results.csv")
    frame.loc[1, "p"] = .001
    frame.to_csv(run / "results.csv", index=False)
    result = load_recovery(run, allow_demo=True)
    assert "complete_case" not in result["coefficients"]
    assert pd.isna(result["results"].loc[1, "p"])
    assert "p" not in result["result_details"]["complete_case"]


@pytest.mark.parametrize("mutation", ["missing_effect", "grid_missing", "grid_duplicate", "count",
                                     "dead_to_alive", "mrs_denominator", "fit", "state_counts"])
def test_trajectory_rejects_inconsistent_saved_evidence(tmp_path, mutation):
    run = make_trajectory(tmp_path)
    if mutation == "missing_effect":
        path = run / "transition_effects.csv"
        frame = pd.read_csv(path).iloc[1:]
    elif mutation in {"grid_missing", "grid_duplicate"}:
        path = run / "state_probabilities_36_48_60.csv"
        frame = pd.read_csv(path)
        frame = frame.iloc[3:] if mutation == "grid_missing" else pd.concat([frame, frame.iloc[:1]])
    elif mutation in {"count", "dead_to_alive"}:
        path = run / "transition_counts.csv"
        frame = pd.read_csv(path)
        frame.loc[0, "transition_rows_n" if mutation == "count" else "death_to_alive_n"] = 1
    elif mutation == "mrs_denominator":
        path = run / "mrs_trajectory_distribution.csv"
        frame = pd.read_csv(path)
        frame.loc[0, "denominator_known"] = 89
    elif mutation == "fit":
        write_json(run / "fit_diagnostics.json", [{"converged": False, "n": 450}] * 2)
    else:
        modify_json(run / "audit.json", "state_counts", {})
    if mutation not in {"fit", "state_counts"}:
        frame.to_csv(path, index=False)
    with pytest.raises(ValueError):
        load_trajectory(run, allow_demo=True)


def test_manifest_is_explicit_relative_and_hash_only(tmp_path):
    run = make_recovery(tmp_path)
    result = source_manifest(run, ["primary/result.json"])
    assert len(result) == 1 and set(result[0]) == {"artifact", "sha256", "bytes"}
    assert result[0]["artifact"] == "primary/result.json" and len(result[0]["sha256"]) == 64
    for bad in ("../outside.txt", str(run / "status.json")):
        with pytest.raises(ValueError):
            source_manifest(run, [bad])


@pytest.mark.parametrize("mutation", ["sequential", "final_label", "final_count", "clinical_count"])
def test_recovery_flow_reconciles_all_stages_and_audit(tmp_path, mutation):
    run = make_recovery(tmp_path)
    path = run / "cohort_flow.csv"
    frame = pd.read_csv(path)
    if mutation == "sequential":
        frame.loc[0, "remaining"] = 100
    elif mutation == "final_label":
        frame.loc[1, "step"] = "other_endpoint"
    elif mutation == "final_count":
        # Keep the subtraction internally coherent, but disagree with the cohort.
        frame.loc[1, ["remaining", "excluded_here"]] = [80, 10]
    else:
        modify_json(run / "audit.json", "clinical_n", 100)
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="flow"):
        load_recovery(run, allow_demo=True)


@pytest.mark.parametrize("mutation", ["quantile", "negative_ci", "excessive_ci"])
def test_primary_probability_grid_matches_known_cohort_and_bounded_intervals(tmp_path, mutation):
    run = make_recovery(tmp_path)
    path = run / "primary" / "standardized_states.csv"
    frame = pd.read_csv(path)
    if mutation == "quantile":
        frame["wmh_ml"] += 1  # Still a complete nine-scenario grid.
    elif mutation == "negative_ci":
        frame.loc[0, "lower"] = -.01
    else:
        frame.loc[0, "upper"] = 1.01
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError):
        load_recovery(run, allow_demo=True)


def test_reader_returns_only_imaging_contrasts_despite_nuisance_ratio_overflow(tmp_path):
    run = make_recovery(tmp_path)
    path = run / "primary" / "coefficients.csv"
    frame = pd.read_csv(path)
    nuisance = {"term": "dead:intercept", "ratio": np.inf,
                "ratio_lower": 0, "ratio_upper": np.inf, "p": .5}
    pd.concat([frame, pd.DataFrame([nuisance])], ignore_index=True).to_csv(path, index=False)
    result = load_recovery(run, allow_demo=True)
    assert result["coefficients"]["primary"].term.tolist() == list(PRIMARY_TERMS)


@pytest.mark.parametrize("name", ["../primary", "..\\primary", "primary/other", "primary\\other",
                                 "stroke_recurrence", "bp", "not_a_prespecified_recovery_analysis"])
def test_analysis_names_reject_cross_platform_paths_and_unrelated_analyses(tmp_path, name):
    run = make_recovery(tmp_path)
    path = run / "results.csv"
    frame = pd.read_csv(path)
    frame.loc[1, "analysis"] = name
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="analysis name"):
        load_recovery(run, allow_demo=True)


@pytest.mark.parametrize("where,value", [("result", -.1), ("result", 1.1), ("result", "missing"),
                                        ("result", .03), ("results", np.nan), ("results", .03)])
def test_result_p_values_are_valid_and_agree_between_saved_sources(tmp_path, where, value):
    run = make_recovery(tmp_path)
    if where == "result":
        modify_json(run / "primary" / "result.json", "p", value)
    else:
        path = run / "results.csv"
        frame = pd.read_csv(path)
        frame.loc[0, "p"] = value
        frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="P value"):
        load_recovery(run, allow_demo=True)


def test_binary_sensitivity_requires_two_imaging_coefficients_and_preserves_or_scale(tmp_path):
    run = make_recovery(tmp_path)
    result = {"analysis": "dependent_or_dead", "status": "ESTIMATED", "n": 81, "p": .04,
              "family": "binary"}
    write_json(run / "dependent_or_dead" / "result.json", result)
    path = run / "results.csv"
    pd.concat([pd.read_csv(path), pd.DataFrame([result])], ignore_index=True).to_csv(path, index=False)
    table = pd.DataFrame([{"term": term, "ratio": 1.2, "ratio_lower": 1.01, "ratio_upper": 1.4,
                           "p": .04, "scale": "odds_ratio"} for term in ("wmh_ml", "gm119_ml")])
    table.to_csv(run / "dependent_or_dead" / "coefficients.csv", index=False)
    coefficients = load_recovery(run, allow_demo=True)["coefficients"]["dependent_or_dead"]
    assert len(coefficients) == 2 and coefficients.scale.eq("odds_ratio").all()
    table.iloc[:1].to_csv(run / "dependent_or_dead" / "coefficients.csv", index=False)
    with pytest.raises(ValueError, match="incomplete imaging"):
        load_recovery(run, allow_demo=True)


@pytest.mark.parametrize("mutation", ["recurrence", "missing_scope_audit", "median_point", "median_volume",
                                     "median_ci"])
def test_trajectory_scope_and_shared_median_are_verified(tmp_path, mutation):
    run = make_trajectory(tmp_path)
    if mutation in {"recurrence", "missing_scope_audit"}:
        modify_json(run / "audit.json", "future_recurrence_fields_requested",
                    ["Y5_STROKE"] if mutation == "recurrence" else None)
    else:
        path = run / "state_probabilities_36_48_60.csv"
        frame = pd.read_csv(path)
        median = frame.scenario.eq("gm_q50")
        if mutation == "median_volume":
            frame.loc[median, "wmh_ml"] += 1
        elif mutation == "median_ci":
            frame.loc[median, "upper"] += .01
        else:
            frame.loc[median & frame.state.eq("independent"), "probability"] -= .01
            frame.loc[median & frame.state.eq("dependent"), "probability"] += .01
        frame.to_csv(path, index=False)
    with pytest.raises(ValueError):
        load_trajectory(run, allow_demo=True)


def test_mrs_state_classification_must_match_visit_audit(tmp_path):
    run = make_trajectory(tmp_path)
    path = run / "mrs_trajectory_distribution.csv"
    frame = pd.read_csv(path)
    first = frame.month.eq(60) & frame.mrs.eq("0")
    dependent = frame.month.eq(60) & frame.mrs.eq("3")
    frame.loc[first, "n"] -= 1
    frame.loc[dependent, "n"] += 1
    known = frame.mrs.ne("unknown")
    frame.loc[known, "proportion"] = frame.loc[known, "n"] / frame.loc[known, "denominator_known"]
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="state audit"):
        load_trajectory(run, allow_demo=True)
