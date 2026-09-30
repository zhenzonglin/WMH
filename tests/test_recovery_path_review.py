"""Saved-run review preserves provenance, privacy, and full-covariance contrasts."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from scipy.stats import norm

from wmh_hcy.studies.recovery_path import recurrence_spec

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/recovery_path_review.py"
MODULE_SPEC = importlib.util.spec_from_file_location("recovery_path_review", SCRIPT)
review = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(review)


def dump(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def make_run(tmp_path: Path, name: str = "20260930T010000000000Z", status: str = "COMPLETED") -> Path:
    root = tmp_path / "01_recovery_path"
    run = root / "runs" / name
    run.mkdir(parents=True)
    dump(run / "status.json", {"study": "recovery_path", "run": name, "path": str(run),
                              "mode": "real", "status": status, "analysis": {"recurrence": "ESTIMATED"}})
    dump(root / "latest_path_attempt.json", {"run": name, "path": str(run), "status": status})
    dump(run / "audit.json", {"clinical_n": 5, "exact_id_intersection": 4, "base_eligible_n": 4,
                             "landmark_n": 3, "events": 1, "chronology": {"unresolved_loss_n": 1}})
    data = pd.DataFrame({"patient_id": [f"PRIVATE_{i}" for i in range(4)], "state60": [0, 1, 2, np.nan],
                         "wmh_ml": [5., 10., 25., 12.], "gm119_ml": [500., 550., 600., 525.],
                         "age": [60., 70., 75., 62.], "icv_ml": [1200., 1300., 1400., 1350.],
                         "lesion_ml": [2., 5., 10., 3.], "nihss": [2., 3., 4., 1.]})
    data.to_csv(run / "base_eligible.csv", index=False)
    data.iloc[:3].to_csv(run / "landmark.csv", index=False)
    pd.DataFrame({"step": ["actual_visit", "before_visit"], "remaining": [4, 3],
                  "excluded_here": [0, 1]}).to_csv(run / "event_flow.csv", index=False)
    path = run / "stroke_recurrence"
    dump(path / "result.json", {"status": "ESTIMATED", "n": 3, "imputations": 2, "p": .05,
                              "primary_terms": "wmh_ml;wmh_ml_rcs"})
    dump(path / "model_definition.json", asdict(recurrence_spec()))
    dump(path / "frozen_design.json", {
        "wmh_ml": {"center": float(np.log1p(10.)), "scale": .5, "knots": [-1., 0., 1.]},
        "gm119_ml": {"scale": 40.}})
    params = np.array([.2, .1, -.1])
    covariance = np.diag([.02, .03, .01])
    covariance[0, 1] = covariance[1, 0] = .01
    np.savez_compressed(path / "pooled_inputs.npz", terms=np.array(["wmh_ml", "wmh_ml_rcs", "gm119_ml"]),
                        params=np.tile(params, (2, 1)), covariance=np.tile(covariance, (2, 1, 1)))
    dump(path / "fit_diagnostics.json", [{"converged": True, "design_condition": 5.,
                                         "numerical": {"max_scaled_score_per_event": 1e-12}}]*2)
    dump(path / "imputation.json", {"m": 2, "iterations": 10, "seed": 123,
                                    "missing": {"education": 1},
                                    "transformed_mean_traces": [{"iteration": 1, "imputation": 0, "education": 2.},
                                                                {"iteration": 10, "imputation": 0, "education": 3.}]})
    dump(run / "mediation_gate.json", {"status": "NOT_ESTIMABLE", "checks": {"unresolved_loss": 1},
                                      "reasons": ["Temporal order unavailable"]})
    dump(run / "mediation.json", {"status": "NOT_ESTIMABLE", "reasons": ["Temporal order unavailable"]})
    dump(run / "prediction_failure.json", {"reason": "Synthetic fixture: no prediction fit"})
    return run


def test_review_is_aggregate_only_and_does_not_modify_saved_run(tmp_path):
    run = make_run(tmp_path)
    before = {p.relative_to(run): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in run.rglob("*") if p.is_file()}
    result = subprocess.run([sys.executable, str(SCRIPT), "--run", str(run), "--no-pause"],
                            capture_output=True, encoding="utf-8", check=True)
    assert "PRIVATE_" not in result.stdout and "patient_id" not in result.stdout
    assert "[1/6]" in result.stdout and "[6/6]" in result.stdout
    assert "base_eligible: 1/1/1/unknown=1" in result.stdout
    assert "landmark: 1/1/1/unknown=0" in result.stdout
    assert "Temporal order unavailable" in result.stdout
    assert "NOT GENERATED" in result.stdout and "READ_ERROR" not in result.stdout
    after = {p.relative_to(run): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in run.rglob("*") if p.is_file()}
    assert before == after


def test_saved_wmh_contrast_uses_covariance_off_diagonal(tmp_path):
    path = make_run(tmp_path) / "stroke_recurrence"
    z = (np.log1p(np.array([25., 10.])) - np.log1p(10.)) / .5
    # Independent three-knot restricted-cubic basis for [-1, 0, 1].
    basis = (np.maximum(z+1, 0)**3 - 2*np.maximum(z, 0)**3 + np.maximum(z-1, 0)**3)/4
    vector = np.array([z[0]-z[1], basis[0]-basis[1], 0.])
    params = np.array([.2, .1, -.1])
    covariance = np.array([[.02, .01, 0.], [.01, .03, 0.], [0., 0., .01]])
    se = np.sqrt(vector @ covariance @ vector)
    expected = np.exp([vector @ params, vector @ params-norm.ppf(.975)*se,
                       vector @ params+norm.ppf(.975)*se])
    actual = review.saved_contrast(path, 25., 10.)
    assert_allclose([actual[k] for k in ("estimate", "lower", "upper")], expected)
    diagonal_only = np.exp(vector @ params+norm.ppf(.975)*np.sqrt(vector @ np.diag(np.diag(covariance)) @ vector))
    assert not np.isclose(actual["upper"], diagonal_only)
    assert review.saved_contrast(path, 10., 10.) == {"estimate": 1., "lower": 1., "upper": 1.}


def test_latest_failed_attempt_does_not_fall_back_to_old_success(tmp_path):
    old = make_run(tmp_path, "20260929T010000000000Z")
    latest = make_run(tmp_path, "20260930T010000000000Z", "FAILED")
    root = latest.parent.parent
    dump(root / "latest_path_results.json", {"run": old.name, "path": str(old), "status": "COMPLETED"})
    selected, state = review.resolve_run(root)
    assert selected == latest and state["status"] == "FAILED"


def test_saved_run_identity_and_mode_are_checked(tmp_path):
    run = make_run(tmp_path)
    state = review.read_json(run / "status.json")
    state["run"] = "wrong_run"
    dump(run / "status.json", state)
    with pytest.raises(ValueError, match="identity/path"):
        review.resolve_run(run.parent.parent, run)
    state.update(run=run.name, mode="synthetic")
    dump(run / "status.json", state)
    with pytest.raises(ValueError, match="mode=real"):
        review.resolve_run(run.parent.parent, run)
