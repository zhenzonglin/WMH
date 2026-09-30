"""Combined result review keeps main/extension provenance separate and never refits."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
PATH_SPEC = importlib.util.spec_from_file_location("recovery_path_review", SCRIPTS / "recovery_path_review.py")
path_review = importlib.util.module_from_spec(PATH_SPEC)
PATH_SPEC.loader.exec_module(path_review)
sys.modules["recovery_path_review"] = path_review
SPEC = importlib.util.spec_from_file_location("recovery_results_review", SCRIPTS / "recovery_results_review.py")
review = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(review)
SCRIPT = SCRIPTS / "recovery_results_review.py"


def dump(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def make_main(root: Path, name: str, status: str = "COMPLETED") -> Path:
    run = root / "01_recovery/runs" / name
    run.mkdir(parents=True)
    dump(run / "status.json", {"study": "recovery", "run": name, "path": str(run), "mode": "real",
                              "status": status, "audit": {"eligible_n": 4, "outcome_observed_n": 3,
                                                           "independent": 1, "dependent": 1, "dead": 1,
                                                           "unknown": 1}})
    data = pd.DataFrame({"patient_id": [f"PRIVATE_{i}" for i in range(4)], "state60": [0, 1, 2, np.nan],
                         "wmh_ml": [5., 10., 25., 12.], "gm119_ml": [500., 550., 600., 525.],
                         "age": [60., 70., 75., 62.], "education": [1., 2., np.nan, 1.],
                         "sex": [1, 2, 1, 2]})
    data.to_csv(run / "eligible.csv", index=False)
    dump(run / "primary/result.json", {"status": "ESTIMATED", "n": 3, "p": .005})
    pd.DataFrame({"term": ["dependent:wmh_ml", "dependent:gm119_ml", "dead:wmh_ml", "dead:gm119_ml"],
                  "ratio": [1.5, .8, 1.4, .7], "ratio_lower": [1.2, .5, 1.1, .4],
                  "ratio_upper": [1.8, 1.1, 1.7, 1.], "p": [.005, .4, .02, .06]}).to_csv(
                      run / "primary/coefficients.csv", index=False)
    pd.DataFrame({"analysis": ["primary", "month12", "month24"],
                  "status": ["ESTIMATED", "NOT_ESTIMABLE", "NOT_ESTIMABLE"],
                  "n": [3, 4, 4], "p": [.005, np.nan, np.nan],
                  "reason": [None, "Multinomial model failed to converge", "Nonfinite model estimates"]}).to_csv(
                      run / "results.csv", index=False)
    dump(run / "primary/frozen_design.json", {"wmh_ml": {"scale": .5}, "gm119_ml": {"scale": 40.}})
    dump(run / "primary/imputation.json", {"m": 2, "iterations": 10, "seed": 123,
                                           "missing": {"education": 1}, "transformed_mean_traces": [
                                               {"iteration": 1, "education": 2.},
                                               {"iteration": 10, "education": 3.}]})
    dump(run / "primary/fit_diagnostics.json", [{"converged": True, "design_condition": 5.}]*2)
    dump(run / "observation_weighted/observation_weights.json", [
        {"minimum_probability": .8, "weight_max": 1.25, "weight_p99": 1.2, "ess": 2.9}])
    pd.DataFrame({"wmh_ml": [10., 10., 10.], "gm119_ml": [550., 550., 550.],
                  "state": ["independent", "dependent", "dead"], "probability": [.7, .1, .2],
                  "lower": [.6, .05, .15], "upper": [.8, .15, .25]}).to_csv(
                      run / "primary/standardized_states.csv", index=False)
    return run


def make_path(root: Path, name: str = "20260929T081847414216Z") -> Path:
    run = root / "01_recovery_path/runs" / name
    dump(run / "status.json", {"study": "recovery_path", "run": name, "path": str(run),
                              "mode": "real", "status": "COMPLETED"})
    dump(run / "audit.json", {"base_eligible_n": 4, "landmark_n": 2, "events": 1})
    pd.DataFrame({"patient_id": ["PRIVATE_PATH_0", "PRIVATE_PATH_1"], "state60": [0, 2]}).to_csv(
        run / "landmark.csv", index=False)
    return run


def hashes(root: Path) -> dict:
    return {p.relative_to(root): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*") if p.is_file()}


def test_explicit_extension_and_latest_main_are_printed_separately_read_only(tmp_path):
    old = make_main(tmp_path, "20260918T115107666042Z")
    main = make_main(tmp_path, "20260928T010000000000Z")
    # An audit-only run must not replace a modeled main run.
    prepared = tmp_path / "01_recovery/runs/20260930T010000000000Z"
    dump(prepared / "status.json", {"status": "PREPARED"})
    dump(main.parent.parent / "latest_results.json", {"run": old.name, "path": str(old)})
    path = make_path(tmp_path)
    newer = make_path(tmp_path, "20260930T010000000000Z")
    dump(path.parent.parent / "latest_path_attempt.json", {"run": newer.name, "path": str(newer)})
    before = hashes(tmp_path)
    result = subprocess.run([sys.executable, str(SCRIPT), "--root", str(tmp_path), "--path-run", str(path),
                             "--no-pause"], capture_output=True, encoding="utf-8", check=True)
    text = result.stdout
    assert "RECOVERY | RUN=" + main.name in text
    assert "RECOVERY_PATH | RUN=" + path.name in text
    assert old.name not in text and newer.name not in text
    assert "[1/12]" in text and "[12/12]" in text
    assert "state60 independent/dependent/dead/unknown=1/1/1/unknown=1" in text
    assert "landmark: 1/0/1/unknown=0" in text
    assert "RRR=1.5 [1.2,1.8]" in text and "70[60,80] | 10[5,15] | 20[15,25]" in text
    assert "Multinomial model failed to converge" in text and "Nonfinite model estimates" in text
    assert "weight_max=1.25" in text
    assert "PRIVATE_" not in text and "patient_id" not in text
    assert "READ_ERROR" not in text
    assert "primary_result.png" in text and "wmh_recurrence_curve.png" in text
    assert hashes(tmp_path) == before


def test_newer_failed_main_does_not_fall_back_to_old_success(tmp_path):
    make_main(tmp_path, "20260918T115107666042Z")
    failed = tmp_path / "01_recovery/runs/20260929T010000000000Z"
    dump(failed / "status.json", {"study": "recovery", "run": failed.name, "path": str(failed),
                                 "mode": "real", "status": "FAILED"})
    assert review.resolve_recovery(tmp_path / "01_recovery")[0] == failed
    pages = review.recovery_pages(failed, path_review.read_json(failed / "status.json"))
    assert "recovery results withheld" in pages[0][1][0]


def test_main_same_run_counts_identity_and_real_mode_are_required(tmp_path):
    main = make_main(tmp_path, "20260928T010000000000Z")
    state = path_review.read_json(main / "status.json")
    state["audit"]["independent"] = 2
    with pytest.raises(ValueError, match="independent differs"):
        review.recovery_cohort(main, state)
    state["mode"] = "synthetic"
    dump(main / "status.json", state)
    with pytest.raises(ValueError, match="mode=real"):
        review.resolve_recovery(main.parent.parent, main)
    state.update(mode="real", run="wrong")
    dump(main / "status.json", state)
    with pytest.raises(ValueError, match="identity/path"):
        review.resolve_recovery(main.parent.parent, main)


def test_primary_cohort_mismatch_withholds_main_but_keeps_named_path(tmp_path):
    main = make_main(tmp_path, "20260928T010000000000Z")
    dump(main / "primary/result.json", {"status": "ESTIMATED", "n": 99})
    path = make_path(tmp_path)
    result = subprocess.run([sys.executable, str(SCRIPT), "--root", str(tmp_path), "--path-run", str(path),
                             "--no-pause"], capture_output=True, encoding="utf-8", check=True)
    assert "Primary model N differs" in result.stdout
    assert "RRR=1.5" not in result.stdout
    assert "RECOVERY_PATH | RUN=" + path.name in result.stdout


def test_missing_named_path_never_uses_another_directory(tmp_path):
    make_main(tmp_path, "20260928T010000000000Z")
    make_path(tmp_path, "20260930T010000000000Z")
    missing = tmp_path / "01_recovery_path/runs/20260929T081847414216Z"
    result = subprocess.run([sys.executable, str(SCRIPT), "--root", str(tmp_path), "--path-run", str(missing),
                             "--no-pause"], capture_output=True, encoding="utf-8", check=False)
    assert result.returncode == 2 and "Stopped" in result.stderr
    assert "RECOVERY_PATH" not in result.stdout
