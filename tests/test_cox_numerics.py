"""Equivalent solver recovery, unidentified fits, and failure-stage provenance."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from statsmodels.duration.hazard_regression import PHReg

from wmh_hcy import absolute_risk, analysis
from wmh_hcy.common import DataError, read_json
from wmh_hcy.models import CoxFit, fit_phreg_checked


def cox_problem():
    rng = np.random.default_rng(1801)
    n = 450
    x = rng.normal(size=(n, 3)) * [1, 100, .001] + [0, 1000, .02]
    entry = rng.uniform(0, 2, n)
    stop = entry + rng.exponential(40, n)
    events = rng.binomial(1, .7, n)
    model = PHReg(stop, x, status=events, entry=entry, ties="efron")
    return model, x, int(events.sum())


def test_scaled_solver_recovers_same_coefficients_covariance_and_residuals(monkeypatch):
    model, x, events = cox_problem()
    # Reference MLE from a well-scaled design, mapped back to original units.
    scale = x.std(axis=0)
    ref = PHReg(model.endog, (x-x.mean(axis=0))/scale, status=model.status,
                entry=model.entry, ties="efron").fit(disp=False, maxiter=150)
    original_fit = PHReg.fit

    def fail_original_newton(self, **kwargs):
        if self is model:
            raise np.linalg.LinAlgError("Singular matrix")
        return original_fit(self, **kwargs)

    monkeypatch.setattr(PHReg, "fit", fail_original_newton)
    result, diag = fit_phreg_checked(model, x, events)
    assert result.params == pytest.approx(ref.params / scale, rel=1e-6, abs=1e-8)
    assert result.cov_params() == pytest.approx(ref.cov_params()/np.outer(scale, scale), rel=1e-6)
    assert model.loglike(result.params) == pytest.approx(ref.model.loglike(ref.params), abs=1e-8)
    assert result.schoenfeld_residuals == pytest.approx(ref.schoenfeld_residuals * scale, nan_ok=True)
    assert diag["optimizer"] == "scaled_bfgs_then_newton"
    assert diag["max_scaled_score_per_event"] < 1e-6
    assert len(diag["failed_attempts"]) == 1


def test_constant_cox_column_is_rejected_without_deletion():
    x = np.ones((30, 1))
    model = PHReg(np.arange(1, 31), x, status=np.ones(30))
    with pytest.raises(DataError, match="no automatic covariate deletion"):
        fit_phreg_checked(model, x, 30)


def test_complete_separation_remains_failure():
    # All events in x=1 occur before any x=0 subjects leave the risk set.
    x = np.repeat([0., 1.], 50).reshape(-1, 1)
    stop = np.r_[np.full(50, 100), np.arange(1, 51)]
    status = x[:, 0].astype(int)
    model = PHReg(stop, x, status=status, ties="efron")
    with pytest.raises(DataError, match="both solvers"):
        fit_phreg_checked(model, x, 50)


@pytest.mark.parametrize("failed_cause,stage", [(1, "ischemic_fit"), (2, "death_fit")])
def test_bootstrap_records_the_actual_failing_cause(monkeypatch, failed_cause, stage):
    data = pd.DataFrame({"hcy": [10., 15.], "wmh_ml": [1., 2.]})
    grid = pd.DataFrame({"hcy": [10., 15.], "wmh_ml": [1., 1.], "wmh_quantile": [.5, .5],
                         "reference": [True, True], "supported": [True, True]})
    monkeypatch.setattr(absolute_risk, "make_grid", lambda *args: grid)
    monkeypatch.setattr(absolute_risk, "standardized_risks", lambda *args: np.array([.1, .2]))

    def fake_fit(_data, _spec, cause):
        if cause == failed_cause:
            raise DataError("Invalid Cox coefficients or covariance")

    monkeypatch.setattr(absolute_risk, "fit_cause", fake_fit)
    cfg = {"analysis": {"horizon": 365, "bootstrap_per_imputation": 3, "seed": 50,
                        "wmh_quantiles": [.5]}}
    _, ci, rd, diagnostics = absolute_risk.risk_analysis(
        [data], SimpleNamespace(hcy="hcy"), [(None, None)], cfg)
    assert ci.empty and rd.empty
    assert diagnostics["status"] == "BOOTSTRAP_UNSTABLE"
    assert diagnostics["valid"] == 0 and diagnostics["requested"] == 3
    assert all(row["stage"] == stage for row in diagnostics["failures"])


def test_death_failure_keeps_ischemic_estimate_and_records_stage(tmp_path, monkeypatch):
    data = pd.DataFrame({"patient_id": ["SYN001", "SYN002"], "entry": [1., 1.],
                         "exit": [30., 60.], "event_type": [1, 2]})

    class Spec:
        columns = ("H_x_W",)

        def fit(self, *args):
            return self

        def to_dict(self):
            return {"columns": self.columns}

    fit = CoxFit(["H_x_W"], np.array([.1]), np.eye(1), np.array([30.]), np.array([.1]), {})
    monkeypatch.setattr(analysis, "impute", lambda *args, **kwargs: ([data, data], {"m": 2}))

    def fake_fit(_data, _spec, cause):
        if cause == 2:
            raise DataError("Cox fit failed with both solvers: Singular matrix")
        return fit

    monkeypatch.setattr(analysis, "fit_cause", fake_fit)
    status = analysis.run_survival(data, {"analysis": {"spline_quantiles": [.1, .5, .9]}},
                                   tmp_path, "H2", spec=Spec(), risk=True)
    assert status["status"] == "ESTIMATED" and status["absolute_risk_status"] == "NOT_ESTIMABLE"
    failure = read_json(tmp_path / "risk_failure.json")
    assert failure["stage"] == "death_fit" and failure["imputation"] == 0
    assert failure["completed_death_fits"] == 0
    assert (tmp_path / "hypothesis_estimate.json").is_file()


def test_check_fit_isolated_command_keeps_existing_results(tmp_path, monkeypatch, capsys):
    from wmh_hcy import cli, fit_check
    from wmh_hcy.common import load_config

    cfg = load_config("config/analysis.yml")
    cfg.update(_out=str(tmp_path), mode="synthetic")
    cfg["analysis"]["imputations"] = 2
    data = pd.DataFrame({"patient_id": ["PRIVATE_SYN001", "PRIVATE_SYN002"], "entry": [1., 1.],
                         "exit": [30., 60.], "event_type": [1, 2]})
    (tmp_path / "prepared").mkdir()
    data.to_csv(tmp_path / "prepared/cohort_month3.csv", index=False)
    (tmp_path / "latest_results.json").write_text('{"path":"unchanged"}')
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    monkeypatch.setattr(cli, "load_workstation", lambda *args: cfg)
    monkeypatch.setattr(fit_check, "load_cohort", lambda *args: data)
    monkeypatch.setattr(fit_check, "impute", lambda *args: ([data, data], {"m": 2}))

    class Spec:
        columns = ("x",)

        def __init__(self, **kwargs):
            pass

        def fit(self, *args):
            return self

        def to_dict(self):
            return {"columns": self.columns}

    monkeypatch.setattr(fit_check, "Design", Spec)

    def fake_fit(_data, _spec, cause):
        if cause == 2:
            raise DataError("Singular matrix")
        diag = {"numerical_fit": {"optimizer": "scaled_bfgs_then_newton",
                                 "max_scaled_score_per_event": 1e-9, "scaled_information_condition": 15}}
        return SimpleNamespace(diagnostics=diag)

    monkeypatch.setattr(fit_check, "fit_cause", fake_fit)
    monkeypatch.setattr("sys.argv", ["wmh-hcy", "check-fit", "--hypothesis", "H3"])
    assert cli.main() == 2
    text = capsys.readouterr().out
    assert "ISCHEMIC: pass=2/2; original=0; recovered=2" in text
    assert "DEATH: pass=0/2" in text and "PRIVATE" not in text
    assert all(path.read_bytes() == content for path, content in before.items())
    summary = read_json(next((tmp_path / "diagnostics").glob("*/summary.json")))
    assert summary["status"] == "FIT_FAILURES" and len(summary["fits"]) == 4
    assert "PRIVATE" not in next((tmp_path / "diagnostics").glob("*/summary.json")).read_text()
