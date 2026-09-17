"""Revision contracts: continuous BP, paired CEC models, and renal-WMH interaction."""
from __future__ import annotations

import copy
import json

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from scipy.special import softmax

from wmh_hcy.common import DataError
from wmh_hcy.studies.analyses import run_one
from wmh_hcy.studies.design import StudyDesign
from wmh_hcy.studies.effects import bp_effects, kidney_effects, standardized_difference
from wmh_hcy.studies.imputation import impute
from wmh_hcy.studies.models import Fit, fit, state_probability_gradients
from wmh_hcy.studies.registry import ModelSpec, primary_spec
from wmh_hcy.studies.runner import read_results, study_root


def kidney_sample(n=3000):
    rng = np.random.default_rng(128)
    d = pd.DataFrame({"patient_id": [str(i) for i in range(n)],
                      "wmh_ml": np.expm1(rng.normal(2, .4, n)), "albuminuria": rng.integers(0, 4, n),
                      "cysc3": np.exp(rng.normal(0, .15, n))})
    spec = ModelSpec("kidney", exposures=("wmh_ml", "albuminuria"), covariates=("cysc3",), splines=(),
                     interactions=primary_spec("kidney").interactions,
                     primary=("dependent:albuminuria_3_x_wmh_ml",))
    design = StudyDesign.freeze(d, spec)
    x = design.transform(d)
    eta = -.4+.2*x.wmh_ml+.4*x.albuminuria_3+.65*x.albuminuria_3_x_wmh_ml
    p = softmax(np.column_stack([np.zeros(n), eta, -.6+.1*x.wmh_ml]), axis=1)
    d["state60"] = [rng.choice(3, p=v) for v in p]
    return d, spec, design


def test_kidney_primary_is_one_prespecified_interaction_and_fixed_reference():
    spec = primary_spec("kidney")
    assert spec.primary == ("dependent:albuminuria_3_x_wmh_ml",)
    assert len(spec.interactions) == 3
    d, simple, _ = kidney_sample(400)
    with pytest.raises(DataError, match="all four"):
        StudyDesign.freeze(d.loc[d.albuminuria.ne(0)], simple)


def test_kidney_interaction_recovery_and_shared_covariance_difference():
    d, _, design = kidney_sample(6000)
    fitted = fit(d, design)
    index = fitted.terms.index("dependent:albuminuria_3_x_wmh_ml")
    assert abs(fitted.params[index]-.65) < .15
    w = d.wmh_ml.median()
    a = design.transform(d.assign(albuminuria=3, wmh_ml=w))
    b = design.transform(d.assign(albuminuria=0, wmh_ml=w))
    q, u = standardized_difference(fitted, a, b)
    gradient = []
    for j in range(len(fitted.params)):
        hi, lo = copy.deepcopy(fitted), copy.deepcopy(fitted)
        hi.params[j] += 1e-5
        lo.params[j] -= 1e-5
        gradient.append((standardized_difference(hi, a, b)[0]-standardized_difference(lo, a, b)[0])/2e-5)
    gradient = np.asarray(gradient)
    assert_allclose(u, gradient @ fitted.covariance @ gradient, rtol=1e-6)
    assert_allclose(q, state_probability_gradients(fitted, a)[1][0]-state_probability_gradients(fitted, b)[1][0])
    curve = kidney_effects([d], design, [fitted])
    e = curve.loc[curve.status.eq("ESTIMATED")]
    assert_allclose(e.difference_estimate, e.persistent_probability-e.both_low_probability)
    assert (e.ratio_lower <= e.relative_probability_ratio).all()


def test_kidney_imputation_includes_categorical_interactions_without_bp_columns():
    d, _, design = kidney_sample(400)
    d.loc[:29, "cysc3"] = np.nan
    completed, info = impute(d, design, {"imputations": 2, "mice_iterations": 1, "mice_threads": 2}, progress=lambda _: None)
    assert "sbp3" not in d
    for i in (1, 2, 3):
        assert f"aux_albuminuria_{i}_x_wmh_ml" in info["predictors"]["cysc3"]
    for frame in completed:
        pd.testing.assert_series_equal(frame.albuminuria, d.albuminuria)
        assert frame.cysc3.notna().all()


def test_continuous_bp_reference_and_no_extrapolation():
    rng = np.random.default_rng(600)
    d = pd.DataFrame({"sbp3": rng.uniform(100, 175, 800), "wmh_ml": np.exp(rng.normal(2, .4, 800))})
    spec = primary_spec("bp").variant("validation", covariates=(), splines=("sbp3", "wmh_ml"))
    design = StudyDesign.freeze(d, spec)
    columns = list(design.transform(d))
    fitted = Fit(columns, np.linspace(.01, .07, len(columns)), np.eye(len(columns))*.002, {})
    continuous = bp_effects([d], design, [fitted], continuous=True)
    reference = continuous.loc[continuous.sbp.eq(140) & continuous.status.eq("ESTIMATED")]
    assert len(reference) == 3
    assert_allclose(reference[["HR", "HR_lower", "HR_upper"]], 1)
    assert len(continuous.sbp.unique()) > 90
    altered = d.assign(sbp3=d.sbp3-90)
    assert not bp_effects([altered], design, [fitted], continuous=True).status.eq("ESTIMATED").any()


def test_cec_pair_reuses_same_completed_data_and_scales(tmp_path, monkeypatch):
    rng = np.random.default_rng(812)
    n = 200
    d = pd.DataFrame({"patient_id": [str(i) for i in range(n)], "cec": rng.normal(size=n), "hdl": rng.normal(size=n)})
    d["gm119_ml"] = 500+5*d.cec+8*d.hdl+rng.normal(size=n)
    d.loc[:24, "hdl"] = np.nan
    spec = ModelSpec("cec", family="ols", exposures=("cec",), covariates=("hdl",),
                     outcome="gm119_ml", primary=("cec",), splines=())
    saved = []
    result, design = run_one(d, spec, tmp_path, {"imputations": 2, "mice_iterations": 1}, capture_completed=saved)
    assert result["status"] == "ESTIMATED"
    monkeypatch.setattr("wmh_hcy.studies.analyses.impute", lambda *args: pytest.fail("Paired model must not reimpute"))
    paired = spec.variant("without_hdl_same_sample", covariates=())
    second, second_design = run_one(d, paired, tmp_path, {}, inherited=design, completed=saved[0])
    assert second["status"] == "ESTIMATED" and second["same_primary_imputations"]
    assert result["n"] == second["n"] == n
    assert second_design.coding["cec"] == design.coding["cec"]
    a = pd.read_csv(tmp_path / "primary/model_membership.csv")
    b = pd.read_csv(tmp_path / "without_hdl_same_sample/model_membership.csv")
    pd.testing.assert_frame_equal(a, b)
    bad, _ = run_one(d.iloc[::-1], paired.variant("bad_pair"), tmp_path, {}, inherited=design, completed=saved[0])
    assert bad["status"] == "NOT_ESTIMABLE" and "row order" in bad["reason"]


def test_old_contract_not_in_current_holm_summary(tmp_path):
    cfg = {"_out": str(tmp_path), "mode": "synthetic"}
    root = study_root(cfg, "kidney")
    run = root / "runs/old"
    run.mkdir(parents=True)
    (run / "primary").mkdir()
    (run / "status.json").write_text(json.dumps({"contract": "imaging_five_studies_20260917_v1"}))
    (run / "primary/result.json").write_text(json.dumps({"status": "ESTIMATED", "p": .001}))
    (root / "latest_results.json").write_text(json.dumps({"path": str(run), "run": "old"}))
    row = read_results(cfg).set_index("study").loc["kidney"]
    assert row.status == "PREVIOUS_VERSION" and np.isnan(row.p) and np.isnan(row.p_holm_five)
