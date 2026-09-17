"""Scientific contracts and independent formula checks for the five-study extension."""
from __future__ import annotations

import copy

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from scipy import stats
from scipy.special import softmax

from wmh_hcy.common import DataError
from wmh_hcy.studies.data import add_bp, build_cohort, derive, functional_state, read_clinical
from wmh_hcy.studies.demo import make_demo
from wmh_hcy.studies.design import StudyDesign, split_time
from wmh_hcy.studies.imputation import impute
from wmh_hcy.studies.models import Fit, fit, state_probabilities, weighted_multinomial
from wmh_hcy.studies.pooling import contrast, joint_test
from wmh_hcy.studies.registry import SOURCES, STUDIES, ModelSpec, primary_spec
from wmh_hcy.studies.runner import prepare, read_results, study_sources


@pytest.fixture
def master(tmp_path):
    cfg = make_demo(tmp_path / "demo", n=700)
    clinical, _ = read_clinical(cfg)
    images = pd.read_csv(cfg["inputs"]["imaging_csv"], dtype={"participant_id": str}).rename(columns={"participant_id": "patient_id"})
    return derive(clinical.merge(images, on="patient_id", validate="one_to_one")), cfg


def test_source_separation_and_no_hcy():
    assert not any("HCY" in s.upper() for s in SOURCES)
    for study in STUDIES:
        assert ("CEC" in study_sources(study)) == (study == "cec")
        assert ("BSL_Cer_16_0" in study_sources(study)) == (study == "ceramide")
        assert ("M03_UACR" in study_sources(study)) == (study == "kidney")


def test_cec_missing_does_not_change_other_cohorts(master):
    d, _ = master
    changed = d.copy()
    changed["cec"] = np.nan
    for study in STUDIES:
        old = build_cohort(d, study)[0]
        new = build_cohort(changed, study)[0]
        if study == "cec":
            assert len(new) == 0
        else:
            assert old.patient_id.tolist() == new.patient_id.tolist()
    changed["wmh_ml"] = np.nan
    assert len(build_cohort(d.assign(wmh_ml=np.nan), "cec")[0]) > 0


def test_functional_death_absorbing_unknown_not_independent():
    d = pd.DataFrame({"mrs3": [1, 1, 1, 1, 1], "mrs12": [6, 1, 1, 1, 1],
                      "mrs60": [np.nan, np.nan, 4, 0, 2], "death24": [1, 2, 2, 1, 2]})
    state, conflicts = functional_state(d, 60)
    assert state[0] == 2
    assert np.isnan(state[1])
    assert state[2] == 1
    assert conflicts[3] and np.isnan(state[3])
    assert state[4] == 0
    state12, _ = functional_state(d, 12)
    assert state12[3] == 0  # future death never propagated backwards


def test_blood_pressure_arm_and_no_component_mixing():
    d = pd.DataFrame({"lsbp3": [130, np.nan, 120, 80], "rsbp3": [140, 130, np.nan, 150],
                      "ldbp3": [90, 99, 75, 90], "rdbp3": [80, 70, 60, 80]})
    result = add_bp(d, 3)
    assert result.sbp3[:3].tolist() == [140, 130, 120]
    assert result.dbp3[:3].tolist() == [80, 70, 75]
    assert result.sbp3_mean[0] == 135
    assert result.bp3_conflict[3] and np.isnan(result.sbp3[3])


def test_bp_landmark_prevents_reverse_time_and_preserves_early_censor(master):
    d, _ = master
    d.loc[0, ["visit3_day", "y5_is_event", "y5_is_day"]] = [100, 1, 99]
    d.loc[1, ["visit3_day", "y5_is_event", "y5_is_day"]] = [100, 0, 200]
    for c in ["is_event", "y2_is_event", "y3_is_event", "y4_is_event", "death_day"]:
        d.loc[[0, 1], c] = np.nan
    cohort, _, _ = build_cohort(d, "bp")
    assert d.patient_id[0] not in set(cohort.patient_id)
    row = cohort.set_index("patient_id").loc[d.patient_id[1]]
    assert row.exit == 200 and row.event_type == 0
    assert (cohort.exit > cohort.entry).all()


def test_uacr_threshold_and_ceramide_direction():
    d = derive(pd.DataFrame({"uacr0": [2.99, 3, 1, 3, np.nan], "uacr3": [2, 2, 3, 3, 9],
                             "cer16": [2, 4, 1, 1, 1], "cer24": [1, 2, 0, -1, np.nan]}))
    assert d.albuminuria[:4].tolist() == [0, 1, 2, 3]
    assert np.isnan(d.albuminuria[4])
    assert_allclose(d.cer_ratio[:2], [1, 1])
    assert d.cer_ratio[2:].isna().all()


def test_split_intervals_no_double_event():
    d = pd.DataFrame({"patient_id": ["01", "02", "03"], "entry": [90., 100, 400],
                      "exit": [365., 700, 1200], "event_type": [1, 1, 0]})
    split = split_time(d)
    assert split.event_type.sum() == d.event_type.sum()
    assert_allclose(split.groupby("patient_id").apply(lambda x: (x.exit-x.entry).sum(), include_groups=False), d.exit-d.entry)
    assert (split.exit > split.entry).all()


def test_rubin_contrast_keeps_offdiagonal_covariance():
    f = Fit(["a", "b"], np.array([1., 2.]), np.array([[2., .7], [.7, 3.]]), {})
    v = np.array([1., -1.])
    result = contrast([f], v)
    assert result["estimate"] == -1
    assert_allclose(result["se"]**2, 2+3-2*.7)


def test_d1_against_independent_formula():
    q = np.array([[.2, .3], [.4, .2], [.3, .5], [.1, .2]])
    u = np.tile(np.array([[.1, .02], [.02, .2]]), (4, 1, 1))
    r = (1+1/4)*np.trace(np.cov(q.T) @ np.linalg.inv(u[0]))/2
    statistic = q.mean(0) @ np.linalg.inv(u[0]) @ q.mean(0) / (2*(1+r))
    df = 4+(6-4)*(1+(1-2/6)/r)**2
    result = joint_test(q, u)
    assert_allclose([result["statistic"], result["df2"], result["p"]], [statistic, df, stats.f.sf(statistic, 2, df)])
    repeated = joint_test(np.tile(q[0], (4, 1)), u)
    assert_allclose(repeated["p"], stats.chi2.sf(q[0] @ np.linalg.inv(u[0]) @ q[0], 2))


def test_multinomial_reference_and_full_covariance_gradient():
    rng = np.random.default_rng(889)
    n = 3000
    d = pd.DataFrame({"x": rng.normal(size=n)})
    logits = np.column_stack([np.zeros(n), -.8+.6*d.x, -1.1-.3*d.x])
    probabilities = softmax(logits, axis=1)
    d["state60"] = [rng.choice(3, p=p) for p in probabilities]
    spec = ModelSpec("validation", exposures=("x",), primary=("dependent:x",), splines=())
    design = StudyDesign.freeze(d, spec)
    fitted = fit(d, design)
    assert abs(fitted.params[fitted.terms.index("dependent:x")]-.6) < .12
    assert abs(fitted.params[fitted.terms.index("dead:x")]+.3) < .12
    x = design.transform(d)
    estimates = state_probabilities(fitted, x)
    assert_allclose(sum(v[0] for v in estimates), 1)
    # Numerical finite differences independently verify all cross-equation covariance blocks.
    gradient = []
    for j in range(len(fitted.params)):
        upper, lower = copy.deepcopy(fitted), copy.deepcopy(fitted)
        upper.params[j] += 1e-5
        lower.params[j] -= 1e-5
        gradient.append((state_probabilities(upper, x)[1][0]-state_probabilities(lower, x)[1][0])/2e-5)
    gradient = np.array(gradient)
    assert_allclose(estimates[1][1], gradient @ fitted.covariance @ gradient, rtol=1e-6)
    unweighted = weighted_multinomial(x.to_numpy(), d.state60.to_numpy(), list(x), np.ones(n), {})
    assert_allclose(unweighted.params, fitted.params, atol=1e-5)


def test_observation_weight_sandwich_against_numeric_stacked_derivatives():
    import statsmodels.api as sm
    from scipy.optimize._numdiff import approx_derivative
    from scipy.special import expit

    rng = np.random.default_rng(555)
    n = 1000
    z = np.column_stack([np.ones(n), rng.normal(size=n)])
    observed = rng.random(n) < expit(.8+.7*z[:, 1])
    p = softmax(np.column_stack([np.zeros(n), -.4+.4*z[:, 1], -.7-.3*z[:, 1]]), axis=1)
    y = np.array([rng.choice(3, p=a) for a in p])
    obsfit = sm.Logit(observed.astype(int), z).fit(disp=False)
    prob = obsfit.predict(z)
    fitted = weighted_multinomial(z[observed], y[observed], ["intercept", "x"], 1/prob[observed], {}, (z, observed, prob))

    def scores(theta):
        beta, alpha = theta[:4], theta[4:]
        pr = expit(z @ alpha)
        out = softmax(np.column_stack([np.zeros(n), z @ beta.reshape(2, 2, order="F")]), axis=1)
        residual = np.eye(3)[y]-out
        beta_scores = np.column_stack([z*residual[:, j, None]*observed[:, None]/pr[:, None] for j in (1, 2)])
        return np.column_stack([beta_scores, z*(observed-pr)[:, None]])

    theta = np.r_[fitted.params, obsfit.params]
    a = -approx_derivative(lambda v: scores(v).sum(0), theta)
    b = np.linalg.inv(a)
    expected = b @ scores(theta).T @ scores(theta) @ b.T
    assert_allclose(fitted.covariance, expected[:4, :4], rtol=1e-5, atol=1e-8)


def test_generic_imputation_preserves_exposures_images_ids_and_observed_values(master):
    d, cfg = master
    cohort = build_cohort(d, "ceramide")[0]
    spec = primary_spec("ceramide")
    design = StudyDesign.freeze(cohort, spec)
    datasets, info = impute(cohort, design, cfg["analysis"], progress=lambda _: None)
    assert info["m"] == 2
    for filled in datasets:
        for c in (*spec.exposures, "patient_id", "wmh_ml", "icv_ml", spec.outcome):
            pd.testing.assert_series_equal(filled[c], cohort[c].reset_index(drop=True))
        for c in spec.covariates:
            assert filled[c].notna().all()
            mask = cohort[c].notna().to_numpy()
            assert_allclose(filled.loc[mask, c], cohort.loc[mask, c])
    with pytest.raises(DataError, match="Entirely"):
        impute(cohort.assign(cysc=np.nan), design, cfg["analysis"], progress=lambda _: None)


def test_no_real_fallback_and_legacy_pointer_preserved(master, tmp_path):
    _, cfg = master
    legacy = tmp_path / "legacy.json"
    legacy.write_text('"unchanged"')
    cfg["mode"] = "real"
    cfg["inputs"]["clinical_csv"] = str(tmp_path / "absent.csv")
    state = prepare(cfg, "cec")
    assert state["status"] == "INPUTS_REQUIRED" and state["mode"] == "real"
    assert legacy.read_text() == '"unchanged"'
    assert read_results(cfg).p.isna().all()


def test_explicit_assay_unit_conflicts_are_not_silently_converted():
    from wmh_hcy.studies.runner import explicit_unit_conflicts

    fields = {"BSL_UACR": ("uacr0", "continuous", None), "BSL_Cer_16_0": ("cer16", "continuous", None)}
    record = {"path": "source.sas7bdat", "field_columns": {s: s for s in fields},
              "labels": {"BSL_UACR": "UACR (mg/g)", "BSL_Cer_16_0": "Ceramide (nmol/L)"}}
    owner = {s: record["path"] for s in fields}
    assert len(explicit_unit_conflicts([record], owner, fields)) == 2
    record["labels"] = {"BSL_UACR": "UACR (mg/mmol)", "BSL_Cer_16_0": "Ceramide (pmol/L)"}
    assert explicit_unit_conflicts([record], owner, fields) == []
