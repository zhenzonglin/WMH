"""Temporal boundary, stochastic-process and no-leakage checks for recovery-path."""
from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose

from wmh_hcy.studies.demo import make_demo
from wmh_hcy.studies.design import StudyDesign
from wmh_hcy.studies.recovery_mediation import feasibility, gformula, vital_state
from wmh_hcy.studies.recovery_path import (
    first_stroke_conflicts,
    landmark,
    observation_after_known_death,
    prepare,
    recurrence_curve,
    source_fields,
)
from wmh_hcy.studies.registry import ModelSpec, primary_spec


def test_recovery_path_fields_do_not_require_retired_studies():
    fields = source_fields()
    assert {"F3_DATE", "Y5_STROKE", "Y5_STROKE_DD", "y2_STROKE", "y4_STROKE_dd"} <= fields.keys()
    assert not {"CEC", "BSL_HCY", "M03_UACR", "IMG_ICAS", "Apo_AI"} & fields.keys()
    spec = primary_spec("recovery")
    assert "y5_stroke_event" not in spec.predictors  # no future event leakage into three-month risk model


def test_landmark_excludes_previsit_event_and_preserves_early_censor():
    frame = pd.DataFrame({"patient_id": ["a", "b", "c", "d"],
                          "visit3_day": [100, 100, 100, 100],
                          "y5_stroke_event": [1, 0, 1, 0],
                          "y5_stroke_day": [90, 200, 150, 100],
                          "state60": [0, 1, 2, np.nan]})
    cohort, flow = landmark(frame)
    assert cohort.patient_id.tolist() == ["b", "c"]
    assert cohort.loc[0, "exit"] == 200 and cohort.loc[0, "event_type"] == 0
    assert cohort.loc[1, "exit"] == 150 and cohort.loc[1, "event_type"] == 1
    assert flow.set_index("step").loc["no_first_recurrence_before_or_on_visit", "excluded_here"] == 1


def test_yearly_first_stroke_conflict_is_detected_without_all_years_required():
    frame = pd.DataFrame({"y5_stroke_event": [1, 0, 1, 1], "y5_stroke_day": [400, 365, 500, 500],
                          "y2_stroke_event": [1, 1, 0, np.nan], "y2_stroke_day": [400, 350, 600, np.nan]})
    assert first_stroke_conflicts(frame).tolist() == [False, True, True, False]


def test_followup_beyond_visit_confirming_death_is_rejected_without_making_up_death_day():
    frame = pd.DataFrame({"y5_stroke_day": [1825, 600, 1400], "state12": [0, 0, 0],
                          "state24": [2, 2, 0], "state36": [2, 2, 2]})
    assert observation_after_known_death(frame).tolist() == [True, False, True]


def test_gformula_zero_and_known_recurrence_path_effect():
    high = np.column_stack([np.ones(50), np.ones(50)])
    low = np.column_stack([np.ones(50), np.zeros(50)])
    rec = np.array([-2., 1., 0., 0., 0., 0.])
    death = np.array([-4., 0., 0., 0., 0., 0., 0.])
    null_dependence = np.array([-2., 0., 0.])
    null_high = gformula(high, high, (rec, death, null_dependence))
    null_low = gformula(high, low, (rec, death, null_dependence))
    assert_allclose(null_high, null_low, atol=1e-12)
    affected = np.array([-2., 0., 1.5])
    high_effect = gformula(high, high, (rec, death, affected))
    shifted = gformula(high, low, (rec, death, affected))
    assert high_effect[1] > shifted[1]
    assert_allclose(high_effect.sum(), 1)


def test_mediation_gate_records_censor_and_ambiguous_order():
    d = pd.DataFrame({"event_type": [1, 0], "exit": [1200., 700.], "entry": [90., 90.],
                      "state12": [0, 0], "state24": [0, 0], "state36": [0, 0],
                      "state48": [0, 0], "state60": [1, 2]})
    gate = feasibility(d, {"death_interval_ambiguous_event_n": 1, "event_after_known_death_n": 0})
    assert gate["status"] == "NOT_ESTIMABLE"
    assert gate["checks"]["early_censor"] == 1
    assert gate["checks"]["ambiguous_event_death_order"] == 1


def test_known_alive_death_flag_does_not_require_interim_mrs():
    frame = pd.DataFrame({"state12": [np.nan, 2], "death12": [2, 1]})
    assert vital_state(frame, 12).tolist() == [0., 2.]


def test_recurrence_curve_does_not_retransform_missing_education(tmp_path):
    data = pd.DataFrame({"wmh_ml": np.linspace(2., 30., 50),
                         "education": [np.nan] + [1, 2] * 24 + [1]})
    spec = ModelSpec(study="recovery", name="stroke_recurrence", family="cox",
                     exposures=("wmh_ml",), covariates=("education",),
                     splines=("wmh_ml",))
    design = StudyDesign.freeze(data, spec)
    completed = data.assign(education=data.education.fillna(1))
    terms = design.transform(completed).columns.tolist()
    params = np.array([.3, -.1, .2])
    covariance = np.diag([.02, .03, .04])
    covariance[0, 1] = covariance[1, 0] = .01
    model = tmp_path / "stroke_recurrence"
    model.mkdir()
    np.savez_compressed(model / "pooled_inputs.npz", terms=np.asarray(terms),
                        params=np.tile(params, (2, 1)),
                        covariance=np.tile(covariance, (2, 1, 1)))

    recurrence_curve(data, design, model, tmp_path)

    curve = pd.read_csv(tmp_path / "wmh_recurrence_curve.csv")
    assert (tmp_path / "wmh_recurrence_curve.png").exists()
    assert_allclose(curve.iloc[40].HR, 1.)
    vector = design.contrast(completed, {"wmh_ml": curve.iloc[0].wmh_ml},
                             {"wmh_ml": curve.iloc[0].reference_wmh_ml})
    expected_log_hr = vector @ params
    expected_se = np.sqrt(vector @ covariance @ vector)
    assert_allclose(curve.iloc[0].HR, np.exp(expected_log_hr))
    assert_allclose(curve.iloc[0].upper, np.exp(expected_log_hr + 1.959963984540054*expected_se))


def test_synthetic_preparation_is_separate_and_audited(tmp_path):
    cfg = make_demo(tmp_path / "fixture", n=300)
    state = prepare(cfg)
    assert state["status"] == "PREPARED"
    assert state["audit"]["landmark_n"] > 0
    assert state["audit"]["events"] > 0
    assert state["audit"]["model_parameters"] > 0
    assert not state["audit"]["covariates_entirely_missing"]
    assert "01_recovery_path" in state["path"]
    assert not (tmp_path / "fixture/outputs/studies/01_recovery/latest_results.json").exists()
