import copy
import json

import numpy as np
import pandas as pd
import pytest
from test_contracts import base_rows
from test_longterm import setup_inputs
from test_statistics import simulation_data

from wmh_hcy.common import DataError, load_config
from wmh_hcy.imputation import impute, pool_scalar
from wmh_hcy.models import CoxFit, fit_cause
from wmh_hcy.recurrence import run_recurrence
from wmh_hcy.recurrence_data import (
    HORIZONS,
    REQUIRED,
    SOURCES,
    build_cohort,
    counts,
    horizon_auxiliaries,
    split_at_90,
    truncate,
)
from wmh_hcy.recurrence_models import (
    TimeDesign,
    clinical_contrasts,
    fit_series,
    freeze_design,
    horizon_table,
    pooled_linear,
)
from wmh_hcy.sas_extract import resolve_owners


def cohort_rows():
    d = base_rows(8).assign(y5_is_event=[0, 1, 1, 1, 1, 0, 0, 0],
                            y5_is_day=[200, 90, 180, 2, 1900, 1825, 1825, 1825],
                            is_event=np.nan, is_day=np.nan)
    d.loc[5, "sample_day"] = 100
    return d


def test_shared_five_year_censoring_and_cutoff_boundaries():
    d, ex, _, _ = build_cohort(cohort_rows())
    assert "0003" not in set(d.patient_id)  # Same day as baseline sample.
    assert ex.exclusion_reason.iloc[3] == "no_recurrence_before_or_on_entry"
    m3 = truncate(d, 3).set_index("patient_id")
    assert "0005" not in m3.index
    assert m3.loc["0001", "event_type"] == 1  # Cutoff day belongs to this horizon.
    assert m3.loc["0002", "event_type"] == 0
    for month in [12, 24, 36, 48, 60]:
        t = truncate(d, month).set_index("patient_id")
        assert t.loc["0000", "exit"] == 200  # Never extend an early non-event.
        assert t.loc["0002", "exit"] == 180 and t.loc["0002", "event_type"] == 1
    assert truncate(d, 60).set_index("patient_id").loc["0004", "event_type"] == 0
    assert counts(d).iloc[0].entry_at_or_after_cutoff == 1


def test_conflicts_do_not_reconstruct_or_substitute_year5():
    d = cohort_rows()
    d.loc[0, ["y2_is_event", "y2_is_day"]] = [1, 35]
    d.loc[1, ["y4_is_event", "y4_is_day"]] = [1, 80]
    c, ex, _, audit = build_cohort(d)
    assert audit["cross_year_conflicts"] == 2
    assert not {"0000", "0001"} & set(c.patient_id)
    assert ex.exclusion_reason.iloc[0] == "consistent_cumulative_first_event"


def test_known_death_requires_exit_at_or_before_death():
    d = cohort_rows()
    d.loc[0, "death_day"] = 100
    c, _, _, audit = build_cohort(d)
    assert "0000" not in set(c.patient_id) and audit["known_death_conflicts"] == 1


def test_year_one_non_event_dd_checked_without_other_visit_dates():
    d = cohort_rows()
    d.loc[1, ["is_event", "is_day", "last_contact_day"]] = [0, 100, np.nan]
    c, _, _, audit = build_cohort(d)
    assert "0001" not in set(c.patient_id) and audit["cross_year_conflicts"] == 1


def test_auxiliary_risk_sets_and_split_intervals():
    d, _, _, _ = build_cohort(cohort_rows())
    a = horizon_auxiliaries(d)
    late_index = d.index[d.patient_id.eq("0005")][0]
    assert a.loc[late_index, "_eligible_m3"] == 0
    assert a.loc[late_index, "_hazard_m3"] == 0
    assert a._event_m3.sum() == 1
    assert len(a.columns) == 21 and np.isfinite(a).all().all()
    s = split_at_90(d)
    assert (s.exit > s.entry).all()
    assert s.event_type.sum() == d.event_type.sum()
    assert s.loc[s.patient_id.eq("0001"), "late"].tolist() == [0]
    for _, p in s.groupby("patient_id"):
        if len(p) == 2:
            p = p.sort_values("entry")
            assert p.exit.iloc[0] == p.entry.iloc[1] == 90


def test_dedicated_fields_skip_unrelated_ambiguity():
    assert set(REQUIRED) <= set(SOURCES)
    assert not {"M03_HCY", "F12_MRS", "m60_mrs", "BSL_Cr"} & set(SOURCES)
    inv = [{"matched": ["code_n", "M03_HCY"], "basename": p, "path": p} for p in ["a", "b"]]
    assert resolve_owners(inv, {}, fields=SOURCES) == ({}, [])
    assert resolve_owners(inv, {})[1]


def test_profile_ignores_bad_unrelated_scores_and_preserves_legacy(tmp_path):
    cfg = setup_inputs(tmp_path)
    raw = pd.read_csv(tmp_path / "clinical.csv", dtype={"code_n": str})
    raw["M03_HCY"] = "unavailable"
    raw["F12_MRS"] = 99
    raw = raw.drop(columns=[c for c in raw if c.lower().startswith("y1_")])
    raw.to_csv(tmp_path / "clinical.csv", index=False)
    previous = tmp_path / "outputs/real/latest_results.json"
    previous.parent.mkdir(parents=True)
    previous.write_text('{"path":"old-one-year"}')
    old = previous.parent / "longterm/latest_results.json"
    old.parent.mkdir()
    old.write_text('{"path":"old-longterm"}')
    original = (tmp_path / "clinical.csv").read_bytes()
    result = run_recurrence(cfg, "prepare")
    assert result["status"] == "PREPARED"
    assert len(result["counts"]) == 7
    assert json.loads(previous.read_text())["path"] == "old-one-year"
    assert json.loads(old.read_text())["path"] == "old-longterm"
    assert (tmp_path / "clinical.csv").read_bytes() == original


def test_same_mi_preserves_endpoints_and_includes_all_horizon_auxiliaries():
    d, _ = simulation_data(891, n=180)
    d = d.reset_index(drop=True)
    d["patient_id"] = [f"{i:05d}" for i in range(len(d))]
    d.loc[:7, "b12"] = np.nan
    cfg = load_config("config/analysis.yml")
    cfg["analysis"].update(imputations=2, mice_iterations=1)
    frames, meta = impute(d, cfg, death_auxiliaries=False, extra_auxiliaries=horizon_auxiliaries(d))
    assert len(frames) == 2
    predictors = meta["predictors"]["b12"]
    assert "_D2" not in predictors and "_NA2" not in predictors
    for month in HORIZONS:
        assert f"_event_m{month}" in predictors and f"_hazard_m{month}" in predictors
    for f in frames:
        pd.testing.assert_frame_equal(f[["hcy", "wmh_ml", "entry", "exit", "event_type"]],
                                      d[["hcy", "wmh_ml", "entry", "exit", "event_type"]])


def test_full_covariance_and_rubin_pooling():
    fits = [CoxFit(["a", "b"], np.array([1., 2.]), np.array([[.3, -.1], [-.1, .2]]),
                   np.array([]), np.array([]), {}),
            CoxFit(["a", "b"], np.array([1.1, 2.2]), np.array([[.4, -.2], [-.2, .5]]),
                   np.array([]), np.array([]), {})]
    r = pooled_linear(fits, [1, 1])
    expected = pool_scalar([3., 3.3], [.3, .5])
    assert r["se"] == pytest.approx(expected["se"])
    assert r["HR"] == pytest.approx(np.exp(3.15))


def test_holm_fixed_six_and_primary_exclusion():
    r = horizon_table([{"month": 3, "p": .01}, {"month": 60, "p": .001}])
    assert r.loc[r.month.eq(3), "p_holm_6"].iloc[0] == pytest.approx(.06)
    assert np.isnan(r.loc[r.month.eq(60), "p_holm_6"].iloc[0])


def test_frozen_scale_and_local_support():
    d, _ = simulation_data(874, n=500)
    spec = freeze_design(d)
    before = copy.deepcopy(spec.to_dict())
    for month in HORIZONS:
        subset = truncate(d, month)
        pd.testing.assert_frame_equal(spec.transform(subset), spec.transform(d).loc[subset.index])
    assert before == spec.to_dict()
    # No fit is needed to detect unsupported contrast; never replace reference values.
    too_high = d.assign(hcy=30.0)
    result = clinical_contrasts(too_high, spec, [])
    assert set(result.status) == {"UNSUPPORTED"}


@pytest.mark.parametrize("beta", [0., .6])
def test_five_year_zero_and_known_interaction(beta):
    d, _ = simulation_data(872, n=2400, beta=beta)
    # Scale the complete survival clock in this synthetic fixture, preserving the true beta.
    d["entry"] *= 5
    d["exit"] *= 5
    spec = freeze_design(d)
    f = fit_cause(d, spec, 1)
    assert abs(f.params[f.columns.index("H_x_W")]-beta) < .2


def test_time_design_recovers_same_likelihood_when_late_coefficients_zero():
    from statsmodels.duration.hazard_regression import PHReg
    d, _ = simulation_data(832, n=400)
    spec = freeze_design(d)
    s = split_at_90(d)
    expanded = TimeDesign(spec)
    x, xs = spec.transform(d), expanded.transform(s)
    beta = np.linspace(-.03, .03, len(x.columns))
    original = PHReg(d.exit, x, status=d.event_type, entry=np.nextafter(d.entry, np.inf), ties="efron")
    split = PHReg(s.exit, xs, status=s.event_type, entry=np.nextafter(s.entry, np.inf), ties="efron")
    assert original.loglike(beta) == pytest.approx(split.loglike(np.r_[beta, np.zeros(5)]), abs=1e-8)


def test_fit_failure_does_not_pool_subset_or_call_death(tmp_path, monkeypatch):
    d, _ = simulation_data(666, n=100)
    d["patient_id"] = [str(i) for i in d.index]
    spec = freeze_design(d)
    calls = []

    def fake(frame, design, cause):
        calls.append(cause)
        raise DataError("intentional singular fit")

    monkeypatch.setattr("wmh_hcy.recurrence_models.fit_cause", fake)
    result, fits = fit_series([d, d], spec, tmp_path, "test")
    assert result["status"] == "NOT_ESTIMABLE" and not fits and calls == [1]
    assert not (tmp_path / "coefficients.csv").exists()


def test_orchestrator_reuses_shared_frames_and_never_legacy_risk(tmp_path, monkeypatch):
    from wmh_hcy.recurrence import analyse_recurrence
    d, _ = simulation_data(378, n=160)
    d = d.reset_index(drop=True).assign(wmh_raw_ml=lambda f: f.wmh_ml*1.02,
                                       lesion_ml=3.0, nihss=4.0, toast=np.nan)
    d["patient_id"] = [str(i) for i in d.index]
    calls, frame_ids = [], []

    def imputer(data, cfg, folder, kind="main"):
        calls.append(kind)
        if kind == "extended":
            raise DataError("TOAST unavailable")
        return [data.copy()]

    def fake_fit(frames, spec, folder, label, month=60, **kwargs):
        if label.startswith("month"):
            frame_ids.append(id(frames))
        return {"analysis": label, "month": month, "status": "NOT_ESTIMABLE", "p": np.nan,
                "reason": "test fixture", "n": len(frames[0]), "events": 1, "parameters": 22}, []

    monkeypatch.setattr("wmh_hcy.recurrence.shared_imputation", imputer)
    monkeypatch.setattr("wmh_hcy.recurrence.fit_series", fake_fit)
    cfg = load_config("config/analysis.yml")
    analyse_recurrence(d, cfg, tmp_path)
    assert len(set(frame_ids)) == 1 and len(frame_ids) == 7
    assert calls == ["main", "main", "extended"]
    assert "bootstrap" not in (tmp_path / "horizon_results.csv").read_text()


def test_failure_only_report_still_explains_every_horizon(tmp_path):
    from wmh_hcy.common import dump_json
    from wmh_hcy.recurrence_reporting import report_recurrence
    dump_json(tmp_path / "status.json", {"mode": "synthetic", "status": "COMPLETED_WITH_MODEL_FAILURES"})
    table = horizon_table([{"month": m, "status": "NOT_ESTIMABLE", "HR": np.nan,
                             "HR_lower": np.nan, "HR_upper": np.nan, "p": np.nan,
                             "reason": "No ischemic events"} for m in HORIZONS])
    table.to_csv(tmp_path / "horizon_results.csv", index=False)
    (tmp_path / "diagnostic_summary.csv").write_text("\n")
    report = report_recurrence(tmp_path).read_text(encoding="utf-8")
    assert "No ischemic events" in report and "五年主模型未成功估计" in report
