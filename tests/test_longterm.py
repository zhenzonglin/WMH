import json
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyreadstat
import pytest
from test_contracts import base_rows
from test_statistics import simulation_data

from wmh_hcy.common import load_config, read_json
from wmh_hcy.fields import LONGTERM_FIELDS
from wmh_hcy.harmonize import harmonize
from wmh_hcy.imputation import impute
from wmh_hcy.longterm import audit_longterm, build_longterm_cohorts, run_longterm, summarize_estimates
from wmh_hcy.models import fit_cause
from wmh_hcy.sas_extract import extract
from wmh_hcy.synthetic import synthetic_frames


def sample(n=8):
    d = base_rows(n).assign(lesion_ml=3.0)
    d["is_event"] = np.nan  # The long-term cohort does not require year-one follow-up.
    d["y2_is_event"] = [0, 0, 1, 1, 1, 0, np.nan, 0][:n]
    d["y2_is_day"] = [730, 400, 40, 2, 800, np.nan, 730, 730][:n]
    d["mrs24"] = [2, 6, np.nan, 1, 0, 3, 2, 0][:n]
    return d


def test_censoring_event_cap_and_early_recurrence():
    cohorts, exclusions, _ = build_longterm_cohorts(sample(), 2)
    main = cohorts["H2"].set_index("patient_id")
    assert main.loc["0001", "exit"] == 400
    assert main.loc["0001", "event_type"] == 0
    assert main.loc["0004", "exit"] == 730  # Event after cutoff is censored, not counted.
    assert main.loc["0004", "event_type"] == 0
    assert main.loc["0002", "event_type"] == 1
    assert "0002" not in set(cohorts["H3"].patient_id)
    assert "0003" not in main.index  # Same-day recurrence cannot be ordered after sampling.
    assert (
        exclusions["H2"].set_index("patient_id").loc["0005", "exclusion_reason"]
        == "observed_event_or_censor_time"
    )
    assert "0006" not in main.index


def test_mrs_independent_of_recurrence_and_no_missing_as_death():
    cohorts, _, _ = build_longterm_cohorts(sample(), 2)
    function = cohorts["H4"].set_index("patient_id")
    assert "0005" in function.index and "0006" in function.index  # Missing IS time or state.
    assert "0002" not in function.index  # Missing mRS is not death or zero.
    assert function.loc["0001", "mrs24"] == 6
    assert np.isnan(function.loc["0001", "death_day"])  # mRS=6 supplies no death date.


def test_later_mrs_death_does_not_backfill_earlier_visits():
    d = sample(1).assign(mrs24=np.nan, mrs60=6.0)
    cohorts, _, _ = build_longterm_cohorts(d, 2)
    assert cohorts["H4"].empty


def test_confirmed_death_and_inconsistent_scores_are_separate():
    d = sample(2)
    d.loc[0, ["death_day", "mrs24", "y2_is_day"]] = [200, np.nan, 200]
    d.loc[1, ["death_day", "mrs24", "y2_is_day"]] = [250, 2, 400]
    cohorts, exclusions, _ = build_longterm_cohorts(d, 2)
    assert cohorts["H4"].mrs24.tolist() == [6]
    assert cohorts["H2"].exit.tolist() == [200]
    assert exclusions["H2"].exclusion_reason.iloc[1] == "not_observed_after_known_death"


def test_cross_year_first_event_conflicts_excluded_and_reported():
    d = sample(2).assign(is_event=1.0, is_day=40.0, y2_is_event=1.0, y2_is_day=40.0)
    d.loc[0, "y2_is_day"] = 50
    d.loc[1, "y2_is_event"] = 0
    assert audit_longterm(d, [2]).cross_year_conflicts.iloc[0] == 2
    cohorts, _, _ = build_longterm_cohorts(d, 2)
    assert cohorts["H2"].empty
    assert len(cohorts["H4"]) == 2  # A recurrence conflict does not remove function records.


def test_earlier_negative_status_contradicts_later_early_event():
    d = sample(1).assign(is_event=0.0, last_contact_day=365.0, y2_is_event=1.0, y2_is_day=30.0)
    assert audit_longterm(d, [2]).cross_year_conflicts.iloc[0] == 1


def test_each_year_uses_its_own_fields():
    d = sample(1).assign(y3_is_event=1.0, y3_is_day=950.0, mrs36=4.0)
    y2, _, _ = build_longterm_cohorts(d, 2)
    y3, _, _ = build_longterm_cohorts(d, 3)
    assert y2["H2"].exit.tolist() == [730]
    assert y3["H2"].exit.tolist() == [950]
    assert y3["H4"].mrs36.tolist() == [4]


def test_missing_year_does_not_borrow_other_year_or_block_function():
    d = sample(1).assign(mrs60=2.0)
    c, _, _ = build_longterm_cohorts(d, 5)
    assert c["H2"].empty and len(c["H4"]) == 1


def test_fixed_holm_four_when_only_one_year_estimated():
    table = summarize_estimates(
        [{"family": "H2", "year": 2, "p": 0.01, "estimate": 0.1, "lower": -0.1, "upper": 0.3}]
    )
    assert table.p_holm_4.iloc[0] == pytest.approx(0.04)
    assert table.ratio.iloc[0] == pytest.approx(np.exp(0.1))


def setup_inputs(tmp_path, n=100):
    (tmp_path / "config").mkdir()
    shutil.copy("config/analysis.yml", tmp_path / "config/analysis.yml")
    clinical, images = synthetic_frames(n, missing=False)
    for year in [2, 3, 4, 5]:
        # Lowercase CSV names exercise SAS-style case matching at harmonization.
        clinical[f"y{year}_is"] = clinical.y1_is
        known_dead = clinical.F12_DEATH.eq(1)
        death_day = (
            pd.to_datetime(clinical.F12_DEATH_D, errors="coerce") - pd.to_datetime(clinical.ONSET_D)
        ).dt.days
        clinical[f"y{year}_is_dd"] = clinical.y1_is_dd.where(
            clinical.y1_is.eq(1),
            pd.Series(float(365 * year), index=clinical.index).where(~known_dead, death_day.astype(float)),
        )
        clinical[f"m{12 * year}_mrs"] = clinical.F12_MRS.fillna(6)
    clinical.to_csv(tmp_path / "clinical.csv", index=False)
    images.to_csv(tmp_path / "imaging.csv", index=False)
    cfg = load_config(tmp_path / "config/analysis.yml")
    cfg["mode"] = "synthetic"
    cfg["inputs"].update(clinical_csv="clinical.csv", imaging_csv="imaging.csv")
    return cfg


def test_harmonization_case_missing_one_year_and_mrs6(tmp_path):
    cfg = setup_inputs(tmp_path)
    raw = pd.read_csv(tmp_path / "clinical.csv").drop(columns=["y1_is", "y1_is_dd"])
    raw.loc[0, "m60_mrs"] = 6
    raw.to_csv(tmp_path / "clinical.csv", index=False)
    d = harmonize(cfg, require_one_year=False)
    assert d.mrs60.iloc[0] == 6
    assert d.y5_is_day.notna().all()
    assert len(LONGTERM_FIELDS) == 28


def test_longterm_prepare_preserves_one_year_results_and_reports(tmp_path):
    cfg = setup_inputs(tmp_path)
    marker = tmp_path / "outputs/real/latest_results.json"
    marker.parent.mkdir(parents=True)
    marker.write_text('{"path":"original-one-year"}')
    result = run_longterm(cfg, through="prepare")
    assert result["status"] == "PREPARED"
    assert len(result["analyses"]) == 16
    assert json.loads(marker.read_text())["path"] == "original-one-year"
    assert len(read_json(marker.parent / "longterm/latest_results.json")["path"]) > 0


def test_risk_model_can_use_longer_observation_times():
    d, spec = simulation_data(32, n=400)
    d["exit"] *= 4  # Synthetic times only; tests never imply actual follow-up extension.
    fit = fit_cause(d, spec, 1)
    assert fit.diagnostics["events"] == int(d.event_type.eq(1).sum())
    assert np.isfinite(fit.params).all()


def test_functional_imputation_uses_correct_outcome_without_survival_fields():
    cfg = load_config("config/analysis.yml")
    cfg["analysis"].update(imputations=2, mice_iterations=1)
    d, _ = simulation_data(122, n=150)
    d = d.reset_index(drop=True).drop(columns=["entry", "exit", "event_type"])
    d = d.assign(nihss=4.0, toast=1.0, lesion_ml=3.0, pre_mrs=0.0, mrs60=np.arange(len(d)) % 7)
    d.loc[0:3, "b12"] = np.nan
    completed, meta = impute(d, cfg, "functional", functional_outcome="mrs60", survival_auxiliaries=False)
    assert "_mrs" in meta["predictors"]["b12"]
    assert "_D1" not in meta["predictors"]["b12"]
    assert completed[0].mrs60.tolist() == d.mrs60.tolist()


def test_longterm_survival_imputation_omits_unmeasured_death_information():
    cfg = load_config("config/analysis.yml")
    cfg["analysis"].update(imputations=2, mice_iterations=1)
    d, _ = simulation_data(111, n=150)
    d.loc[d.index[:3], "b12"] = np.nan
    _, meta = impute(d, cfg, death_auxiliaries=False)
    assert "_NA1" in meta["predictors"]["b12"]
    assert "_D2" not in meta["predictors"]["b12"] and "_NA2" not in meta["predictors"]["b12"]


def test_ordinal_model_uses_target_year():
    # Run in a fresh analysis process, as the CLI does. In a combined Windows
    # pytest process, nibabel's longdouble finfo initialization can change
    # statsmodels' cached float epsilon into a NumPy longdouble (linalg rejects
    # that dtype on Windows). Do not patch dependencies or change the estimator.
    code = '''
import numpy as np
import runpy
from wmh_hcy.design import Design
from wmh_hcy.models import ordinal
simulate = runpy.run_path("tests/test_statistics.py")["simulation_data"]
d, _ = simulate(891, n=350)
rng = np.random.default_rng(800)
d = d.assign(nihss=rng.poisson(4, len(d)), toast=rng.integers(1, 6, len(d)),
    lesion_ml=rng.exponential(3, len(d)), pre_mrs=rng.integers(0, 4, len(d)),
    mrs60=rng.integers(0, 7, len(d)), mrs12=0.)
spec = Design(expanded=True, functional=True).fit(d)
table, diag, _ = ordinal(d, spec, outcome="mrs60")
assert diag["converged"] and diag["outcome"] == "mrs60"
assert len(diag["mrs_distribution"]) == 7 and "H_x_W" in set(table.term)
'''
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=90, check=False)
    assert result.returncode == 0, result.stderr


def test_longterm_sas_directory_extracts_fields_from_multiple_sources(tmp_path, monkeypatch):
    cfg = setup_inputs(tmp_path)
    raw = pd.read_csv(tmp_path / "clinical.csv", dtype=str).rename(columns={"code_n": "CODE_N"})
    raw["CODE_N"] = [f"{i:05d}" for i in range(len(raw))]
    tables = {}
    columns = [c for c in raw if c != "CODE_N"]
    for i in range(3):
        path = tmp_path / f"part{i}.sas7bdat"
        path.touch()
        tables[str(path)] = raw[["CODE_N", *columns[i::3]]]

    def metadata(path, **kwargs):
        frame = tables[str(path)]
        return None, SimpleNamespace(
            column_names=list(frame),
            number_rows=len(frame),
            column_names_to_labels={},
            original_variable_types={},
            file_encoding="UTF-8",
        )

    def chunks(reader, path, **kwargs):
        frame = tables[str(path)][kwargs["usecols"]]
        yield frame.copy(), SimpleNamespace(missing_user_values={}, readstat_variable_types={})

    monkeypatch.setattr(pyreadstat, "read_sas7bdat", metadata)
    monkeypatch.setattr(pyreadstat, "read_file_in_chunks", chunks)
    cfg["inputs"].update(clinical_csv="", sas_globs=[str(tmp_path / "**/*.sas7bdat")])
    result = pd.read_csv(extract(cfg), dtype=str)
    assert {"y2_IS", "y3_IS_dd", "Y5_IS_DD", "m60_mrs"} <= set(result)
    assert result.code_n.iloc[0] == "00000" and len(result) == len(raw)
    assert result.Y5_IS_DD.notna().all()
