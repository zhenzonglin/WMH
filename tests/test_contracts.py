
import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from wmh_hcy.cohorts import build_cohorts
from wmh_hcy.common import DataError, dump_json, unique_ids
from wmh_hcy.harmonize import sas_date
from wmh_hcy.imaging import mask_volume_ml, read_subject
from wmh_hcy.imputation import nelson_aalen_increment, pool_scalar


def base_rows(n=1):
    return pd.DataFrame({
        "patient_id": [f"{i:04d}" for i in range(n)], "age": 60, "diagnosis": 1,
        "image_valid": True, "hcy": 14., "sample_day": 2.,
        "is_event": 0., "is_day": np.nan, "last_contact_day": 365.,
        "death_day": np.nan, "death_date_missing": False, "sample3_day": 92.,
        "hcy3": 12., "mrs12": 2., "gm119_ml": 600., "t1_qc": "pass",
    })


def test_sas_date_and_datetime_are_distinct():
    x = pd.Series(["0", "1", ""], name="date")
    result = sas_date(x, "DATE9")
    assert result.iloc[0] == pd.Timestamp("1960-01-01")
    assert result.iloc[1] == pd.Timestamp("1960-01-02")
    t = sas_date(pd.Series(["86400"]), "DATETIME20")
    assert t.iloc[0] == pd.Timestamp("1960-01-02")
    with pytest.raises(DataError):
        sas_date(pd.Series(["23000"]), "BEST")


def test_empty_and_duplicate_ids_rejected():
    unique_ids(pd.DataFrame({"id": ["0001", "0010"]}), "id", "test")
    with pytest.raises(DataError):
        unique_ids(pd.DataFrame({"id": ["01", "01"]}), "id", "test")


def test_early_recurrence_is_in_main_not_month3():
    d = base_rows()
    d.loc[0, ["is_event", "is_day"]] = [1, 40]
    c, _, _ = build_cohorts(d)
    assert len(c["main"]) == 1
    assert c["main"].event_type.iloc[0] == 1
    assert c["month3"].empty


@pytest.mark.parametrize("event_day", [1, 2])
def test_events_before_or_on_entry_excluded(event_day):
    d = base_rows()
    d.loc[0, ["is_event", "is_day"]] = [1, event_day]
    c, audit, _ = build_cohorts(d)
    assert c["main"].empty
    assert audit.exclusion_reason.iloc[0] == "no_event_before_or_on_entry"


def test_entry_is_blood_draw_and_next_day_recurrence_is_retained():
    d = base_rows()
    d.loc[0, ["is_event", "is_day"]] = [1, 3]
    c, _, _ = build_cohorts(d)
    assert c["main"].entry.tolist() == [2]
    assert c["main"].exit.tolist() == [3]
    assert c["main"].event_type.tolist() == [1]
    assert c["month3"].empty


@pytest.mark.parametrize("sample_day", [np.nan, -1, np.inf])
def test_blood_draw_time_is_still_required(sample_day):
    d = base_rows()
    d.loc[0, "sample_day"] = sample_day
    cohorts, audit, _ = build_cohorts(d)
    assert cohorts["main"].empty
    assert audit.exclusion_reason.iloc[0] == "known_baseline_measurement_time"


def test_death_is_competing_event_not_ischemic_recurrence():
    d = base_rows()
    d.loc[0, "death_day"] = 30
    c, _, _ = build_cohorts(d)
    assert c["main"].event_type.iloc[0] == 2
    assert c["functional"].mrs12.iloc[0] == 6


def test_confirmed_ischemia_same_day_as_death_has_event_priority():
    d = base_rows()
    d.loc[0, ["is_event", "is_day", "death_day"]] = [1, 30, 30]
    c, _, _ = build_cohorts(d)
    assert c["main"].event_type.iloc[0] == 1


def test_event_after_death_rejected():
    d = base_rows()
    d.loc[0, ["is_event", "is_day", "death_day"]] = [1, 50, 30]
    c, _, _ = build_cohorts(d)
    assert c["main"].empty


def test_missing_event_date_not_imputed():
    d = base_rows()
    d.loc[0, "is_event"] = 1
    c, a, _ = build_cohorts(d)
    assert c["main"].empty
    assert a.exclusion_reason.iloc[0] == "known_event_date_if_event"


def test_loss_is_censoring_not_death():
    d = base_rows()
    d.loc[0, "last_contact_day"] = 180
    c, _, _ = build_cohorts(d)
    assert c["main"].exit.iloc[0] == 180
    assert c["main"].event_type.iloc[0] == 0


def test_actual_month3_date_not_nominal_90():
    d = base_rows()
    d.loc[0, "sample3_day"] = 103
    c, _, _ = build_cohorts(d)
    assert c["month3"].entry.iloc[0] == 103


def test_missing_month3_does_not_reduce_main():
    d = base_rows()
    d.loc[0, "sample3_day"] = np.nan
    c, _, _ = build_cohorts(d)
    assert len(c["main"]) == 1 and c["month3"].empty


def test_missing_followup_not_assumed_one_year_event_free():
    d = base_rows()
    d.loc[0, "last_contact_day"] = np.nan
    c, _, _ = build_cohorts(d)
    assert c["main"].empty


def test_mask_volume_physical_units_and_affine(tmp_path):
    img = nib.Nifti1Image(np.ones((10, 10, 10), np.uint8), np.diag([2, 3, 4, 1]))
    img.header.set_xyzt_units("mm")
    path = tmp_path / "mask.nii.gz"
    nib.save(img, path)
    assert mask_volume_ml(path) == pytest.approx(24)


def test_unknown_mask_units_and_probabilities_rejected(tmp_path):
    for value, unit in [(1, "unknown"), (0.5, "mm")]:
        img = nib.Nifti1Image(np.full((2, 2, 2), value, dtype=float), np.eye(4))
        img.header.set_xyzt_units(unit)
        path = tmp_path / f"mask_{unit}.nii.gz"
        nib.save(img, path)
        with pytest.raises(DataError):
            mask_volume_ml(path)


def test_original_product_total_and_raw_icv_without_normalization(tmp_path):
    subject = tmp_path / "sub-0001"
    dump_json(subject / "wmh/wmh_features.json", {
        "participant_id": "0001", "raw_ml": {str(i): 1 for i in range(20)},
        "contralateral_correction": {"wmh_volume_after_correction_ml": 35,
                                      "wmh_volume_before_correction_ml": 38}})
    raw = subject / "t1/nichart_tool_output"
    raw.mkdir(parents=True)
    pd.DataFrame({"702": [1500000], "701": [1200000]}).to_csv(raw / "DLMUSE_Volumes.csv", index=False)
    cfg = {"inputs": {}, "imaging": {"path_prefix_map": {}}}
    row = read_subject(cfg, subject)
    assert row["wmh_ml"] == 35  # Not the 20-region sum.
    assert row["icv_ml"] == 1500  # Not other tissue volumes; no t1_features required.
    assert row["participant_id"] == "0001"


def test_nelson_aalen_respects_delayed_entry():
    d = pd.DataFrame({"entry": [0., 5.], "exit": [3., 8.], "event_type": [1, 0]})
    assert nelson_aalen_increment(d, 1).tolist() == [1, 0]


def test_rubin_pooling_accounts_for_between_imputation_variance():
    r = pool_scalar([0.1, 0.3], [0.04, 0.04])
    assert r["estimate"] == pytest.approx(0.2)
    assert r["se"]**2 == pytest.approx(0.04 + 1.5*0.02)
