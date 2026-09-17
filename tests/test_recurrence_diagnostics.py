import warnings

import numpy as np
import pandas as pd
from test_recurrence_v3 import cohort_rows

from wmh_hcy.common import dump_json
from wmh_hcy.recurrence_data import build_cohort
from wmh_hcy.recurrence_diagnostics import timing_summary


def test_missing_image_flag_has_same_eligibility_without_futurewarning():
    data = cohort_rows()
    data["image_valid"] = data.image_valid.astype(object)
    data.loc[0, "image_valid"] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        cohort, exclusions, _, _ = build_cohort(data)
    assert "0000" not in set(cohort.patient_id)
    assert exclusions.exclusion_reason.iloc[0] == "available_wmh_icv"
    assert "0001" in set(cohort.patient_id)


def test_timing_summary_is_read_only_aggregate_and_flags_mismatch(tmp_path):
    data = pd.DataFrame({"patient_id": ["SECRET_ID_A", "SECRET_ID_B", "SECRET_ID_C"],
                         "entry": [3, 90, 400], "sample_day": [3, 90, 400],
                         "exit": [1825, 110, 1825], "event_type": [0, 1, 0],
                         "onset_date": ["2020-01-01"]*3,
                         "sample_date": ["2020-01-04", "2020-03-31", "2021-02-05"]})
    data.to_csv(tmp_path / "cohort.csv", index=False)
    dump_json(tmp_path / "config_snapshot.json", {"sas": {"date_formats": {}}})
    before = {str(p): p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    text = "\n".join(timing_summary(tmp_path))
    assert "N=3; IS=1; early_censor=0; event_free_at_day1825=2" in text
    assert "checked=3; mismatches=1" in text
    assert "sample_day vs entry mismatches=0" in text
    assert "90-179" in text and ">=365" in text
    assert "SECRET_ID" not in text and "2020-01" not in text and "2021-02" not in text
    assert before == {str(p): p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
