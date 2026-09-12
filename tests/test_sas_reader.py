from pathlib import Path

import pandas as pd
import pyreadstat
import pytest

from wmh_hcy.sas_extract import sas_to_csv


@pytest.mark.parametrize("name", ["basic__sample", "basic__dates", "basic__sample_bincompressed",
                                  "missing_data__missing_test"])
def test_actual_sas_chunk_csv_preserves_raw_values(tmp_path, name):
    source = Path("tests/fixtures/sas") / f"test_data__{name}.sas7bdat"
    full, meta = pyreadstat.read_sas7bdat(str(source), disable_datetime_conversion=True,
                                       user_missing=True)
    target = tmp_path / "中文路径_提取.csv"
    count, missing = sas_to_csv(source, target, meta.column_names, {"chunksize": 2})
    expected = tmp_path / "reference.csv"
    full.to_csv(expected, index=False, encoding="utf-8")
    assert count == len(full)
    pd.testing.assert_frame_equal(pd.read_csv(target, dtype=str, keep_default_na=False),
                                  pd.read_csv(expected, dtype=str, keep_default_na=False))
    assert missing == (meta.missing_user_values or {})
    if name == "basic__sample":
        assert "DATE" in meta.original_variable_types["mydate"] or "YYMMDD" in meta.original_variable_types["mydate"]
        assert pd.api.types.is_numeric_dtype(full.mydate)
        assert pd.api.types.is_numeric_dtype(full.dtime)


def test_optional_lesion_failure_preserves_whole_wmh(tmp_path):
    import json

    import nibabel as nib
    import numpy as np

    from wmh_hcy.imaging import read_subject

    for folder in ["wmh", "t1", "lesion"]:
        (tmp_path / folder).mkdir()
    (tmp_path / "wmh/wmh_features.json").write_text(json.dumps({
        "contralateral_correction": {"wmh_volume_after_correction_ml": 12}}))
    (tmp_path / "t1/t1_features.json").write_text(json.dumps({"dlicv_icv": {"icv_label702_ml": 1400}}))
    image = nib.Nifti1Image(np.full((3, 3, 3), .7), np.eye(4))
    image.header.set_xyzt_units("mm")
    nib.save(image, tmp_path / "lesion/lesion_space-FLAIR.nii.gz")
    row = read_subject({"imaging": {}, "inputs": {}}, tmp_path)
    assert row["wmh_ml"] == 12
    assert row["icv_ml"] == 1400
    assert np.isnan(row["lesion_ml"])
    assert "binary mask" in row["lesion_issue"]
