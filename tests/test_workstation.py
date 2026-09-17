import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyreadstat
import pytest
import yaml

from wmh_hcy.audit import clinical_audit
from wmh_hcy.common import DataError, load_config, read_json
from wmh_hcy.harmonize import sas_date
from wmh_hcy.sas_extract import inventory, iter_sas_chunks, resolve_owners
from wmh_hcy.synthetic import synthetic_frames
from wmh_hcy.workstation import configure, discover_sustain, image_audit, run_pipeline


@pytest.fixture
def setup(tmp_path):
    (tmp_path / "config").mkdir()
    shutil.copy("config/analysis.yml", tmp_path / "config/analysis.yml")
    shutil.copy("config/gm119_labels.json", tmp_path / "config/gm119_labels.json")
    (tmp_path / "pyproject.toml").write_text("[project]\nname='fixture'\n")
    clinical, images = synthetic_frames(n=100, missing=False)
    clinical["code_n"] = [f"{i:05d}" for i in range(100)]
    images["participant_id"] = clinical.code_n
    clinical.to_csv(tmp_path / "clinical.csv", index=False)
    images.to_csv(tmp_path / "imaging.csv", index=False)
    cfg = load_config(tmp_path / "config/analysis.yml")
    cfg["mode"] = "synthetic"
    cfg["inputs"].update(clinical_csv="clinical.csv", imaging_csv="imaging.csv")
    return tmp_path, cfg, clinical, images


def test_audit_preserves_leading_zeros_counts_people_and_needs_no_images(setup):
    root, cfg, _, _ = setup
    cfg["inputs"]["imaging_csv"] = "does_not_exist.csv"
    result = clinical_audit(cfg)
    assert result["status"] == "READY_FOR_EXTRACTION"
    assert result["patients_in_selected_source_union"] == 100
    assert result["whitelist_fields_present"] == 38
    assert (root / "outputs/real/audit/report.html").is_file()


def test_audit_field_exists_but_all_missing_is_reported(setup):
    root, cfg, clinical, _ = setup
    clinical["BSL_B12"] = np.nan
    clinical.to_csv(root / "clinical.csv", index=False)
    result = clinical_audit(cfg)
    assert result["status"] == "REVIEW_REQUIRED"
    assert "BSL_B12" in result["required_all_missing_or_invalid"]
    assert "BSL_B12" not in result["required_missing"]


def test_audit_no_events_does_not_require_imputed_event_dates(setup):
    root, cfg, clinical, _ = setup
    clinical["y1_is"] = 0
    clinical["y1_is_dd"] = np.nan
    clinical.to_csv(root / "clinical.csv", index=False)
    result = clinical_audit(cfg)
    assert result["status"] == "READY_FOR_EXTRACTION"
    assert result["profiles"][1]["n_observed_all_listed_values"] == 100


def test_duplicate_ids_are_not_counted_as_eligible_patients(setup):
    root, cfg, clinical, _ = setup
    pd.concat([clinical, clinical.iloc[[0]]]).to_csv(root / "clinical.csv", index=False)
    result = clinical_audit(cfg)
    assert result["status"] == "REVIEW_REQUIRED"
    assert result["patients_in_selected_source_union"] == 0
    assert result["unique_ids_in_any_scanned_file"] == 100
    report = pd.read_csv(root / "outputs/real/audit/files.csv")
    assert report.rows_scanned.iloc[0] == 101
    assert report.duplicate_rows.iloc[0] == 1


def test_invalid_core_code_blocks_and_unknown_98_is_missing(setup):
    root, cfg, clinical, _ = setup
    clinical.loc[0, "GENDER"] = 98
    clinical.loc[1, "H_DIAB"] = 98
    clinical.to_csv(root / "clinical.csv", index=False)
    result = clinical_audit(cfg)
    assert result["invalid_core_values"] == {"GENDER": 1}
    assert result["status"] == "REVIEW_REQUIRED"


def test_sas_name_case_and_selected_source_are_explicit(setup, monkeypatch):
    root, cfg, _, _ = setup
    source = root / "test.SAS7BDAT"
    source.touch()
    meta = SimpleNamespace(column_names=["CODE_N", "BSL_HCY"], number_rows=3,
                           column_names_to_labels={}, original_variable_types={}, file_encoding="UTF-8")
    monkeypatch.setattr(pyreadstat, "read_sas7bdat", lambda *a, **k: (None, meta))
    cfg["inputs"]["sas_globs"] = [str(root / "**/*")]
    inv = inventory(cfg)
    assert inv[0]["field_columns"]["code_n"] == "CODE_N"
    other = {**inv[0], "path": str(root / "other.sas7bdat"), "basename": "other.sas7bdat"}
    owner, issues = resolve_owners([*inv, other], cfg)
    assert "BSL_HCY" not in owner and len(issues) == 1
    cfg["variable_sources"]["BSL_HCY"] = str(source)
    owner, issues = resolve_owners([*inv, other], cfg)
    assert not issues and owner["BSL_HCY"] == str(source)


def test_numeric_sas_id_missing_not_padded_to_fake_patient(monkeypatch):
    meta = SimpleNamespace(readstat_variable_types={"CODE_N": "double"})
    monkeypatch.setattr(pyreadstat, "read_file_in_chunks", lambda *a, **k:
                        iter([(pd.DataFrame({"CODE_N": [1.0, "A", np.nan]}), meta)]))
    chunks = list(iter_sas_chunks("fake", ["code_n"], {"chunksize": 2, "id_width": 4}, {"code_n": "CODE_N"}))
    assert chunks[0][0].code_n.tolist() == ["0001", "", ""]


def test_real_public_sas_inventory_has_no_fabricated_cnsr_fields(setup):
    _, cfg, _, _ = setup
    cfg["inputs"].update(clinical_csv="", sas_globs=[str(Path("tests/fixtures/sas").resolve() / "*.sas7bdat")])
    result = clinical_audit(cfg)
    assert result["sas_or_csv_files"] == 4
    assert not result["read_errors"]
    assert result["status"] == "REVIEW_REQUIRED"
    assert result["patients_in_selected_source_union"] == 0


def test_configure_in_two_stages_preserves_clinical_inputs(setup):
    root, _, _, _ = setup
    target = root / "config/workstation.local.yml"
    sas = root / "sas"
    sas.mkdir()
    result = configure(str(target), sas_dir=str(sas))
    assert result["sas_globs"] == [str(sas / "**/*")]
    sustain = root / "SuStaIn"
    derivative = sustain / "derivatives/sub-00001/wmh"
    derivative.mkdir(parents=True)
    configure(str(target), sustain_dir=str(sustain))
    saved = yaml.safe_load(target.read_text())
    assert saved["inputs"]["sas_globs"] == result["sas_globs"]
    assert saved["inputs"]["derivatives_root"] == str(sustain / "derivatives")
    assert list(target.parent.glob("*.local.yml.*.bak"))
    with pytest.raises(DataError, match="local.yml"):
        configure(str(root / "config/analysis.yml"), sas_dir=str(sas))


def test_sustain_ambiguity_requires_exact_derivatives_path(tmp_path):
    for name in ["first", "second"]:
        (tmp_path / "derivatives" / name / "sub-01/wmh").mkdir(parents=True)
    with pytest.raises(DataError, match="found 2"):
        discover_sustain(tmp_path)
    assert discover_sustain(tmp_path / "derivatives/first")["subject_directories"] == 1


def test_image_join_uses_exact_ids_without_review_gate(setup):
    root, cfg, _, images = setup
    images.loc[0, "participant_id"] = "0"  # Must not match clinical 00000.
    images.loc[1, "wmh_qc"] = "unreviewed"
    images.to_csv(root / "imaging.csv", index=False)
    result = image_audit(cfg)
    assert result["clinical_image_id_intersection"] == 99
    assert result["wmh_icv_eligible_matched"] == 99
    assert result["manual_image_review_required"] is False


@pytest.mark.parametrize("state", [None, "unreviewed", "fail", "stale"])
def test_prepare_ignores_review_labels_and_legacy_config(setup, state):
    root, cfg, _, images = setup
    if state is not None:
        for col in ["wmh_qc", "icv_qc", "t1_qc", "lesion_qc"]:
            images[col] = state
    images.to_csv(root / "imaging.csv", index=False)
    cfg["imaging"]["require_qc"] = True  # Old workstation settings must not re-enable the gate.
    cfg["inputs"]["qc_csv"] = "missing_legacy_reviews.csv"
    result = run_pipeline(cfg, "prepare")
    assert result["status"] == "COMPLETED"
    assert result["stages"]["image_audit"]["wmh_icv_eligible_matched"] == 100
    assert result["stages"]["image_audit"]["t1_available_matched"] == 100
    assert result["stages"]["image_audit"]["acute_lesion_available_matched"] == 100
    cohorts = result["stages"]["prepare"]
    assert cohorts["functional_t1"]["n"] == cohorts["functional"]["n"] > 0


def test_image_availability_still_rejects_invalid_numeric_values(setup):
    root, cfg, _, images = setup
    images.loc[0, "wmh_ml"] = -1
    images.loc[1, "icv_ml"] = 0
    images.loc[2, "wmh_ml"] = np.inf
    images.loc[3, "icv_ml"] = np.nan
    images.loc[4, "wmh_ml"] = images.loc[4, "icv_ml"] + 1
    images.loc[5, "wmh_ml"] = 0  # Zero burden is a valid observed value.
    images.loc[6, "gm119_ml"] = np.inf
    images.loc[7, "lesion_ml"] = -1
    images.loc[8, "lesion_ml"] = 0
    images.to_csv(root / "imaging.csv", index=False)
    result = image_audit(cfg)
    assert result["wmh_icv_eligible_matched"] == 95
    assert result["t1_available_matched"] == 99
    assert result["acute_lesion_available_matched"] == 99


def test_automatic_review_table_is_not_read(setup):
    root, cfg, _, _ = setup
    from wmh_hcy.common import dump_json

    derivative = root / "derivatives"
    subject = derivative / "sub-00001"
    dump_json(subject / "wmh/wmh_features.json", {
        "participant_id": "00001", "contralateral_correction": {"wmh_volume_after_correction_ml": 10}})
    dump_json(subject / "t1/t1_features.json", {"dlicv_icv": {"icv_label702_ml": 1500}})
    (derivative / "tables").mkdir()
    (derivative / "tables/qc_reviews.tsv").write_text("malformed obsolete review table\n")
    cfg["inputs"].update(imaging_csv="", derivatives_root="derivatives")
    result = image_audit(cfg)
    assert result["wmh_icv_eligible_matched"] == 1
    assert result["status"] == "READY_FOR_PREPARE"


def test_pipeline_prepare_never_starts_statistics_and_audit_never_needs_images(setup, monkeypatch):
    _, cfg, _, _ = setup
    import wmh_hcy.analysis
    monkeypatch.setattr(wmh_hcy.analysis, "analyse", lambda *a, **kw: pytest.fail("Unexpected model run"))
    assert run_pipeline(cfg, "prepare")["status"] == "COMPLETED"
    cfg["inputs"]["imaging_csv"] = "missing.csv"
    result = run_pipeline(cfg, "audit")
    assert list(result["stages"]) == ["audit"]


def test_pipeline_stops_at_failed_audit_and_records_status(setup):
    root, cfg, clinical, _ = setup
    clinical.drop(columns="BSL_HCY").to_csv(root / "clinical.csv", index=False)
    result = run_pipeline(cfg, "prepare")
    assert result["status"] == "REVIEW_REQUIRED"
    assert list(result["stages"]) == ["audit"]
    assert read_json(root / "outputs/real/pipeline_status.json")["status"] == "REVIEW_REQUIRED"


def test_date_special_missing_remains_missing():
    parsed = sas_date(pd.Series(["0", "A", ".Z", "_", ""], name="ONSET_D"), "DATE9")
    assert parsed.notna().tolist() == [True, False, False, False, False]


def test_full_clinical_audit_and_prepare_without_mri_interval(setup):
    root, cfg, clinical, _ = setup
    from wmh_hcy.fields import FIELDS
    from wmh_hcy.harmonize import harmonize

    assert "IMG_ONSET_TO_MRI_D" not in FIELDS
    assert "IMG_ONSET_TO_MRI_D" not in clinical
    result = clinical_audit(cfg)
    assert result["status"] == "READY_FOR_EXTRACTION"
    assert result["whitelist_fields_expected"] == 66  # 38 original + 28 optional long-term fields.
    assert result["required_missing"] == []
    assert result["profiles"][0]["n_observed_all_listed_values"] == 100
    assert result["profiles"][1]["n_observed_all_listed_values"] == 100
    standardized = harmonize(cfg)
    assert "mri_day" not in standardized
    assert run_pipeline(cfg, "prepare")["status"] == "COMPLETED"
    cohort = pd.read_csv(root / "outputs/real/prepared/cohort_main.csv")
    assert len(cohort) > 0
    assert cohort.entry.equals(cohort.sample_day)


def test_old_prepared_cohorts_require_regeneration(setup):
    root, cfg, _, _ = setup
    from wmh_hcy.analysis import load_cohort

    folder = root / "outputs/real/prepared"
    folder.mkdir(parents=True)
    pd.DataFrame({"patient_id": ["00001"], "entry": [10]}).to_csv(folder / "cohort_main.csv", index=False)
    with pytest.raises(DataError, match="Run prepare again"):
        load_cohort(cfg, "main")
    # Cohorts from the previous blood-draw revision still have the old image gate.
    from wmh_hcy.cohorts import COHORT_ENTRY_RULE
    from wmh_hcy.common import dump_json
    dump_json(folder / "cohort_contract.json", {"entry_rule": COHORT_ENTRY_RULE})
    with pytest.raises(DataError, match="Run prepare again"):
        load_cohort(cfg, "main")
    assert run_pipeline(cfg, "prepare")["status"] == "COMPLETED"
    assert len(load_cohort(cfg, "main")) > 0


def test_four_sas_sources_extract_and_prepare_with_38_fields(setup, monkeypatch):
    root, cfg, clinical, _ = setup
    columns = [c for c in clinical if c != "code_n"]
    tables = {}
    for i in range(4):
        path = root / f"source_{i}.sas7bdat"
        path.touch()
        tables[str(path)] = clinical[["code_n", *columns[i::4]]].copy()

    def metadata(path, **kwargs):
        assert kwargs.get("metadataonly")
        frame = tables[str(path)]
        meta = SimpleNamespace(column_names=list(frame), number_rows=len(frame),
                               column_names_to_labels={}, original_variable_types={}, file_encoding="UTF-8")
        return None, meta

    def chunks(reader, path, **kwargs):
        frame = tables[str(path)][kwargs["usecols"]]
        meta = SimpleNamespace(missing_user_values={}, readstat_variable_types={})
        yield frame.iloc[:50].copy(), meta
        yield frame.iloc[50:].copy(), meta

    monkeypatch.setattr(pyreadstat, "read_sas7bdat", metadata)
    monkeypatch.setattr(pyreadstat, "read_file_in_chunks", chunks)
    cfg["inputs"].update(clinical_csv="", sas_globs=[str(root / "*.sas7bdat")])
    result = run_pipeline(cfg, "prepare")
    census = result["stages"]["audit"]
    assert result["status"] == "COMPLETED"
    assert census["sas_or_csv_files"] == 4
    assert census["whitelist_fields_present"] == 38
    assert census["whitelist_fields_expected"] == 66
    assert census["patients_in_selected_source_union"] == 100
    extracted = pd.read_csv(root / "outputs/real/extracted/clinical_raw.csv", dtype=str)
    assert len(extracted) == 100 and len(extracted.columns) == 38
    assert "IMG_ONSET_TO_MRI_D" not in extracted
