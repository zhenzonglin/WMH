"""Screenshot diagnostics must locate missing sources without mutating patient data."""
import copy
from pathlib import Path

import pandas as pd
import pytest

from wmh_hcy.common import DataError, dump_json, sha256
from wmh_hcy.studies.data import cohort_audit, read_clinical
from wmh_hcy.studies.demo import make_demo
from wmh_hcy.studies.diagnostics import candidate_page, diagnose, numeric_parts, raw_columns
from wmh_hcy.studies.registry import COVARIATES, FOLDERS
from wmh_hcy.studies.reporting import print_audit
from wmh_hcy.studies.runner import prepare


def test_empty_cohort_is_not_evidence_all_covariates_absent(capsys):
    empty = pd.DataFrame(columns=COVARIATES['cec'])
    audit = cohort_audit(empty, 'cec')
    assert audit['covariates_entirely_missing'] == []
    assert audit['covariates_assessed'] is False
    print_audit({'study': 'cec', 'audit': audit})
    shown = capsys.readouterr().out
    assert 'NOT ASSESSED' in shown and 'Entirely missing covariates:' not in shown


def test_saved_metadata_candidates_do_not_map_names(tmp_path):
    dump_json(tmp_path/'source_inventory.json', [{
        'path': '/private/chemistry.sas7bdat', 'mtime_ns': 1, 'bytes': 500,
        'columns': ['other_name', 'CEC', 'Apo_AI'],
        'labels': {'other_name': 'cholesterol efflux'},
    }])
    before = sha256(tmp_path/'source_inventory.json')
    text = '\n'.join(candidate_page({'cec': tmp_path}))
    assert 'CEC: candidates=2' in text and 'other_name' in text
    assert 'ICAS' not in text and 'Apo_AI' not in text and 'APO_AI' not in text
    assert 'NOT mapped' in text and '/private/' not in text
    assert sha256(tmp_path/'source_inventory.json') == before


@pytest.mark.parametrize('source, expected', [('M03_CYSC', 5), ('CEC', 4)])
def test_numeric_reason_counts_agree_with_preparation(tmp_path, source, expected):
    cfg = make_demo(tmp_path/'demo', n=80)
    path = cfg['inputs']['clinical_csv']
    raw = pd.read_csv(path, dtype=str, keep_default_na=False)
    values = ['', '.', '.A', 'A', '_', '.Z', '0', '-5', 'inf', '-inf', 'broken', '1.1']
    raw[source] = '1'
    raw.loc[:len(values)-1, source] = values
    raw.to_csv(path, index=False)
    _, fields = read_clinical(cfg)
    _, parts = numeric_parts(raw[source])
    bad = parts['negative'] | parts['nonfinite'] | parts['unparseable']
    if source == 'M03_CYSC':
        bad |= parts['zero']
    assert int(bad.sum()) == int(fields.set_index('source').loc[source, 'invalid']) == expected
    assert parts['sas_or_blank_missing'].sum() == 6
    assert sum(int(v.sum()) for v in parts.values()) == len(raw)


def test_all_three_pages_are_readonly_and_never_show_ids(tmp_path):
    cfg = make_demo(tmp_path/'demo', n=150)
    states = {}
    for study in ('bp', 'kidney', 'cec'):
        states[study] = prepare(cfg, study)
    # Use audited source hashes for the CSV path; no new SAS scan is required.
    before = {str(p): sha256(p) for p in tmp_path.rglob('*') if p.is_file()}
    source = pd.read_csv(cfg['inputs']['clinical_csv'], dtype=str)
    text = '\n'.join(line for page in (1, 2, 3) for line in diagnose(cfg, page))
    after = {str(p): sha256(p) for p in tmp_path.rglob('*') if p.is_file()}
    assert before == after
    assert source.code_n.iloc[0] not in text
    assert 'M03_CYSC all_extracted:' in text and 'H_CHD eligible N=' in text
    assert 'CEC: PRESENT' in text
    assert 'CEC all_extracted:' in text and 'CEC eligible:' in text
    assert 'No saved SAS inventory' in text
    # Exact source absence is distinct from zero coverage.
    path = Path(states['cec']['path'])/'field_audit.csv'
    table = pd.read_csv(path)
    table.loc[table.source.eq('CEC'), ['present', 'observed']] = [False, 0]
    table.to_csv(path, index=False)
    assert 'CEC: ABSENT_COLUMN' in '\n'.join(diagnose(cfg, 1))


def test_external_csv_change_and_no_saved_audit_are_explicit(tmp_path):
    cfg = make_demo(tmp_path/'demo', n=80)
    with pytest.raises(DataError, match='No saved study audits'):
        diagnose(cfg)
    state = prepare(cfg, 'kidney')
    root = Path(state['path'])
    assert 'M03_CYSC' in raw_columns(root, ['M03_CYSC'])
    source = Path(cfg['inputs']['clinical_csv'])
    source.write_text(source.read_text()+'\n')
    with pytest.raises(DataError, match='changed since saved audit'):
        raw_columns(root, ['M03_CYSC'])
    changed = copy.deepcopy(cfg)
    changed['_out'] = str(tmp_path/'missing_outputs')
    with pytest.raises(DataError):
        diagnose(changed)
    assert not (tmp_path/'missing_outputs'/FOLDERS['kidney']).exists()


@pytest.mark.parametrize('in_cohort', [False, True])
def test_invalid_month3_cysc_only_blocks_its_kidney_cohort(tmp_path, in_cohort):
    cfg = make_demo(tmp_path/'demo', n=150)
    source = cfg['inputs']['clinical_csv']
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    raw['M03_CYSC'] = raw.M03_CYSC.replace('', '1.1')
    raw.loc[0, ['M03_CYSC', 'D_DIAG', 'F3_MRS', 'F3_DEATH', 'BSL_UACR', 'M03_UACR']] = [
        '0', '1' if in_cohort else '2', '1', '2', '1', '1']
    # Missing final outcome still belongs to the starting cohort and must not bypass review.
    raw.loc[0, 'm60_mrs'] = ''
    for column in ('D_DEATH', 'F6_DEATH', 'F12_DEATH', 'F2Y_DEATH', 'F3Y_DEATH', 'F4Y_DEATH', 'F5Y_DEATH'):
        raw.loc[0, column] = '2'
    for column in ('F12_MRS', 'm24_mrs', 'm36_mrs', 'm48_mrs'):
        raw.loc[0, column] = '1'
    raw.to_csv(source, index=False)
    before = sha256(Path(source))
    state = prepare(cfg, 'kidney')
    assert sha256(Path(source)) == before
    assert state['status'] == ('REVIEW_REQUIRED' if in_cohort else 'PREPARED')
    audit = state['audit']
    assert audit['invalid_fields'] == [{'source': 'M03_CYSC', 'invalid': 1}]
    assert audit['invalid_primary_covariates_require_review'] == (['M03_CYSC'] if in_cohort else [])
    fields = pd.read_csv(Path(state['path'])/'field_audit.csv').set_index('source')
    assert fields.loc['M03_CYSC', 'invalid_eligible'] == int(in_cohort)
    eligible = pd.read_csv(Path(state['path'])/'eligible.csv', dtype={'patient_id': str})
    row = eligible.loc[eligible.patient_id.eq(raw.loc[0, 'code_n'])]
    assert len(row) == int(in_cohort)
    if in_cohort:
        assert row.state60.isna().all() and row.cysc3.isna().all()


def test_invalid_cec_exposure_is_excluded_without_blocking(tmp_path):
    cfg = make_demo(tmp_path/'demo', n=150)
    source = cfg['inputs']['clinical_csv']
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    raw['CEC'] = raw.CEC.replace('', '12')
    raw.loc[:3, 'CEC'] = ['-1', 'inf', 'bad', '0']
    raw.to_csv(source, index=False)
    state = prepare(cfg, 'cec')
    assert state['status'] == 'PREPARED'
    assert state['audit']['invalid_primary_covariates_require_review'] == []
    assert state['audit']['cec_invalid_policy']['all_extracted'] == 3
    assert state['audit']['cec_invalid_policy']['remaining_eligible'] == 0
    assert state['audit']['invalid_fields_eligible'] == []
    lines = '\n'.join(diagnose(cfg, 3))
    assert 'unparseable=1; nonfinite=1; zero=1; negative=1; positive=146' in lines
    assert 'CEC zero is allowed' in lines
    assert raw.code_n.iloc[0] not in lines
