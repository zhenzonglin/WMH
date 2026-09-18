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
        'labels': {'other_name': 'intracranial stenosis'},
    }])
    before = sha256(tmp_path/'source_inventory.json')
    text = '\n'.join(candidate_page({'cec': tmp_path}))
    assert 'CEC: candidates=1' in text and 'APO_AI: candidates=1' in text
    assert 'ICAS: candidates=1' in text and 'other_name' in text
    assert 'NOT mapped' in text and '/private/' not in text
    assert sha256(tmp_path/'source_inventory.json') == before


def test_numeric_reason_counts_agree_with_preparation(tmp_path):
    cfg = make_demo(tmp_path/'demo', n=80)
    source = cfg['inputs']['clinical_csv']
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    values = ['', '.', '.A', 'A', '_', '.Z', '0', '-5', 'inf', '-inf', 'broken', '1.1']
    raw['M03_CYSC'] = '1'
    raw.loc[:len(values)-1, 'M03_CYSC'] = values
    raw.to_csv(source, index=False)
    _, fields = read_clinical(cfg)
    _, parts = numeric_parts(raw.M03_CYSC)
    bad = parts['zero'] | parts['negative'] | parts['nonfinite'] | parts['unparseable']
    assert int(bad.sum()) == int(fields.set_index('source').loc['M03_CYSC', 'invalid']) == 5
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
