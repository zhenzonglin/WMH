"""One-command execution and researcher-authorized CEC exposure exclusions."""
from pathlib import Path

import pandas as pd
import pytest

from wmh_hcy.common import DataError, dump_json, sha256
from wmh_hcy.studies import CONTRACT, cli, runner
from wmh_hcy.studies.demo import make_demo
from wmh_hcy.studies.registry import STUDIES


def test_cec_invalid_exclusion_preserves_other_cohorts_and_raw_source(tmp_path):
    cfg = make_demo(tmp_path/'demo', n=300)
    source = Path(cfg['inputs']['clinical_csv'])
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    raw['CEC'] = '12'
    raw.to_csv(source, index=False)
    before = {s: runner.prepare(cfg, s) for s in STUDIES}
    original = pd.read_csv(Path(before['cec']['path'])/'eligible.csv', dtype={'patient_id': str})
    ids = original.patient_id.iloc[:8].tolist()
    row_indices = [raw.index[raw.code_n.eq(pid)][0] for pid in ids]
    raw.loc[row_indices, 'CEC'] = ['-1', 'inf', '-inf', 'bad', '0', '', '.A', '999']
    raw.to_csv(source, index=False)
    checksum = sha256(source)
    for study in STUDIES:
        state = runner.prepare(cfg, study)
        assert state['status'] == 'PREPARED'
        current = pd.read_csv(Path(state['path'])/'eligible.csv', dtype={'patient_id': str})
        if study != 'cec':
            old = pd.read_csv(Path(before[study]['path'])/'eligible.csv', dtype={'patient_id': str})
            assert current.patient_id.tolist() == old.patient_id.tolist()
        else:
            assert not set(ids[:4]+ids[5:7]) & set(current.patient_id)
            assert set(ids[4:5]+ids[7:8]) <= set(current.patient_id)
            assert current.cec.notna().all() and current.cec.ge(0).all()
            assert state['audit']['cec_invalid_policy']['excluded_at_cec_step'] == 4
            assert state['audit']['flow_excluded']['observed_baseline_cec'] == 2
            exclusions = pd.read_csv(Path(state['path'])/'exclusions.csv', dtype={'patient_id': str})
            assert exclusions.set_index('patient_id').loc[ids[:4], 'exclusion_reason'].eq('no_invalid_baseline_cec').all()
    assert sha256(source) == checksum


@pytest.mark.parametrize('case', ['all_invalid', 'other_covariate_invalid'])
def test_cec_exclusion_does_not_bypass_other_errors(tmp_path, case):
    cfg = make_demo(tmp_path/'demo', n=300)
    source = cfg['inputs']['clinical_csv']
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    raw['CEC'] = '-1' if case == 'all_invalid' else '12'
    if case == 'other_covariate_invalid':
        raw.loc[0, 'BSL_HDL'] = '-1'
    raw.to_csv(source, index=False)
    state = runner.prepare(cfg, 'cec')
    assert state['status'] == 'REVIEW_REQUIRED'
    if case == 'all_invalid':
        assert state['audit']['eligible_n'] == 0
    else:
        assert state['audit']['invalid_primary_covariates_require_review'] == ['BSL_HDL']


@pytest.mark.parametrize('failure', [None, 'blocked', 'exception'])
def test_start_runs_all_to_report_without_prompt(monkeypatch, capsys, failure):
    calls = []
    monkeypatch.setattr('sys.argv', ['wmh-study', 'start'])
    monkeypatch.setattr(cli, 'load_workstation', lambda *a: {'mode': 'real'})
    monkeypatch.setattr('builtins.input', lambda *a: pytest.fail('start must not prompt'))
    monkeypatch.setattr('wmh_hcy.studies.reporting.print_audit', lambda state: None)
    monkeypatch.setattr('wmh_hcy.studies.reporting.summary', lambda cfg, **kwargs: calls.append('summary'))

    def fake_run(cfg, study, through, *, fresh):
        calls.append((study, through, fresh))
        if study == 'bp' and failure == 'exception':
            raise DataError('test missing input')
        return {'study': study, 'status': 'REVIEW_REQUIRED' if study == 'bp' and failure else 'COMPLETED'}

    monkeypatch.setattr(runner, 'run', fake_run)
    if failure:
        with pytest.raises(SystemExit) as stopped:
            cli.main()
        assert stopped.value.code == 2
    else:
        cli.main()
    assert calls == [(s, 'report', True) for s in STUDIES]+['summary']
    assert '[4/4] kidney' in capsys.readouterr().out


def test_fresh_report_prepares_then_fits_without_loading_stale_snapshot(tmp_path, monkeypatch):
    cfg = make_demo(tmp_path/'demo', n=300)
    original_prepare = runner.prepare
    events = []

    def prepare(cfg, study):
        events.append('prepare')
        return original_prepare(cfg, study)

    def fit(*args):
        events.append('fit')
        return pd.DataFrame([{'status': 'ESTIMATED'}])

    monkeypatch.setattr(runner, 'prepare', prepare)
    monkeypatch.setattr(runner, 'load_prepared', lambda *a: pytest.fail('fresh start cannot reuse old input'))
    monkeypatch.setattr(runner, 'run_study', fit)
    monkeypatch.setattr('wmh_hcy.studies.reporting.report', lambda *a: events.append('report'))
    state = runner.run(cfg, 'cec', 'report', fresh=True)
    assert events == ['prepare', 'fit', 'report']
    assert state['status'] == 'COMPLETED'


def test_failed_start_does_not_report_older_success_as_current(tmp_path):
    cfg = make_demo(tmp_path/'demo', n=80)
    root = runner.study_root(cfg, 'cec')
    old = root/'runs/old'
    (old/'primary').mkdir(parents=True)
    dump_json(old/'status.json', {'contract': CONTRACT})
    dump_json(old/'primary/result.json', {'status': 'ESTIMATED', 'p': .001})
    dump_json(root/'latest_results.json', {'path': str(old), 'run': 'old'})
    table = runner.read_results(cfg, attempts={'cec': {'status': 'FAILED', 'error': 'input changed'}})
    row = table.set_index('study').loc['cec']
    assert row.status == 'FAILED' and pd.isna(row.p) and pd.isna(row.p_holm_four)
    assert (old/'primary/result.json').is_file()
