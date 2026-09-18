"""Four-study contract, questionnaire derivation and explicit contradiction handling."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from wmh_hcy.common import DataError, dump_json
from wmh_hcy.studies import CONTRACT
from wmh_hcy.studies.analyses import run_one, run_study
from wmh_hcy.studies.data import build_cohort, chd_rule_summary, read_clinical, reconcile_chd
from wmh_hcy.studies.demo import make_demo
from wmh_hcy.studies.registry import FOLDERS, SOURCES, STUDIES, primary_spec
from wmh_hcy.studies.runner import prepare, read_results, study_sources


def test_chd_skip_rules_preserve_recorded_and_expose_contradictions():
    raw = pd.DataFrame({'chd': [np.nan, np.nan, 0, 1, np.nan, 1, 0, np.nan, 0, 1],
                        'heart_disease_gate': [1, 2, 1, 2, 1, 1, 2, 2, 2, 2],
                        'chd_type_present': [np.nan, 1, np.nan, 1, 1, np.nan, 1, np.nan, np.nan, np.nan]})
    d = reconcile_chd(raw)
    assert_allclose(d.chd, [0, 1, 0, 1, np.nan, np.nan, np.nan, np.nan, 0, 1], equal_nan=True)
    assert_allclose(d.chd_recorded, raw.chd, equal_nan=True)
    assert d.chd_rule_conflict.tolist() == [False]*4+[True]*3+[False]*3
    pd.testing.assert_frame_equal(reconcile_chd(d), d)
    counts = chd_rule_summary(d)
    assert counts['filled_from_H_HD'] == counts['filled_from_H_CHD_TP'] == 1
    assert counts['both_rules'] == counts['HD1_vs_recorded1'] == counts['TP_present_vs_recorded0'] == 1
    assert counts['conflicts'] == 3 and counts['unresolved_missing'] == 1


def test_chd_type_nonmissing_rule_and_sas_missing_tokens(tmp_path):
    cfg = make_demo(tmp_path/'demo', n=80)
    source = cfg['inputs']['clinical_csv']
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    raw['H_HD'], raw['H_CHD'], raw['H_CHD_TP'] = '2', '', ''
    tokens = ['', ' ', '.', '.A', 'A', '_', None, '0', '1', 'ACS']
    raw.loc[:len(tokens)-1, 'H_CHD_TP'] = tokens
    raw.to_csv(source, index=False)
    d, audit = read_clinical(cfg)
    assert d.chd[:7].isna().all()
    assert d.chd[7:10].eq(1).all()  # The user specified any nonmissing type, including numeric zero.
    assert not d.chd_rule_conflict.any()
    assert audit.set_index('source').loc['H_CHD_TP', 'observed'] == 3
    raw = raw.drop(columns=['H_HD', 'H_CHD_TP'])
    raw.loc[0, 'H_CHD'] = '1'
    raw.to_csv(source, index=False)
    d, _ = read_clinical(cfg)
    assert d.chd.iloc[0] == 1 and d.chd.iloc[1:].isna().all()


def test_four_active_studies_no_ceramides_no_visual_score_dependency():
    assert STUDIES == ('recovery', 'bp', 'cec', 'kidney')
    assert not any('cer' in s.lower() or 'WMH_Score' in s for s in SOURCES)
    assert {'H_HD', 'H_CHD_TP', 'H_CHD'} <= study_sources('bp').keys()
    assert 'Apo_AI' in study_sources('cec')
    for study in ('recovery', 'cec', 'kidney'):
        assert not {'H_HD', 'H_CHD_TP', 'H_CHD'} & study_sources(study).keys()
    with pytest.raises(DataError, match='not active'):
        build_cohort(pd.DataFrame(), 'ceramide')


def test_chd_conflict_blocks_bp_without_excluding_patients_or_other_studies(tmp_path, monkeypatch):
    cfg = make_demo(tmp_path/'demo', n=500)
    first = prepare(cfg, 'bp')
    assert first['status'] == 'PREPARED'
    cohort = pd.read_csv(Path(first['path'])/'eligible.csv', dtype={'patient_id': str})
    source = cfg['inputs']['clinical_csv']
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    row = raw.code_n.eq(cohort.patient_id.iloc[0])
    raw.loc[row, ['H_CHD', 'H_HD', 'H_CHD_TP']] = ['0', '2', '1']
    raw.to_csv(source, index=False)
    state = prepare(cfg, 'bp')
    assert state['status'] == 'REVIEW_REQUIRED'
    assert state['audit']['eligible_n'] == len(cohort)
    assert state['audit']['chd_rules']['eligible']['conflicts'] == 1
    conflict = pd.read_csv(Path(state['path'])/'chd_conflicts.csv')
    assert len(conflict) == 1 and conflict.chd_recorded.iloc[0] == 0
    assert prepare(cfg, 'cec')['status'] == 'PREPARED'
    d = pd.read_csv(Path(state['path'])/'eligible.csv', dtype={'patient_id': str})
    monkeypatch.setattr('wmh_hcy.studies.analyses.impute', lambda *a, **k: pytest.fail('Must not impute contradictions'))
    result, _ = run_one(d, primary_spec('bp'), tmp_path/'analysis', cfg['analysis'], complete_case=True)
    assert result['status'] == 'NOT_ESTIMABLE' and 'Conflicting H_CHD' in result['reason']


def test_bp_no_longer_schedules_visual_severe_wmh_analysis(tmp_path, monkeypatch):
    seen = []
    def fake_run(data, spec, *args, **kwargs):
        seen.append(spec.name)
        return {'analysis': spec.name, 'study': spec.study, 'tier': spec.tier, 'status': 'ESTIMATED', 'p': .5}, None
    monkeypatch.setattr('wmh_hcy.studies.analyses.run_one', fake_run)
    monkeypatch.setattr('wmh_hcy.studies.analyses.build_cohort', lambda *a, **k: (pd.DataFrame(), None, None))
    run_study(pd.DataFrame({'sbp3_mean': [130]}), pd.DataFrame(), 'bp', tmp_path, {})
    assert seen == ['primary', 'complete_case', 'mean_arms', 'month12_landmark', 'time_varying']


def test_four_test_holm_excludes_retired_and_previous_contracts(tmp_path):
    cfg = {'_out': str(tmp_path)}
    pvalues = [.01, .02, .2, .4]
    for study, p in zip(STUDIES, pvalues, strict=True):
        root = tmp_path/'studies'/FOLDERS[study]
        run = root/'runs'/'new'
        dump_json(root/'latest_results.json', {'path': str(run), 'run': 'new'})
        dump_json(run/'status.json', {'contract': CONTRACT})
        dump_json(run/'primary/result.json', {'status': 'ESTIMATED', 'p': p})
    retired = tmp_path/'studies/04_ceramide/latest_results.json'
    dump_json(retired, {'path': 'must_not_be_read'})
    result = read_results(cfg)
    assert_allclose(result.p_holm_four, [.04, .06, .4, .4])
    assert len(result) == 4 and 'p_holm_five' not in result
    old = tmp_path/'studies'/FOLDERS['bp']/'runs/new/status.json'
    dump_json(old, {'contract': 'imaging_five_studies_20260917_v2'})
    result = read_results(cfg).set_index('study')
    assert result.loc['bp', 'status'] == 'PREVIOUS_VERSION' and pd.isna(result.loc['bp', 'p_holm_four'])
    assert retired.exists()
