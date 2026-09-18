"""Retired fields cannot change active cohorts, models, MI inputs or diagnostics."""
from pathlib import Path

import pandas as pd

from wmh_hcy.studies import CONTRACT
from wmh_hcy.studies.analyses import run_study
from wmh_hcy.studies.demo import make_demo
from wmh_hcy.studies.diagnostics import diagnose
from wmh_hcy.studies.registry import SOURCES, STUDIES, primary_spec, roles
from wmh_hcy.studies.runner import prepare, study_sources


def test_v4_removes_retired_fields_from_every_study_but_retains_hdl():
    assert CONTRACT == 'imaging_four_studies_20260918_v5'
    retired = {'IMG_ICAS', 'Apo_AI'}
    assert not retired & SOURCES.keys()
    for study in STUDIES:
        assert not retired & study_sources(study).keys()
        assert not {'icas', 'apo_ai'} & set(primary_spec(study).predictors)
        assert not {'icas', 'apo_ai'} & {r['variable'] for r in roles(study)}
    assert 'BSL_HDL' in study_sources('cec')
    assert 'hdl' in primary_spec('cec').covariates


def test_retired_fields_absent_or_invalid_do_not_change_preparation(tmp_path):
    cfg = make_demo(tmp_path/'demo', n=300)
    before = {s: prepare(cfg, s) for s in ('bp', 'cec')}
    raw = pd.read_csv(cfg['inputs']['clinical_csv'], dtype=str, keep_default_na=False)
    raw['IMG_ICAS'], raw['Apo_AI'] = 'invalid-unused', '-999'
    raw.to_csv(cfg['inputs']['clinical_csv'], index=False)
    for study, original in before.items():
        after = prepare(cfg, study)
        assert original['status'] == after['status'] == 'PREPARED'
        assert original['audit']['model_parameters'] == after['audit']['model_parameters']
        a = pd.read_csv(Path(original['path'])/'eligible.csv')
        b = pd.read_csv(Path(after['path'])/'eligible.csv')
        pd.testing.assert_frame_equal(a, b)
        fields = pd.read_csv(Path(after['path'])/'field_audit.csv')
        assert not fields.source.isin(['IMG_ICAS', 'Apo_AI']).any()
        assert not {'icas', 'apo_ai'} & set(b.columns)
    text = '\n'.join(line for page in (1, 2, 3) for line in diagnose(cfg, page))
    assert 'IMG_ICAS' not in text and 'Apo_AI' not in text


def test_cec_no_longer_schedules_apo_ai_extension(tmp_path, monkeypatch):
    cfg = make_demo(tmp_path/'demo', n=300)
    state = prepare(cfg, 'cec')
    root = Path(state['path'])
    data = pd.read_csv(root/'eligible.csv')
    master = pd.read_csv(root/'master.csv')
    seen = []

    def fake_run(d, spec, *args, **kwargs):
        seen.append(spec)
        if kwargs.get('capture_completed') is not None:
            kwargs['capture_completed'].append([d])
        return {'analysis': spec.name, 'study': spec.study, 'tier': spec.tier,
                'status': 'ESTIMATED', 'p': .5}, None

    monkeypatch.setattr('wmh_hcy.studies.analyses.run_one', fake_run)
    results = run_study(data, master, 'cec', tmp_path, {})
    assert results.analysis.tolist() == [
        'primary', 'complete_case', 'without_hdl_same_sample', 'white_matter',
        'function60', 'function60_with_gm', 'function60_observation_weighted', 'cec_spline']
    assert not any('apo_ai' in spec.predictors for spec in seen)
    assert results.loc[results.tier.eq('secondary'), 'analysis'].tolist() == [
        'white_matter', 'function60', 'function60_with_gm']
