"""Synthetic contracts for the isolated longitudinal recovery trajectory patch."""
from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose

from wmh_hcy.studies.design import StudyDesign
from wmh_hcy.studies.recovery_trajectory import (
    _transition_y,
    fit_transition,
    root,
    source_fields,
    state_probabilities,
    trajectory_spec,
    transition_effects,
    transition_panel,
)
from wmh_hcy.studies.registry import COVARIATES, ModelSpec


def test_trajectory_fields_are_recovery_only():
    fields = source_fields()
    recurrence = {name for name in fields if "STROKE" in name.upper()} - {"H_STROKE"}
    assert not recurrence
    assert not {"CEC", "BSL_HCY", "M03_UACR", "M03_CYSC"} & set(fields)
    assert {"F3_MRS", "F12_MRS", "m36_mrs", "m60_mrs"} <= set(fields)


def test_transition_panel_does_not_forward_fill_and_death_is_absorbing():
    frame = pd.DataFrame({
        "patient_id": ["a", "b", "c", "d"],
        "state3": [0, 0, 0, 0], "state3_conflict": [False] * 4,
        "state12": [1, 2, np.nan, 2], "state12_conflict": [False] * 4,
        "state24": [np.nan, 2, 0, 0], "state24_conflict": [False] * 4,
        "state36": [0, 2, 0, 0], "state36_conflict": [False] * 4,
        "state48": [0, 2, 0, 0], "state48_conflict": [False] * 4,
        "state60": [0, 2, 0, 0], "state60_conflict": [False] * 4,
    })
    panel, audit = transition_panel(frame)
    assert set(zip(panel.patient_id, panel.from_month, panel.to_month)) == {
        ("a", 3, 12), ("a", 36, 48), ("a", 48, 60),
        ("b", 3, 12),
        ("c", 24, 36), ("c", 36, 48), ("c", 48, 60),
        ("d", 3, 12),
    }
    assert audit.loc[audit.to_month.eq(24), "death_to_alive_n"].iloc[0] == 1
    assert not audit.death_absorbing_check.all()
    assert not (panel.origin == 2).any()


def test_transition_y_represents_recovery_as_a_switch():
    panel = pd.DataFrame({"origin": [0, 0, 1, 1, 1], "destination": [0, 1, 0, 1, 2]})
    assert _transition_y(panel).tolist() == [0, 1, 1, 0, 2]


def _small_design_data(seed: int = 11):
    rng = np.random.default_rng(seed)
    n = 500
    base = pd.DataFrame({
        "patient_id": [f"P{i:04d}" for i in range(n)],
        "wmh_ml": rng.lognormal(2.5, .4, n),
        "gm119_ml": rng.normal(560, 35, n),
    })
    rows = []
    for i in range(n):
        for interval in range(5):
            origin = int((i + interval) % 7 == 0)
            score = .15 * (base.wmh_ml.iloc[i] - base.wmh_ml.mean()) / base.wmh_ml.std()
            u = rng.random()
            if origin == 0:
                destination = 1 if u < .25 + .08 * (score > 0) else 2 if u < .33 else 0
            else:
                destination = 0 if u < .30 else 2 if u < .40 else 1
            rows.append({"row_id": i, "patient_id": base.patient_id.iloc[i],
                         "interval_index": interval, "from_month": [3, 12, 24, 36, 48][interval],
                         "to_month": [12, 24, 36, 48, 60][interval],
                         "origin": origin, "destination": destination})
    return base, pd.DataFrame(rows)


def test_transition_fit_effects_and_state_probabilities_are_complete():
    base, panel = _small_design_data()
    spec = ModelSpec("trajectory_test", exposures=("wmh_ml", "gm119_ml"), covariates=(),
                     primary=(), splines=())
    design = StudyDesign.freeze(base, spec)
    fitted = fit_transition(base, panel, design)
    effects, joint = transition_effects([fitted])
    assert set(effects.transition) == {
        "independent_to_dependent", "dependent_to_independent",
        "independent_to_dead", "dependent_to_dead",
    }
    assert set(joint) == {"wmh_ml", "gm119_ml"}
    probabilities = state_probabilities([fitted], [base], design, draws=2)
    assert set(probabilities.month) == {36, 48, 60}
    totals = probabilities.groupby(["scenario", "month"]).probability.sum()
    assert_allclose(totals.to_numpy(), 1.0, atol=1e-10)


def test_trajectory_spec_preserves_recovery_exposures_and_covariates():
    spec = trajectory_spec()
    assert spec.exposures == ("wmh_ml", "gm119_ml")
    assert tuple(spec.covariates) == tuple(COVARIATES["recovery"])


def test_trajectory_writes_to_a_new_root_only(tmp_path):
    configured = {"_out": str(tmp_path)}
    trajectory_root = root(configured)
    assert trajectory_root.name == "01_recovery_trajectory"
    assert trajectory_root.parent.name == "studies"
    assert not (trajectory_root.parent / "01_recovery").exists()
    assert not (trajectory_root.parent / "01_recovery_path").exists()
