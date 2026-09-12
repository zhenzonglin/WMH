import copy

import numpy as np
import pandas as pd
import pytest
from test_statistics import simulation_data

from wmh_hcy.adjustment import adjustment_columns, select_background, selection_registry
from wmh_hcy.design import Design
from wmh_hcy.inference import contrast_vector, decide_hypotheses, pooled_contrast
from wmh_hcy.models import CoxFit


def test_selection_applies_disjunctive_rule_not_variable_count():
    r = selection_registry()[0]
    cases = []
    for name, changes in [
        ("exposure_only", {"outcome_parent": False}),
        ("outcome_only", {"exposure_parent": False}),
        ("proxy", {"exposure_parent": False, "outcome_parent": False, "shared_cause_proxy": True}),
        ("instrument", {"known_instrument": True}),
        ("descendant", {"structural_descendant": True}),
        ("unrelated", {"exposure_parent": False, "outcome_parent": False}),
    ]:
        cases.append({**r, "variable": name, **changes})
    assert select_background(cases) == ["exposure_only", "outcome_only", "proxy"]
    assert len(select_background()) == 10
    assert "prior_stroke" in adjustment_columns()
    assert "nihss" not in adjustment_columns()
    assert {"nihss", "pre_mrs", "lesion_ml", "toast"} <= set(adjustment_columns("functional"))


def test_month3_keeps_baseline_nonlinearity_and_interaction():
    d, _ = simulation_data(821, n=300)
    d = d.assign(hcy3=d.hcy*.9+2, b123=d.b12, folate3=d.folate, cysc3=d.cysc,
                 sample3_day=np.linspace(85, 110, len(d)))
    spec = Design(month3=True).fit(d)
    assert {"H", "H_rcs", "baseline_log2_hcy", "baseline_hcy_rcs", "baseline_H_x_W"} <= set(spec.columns)
    assert "H_x_W" not in spec.columns
    assert "log2_b123" in spec.columns
    assert "log2_b12" not in spec.columns
    m0 = copy.deepcopy(spec)
    m0.columns = [c for c in spec.columns if c not in ["H", "H_rcs"]]
    assert spec.transform(d).shape[1]-m0.transform(d).shape[1] == 2


def test_contrast_uses_off_diagonal_covariance():
    d, spec = simulation_data(123, n=300)
    columns = ["H", "H_rcs"]
    cov = np.array([[.4, -.12], [-.12, .3]])
    fit = CoxFit(columns, np.array([.2, .1]), cov, np.array([]), np.array([]), {})
    c = contrast_vector(d, spec, columns)
    result = pooled_contrast(d, spec, [fit])
    assert result["estimate"] == pytest.approx(c @ fit.params)
    assert result["se"]**2 == pytest.approx(c @ cov @ c)
    assert result["se"]**2 != pytest.approx(np.sum(c*c*np.diag(cov)))


def test_holm_family_stays_three_even_when_h4_unavailable():
    rows = [{"hypothesis": h, "estimate": e, "p": p} for h, e, p in
            [("H1", .2, .02), ("H2", .3, .03), ("H3", -.2, .01), ("H4", np.nan, np.nan)]]
    result = decide_hypotheses(rows).set_index("hypothesis")
    assert result.loc["H1", "p_decision"] == pytest.approx(.04)
    assert result.loc["H2", "decision"] == "SUPPORTED_POSITIVE"
    assert result.loc["H3", "decision"] == "SUPPORTED_OPPOSITE"
    assert result.loc["H4", "decision"] == "NOT_ESTIMABLE"
    single = pd.DataFrame(rows).iloc[:1].to_dict("records")
    assert decide_hypotheses(single).iloc[0].p_decision == pytest.approx(.06)
