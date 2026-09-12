import numpy as np
import pandas as pd
import pytest

from wmh_hcy.absolute_risk import cumulative_incidence
from wmh_hcy.common import DataError, load_config
from wmh_hcy.design import Design, rcs_nonlinear
from wmh_hcy.imputation import impute
from wmh_hcy.models import CoxFit, fit_cause


class OneColumn:
    def transform(self, data):
        return pd.DataFrame({"x": np.zeros(len(data))})


def test_competing_risk_matches_constant_hazard_solution():
    t = np.arange(1, 101, dtype=float)
    fit1 = CoxFit(["x"], np.array([0.]), np.eye(1), t, np.full(100, .01), {})
    fit2 = CoxFit(["x"], np.array([0.]), np.eye(1), t, np.full(100, .02), {})
    d = pd.DataFrame({"entry": [0., 50.]})
    risks = cumulative_incidence(d, OneColumn(), fit1, fit2, 100)
    expected = (1/3)*(1-np.exp(-.03*np.array([100, 50])))
    assert risks == pytest.approx(expected, abs=1e-12)
    assert risks[0] < 1-np.exp(-1)  # Competing risk is not 1-KM for ischemic stroke.


def test_rcs_has_linear_tails():
    y = rcs_nonlinear(np.array([4., 5., 6.]), [0., 1., 2.])
    assert np.diff(y, 2)[0] == pytest.approx(0)


def simulation_data(seed, n=1400, beta=0.0):
    rng = np.random.default_rng(seed)
    d = pd.DataFrame({"age": rng.normal(60, 9, n), "hcy": np.exp(rng.normal(2.6, .4, n)),
                      "wmh_ml": np.exp(rng.normal(2., .7, n)), "icv_ml": rng.normal(1450, 90, n),
                      "b12": np.exp(rng.normal(5.7, .3, n)), "folate": np.exp(rng.normal(2.8, .3, n)),
                      "cysc": np.exp(rng.normal(0, .2, n)), "sample_day": rng.integers(1, 4, n),
                      "sex": rng.integers(1, 3, n), "smoking": rng.integers(1, 5, n),
                      "drinking": rng.integers(1, 5, n), "hypertension": rng.integers(1, 3, n),
                      "diabetes": rng.integers(1, 3, n), "prior_stroke": rng.integers(1, 3, n)})
    spec = Design().fit(d)
    x = spec.transform(d)
    event = rng.exponential(1/(0.003*np.exp(.2*x.H + .15*x.W + beta*x.H_x_W)))
    d["entry"] = np.maximum(d.sample_day, 3)
    d["exit"] = np.minimum(event, 365)
    d["event_type"] = (event <= 365).astype(int)
    d = d.loc[d.exit.gt(d.entry)].copy()
    return d, spec


@pytest.mark.parametrize("beta", [0., 0.65])
def test_known_and_zero_interaction_recovery(beta):
    d, spec = simulation_data(872, n=2400, beta=beta)
    fit = fit_cause(d, spec, 1)
    j = fit.columns.index("H_x_W")
    assert abs(fit.params[j]-beta) < .2
    assert fit.diagnostics["max_abs_gradient"] < 1e-4


def test_no_supplement_or_center_in_model():
    _d, spec = simulation_data(92, n=300)
    assert not any("site" in c or "treat" in c for c in spec.columns)
    assert {"log2_b12", "log2_folate", "log2_cysc", "H_x_W"} <= set(spec.columns)


def test_imputation_preserves_observed_exposures_and_categories():
    cfg = load_config("config/analysis.yml")
    cfg["analysis"].update({"imputations": 2, "mice_iterations": 2})
    d, _ = simulation_data(91, n=250)
    d = d.reset_index(drop=True)
    d.loc[0:4, "b12"] = np.nan
    d.loc[8:12, "smoking"] = np.nan
    frames, _meta = impute(d, cfg)
    assert len(frames) == 2
    for f in frames:
        assert f.b12.notna().all()
        assert set(f.smoking.unique()) <= {1, 2, 3, 4}
        assert np.array_equal(f.hcy, d.hcy)
        assert np.array_equal(f.wmh_ml, d.wmh_ml)
        assert np.array_equal(f.exit, d.exit)
        observed = d.b12.notna()
        assert np.array_equal(f.loc[observed, "b12"], d.loc[observed, "b12"])


def test_whole_unmeasured_b12_not_imputed():
    cfg = load_config("config/analysis.yml")
    d, _ = simulation_data(90, n=200)
    d["b12"] = np.nan
    with pytest.raises(DataError, match="entirely unavailable"):
        impute(d, cfg)
