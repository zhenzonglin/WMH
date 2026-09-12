"""Small MAR stress check with a prognostic, Hcy-correlated missing covariate."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from test_statistics import simulation_data

from wmh_hcy.common import load_config
from wmh_hcy.design import Design
from wmh_hcy.imputation import impute, pool_scalar
from wmh_hcy.models import fit_cause

root = Path(__file__).resolve().parents[1]
cfg = load_config(root / "config/analysis.yml")
cfg["analysis"].update(imputations=5, mice_iterations=5)
records = []
for beta in [0., .65]:
    for replicate in range(5):
        rng = np.random.default_rng(7230 + replicate)
        full, _ = simulation_data(7230 + replicate, n=1400)
        full = full.reset_index(drop=True)
        full["hcy"] *= (full.b12 / full.b12.median()) ** -.5
        spec = Design().fit(full)
        x = spec.transform(full)
        eta = .2*x.H + .15*x.W + beta*x.H_x_W + .45*x.log2_b12 + .3*x.log2_cysc
        eta = eta - eta.mean()
        event = rng.exponential(1/(.003*np.exp(eta)))
        death = rng.exponential(1/(.0005*np.exp(.3*x.log2_cysc)))
        full["exit"] = np.minimum(np.minimum(event, death), 365)
        full["event_type"] = np.where((event <= death) & (event <= 365), 1,
                                      np.where(death <= 365, 2, 0))
        full = full.loc[full.exit.gt(full.entry)].reset_index(drop=True)
        complete_fit = fit_cause(full, spec, 1)
        masked = full.copy()
        probability = expit(-1.5 + .4*np.log2(masked.hcy) - .4*3.8
                            + .5*masked.event_type.eq(1))
        missing = rng.uniform(size=len(masked)) < probability
        masked.loc[missing, "b12"] = np.nan
        cfg["analysis"]["seed"] = 900 + replicate
        completed, _ = impute(masked, cfg)
        estimates, variances = [], []
        for d in completed:
            fit = fit_cause(d, spec, 1)
            j = fit.columns.index("H_x_W")
            estimates.append(fit.params[j]); variances.append(fit.covariance[j, j])
        pooled = pool_scalar(estimates, variances)
        cc = fit_cause(masked.dropna(subset=["b12"]), spec, 1)
        j = complete_fit.columns.index("H_x_W")
        records.append({"true_beta": beta, "replicate": replicate, "n": len(full),
                        "fraction_b12_missing": missing.mean(),
                        "full_data": complete_fit.params[j], "mi": pooled["estimate"],
                        "mi_lower": pooled["lower"], "mi_upper": pooled["upper"],
                        "complete_case": cc.params[j]})
    print(f"MAR checks completed: beta={beta}", flush=True)
frame = pd.DataFrame(records)
out = root / "outputs/validation"
out.mkdir(exist_ok=True, parents=True)
frame.to_csv(out / "mi_stress_replicates.csv", index=False)
summary = frame.assign(mi_minus_full=frame.mi-frame.full_data).groupby("true_beta").agg(
    replicates=("replicate", "count"), full_mean=("full_data", "mean"), mi_mean=("mi", "mean"),
    cc_mean=("complete_case", "mean"), mean_mi_minus_full=("mi_minus_full", "mean"),
    max_abs_mi_minus_full=("mi_minus_full", lambda x: x.abs().max())).reset_index()
(out / "mi_stress_summary.json").write_text(json.dumps({"synthetic_only": True,
    "scope": "10 small MAR checks, prognostic B12 correlated with Hcy, 5 imputations each. Not a clinical MI validation or MNAR sensitivity analysis.",
    "results": summary.to_dict("records")}, indent=2))
print(summary.to_string(index=False))
