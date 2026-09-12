"""Limited Monte Carlo software check, not a power study for CNSR-III."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from test_statistics import simulation_data

from wmh_hcy.models import fit_cause

root = Path(__file__).resolve().parents[1]
records = []
for beta in [0.0, 0.65]:
    for replicate in range(20):
        data, spec = simulation_data(62000 + replicate, n=1000, beta=beta)
        fit = fit_cause(data, spec, 1)
        j = fit.columns.index("H_x_W")
        estimate = float(fit.params[j])
        se = float(np.sqrt(fit.covariance[j, j]))
        records.append({"true_beta": beta, "replicate": replicate, "estimate": estimate,
                        "se": se, "covered": abs(estimate - beta) <= 1.96 * se,
                        "events": int(data.event_type.sum()), "n": len(data)})
    print(f"Completed 20 simulations: true interaction {beta}", flush=True)
frame = pd.DataFrame(records)
folder = root / "outputs/validation"
folder.mkdir(parents=True, exist_ok=True)
frame.to_csv(folder / "simulation_replicates.csv", index=False)
summary = []
for beta, group in frame.groupby("true_beta"):
    summary.append({"true_beta": beta, "replicates": len(group),
                    "mean_estimate": group.estimate.mean(), "bias": group.estimate.mean() - beta,
                    "empirical_sd": group.estimate.std(), "mean_se": group.se.mean(),
                    "coverage": group.covered.mean(),
                    "bias_mcse": group.estimate.std() / np.sqrt(len(group))})
(folder / "simulation_summary.json").write_text(json.dumps({"synthetic_only": True,
    "scope": "40 finite-sample numerical checks; not a clinical power study or proof of MI validity",
    "scenarios": summary}, indent=2))
print(json.dumps(summary, indent=2))
