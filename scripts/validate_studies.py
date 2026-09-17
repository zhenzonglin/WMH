"""Monte Carlo verification, synthetic data only; sampling uncertainty is reported."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, softmax
from statsmodels.stats.proportion import proportion_confint

from wmh_hcy.studies.design import StudyDesign
from wmh_hcy.studies.imputation import impute
from wmh_hcy.studies.models import fit
from wmh_hcy.studies.pooling import coefficients, terms_test
from wmh_hcy.studies.registry import ModelSpec


def simulate(repetitions=60, n=900, seed=20260917, families=None):
    rng = np.random.default_rng(seed)
    rows = []
    for family in (families or ("ols", "multinomial", "cox", "ols_mi", "cox_interaction")):
        for alternative in (False, True):
            target = .4 if alternative else 0.0
            estimates, covered, rejected, failures = [], [], [], []
            for rep in range(repetitions):
                try:
                    d = pd.DataFrame({"patient_id": [str(i) for i in range(n)], "x": rng.normal(size=n), "z": rng.normal(size=n)})
                    model_family = ("ols" if family == "ols_mi" else "cox" if family == "cox_interaction"
                                    else "multinomial" if family == "kidney_interaction" else family)
                    spec = ModelSpec("simulation", family=model_family, outcome="y", exposures=("x",),
                                     covariates=("z",), splines=(), primary=(("dependent:x",) if family == "multinomial" else ("x",)))
                    d.z += .4*d.x
                    if model_family == "ols":
                        d["y"] = target*d.x+.5*d.z+rng.normal(size=n)*(1+.2*np.abs(d.x))
                    elif family == "kidney_interaction":
                        d["wmh_ml"] = np.exp(d.x)
                        d["albuminuria"] = rng.integers(0, 4, n)
                        spec = ModelSpec("simulation", family="multinomial", outcome="y",
                                         exposures=("wmh_ml", "albuminuria"), covariates=("z",), splines=(),
                                         interactions=tuple((f"albuminuria_{i}", "wmh_ml") for i in (1, 2, 3)),
                                         primary=("dependent:albuminuria_3_x_wmh_ml",))
                        fixed = StudyDesign.freeze(d, spec)
                        x = fixed.transform(d)
                        eta = -.4+.2*x.wmh_ml+.3*x.albuminuria_3+.2*x.z+target*x.albuminuria_3_x_wmh_ml
                        p = softmax(np.column_stack([np.zeros(n), eta, -.6+.1*x.wmh_ml]), axis=1)
                        d["y"] = [rng.choice(3, p=v) for v in p]
                    elif family == "multinomial":
                        p = softmax(np.column_stack([np.zeros(n), -.5+target*d.x+.4*d.z, -.7+.2*d.z]), axis=1)
                        d["y"] = [rng.choice(3, p=v) for v in p]
                    else:
                        if family == "cox_interaction":
                            spec = ModelSpec("simulation", family="cox", outcome="y", exposures=("x", "z"),
                                             splines=("x", "z"), interactions=(("x", "z"), ("x_rcs", "z")),
                                             primary=("x_x_z", "x_rcs_x_z"))
                            fixed = StudyDesign.freeze(d, spec)
                            matrix = fixed.transform(d)
                            eta = .2*matrix.x+.2*matrix.z+target*matrix.x_x_z-(.1 if alternative else 0)*matrix.x_rcs_x_z
                        else:
                            eta = target*d.x+.3*d.z
                        time = rng.exponential(3000*np.exp(-eta))
                        censor = rng.uniform(1000, 2000, n)
                        d["entry"] = rng.uniform(30, 150, n)
                        d["exit"] = np.minimum(time, censor)
                        d["y"] = d["event_type"] = (time <= censor).astype(int)
                        d = d.loc[d.exit.gt(d.entry)].reset_index(drop=True)
                    if family == "ols_mi":
                        d.loc[rng.random(len(d)) < expit(-1.6+.5*d.x), "z"] = np.nan
                    design = fixed if family in {"cox_interaction", "kidney_interaction"} else StudyDesign.freeze(d, spec)
                    completed = impute(d, design, {"imputations": 5, "mice_iterations": 3,
                                                   "mice_threads": 2, "seed": seed+rep}, progress=lambda _: None)[0] if family == "ols_mi" else [d]
                    fits = [fit(a, design) for a in completed]
                    row = coefficients(fits, model_family).set_index("term").loc[spec.primary[0]]
                    estimates.append(row.estimate)
                    covered.append(row.lower <= target <= row.upper)
                    rejected.append(terms_test(fits, spec.primary)["p"] < .05)
                except (ValueError, np.linalg.LinAlgError) as exc:
                    failures.append(str(exc))
            valid = len(estimates)
            coverage_ci = proportion_confint(sum(covered), valid, method="wilson") if valid else [np.nan]*2
            rejection_ci = proportion_confint(sum(rejected), valid, method="wilson") if valid else [np.nan]*2
            rows.append({"family": family, "scenario": "known_effect" if alternative else "null", "true_effect": target,
                         "requested": repetitions, "valid": valid, "failed": len(failures),
                         "bias": np.mean(estimates)-target, "bias_mcse": np.std(estimates, ddof=1)/np.sqrt(valid),
                         "coverage": np.mean(covered), "coverage_mc_lower": coverage_ci[0], "coverage_mc_upper": coverage_ci[1],
                         "rejection_rate": np.mean(rejected), "rejection_mc_lower": rejection_ci[0], "rejection_mc_upper": rejection_ci[1],
                         "failures": "; ".join(set(failures))})
            print(f"{family} / {'known' if alternative else 'null'}: {valid}/{repetitions}", flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repetitions", type=int, default=60)
    parser.add_argument("--n", type=int, default=900)
    parser.add_argument("--output", default="outputs/synthetic/studies_validation")
    parser.add_argument("--family", choices=["ols", "multinomial", "cox", "ols_mi", "cox_interaction", "kidney_interaction"], action="append")
    args = parser.parse_args()
    folder = Path(args.output)
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / "monte_carlo.csv").exists():
        parser.error("Validation destination already contains results; choose a new directory")
    result = simulate(args.repetitions, args.n, families=args.family)
    result.to_csv(folder / "monte_carlo.csv", index=False)
    (folder / "manifest.json").write_text(json.dumps({"mode": "synthetic", "seed": 20260917,
        "repetitions_per_scenario": args.repetitions, "n": args.n,
        "note": "Model/MI verification under specified DGPs, not proof of MAR or adequacy for patient data; Wilson intervals quantify Monte Carlo uncertainty."}, indent=2), encoding="utf-8")
    print(result.to_string(index=False))
