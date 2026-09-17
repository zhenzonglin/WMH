"""Isolated solver checks on prepared cohorts; no risk bootstrap or result-pointer changes."""
from __future__ import annotations

import warnings
from collections import Counter
from datetime import UTC, datetime

import numpy as np
import pandas as pd

from .analysis import load_cohort
from .common import DataError, dump_json, outdir, sha256
from .design import Design
from .imputation import impute
from .models import fit_cause


def check_fit(cfg: dict, hypothesis: str = "H3") -> dict:
    if hypothesis not in {"H2", "H3"}:
        raise DataError("check-fit supports H2 or H3 only")
    month3 = hypothesis == "H3"
    name, kind = ("month3", "month3") if month3 else ("main", "main")
    data = load_cohort(cfg, name)
    if data.empty:
        raise DataError("Empty prepared cohort")
    base = outdir(cfg)
    folder = base / "diagnostics" / ("fit_check_" + datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ"))
    folder.mkdir(parents=True, exist_ok=False)
    summary = {"mode": cfg["mode"], "hypothesis": hypothesis, "n": len(data),
               "ischemic_events": int(data.event_type.eq(1).sum()),
               "death_events": int(data.event_type.eq(2).sum()), "result_dir": str(folder),
               "prepared_sha256": sha256(base / f"prepared/cohort_{name}.csv"),
               "configuration": {key: cfg["analysis"][key] for key in
                                 ["seed", "imputations", "mice_iterations", "spline_quantiles"]},
               "status": "STARTED", "fits": []}
    dump_json(folder / "summary.json", summary)
    print(f"FIT CHECK | {hypothesis} | mode={cfg['mode']} | n={len(data)} | "
          f"IS={summary['ischemic_events']} | death={summary['death_events']}", flush=True)
    print("Prepared cohort only. Recomputes configured MI; no bootstrap, no result-pointer changes.", flush=True)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning, module=r"miceforest\..*")
        completed, mi = impute(data, cfg, kind)
    dump_json(folder / "imputation.json", mi)
    spec = Design(month3=month3).fit(completed[0], cfg["analysis"]["spline_quantiles"])
    dump_json(folder / "design.json", spec.to_dict())
    print(f"MI complete: {len(completed)} datasets; {len(spec.columns)} model parameters.", flush=True)
    for index, frame in enumerate(completed):
        for cause in (1, 2):
            row = {"imputation": index, "cause": cause}
            try:
                fit = fit_cause(frame, spec, cause)
                row.update(status="PASS", diagnostics=fit.diagnostics)
            except (DataError, ValueError, np.linalg.LinAlgError) as exc:
                row.update(status="FAIL", reason=str(exc))
            summary["fits"].append(row)
        if (index + 1) % 10 == 0:
            print(f"Checked {index + 1}/{len(completed)} imputed datasets.", flush=True)
    for cause, label in [(1, "ISCHEMIC"), (2, "DEATH")]:
        rows = [r for r in summary["fits"] if r["cause"] == cause]
        passed = [r for r in rows if r["status"] == "PASS"]
        solvers = Counter(r["diagnostics"]["numerical_fit"]["optimizer"] for r in passed)
        print(f"{label}: pass={len(passed)}/{len(rows)}; "
              f"original={solvers['newton']}; recovered={solvers['scaled_bfgs_then_newton']}")
        if passed:
            grad = max(r["diagnostics"]["numerical_fit"]["max_scaled_score_per_event"] for r in passed)
            cond = max(r["diagnostics"]["numerical_fit"]["scaled_information_condition"] for r in passed)
            print(f"  max scaled score/event={grad:.3g}; max information condition={cond:.3g}")
        failures = [r for r in rows if r["status"] == "FAIL"]
        if failures:
            groups = Counter("singular_or_ill_conditioned" if any(
                word in r["reason"].lower() for word in ("singular", "rank-deficient", "ill-conditioned"))
                else "convergence_or_invalid_estimate" for r in failures)
            print("  failures: " + ", ".join(f"{key}={count}" for key, count in groups.items()))
    summary["status"] = "PASS" if all(r["status"] == "PASS" for r in summary["fits"]) else "FIT_FAILURES"
    dump_json(folder / "summary.json", summary)
    print("No variables removed, no penalization; passing checks does not establish model adequacy.")
    print("Saved diagnostic summary:", folder / "summary.json")
    return summary
