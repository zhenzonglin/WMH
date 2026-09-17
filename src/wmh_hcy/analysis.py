"""Frozen H1–H4 orchestration; every hypothesis has its own estimand and status."""
from __future__ import annotations

import copy
from datetime import UTC, datetime

import numpy as np
import pandas as pd

from .absolute_risk import cumulative_incidence, risk_analysis
from .adjustment import adjustment_columns, select_background, selection_registry
from .cohorts import COHORT_ENTRY_RULE
from .common import DataError, dump_json, outdir, read_json, record_run
from .design import Design
from .imaging import IMAGE_ELIGIBILITY_RULE, available_volume
from .imputation import impute, required_covariates
from .inference import decide_hypotheses, pooled_contrast, structural_curve
from .models import cross_sectional, fit_cause, ordinal, pool_coefficients, pool_regression


def load_cohort(cfg: dict, name: str) -> pd.DataFrame:
    contract = read_json(outdir(cfg) / "prepared/cohort_contract.json")
    if (contract.get("entry_rule") != COHORT_ENTRY_RULE
            or contract.get("image_eligibility_rule") != IMAGE_ELIGIBILITY_RULE):
        raise DataError("Prepared cohorts use older or unknown eligibility rules. Run prepare again after updating.")
    path = outdir(cfg) / f"prepared/cohort_{name}.csv"
    if not path.is_file():
        raise DataError("Prepared cohort absent. Run prepare first")
    return pd.read_csv(path, dtype={"patient_id": str, "participant_id": str})


def run_survival(data, cfg, folder, label, spec=None, kind="main", risk=False, complete_case=False,
                 death_auxiliaries=True):
    if data.empty:
        raise DataError("Empty analysis cohort")
    folder.mkdir(parents=True, exist_ok=True)
    data[["patient_id", "entry", "exit", "event_type"]].to_csv(folder / "analysis_participants.csv", index=False)
    completed, mi = ([data.reset_index(drop=True)], {"method": "complete_case", "m": 1}) if complete_case \
        else impute(data, cfg, kind, death_auxiliaries=death_auxiliaries)
    dump_json(folder / "imputation.json", mi)
    spec = (spec or Design()).fit(completed[0], cfg["analysis"]["spline_quantiles"])
    dump_json(folder / "design.json", spec.to_dict())
    fits = [fit_cause(d, spec, 1) for d in completed]
    table = pool_coefficients(fits, label)
    table.to_csv(folder / "coefficients.csv", index=False)
    dump_json(folder / "diagnostics.json", [f.diagnostics for f in fits])
    if kind == "month3":
        primary = pooled_contrast(completed[0], spec, fits)
        primary["estimand"] = "log_HR_Hcy3_15_vs_10_conditional_on_baseline"
        base = copy.deepcopy(spec)
        base.columns = [c for c in spec.columns if c not in ["H", "H_rcs"]]
        base_fits = [fit_cause(d, base, 1) for d in completed]
        pool_coefficients(base_fits, "H3_M0").to_csv(folder / "M0_coefficients.csv", index=False)
        dump_json(folder / "M0_design.json", base.to_dict())
    else:
        primary = table.set_index("term").loc["H_x_W"].to_dict()
        primary["estimand"] = "log_ratio_of_HRs_per_Hcy_doubling_and_WMH_log_SD"
    dump_json(folder / "hypothesis_estimate.json", primary)
    status = {"status": "ESTIMATED", "n": len(data), "ischemic_events": int(data.event_type.eq(1).sum()),
              "parameters": len(spec.columns), "imputations": len(completed), "estimand": primary["estimand"],
              "primary_p": primary["p"], "estimate": primary["estimate"]}
    if risk:
        risk_stage = "death_fit"
        risk_imputation = None
        death_diagnostics = []
        try:
            pairs = []
            for risk_imputation, (d, s) in enumerate(zip(completed, fits, strict=True)):
                death = fit_cause(d, spec, 2)
                pairs.append((s, death))
                death_diagnostics.append({"imputation": risk_imputation, **death.diagnostics})
            dump_json(folder / "death_diagnostics.json", death_diagnostics)
            risk_stage = "death_pooling"
            risk_imputation = None
            pool_coefficients([v[1] for v in pairs], label + "_death").to_csv(folder / "death_coefficients.csv", index=False)
            risk_stage = "risk_bootstrap"
            curves, intervals, differences, diag = risk_analysis(completed, spec, pairs, cfg)
            curves.to_csv(folder / "risk_curves.csv", index=False)
            intervals.to_csv(folder / "risk_intervals.csv", index=False)
            differences.to_csv(folder / "risk_differences.csv", index=False)
            dump_json(folder / "risk_diagnostics.json", diag)
            status["absolute_risk_status"] = diag["status"]
            if kind == "month3":
                risk_stage = "month3_M0_M1_risk_update"
                horizon = cfg["analysis"]["horizon"]
                m1 = [cumulative_incidence(d, spec, *pair, horizon) for d, pair in zip(completed, pairs, strict=True)]
                m0 = [cumulative_incidence(d, base, s, fit_cause(d, base, 2), horizon)
                      for d, s in zip(completed, base_fits, strict=True)]
                update = completed[0][["patient_id", "entry"]].copy()
                update["M0_fitted_risk"] = np.mean(m0, axis=0)
                update["M1_fitted_risk"] = np.mean(m1, axis=0)
                update.to_csv(folder / "risk_update_M0_M1.csv", index=False)
        except (DataError, ValueError, np.linalg.LinAlgError) as exc:
            status.update(absolute_risk_status="NOT_ESTIMABLE", absolute_risk_reason=str(exc))
            dump_json(folder / "death_diagnostics.json", death_diagnostics)
            dump_json(folder / "risk_failure.json", {"stage": risk_stage, "imputation": risk_imputation,
                      "reason": str(exc), "completed_death_fits": len(death_diagnostics)})
    return status


def run_structural(data, cfg, folder):
    if data.empty:
        raise DataError("Empty structural cohort")
    folder.mkdir(parents=True, exist_ok=True)
    frames, mi = impute(data, cfg, "cross_sectional")
    spec = Design().fit(frames[0], cfg["analysis"]["spline_quantiles"])
    tables = [cross_sectional(d, spec) for d in frames]
    data[["patient_id"]].to_csv(folder / "analysis_participants.csv", index=False)
    pool_regression(tables, "H1").to_csv(folder / "coefficients.csv", index=False)
    structural_curve(frames, spec, tables).to_csv(folder / "structural_curve.csv", index=False)
    primary = pooled_contrast(frames[0], spec, tables, kind="ols")
    primary["estimand"] = "difference_in_log1p_WMH_Hcy15_vs_10"
    dump_json(folder / "hypothesis_estimate.json", primary)
    dump_json(folder / "imputation.json", mi)
    dump_json(folder / "design.json", spec.to_dict())
    return {"status": "ESTIMATED", "n": len(data), "estimand": primary["estimand"],
            "primary_p": primary["p"], "estimate": primary["estimate"]}


def run_functional(data, cfg, folder, kind="functional", outcome="mrs12", survival_auxiliaries=True):
    data = data.loc[available_volume(data.lesion_ml, allow_zero=True)].copy()
    if data.empty:
        raise DataError("No patients with observed function and eligible acute lesion imaging")
    folder.mkdir(parents=True, exist_ok=True)
    frames, mi = impute(data, cfg, kind, functional_outcome=outcome,
                        survival_auxiliaries=survival_auxiliaries)
    spec = Design(expanded=True, functional=True, t1=kind == "functional_t1").fit(
        frames[0], cfg["analysis"]["spline_quantiles"])
    results = [ordinal(d, spec, outcome=outcome) for d in frames]
    data[["patient_id", outcome]].to_csv(folder / "analysis_participants.csv", index=False)
    table = pool_regression([r[0] for r in results], kind, exponentiate=True)
    table.to_csv(folder / "coefficients.csv", index=False)
    primary = table.set_index("term").loc["H_x_W"].to_dict()
    primary["estimand"] = "log_ratio_of_common_ORs_Hcy_by_WMH"
    dump_json(folder / "hypothesis_estimate.json", primary)
    pd.concat([r[2].assign(imputation=i) for i, r in enumerate(results)]).to_csv(
        folder / "threshold_diagnostics.csv", index=False)
    dump_json(folder / "diagnostics.json", [r[1] for r in results])
    dump_json(folder / "imputation.json", mi)
    dump_json(folder / "design.json", spec.to_dict())
    return {"status": "ESTIMATED", "n": len(data), "parameters": len(spec.columns),
            "estimand": primary["estimand"], "primary_p": primary["p"], "estimate": primary["estimate"]}


def analyse(cfg: dict, only: str | None = None) -> dict:
    from .hypotheses import ORDER, h01_structure, h02_recurrence, h03_month3_update, h04_function

    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    root = outdir(cfg) / "results" / stamp
    root.mkdir(parents=True, exist_ok=True)
    statuses = {"protocol_version": 2, "mode": cfg["mode"], "result_dir": str(root),
                "hypothesis_order": [x[0] for x in ORDER], "analyses": {}}
    dump_json(root / "config_snapshot.json", cfg)
    dump_json(outdir(cfg) / "latest_results.json", {"path": str(root)})
    dump_json(root / "adjustment_certificate.json", {"rule": "modified_disjunctive_cause_criterion",
              "registry": selection_registry(), "C0_background": select_background(),
              "D0_measurement": ["icv_ml", "sample_day"],
              "required_covariates_by_model": {"H1": adjustment_columns(), "H2": adjustment_columns(),
                 "H3": adjustment_columns("month3"), "H4": adjustment_columns("functional")},
              "H3_additional_conditioning": ["baseline_log2_hcy", "baseline_hcy_rcs", "baseline_H_x_W", "sample3_day_30"]})
    main = load_cohort(cfg, "main")
    pd.DataFrame({"variable": main.columns, "missing": main.isna().sum().to_numpy(),
                  "fraction": main.isna().mean().to_numpy()}).to_csv(root / "missingness.csv", index=False)

    def execute(label, function):
        print(f"Analysis: {label}", flush=True)
        try:
            statuses["analyses"][label] = function()
        except (DataError, ValueError, np.linalg.LinAlgError) as exc:
            statuses["analyses"][label] = {"status": "NOT_ESTIMABLE", "reason": str(exc)}
        dump_json(root / "status.json", statuses)

    sensitivity = cfg["analysis"]["run_sensitivity"]
    if only in [None, "H1"]:
        structural = load_cohort(cfg, "cross_sectional")
        execute("01_structure", lambda: h01_structure.run(structural, cfg, root / "01_structure"))
        if sensitivity:
            raw = structural.loc[structural.wmh_raw_ml.notna() & structural.wmh_raw_ml.ge(0)].copy()
            raw["wmh_ml"] = raw.wmh_raw_ml
            execute("01_structure/sensitivity_raw_wmh", lambda: run_structural(raw, cfg, root / "01_structure/sensitivity_raw_wmh"))
    if only in [None, "H2"]:
        execute("02_recurrence", lambda: h02_recurrence.run(main, cfg, root / "02_recurrence"))
        if sensitivity:
            base = root / "02_recurrence"
            cc = main.dropna(subset=required_covariates())
            execute("02_recurrence/complete_case", lambda: run_survival(cc, cfg, base / "complete_case", "H2_CC", complete_case=True))
            raw = main.loc[main.wmh_raw_ml.notna() & main.wmh_raw_ml.ge(0)].copy()
            raw["wmh_ml"] = raw.wmh_raw_ml
            execute("02_recurrence/raw_wmh", lambda: run_survival(raw, cfg, base / "raw_wmh", "H2_raw"))
            execute("02_recurrence/creatinine", lambda: run_survival(main, cfg, base / "creatinine", "H2_Cr",
                    spec=Design(renal="creatinine"), kind="creatinine"))
            extended = main.loc[available_volume(main.lesion_ml, allow_zero=True)].copy()
            execute("02_recurrence/acute_adjusted", lambda: run_survival(extended, cfg, base / "acute_adjusted",
                    "H2_acute", spec=Design(expanded=True), kind="extended"))
            early = main.loc[main.entry.lt(90)].copy()
            early.loc[early.exit.gt(90), "event_type"] = 0
            early["exit"] = early.exit.clip(upper=90)
            late = main.loc[main.exit.gt(90)].copy()
            late["entry"] = late.entry.clip(lower=90)
            execute("02_recurrence/early_90d", lambda: run_survival(early, cfg, base / "early_90d", "H2_early"))
            execute("02_recurrence/late_90_365d", lambda: run_survival(late, cfg, base / "late_90_365d", "H2_late"))
    if only in [None, "H3"]:
        execute("03_month3_update", lambda: h03_month3_update.run(load_cohort(cfg, "month3"), cfg, root / "03_month3_update"))
    if only in [None, "H4"]:
        execute("04_function", lambda: h04_function.run(load_cohort(cfg, "functional"), cfg, root / "04_function"))
        if cfg["analysis"].get("run_t1_extension", True):
            execute("04_function/t1_extension", lambda: h04_function.run(load_cohort(cfg, "functional_t1"), cfg,
                    root / "04_function/t1_extension", t1=True))
    rows = []
    for hypothesis, folder, statement in ORDER:
        estimate = read_json(root / folder / "hypothesis_estimate.json") or {}
        rows.append({"hypothesis": hypothesis, "statement": statement,
                     "estimate": estimate.get("estimate", np.nan), "lower": estimate.get("lower", np.nan),
                     "upper": estimate.get("upper", np.nan), "p": estimate.get("p", np.nan),
                     "estimand": estimate.get("estimand", "not estimated")})
    decide_hypotheses(rows).to_csv(root / "hypothesis_summary.csv", index=False)
    statuses["status"] = "COMPLETED_WITH_HYPOTHESIS_STATUS"
    dump_json(root / "status.json", statuses)
    record_run(cfg, "analyse", statuses)
    return statuses
