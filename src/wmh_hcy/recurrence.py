"""Independent five-year recurrence workflow; legacy outputs are never touched."""
from __future__ import annotations

import copy
from datetime import UTC, datetime

import numpy as np
import pandas as pd

from .adjustment import adjustment_columns, selection_registry
from .common import DataError, dump_json, outdir, record_run
from .harmonize import harmonize
from .imaging import available_volume, read_imaging
from .imputation import impute
from .recurrence_data import (
    CONTRACT,
    REQUIRED,
    SOURCES,
    build_cohort,
    counts,
    field_audit,
    horizon_auxiliaries,
    validate_covariates,
)
from .recurrence_models import (
    TimeDesign,
    clinical_contrasts,
    fit_series,
    freeze_design,
    horizon_table,
)
from .sas_extract import extract


def shared_imputation(data, cfg, folder, kind="main"):
    frames, meta = impute(data, cfg, kind, death_auxiliaries=False,
                          extra_auxiliaries=horizon_auxiliaries(data),
                          progress=lambda text: print(f"  {folder.name}: {text}", flush=True))
    dump_json(folder / "imputation.json", meta)
    return frames


def analyse_recurrence(data, cfg, root):
    validate_covariates(data)
    print(f"Recurrence v3: shared imputation, n={len(data)}, requested m={cfg['analysis']['imputations']}",
          flush=True)
    completed = shared_imputation(data, cfg, root)
    spec = freeze_design(completed[0])
    dump_json(root / "frozen_design.json", spec.to_dict())
    rows, diagnostic_rows, robustness = [], [], []

    def fit(frames, design, folder, label, month=60, **kwargs):
        r, fits = fit_series(frames, design, root / folder, label, month, **kwargs)
        for i, f in enumerate(fits):
            d = f.diagnostics
            diagnostic_rows.append({
                "analysis": label, "imputation": i, "patients": r["n"], "events": r["events"],
                "parameters": r["parameters"], "rank": d["rank"],
                "information_condition": d["numerical_fit"]["scaled_information_condition"],
                "max_scaled_score_per_event": d["numerical_fit"]["max_scaled_score_per_event"],
                "optimizer": d["numerical_fit"]["optimizer"], "warnings": "; ".join(d["warnings"]),
                **{f"PH_{c}_p_descriptive": d["schoenfeld_time_correlations"].get(c, {}).get("p_descriptive")
                   for c in ["H", "W", "H_x_W"]},
            })
        return r, fits

    # Fit the revised primary question first; all months use these exact completed frames/design.
    for month in [60, 3, 6, 12, 24, 36, 48]:
        print(f"Recurrence v3: month {month}", flush=True)
        result, fits = fit(completed, spec, f"models/month{month:02d}", f"month{month:02d}", month)
        rows.append(result)
        if month == 60 and fits:
            clinical_contrasts(completed[0], spec, fits).to_csv(root / "clinical_contrasts.csv", index=False)
            clinical_contrasts(completed[0], spec, fits, np.linspace(.05, .95, 37)).to_csv(
                root / "clinical_contrast_curve.csv", index=False)
    horizon_table(rows).to_csv(root / "horizon_results.csv", index=False)

    def failed(label, exc, n=None):
        return {"analysis": label, "month": 60, "status": "NOT_ESTIMABLE", "reason": str(exc),
                "n": n, "estimate": np.nan, "p": np.nan, "HR": np.nan,
                "HR_lower": np.nan, "HR_upper": np.nan}

    cc = data.dropna(subset=adjustment_columns()).reset_index(drop=True)
    result, _ = fit([cc], spec, "robustness/complete_case", "complete_case")
    robustness.append(result)

    # Raw WMH requires a compatible imputation model using the replacement exposure.
    raw = data.loc[available_volume(data.wmh_raw_ml, allow_zero=True)
                   & data.wmh_raw_ml.le(data.icv_ml)].copy().reset_index(drop=True)
    raw["wmh_ml"] = raw.wmh_raw_ml
    try:
        raw_frames = shared_imputation(raw, cfg, root / "robustness/raw_wmh")
        raw_spec = freeze_design(raw_frames[0])
        # Keep Hcy/age scales fixed; raw WMH has its own prespecified transformation.
        for name in ["H", "H0", "age"]:
            raw_spec.centers[name] = spec.centers[name]
        for name in ["H", "age10"]:
            raw_spec.knots[name] = spec.knots[name]
        result, _ = fit(raw_frames, raw_spec, "robustness/raw_wmh", "raw_wmh")
    except (DataError, ValueError, np.linalg.LinAlgError) as exc:
        result = failed("raw_wmh", exc, len(raw))
    robustness.append(result)

    # The lesion is an observed imaging requirement only for this sensitivity subset.
    acute = data.loc[available_volume(data.lesion_ml, allow_zero=True)].copy().reset_index(drop=True)
    try:
        acute_frames = shared_imputation(acute, cfg, root / "robustness/clinical_extended", "extended")
        extended_spec = freeze_design(acute_frames[0], parent=spec, expanded=True)
        for design, label in [(spec, "clinical_same_subset_core"), (extended_spec, "clinical_extended")]:
            r, _ = fit(acute_frames, design, f"robustness/{label}", label)
            robustness.append(r)
    except (DataError, ValueError, np.linalg.LinAlgError) as exc:
        for label in ["clinical_same_subset_core", "clinical_extended"]:
            robustness.append(failed(label, exc, len(acute)))

    result, _ = fit(completed, TimeDesign(spec), "robustness/time_varying", "time_varying",
                    time_varying=True)
    robustness.append(result)
    pd.DataFrame(robustness).to_csv(root / "robustness_results.csv", index=False)
    pd.DataFrame(diagnostic_rows).to_csv(root / "diagnostic_summary.csv", index=False)
    return {"horizons": rows, "robustness": robustness, "shared_imputations": len(completed)}


def run_recurrence(cfg, through="prepare"):
    if through not in {"audit", "prepare", "analyse", "report"}:
        raise DataError("Unknown recurrence stage")
    for name in ["imputations", "mice_iterations"]:
        value = cfg["analysis"][name]
        if not isinstance(value, int) or value < 1:
            raise DataError(f"analysis.{name} must be a positive integer")
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    target = outdir(cfg) / "recurrence_v3"
    root = target / "runs" / stamp
    work = copy.deepcopy(cfg)
    work["_out"] = str(root / "inputs_derived")
    status = {"status": "RUNNING", "mode": cfg["mode"], "contract": CONTRACT,
              "result_dir": str(root), "requested_through": through,
              "revision": "After examining year-one results; primary horizon revised to 60 months",
              "death_regression": False, "absolute_risk": False, "bootstrap": False}
    dump_json(root / "status.json", status)
    dump_json(root / "config_snapshot.json", cfg)
    dump_json(target / "latest_run.json", {"path": str(root)})
    try:
        if not work["inputs"].get("clinical_csv"):
            print("Recurrence v3: extracting its dedicated clinical whitelist", flush=True)
            extract(work, fields=SOURCES)
        clinical = harmonize(work, require_one_year=False, source_fields=SOURCES, required_sources=REQUIRED)
        audit = field_audit(clinical)
        audit.to_csv(root / "field_audit.csv", index=False)
        if through == "audit":
            status["status"] = "AUDITED"
        else:
            images = read_imaging(work, include_gm=False)
            master = clinical.merge(images, on="patient_id", how="left", validate="one_to_one")
            cohort, exclusions, flow, audit = build_cohort(master)
            cohort.to_csv(root / "cohort.csv", index=False)
            exclusions.to_csv(root / "exclusions.csv", index=False)
            flow.to_csv(root / "cohort_flow.csv", index=False)
            dump_json(root / "endpoint_audit.json", audit)
            table = counts(cohort)
            table.to_csv(root / "cohort_counts.csv", index=False)
            print(table.to_string(index=False), flush=True)
            status.update(status="PREPARED", counts=table.to_dict("records"), endpoint_audit=audit)
            validate_covariates(cohort)
            selected = set(adjustment_columns("extended"))
            pd.DataFrame([r for r in selection_registry() if r["variable"] in selected]).to_csv(
                root / "covariate_roles.csv", index=False)
            if through in {"analyse", "report"}:
                result = analyse_recurrence(cohort, work, root)
                status.update(result)
                failed = [r for r in result["horizons"]+result["robustness"] if r["status"] != "ESTIMATED"]
                status["status"] = "COMPLETED_WITH_MODEL_FAILURES" if failed else "COMPLETED"
                dump_json(target / "latest_results.json", {"path": str(root)})
        dump_json(root / "status.json", status)
        if through == "report":
            from .recurrence_reporting import report_recurrence
            status["report"] = str(report_recurrence(root))
        record_run(work, "recurrence", {"contract": CONTRACT, "status": status["status"]})
    except (DataError, OSError, ValueError, np.linalg.LinAlgError) as exc:
        status.update(status="FAILED", reason=str(exc))
        dump_json(root / "status.json", status)
        raise
    dump_json(root / "status.json", status)
    return status
