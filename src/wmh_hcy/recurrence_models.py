"""Shared-scale recurrence Cox fits and covariance-based conditional contrasts."""
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .common import DataError, dump_json
from .design import Design
from .imputation import pool_scalar
from .models import fit_cause, pool_coefficients
from .recurrence_data import HORIZONS, split_at_90, truncate

TIME_TERMS = ("H", "H_rcs", "W", "W_rcs", "H_x_W")


def freeze_design(data, parent=None, *, expanded=False):
    spec = Design(expanded=expanded).fit(data, (0.1, 0.5, 0.9))
    if spec.constants:
        raise DataError(f"Constant model columns: {spec.constants}; no automatic removal in v3")
    if parent is not None:
        spec.centers, spec.knots = copy.deepcopy(parent.centers), copy.deepcopy(parent.knots)
        spec.transform(data)
    return spec


@dataclass
class TimeDesign:
    base: Design

    @property
    def columns(self):
        return self.base.columns + [f"{c}_late" for c in TIME_TERMS]

    def transform(self, data):
        x = self.base.transform(data)
        return pd.concat([x, pd.DataFrame({f"{c}_late": x[c]*data.late for c in TIME_TERMS},
                                          index=data.index)], axis=1)

    def to_dict(self):
        return {"base": self.base.to_dict(), "columns": self.columns,
                "split_day_from_onset": 90, "late_terms": list(TIME_TERMS)}


def pooled_linear(fits, vector):
    c = np.asarray(vector, dtype=float)
    result = pool_scalar([float(c @ f.params) for f in fits],
                         [float(c @ f.covariance @ c) for f in fits])
    result.update(HR=float(np.exp(result["estimate"])), HR_lower=float(np.exp(result["lower"])),
                  HR_upper=float(np.exp(result["upper"])))
    return result


def fit_series(frames, spec, folder, label, month=60, *, time_varying=False):
    """Fit every completed frame. Never discard a failed imputation when pooling."""
    folder.mkdir(parents=True, exist_ok=True)
    first = truncate(frames[0], month)
    result = {"analysis": label, "month": month, "n": len(first),
              "events": int(first.event_type.eq(1).sum()), "imputations": len(frames),
              "parameters": len(spec.columns), "status": "NOT_ESTIMABLE",
              "estimate": np.nan, "lower": np.nan, "upper": np.nan, "p": np.nan,
              "HR": np.nan, "HR_lower": np.nan, "HR_upper": np.nan}
    fits = []
    dump_json(folder / "design.json", spec.to_dict())
    first[["patient_id", "entry", "exit", "event_type"]].to_csv(folder / "participants.csv", index=False)
    try:
        for index, frame in enumerate(frames):
            d = truncate(frame, month)
            if not d.patient_id.equals(first.patient_id):
                raise DataError("Imputed frames do not share identical membership/order")
            if time_varying:
                d = split_at_90(d)
            fit = fit_cause(d, spec, 1)
            if time_varying:
                fit.diagnostics["n_patients"] = len(first)
                fit.diagnostics["influence_unit"] = "interval row, descriptive; not patient deletion"
            fits.append(fit)
            if (index+1) % 10 == 0 or index+1 == len(frames):
                print(f"  {label}: {index+1}/{len(frames)} Cox fits", flush=True)
        table = pool_coefficients(fits, label)
        table.to_csv(folder / "coefficients.csv", index=False)
        primary_term = "H_x_W_late" if time_varying else "H_x_W"
        result.update(table.set_index("term").loc[primary_term].to_dict(), status="ESTIMATED")
        np.savez_compressed(folder / "fit_parameters_covariance.npz",
                            parameters=np.array([f.params for f in fits]),
                            covariance=np.array([f.covariance for f in fits]),
                            columns=np.asarray(spec.columns, dtype=str))
        if time_varying:
            rows = []
            for period in ["early_to_day90", "late_day91_to1825", "late_minus_early"]:
                c = np.zeros(len(spec.columns))
                if period != "late_minus_early":
                    c[spec.columns.index("H_x_W")] = 1
                if period != "early_to_day90":
                    c[spec.columns.index("H_x_W_late")] = 1
                rows.append({"period": period, **pooled_linear(fits, c)})
            pd.DataFrame(rows).to_csv(folder / "time_interactions.csv", index=False)
            result["estimand"] = "late_vs_early_ratio_of_interaction_HR_ratios"
        else:
            result["estimand"] = "HR_ratio_per_Hcy_doubling_and_WMH_log_SD"
    except (DataError, ValueError, np.linalg.LinAlgError) as exc:
        result.update(reason=str(exc), failed_imputation=len(fits))
    dump_json(folder / "diagnostics.json", [f.diagnostics for f in fits])
    dump_json(folder / "estimate.json", result)
    return result, fits if result["status"] == "ESTIMATED" else []


def clinical_contrasts(data, spec, fits, quantiles=None):
    """Hcy 15 vs 10, full covariance; local observed support assessed before plotting."""
    quantiles = [0.25, 0.5, 0.75] if quantiles is None else quantiles
    rows = []
    rank = data.wmh_ml.rank(pct=True)
    for q in quantiles:
        wmh = float(data.wmh_ml.quantile(q))
        local = data.loc[rank.between(max(0, q-.10), min(1, q+.10)), "hcy"]
        lo, hi = local.quantile([.05, .95]) if len(local) else (np.nan, np.nan)
        r = {"wmh_quantile": q, "wmh_ml": wmh, "local_n": len(local),
             "hcy_local_p05": lo, "hcy_local_p95": hi, "low_hcy": 10, "high_hcy": 15,
             "status": "UNSUPPORTED", "HR": np.nan, "HR_lower": np.nan, "HR_upper": np.nan}
        if np.isfinite(lo) and lo <= 10 and hi >= 15:
            low = data.iloc[:1].copy().assign(hcy=10.0, wmh_ml=wmh)
            high = low.assign(hcy=15.0)
            c = (spec.transform(high)-spec.transform(low)).iloc[0].to_numpy()
            r.update(pooled_linear(fits, c), status="ESTIMATED")
        rows.append(r)
    return pd.DataFrame(rows)


def horizon_table(rows):
    table = pd.DataFrame(rows).set_index("month").reindex(HORIZONS).reset_index()
    p = pd.to_numeric(table.p, errors="coerce")
    sensitivity = table.month.ne(60)
    adjusted = multipletests(p.loc[sensitivity].fillna(1), method="holm")[1]
    table["p_holm_6"] = np.nan
    table.loc[sensitivity, "p_holm_6"] = np.where(p.loc[sensitivity].notna(), adjusted, np.nan)
    table["role"] = np.where(table.month.eq(60), "primary_after_protocol_revision", "sensitivity")
    return table
