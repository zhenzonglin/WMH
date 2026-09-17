"""Python MICE, with observed exposures and delayed-entry cause-specific hazard auxiliaries."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .adjustment import adjustment_columns
from .common import DataError
from .fields import CATEGORIES


def nelson_aalen_increment(data: pd.DataFrame, cause: int) -> np.ndarray:
    times, counts = np.unique(data.loc[data.event_type.eq(cause), "exit"], return_counts=True)
    entry = data.entry.to_numpy(float)
    stop = data.exit.to_numpy(float)
    total = np.zeros(len(data))
    for time, deaths in zip(times, counts, strict=True):
        risk = (entry < time) & (stop >= time)
        if risk.sum() == 0:
            raise DataError("Empty risk set at an observed event")
        total += ((entry < time) & (stop >= time)) * deaths / risk.sum()
    return total


def required_covariates(kind: str = "main") -> list[str]:
    return adjustment_columns(kind)


def impute(data: pd.DataFrame, cfg: dict, kind: str = "main", *,
           survival_auxiliaries: bool = True, death_auxiliaries: bool = True,
           functional_outcome: str = "mrs12",
           extra_auxiliaries: pd.DataFrame | None = None,
           progress=None) -> tuple[list[pd.DataFrame], dict]:
    cols = required_covariates(kind)
    for col in cols:
        if col not in data or data[col].notna().sum() == 0:
            raise DataError(f"{kind}: {col} entirely unavailable; no whole-substudy imputation")
    working = data[cols].copy()
    for c in cols:
        if c in CATEGORIES:
            working[c] = pd.Categorical(working[c], categories=CATEGORIES[c])
        elif c in {"b12", "folate", "cysc", "creatinine", "b123", "folate3", "cysc3"}:
            if (working[c].dropna() <= 0).any():
                raise DataError(f"Nonpositive {c}; review units/values before imputation")
            working[c] = np.log2(working[c])
    hcol = "hcy3" if kind == "month3" else "hcy"
    working["_H"] = np.log2(data[hcol])
    working["_W"] = np.log1p(data.wmh_ml)
    working["_HW"] = (working._H - working._H.mean()) * (working._W - working._W.mean())
    working["_H0"] = np.log2(data.hcy)
    if kind == "month3":
        working["_H0W"] = (working._H0-working._H0.mean())*(working._W-working._W.mean())
        working["_sample3_day"] = data.sample3_day
    survival = kind != "cross_sectional" and survival_auxiliaries
    if survival:
        working["_D1"] = data.event_type.eq(1).astype(int)
        if death_auxiliaries:
            working["_D2"] = data.event_type.eq(2).astype(int)
        working["_NA1"] = nelson_aalen_increment(data, 1)
        if death_auxiliaries:
            working["_NA2"] = nelson_aalen_increment(data, 2)
        working["_entry"] = data.entry
        working["_stop"] = data.exit
    if kind.startswith("functional"):
        working["_mrs"] = data[functional_outcome]
    if extra_auxiliaries is not None:
        if not extra_auxiliaries.index.equals(data.index):
            raise DataError("Imputation auxiliaries must have the exact cohort index/order")
        if set(extra_auxiliaries) & set(working):
            raise DataError("Duplicate imputation auxiliary names")
        if not np.isfinite(extra_auxiliaries.to_numpy(float)).all():
            raise DataError("Non-finite imputation auxiliary values")
        working = pd.concat([working, extra_auxiliaries], axis=1)
    working = working.reset_index(drop=True)
    missing = {c: int(data[c].isna().sum()) for c in cols}
    if not any(missing.values()):
        return [data.reset_index(drop=True).copy()], {"method": "no_missing_covariates", "m": 1}
    import miceforest as mf
    targets = {c: [x for x in working if x != c] for c in cols if working[c].isna().any()}
    kernel = mf.ImputationKernel(working, num_datasets=cfg["analysis"]["imputations"],
                                 variable_schema=targets, mean_match_candidates=5,
                                 random_state=cfg["analysis"]["seed"])
    traces = []
    for iteration in range(cfg["analysis"]["mice_iterations"]):
        kernel.mice(1, num_threads=cfg["analysis"]["mice_threads"], num_iterations=60, verbosity=-1)
        if progress is not None:
            progress(f"MI iteration {iteration+1}/{cfg['analysis']['mice_iterations']}; "
                     f"{cfg['analysis']['imputations']} completed datasets")
        for index in range(cfg["analysis"]["imputations"]):
            current = kernel.complete_data(dataset=index)
            traces.append({"iteration": iteration + 1, "imputation": index,
                           **{c: float(pd.to_numeric(current.loc[working[c].isna(), c]).mean())
                              for c in targets}})
    completed = []
    summaries = []
    for index in range(cfg["analysis"]["imputations"]):
        filled = kernel.complete_data(dataset=index)
        d = data.reset_index(drop=True).copy()
        for col in cols:
            values = pd.to_numeric(filled[col], errors="raise").astype(float)
            if col in {"b12", "folate", "cysc", "creatinine", "b123", "folate3", "cysc3"}:
                values = 2**values
            # Never change an observed value, even by round-off.
            d[col] = d[col].where(d[col].notna(), values)
        completed.append(d)
        summaries.append({"imputation": index, **{c: float(d[c].mean()) for c in cols}})
    return completed, {"method": "miceforest_FCS_approximate", "m": len(completed),
                       "missing": missing, "predictors": targets, "completed_means": summaries,
                       "imputed_mean_trace_transformed_scale": traces,
                       "assumption": "MAR; approximate compatibility with Cox and nonlinear terms"}


def pool_scalar(estimates: list[float], variances: list[float]) -> dict:
    q = np.asarray(estimates, float)
    u = np.asarray(variances, float)
    if not np.isfinite(q).all() or not np.isfinite(u).all():
        raise DataError("Cannot pool failed or non-finite estimates")
    m = len(q)
    between = float(np.var(q, ddof=1)) if m > 1 else 0.0
    within = float(u.mean())
    total = within + (1 + 1/m) * between
    se = np.sqrt(total)
    df = (m-1) * (1 + within / ((1+1/m)*between))**2 if between > 0 and m > 1 else np.inf
    critical = stats.t.ppf(0.975, df)
    p = 2 * stats.t.sf(abs(q.mean()/se), df) if se > 0 else np.nan
    return {"estimate": float(q.mean()), "se": float(se), "lower": float(q.mean()-critical*se),
            "upper": float(q.mean()+critical*se), "p": float(p), "df": float(df),
            "m": m, "within_variance": within, "between_variance": between}
