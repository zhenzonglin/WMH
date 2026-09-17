"""Frozen coding, transformations, spline knots and full-covariance contrasts."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..common import DataError
from ..design import rcs_nonlinear
from .registry import CATEGORIES, LOG_COLUMNS, ModelSpec


def raw_transform(values, name):
    v = pd.to_numeric(values, errors="coerce").astype(float)
    if name in {"wmh_ml", "wmh_raw_ml", "lesion_ml"}:
        return np.log1p(v)
    if name in LOG_COLUMNS:
        return np.log2(v)
    return v


@dataclass
class StudyDesign:
    spec: ModelSpec
    coding: dict

    @classmethod
    def freeze(cls, data, spec, inherited=None):
        coding = dict(inherited.coding) if inherited else {}
        for name in spec.predictors:
            if name in coding:
                continue
            observed = pd.to_numeric(data[name], errors="coerce").dropna()
            if observed.empty:
                raise DataError(f"Entirely unmeasured predictor: {name}; no surrogate or automatic deletion")
            if name in CATEGORIES:
                levels = sorted(observed.unique().tolist())
                if not set(levels).issubset(CATEGORIES[name]):
                    raise DataError(f"Invalid categorical codes: {name}")
                if len(levels) < 2:
                    raise DataError(f"Constant covariate: {name}; model is not estimable")
                coding[name] = {"kind": "category", "levels": levels, "reference": levels[0]}
            else:
                values = raw_transform(observed, name)
                # These are effect scales, not learned variable selection.
                standardized = name in {"wmh_ml", "wmh_raw_ml", "gm119_ml", "cec"}
                scale = float(values.std(ddof=1)) if standardized else (
                    100.0 if name == "icv_ml" else 10.0 if name in {"age", "sbp3", "sbp0"} else 1.0)
                if not np.isfinite(scale) or scale <= 0 or values.nunique() < 2:
                    raise DataError(f"Constant/nonfinite predictor: {name}")
                center = float(values.mean())
                z = (values - center) / scale
                coding[name] = {"kind": "continuous", "center": center, "scale": scale,
                                "observed_range": [float(observed.min()), float(observed.max())],
                                "knots": np.quantile(z, [.1, .5, .9]).tolist()}
            if name in spec.splines and len(set(coding[name]["knots"])) < 3:
                raise DataError(f"Three distinct spline knots unavailable: {name}")
        return cls(spec, coding)

    def transform(self, data):
        columns = {}
        if self.spec.family != "cox":
            columns["intercept"] = np.ones(len(data))
        for name in self.spec.predictors:
            code = self.coding[name]
            v = pd.to_numeric(data[name], errors="coerce")
            if v.isna().any():
                raise DataError(f"Unresolved missing predictor {name}")
            if code["kind"] == "category":
                if not v.isin(code["levels"]).all():
                    raise DataError(f"Category outside frozen design: {name}")
                for level in code["levels"][1:]:
                    columns[f"{name}_{level:g}"] = v.eq(level).astype(float).to_numpy()
            else:
                z = (raw_transform(v, name) - code["center"]) / code["scale"]
                columns[name] = z.to_numpy()
                if name in self.spec.splines:
                    columns[name + "_rcs"] = rcs_nonlinear(z.to_numpy(), code["knots"])
        for left, right in self.spec.interactions:
            columns[left + "_x_" + right] = columns[left] * columns[right]
        # A counting-process sensitivity permits only the prespecified terms to vary.
        if "late_period" in data:
            for term in ("sbp3", "sbp3_rcs", "wmh_ml", "wmh_ml_rcs",
                         "sbp3_x_wmh_ml", "sbp3_rcs_x_wmh_ml"):
                if term in columns:
                    columns[term + "_late"] = columns[term] * data.late_period.to_numpy()
        x = pd.DataFrame(columns, index=data.index, dtype=float)
        if not np.isfinite(x.to_numpy()).all():
            raise DataError("Nonfinite design matrix")
        return x

    def contrast(self, data, changes, reference):
        a, b = data.iloc[:1].copy(), data.iloc[:1].copy()
        for name, value in changes.items():
            a[name] = value
        for name, value in reference.items():
            b[name] = value
        return (self.transform(a) - self.transform(b)).iloc[0].to_numpy()


def split_time(data, cutoff=365):
    """Attained-time intervals (entry,exit]; never duplicate the event."""
    early = data.loc[data.entry.lt(cutoff)].copy()
    early["exit"] = early.exit.clip(upper=cutoff)
    early["event_type"] = early.event_type.where(data.loc[early.index, "exit"].le(cutoff), 0)
    early["late_period"] = 0
    late = data.loc[data.exit.gt(cutoff)].copy()
    late["entry"] = late.entry.clip(lower=cutoff)
    late["late_period"] = 1
    return pd.concat([early, late], ignore_index=True)
