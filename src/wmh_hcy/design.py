"""Outcome-independent, frozen transformations and model matrices."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .common import DataError
from .fields import CATEGORIES


def rcs_nonlinear(x: np.ndarray, knots: list[float]) -> np.ndarray:
    a, b, c = knots
    if not a < b < c:
        raise DataError("Spline requires three distinct knots")
    return (np.maximum(x-a, 0)**3 - (c-a)/(c-b)*np.maximum(x-b, 0)**3
            + (b-a)/(c-b)*np.maximum(x-c, 0)**3) / (c-a)**2


@dataclass
class Design:
    month3: bool = False
    expanded: bool = False
    renal: str = "cysc"
    functional: bool = False
    t1: bool = False
    wmh_column: str = "wmh_ml"
    centers: dict = field(default_factory=dict)
    knots: dict = field(default_factory=dict)
    columns: list = field(default_factory=list)
    constants: list = field(default_factory=list)

    @property
    def hcy(self) -> str:
        return "hcy3" if self.month3 else "hcy"

    @property
    def nutrients(self) -> list[str]:
        return ["b123", "folate3", "cysc3"] if self.month3 else ["b12", "folate", self.renal]

    def fit(self, data: pd.DataFrame, quantiles=(0.1, 0.5, 0.9)) -> Design:
        self.centers = {
            "H": float(np.log2(data[self.hcy]).mean()),
            "W": float(np.log1p(data[self.wmh_column]).mean()),
            "W_sd": float(np.log1p(data[self.wmh_column]).std(ddof=1)),
            "age": float(data.age.mean()), "H0": float(np.log2(data.hcy).mean()),
        }
        if self.centers["W_sd"] <= 0:
            raise DataError("WMH has no variation")
        basis = self.basic(data)
        for name in ["H", "W", "age10"]:
            self.knots[name] = np.quantile(basis[name], quantiles).tolist()
        if self.month3:
            self.knots["H0"] = np.quantile(np.log2(data.hcy)-self.centers["H0"], quantiles).tolist()
        matrix = self._matrix(data)
        self.constants = [c for c in matrix if matrix[c].std() < 1e-12]
        self.columns = [c for c in matrix if c not in self.constants]
        self.transform(data)
        return self

    def basic(self, data: pd.DataFrame) -> pd.DataFrame:
        return pd.DataFrame({
            "H": np.log2(data[self.hcy]) - self.centers["H"],
            "W": (np.log1p(data[self.wmh_column]) - self.centers["W"]) / self.centers["W_sd"],
            "age10": (data.age - self.centers["age"]) / 10,
        }, index=data.index)

    def _matrix(self, data: pd.DataFrame) -> pd.DataFrame:
        x = self.basic(data)
        for name in ["H", "W", "age10"]:
            x[f"{name}_rcs"] = rcs_nonlinear(x[name].to_numpy(), self.knots[name])
        if not self.month3:
            x["H_x_W"] = x.H * x.W
        for name in self.nutrients:
            x[f"log2_{name}"] = np.log2(data[name])
        x["icv_100ml"] = data.icv_ml / 100
        x["sample_day"] = data.sample_day
        categories = ["sex", "smoking", "drinking", "hypertension", "diabetes", "prior_stroke"]
        if self.expanded:
            categories += ["toast"]
            x["nihss_5"] = data.nihss / 5
            x["log1p_lesion"] = np.log1p(data.lesion_ml)
        if self.month3:
            x["baseline_log2_hcy"] = np.log2(data.hcy) - self.centers["H0"]
            x["baseline_hcy_rcs"] = rcs_nonlinear(x.baseline_log2_hcy.to_numpy(), self.knots["H0"])
            x["baseline_H_x_W"] = x.baseline_log2_hcy*x.W
            x["sample3_day_30"] = data.sample3_day / 30
        if self.functional:
            x["pre_mrs"] = data.pre_mrs
            if not self.expanded:
                x["nihss_5"] = data.nihss / 5
        if self.t1:
            x["gm119_100ml"] = data.gm119_ml / 100
        for name in categories:
            for level in CATEGORIES[name][1:]:
                x[f"{name}_{level}"] = data[name].eq(level).astype(float)
            if data[name].isna().any():
                raise DataError(f"Missing {name} before design construction")
        return x

    def transform(self, data: pd.DataFrame) -> pd.DataFrame:
        x = self._matrix(data)[self.columns]
        if not np.isfinite(x.to_numpy()).all():
            raise DataError("Non-finite model matrix; required covariates were not properly imputed")
        return x

    def to_dict(self) -> dict:
        return dict(self.__dict__)
