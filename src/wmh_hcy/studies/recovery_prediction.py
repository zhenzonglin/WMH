"""Internal validation of three-month five-year state predictions; no action threshold."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.special import logit, softmax

from ..common import DataError, dump_json
from .design import StudyDesign
from .imputation import impute
from .models import fit
from .registry import ModelSpec, primary_spec


def probabilities(fitted, x: np.ndarray) -> np.ndarray:
    k = x.shape[1]
    if len(fitted.params) != 2*k:
        raise DataError("Multinomial prediction dimension mismatch")
    return softmax(np.column_stack([np.zeros(len(x)), x @ fitted.params.reshape(k, 2, order="F")]), axis=1)


def metrics(y: np.ndarray, probability: np.ndarray) -> dict:
    positive = y == 1
    if positive.min() == positive.max():
        raise DataError("Validation fold contains only one outcome class")
    ranks = stats.rankdata(probability)
    n1, n0 = positive.sum(), (~positive).sum()
    auc = (ranks[positive].sum()-n1*(n1+1)/2)/(n1*n0)
    linear = logit(np.clip(probability, 1e-6, 1-1e-6))
    calibration = sm.GLM(positive.astype(int), np.column_stack([np.ones(len(y)), linear]),
                         family=sm.families.Binomial()).fit()
    if not calibration.converged or not np.isfinite(calibration.params).all():
        raise DataError("Calibration model failed")
    return {"n": len(y), "events": int(n1), "auc": float(auc),
            "brier": float(np.mean((positive-probability)**2)),
            "calibration_intercept": float(calibration.params[0]),
            "calibration_slope": float(calibration.params[1]),
            "observed_risk": float(positive.mean()), "mean_predicted_risk": float(probability.mean())}


def fit_validation(datasets: list[pd.DataFrame], spec: ModelSpec, settings: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    design = StudyDesign.freeze(datasets[0], spec)
    predicted = []
    for completed in datasets:
        model = fit(completed, design)
        predicted.append(probabilities(model, design.transform(completed).to_numpy()))
    all_prob = np.mean(predicted, axis=0)
    # Bootstrap optimism is conditional on one imputation; the primary fit uses all M.
    first = datasets[0]
    x = design.transform(first).to_numpy()
    y = first.state60.to_numpy(int)
    rng = np.random.default_rng(int(settings.get("seed", 20260917))+907)
    nboot = int(settings.get("prediction_bootstrap", 100))
    if nboot < 20:
        raise DataError("Prediction validation requires at least 20 bootstrap replicates")
    optimism = {"dependent": [], "dependent_or_dead": []}
    failures = []
    for index in range(nboot):
        sampled = rng.integers(0, len(first), len(first))
        try:
            fitted = fit(first.iloc[sampled].reset_index(drop=True), design)
            train_p = probabilities(fitted, x[sampled])
            original_p = probabilities(fitted, x)
            for label, target in (("dependent", 1), ("dependent_or_dead", -1)):
                train_y = (y[sampled] == target) if target != -1 else (y[sampled] != 0)
                original_y = (y == target) if target != -1 else (y != 0)
                train_r = train_p[:, 1] if target != -1 else 1-train_p[:, 0]
                original_r = original_p[:, 1] if target != -1 else 1-original_p[:, 0]
                train_m = metrics(train_y.astype(int), train_r)
                test_m = metrics(original_y.astype(int), original_r)
                optimism[label].append({m: train_m[m]-test_m[m] for m in
                                        ("auc", "brier", "calibration_intercept", "calibration_slope")})
        except (DataError, ValueError, np.linalg.LinAlgError) as exc:
            failures.append({"replicate": index, "reason": str(exc)})
    if any(len(v) < max(20, nboot//2) for v in optimism.values()):
        raise DataError(f"Bootstrap validation failed too often: {len(failures)}/{nboot}")
    result = {"m": len(datasets), "bootstrap_attempts": nboot,
              "bootstrap_success": len(optimism["dependent"]), "bootstrap_failures": failures,
              "bootstrap_scope": "first completed dataset only; imputation uncertainty not included in optimism"}
    for label, target in (("dependent", 1), ("dependent_or_dead", -1)):
        binary = (y == target) if target != -1 else (y != 0)
        risk = all_prob[:, 1] if target != -1 else 1-all_prob[:, 0]
        apparent = metrics(binary.astype(int), risk)
        names = ("auc", "brier", "calibration_intercept", "calibration_slope")
        mean_optimism = {m: float(np.mean([v[m] for v in optimism[label]])) for m in names}
        result[label] = {"apparent": apparent, "mean_bootstrap_optimism": mean_optimism,
                         "optimism_corrected_auc": apparent["auc"]-mean_optimism["auc"],
                         "optimism_corrected_brier": apparent["brier"]-mean_optimism["brier"],
                         "optimism_corrected_calibration_intercept": apparent["calibration_intercept"]-mean_optimism["calibration_intercept"],
                         "optimism_corrected_calibration_slope": apparent["calibration_slope"]-mean_optimism["calibration_slope"],
                         "optimism_quantiles": {m: np.quantile([v[m] for v in optimism[label]], [.025, .5, .975]).tolist()
                                                for m in names}}
    return all_prob, y, result


def validate(data: pd.DataFrame, directory, settings: dict) -> dict:
    """Death remains a separate multinomial state for the dependence target."""
    directory.mkdir(parents=True, exist_ok=True)
    known = data.loc[data.state60.notna()].reset_index(drop=True)
    if len(known) < 100 or known.state60.nunique() != 3:
        raise DataError("Three five-year states with adequate observed sample required")
    full = primary_spec("recovery").variant("risk_full", primary=("dependent:wmh_ml",))
    simple = ModelSpec("recovery", name="risk_wmh_only", exposures=("wmh_ml",),
                       primary=("dependent:wmh_ml",), splines=())
    output = {}
    for name, spec in (("full", full), ("wmh_only", simple)):
        design = StudyDesign.freeze(known, spec)
        completed, mi = impute(known, design, settings)
        dump_json(directory / f"{name}_imputation.json", mi)
        probability, y, summary = fit_validation(completed, spec, settings)
        output[name] = summary
        prediction = pd.DataFrame({"patient_id": known.patient_id, "state60": y,
                                   "p_independent": probability[:, 0],
                                   "p_dependent": probability[:, 1], "p_dead": probability[:, 2]})
        prediction.to_csv(directory / f"{name}_predictions.csv", index=False)
        if name == "full":
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            for ax, label, risk, observed in ((axes[0], "Dependent; death separate", probability[:, 1], y == 1),
                                               (axes[1], "Dependent or dead", 1-probability[:, 0], y != 0)):
                table = pd.DataFrame({"risk": risk, "observed": observed})
                table["bin"] = pd.qcut(table.risk.rank(method="first"), q=10, labels=False)
                points = table.groupby("bin").agg(predicted=("risk", "mean"), actual=("observed", "mean"))
                ax.plot(points.predicted, points.actual, "o-")
                ax.plot([0, 1], [0, 1], "--", color="gray")
                ax.set(xlabel="Mean predicted probability", ylabel="Observed proportion", title=label,
                       xlim=(0, 1), ylim=(0, 1))
            fig.suptitle("Apparent calibration; internal bootstrap metrics in JSON")
            fig.tight_layout()
            fig.savefig(directory / "apparent_calibration.png", dpi=180)
            plt.close(fig)
    output["limits"] = ("Exploratory development/internal validation only; no external validation, "
                        "action threshold or WMH cutoff; patients with unknown state are not assigned a class")
    dump_json(directory / "validation.json", output)
    return output
