"""Isolated longitudinal functional-state trajectory analysis.

This module is intentionally independent from the historical recovery model and
the recovery-path recurrence extension.  It reads only the recovery baseline,
imaging, repeated mRS/death fields and writes to its own output root.
"""
from __future__ import annotations

import copy
import traceback
from dataclasses import dataclass
from datetime import UTC, datetime
from itertools import pairwise
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import softmax

from ..common import DataError, dump_json, outdir, read_csv, resolve
from ..imaging import read_imaging
from ..sas_extract import extract, inventory, resolve_owners
from .data import build_cohort, cohort_audit, derive, read_clinical
from .design import StudyDesign
from .imputation import impute
from .pooling import joint_test
from .registry import ModelSpec, primary_spec
from .runner import code_hashes, settings, study_sources

MONTHS = (3, 12, 24, 36, 48, 60)
INTERVALS = tuple(pairwise(MONTHS))
REPORT_MONTHS = (36, 48, 60)
ROOT_NAME = "01_recovery_trajectory"
CONTRACT = "recovery_trajectory_discrete_multistate_20261012_v1"


def root(cfg: dict) -> Path:
    return outdir(cfg) / "studies" / ROOT_NAME


def trajectory_spec() -> ModelSpec:
    """Freeze the existing recovery adjustment set and exposure scales."""
    return primary_spec("recovery").variant(
        "functional_trajectory",
        family="multinomial",
        outcome="state60",
        tier="secondary",
    )


def source_fields() -> dict:
    """Recovery fields only; no recurrence, Hcy, BP, CEC or kidney fields."""
    return dict(study_sources("recovery"))


def _number(data: pd.DataFrame, name: str) -> pd.Series:
    return pd.to_numeric(data.get(name, pd.Series(np.nan, index=data.index)), errors="coerce")


def transition_panel(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build adjacent observed state pairs without forward filling outcomes."""
    rows: list[dict] = []
    audits: list[dict] = []
    for interval_index, (start, end) in enumerate(INTERVALS):
        start_state = _number(data, f"state{start}")
        end_state = _number(data, f"state{end}")
        start_conflict = data.get(f"state{start}_conflict", pd.Series(False, index=data.index)).fillna(False)
        end_conflict = data.get(f"state{end}_conflict", pd.Series(False, index=data.index)).fillna(False)
        known = start_state.isin((0, 1, 2)) & end_state.isin((0, 1, 2)) & ~start_conflict & ~end_conflict
        prior_death = pd.Series(False, index=data.index)
        for prior in MONTHS[:interval_index]:
            prior_death |= _number(data, f"state{prior}").eq(2)
        death_to_alive = known & (start_state.eq(2) | prior_death) & end_state.isin((0, 1))
        at_risk = known & start_state.isin((0, 1)) & ~prior_death
        selected = data.index[at_risk]
        for row_id in selected:
            rows.append({
                "row_id": int(row_id),
                "patient_id": str(data.loc[row_id, "patient_id"]),
                "interval_index": interval_index,
                "from_month": start,
                "to_month": end,
                "origin": int(start_state.loc[row_id]),
                "destination": int(end_state.loc[row_id]),
            })
        audits.append({
            "interval_index": interval_index,
            "from_month": start,
            "to_month": end,
            "eligible_start_n": int(start_state.isin((0, 1)).sum()),
            "known_adjacent_state_n": int(known.sum()),
            "transition_rows_n": int(at_risk.sum()),
            "unknown_start_n": int((~start_state.isin((0, 1, 2))).sum()),
            "unknown_end_n": int((~end_state.isin((0, 1, 2))).sum()),
            "conflict_n": int((start_conflict | end_conflict).sum()),
            "death_start_n": int(start_state.eq(2).sum()),
            "death_to_alive_n": int(death_to_alive.sum()),
            "death_absorbing_check": bool(not death_to_alive.any()),
        })
        transition_names = {
            (0, 0): "independent_to_independent_n",
            (0, 1): "independent_to_dependent_n",
            (0, 2): "independent_to_dead_n",
            (1, 0): "dependent_to_independent_n",
            (1, 1): "dependent_to_dependent_n",
            (1, 2): "dependent_to_dead_n",
        }
        for (origin, destination), name in transition_names.items():
            audits[-1][name] = int(
                (at_risk & start_state.eq(origin) & end_state.eq(destination)).sum()
            )
    panel = pd.DataFrame(rows, columns=[
        "row_id", "patient_id", "interval_index", "from_month", "to_month", "origin", "destination"
    ])
    return panel, pd.DataFrame(audits)


@dataclass
class TransitionFit:
    terms: list[str]
    params: np.ndarray
    covariance: np.ndarray
    diagnostics: dict


def _transition_y(panel: pd.DataFrame) -> np.ndarray:
    """Encode stay=0, switch to other living state=1, death=2."""
    origin = panel.origin.to_numpy(int)
    destination = panel.destination.to_numpy(int)
    y = np.where(destination == 2, 2, np.where(destination == origin, 0, 1))
    return y.astype(int)


def _design_matrix(panel: pd.DataFrame, completed: pd.DataFrame, design: StudyDesign) -> tuple[np.ndarray, list[str]]:
    row = completed.iloc[panel.row_id.to_numpy(int)].reset_index(drop=True)
    base = design.transform(row)
    columns = {name: base[name].to_numpy(float) for name in base}
    interval = panel.interval_index.to_numpy(int)
    origin_one = panel.origin.to_numpy(dtype=int) == 1
    for k in range(1, len(INTERVALS)):
        name = f"interval_{k}"
        value = (interval == k).astype(float)
        columns[name] = value
        columns[f"origin_1_x_{name}"] = value * origin_one
    for exposure in ("wmh_ml", "gm119_ml"):
        if exposure not in columns:
            raise DataError(f"Transition design lacks {exposure}")
        columns[f"origin_1_x_{exposure}"] = columns[exposure] * origin_one
    frame = pd.DataFrame(columns, index=panel.index, dtype=float)
    if not np.isfinite(frame.to_numpy()).all():
        raise DataError("Nonfinite transition design matrix")
    return frame.to_numpy(), list(frame.columns)


def fit_transition(completed: pd.DataFrame, panel: pd.DataFrame, design: StudyDesign) -> TransitionFit:
    """Fit the stacked discrete-time transition model with patient clusters."""
    if panel.empty:
        raise DataError("No observed adjacent state transitions")
    x, base_terms = _design_matrix(panel, completed, design)
    y = _transition_y(panel)
    counts = np.bincount(y, minlength=3)
    if np.any(counts == 0):
        raise DataError(f"Transition category has zero rows: stay/switch/death={counts.tolist()}")
    if np.linalg.matrix_rank(x) < x.shape[1]:
        raise DataError("Rank deficient transition design; no automatic covariate deletion")
    variable = x.std(axis=0) > 0
    scaled = x[:, variable]
    scaled = (scaled - scaled.mean(axis=0)) / np.where(scaled.std(axis=0) > 0, scaled.std(axis=0), 1)
    condition = float(np.linalg.cond(scaled))
    groups = panel.patient_id.to_numpy()
    with np.errstate(all="ignore"):
        model = sm.MNLogit(y, x, missing="raise")
        result = model.fit(
            method="newton",
            maxiter=150,
            disp=False,
            cov_type="cluster",
            cov_kwds={"groups": groups, "use_correction": True},
        )
    if not result.mle_retvals.get("converged", False):
        raise DataError("Stacked transition multinomial model failed to converge")
    params_matrix = np.asarray(result.params, dtype=float)
    params = params_matrix.reshape(-1, order="F")
    covariance = np.asarray(result.cov_params(), dtype=float)
    if not np.isfinite(params).all() or not np.isfinite(covariance).all():
        raise DataError("Nonfinite transition estimates")
    if np.linalg.eigvalsh(covariance).min() <= 0:
        raise DataError("Non-positive-definite transition covariance")
    score = np.asarray(result.model.score(params), dtype=float)
    max_score = float(np.max(np.abs(score)) / len(y))
    if max_score > 1e-5:
        raise DataError(f"Transition score check failed: {max_score:g}")
    terms = [f"other:{term}" for term in base_terms] + [f"dead:{term}" for term in base_terms]
    transition_labels = {
        (0, 0): "independent_to_independent",
        (0, 1): "independent_to_dependent",
        (0, 2): "independent_to_dead",
        (1, 0): "dependent_to_independent",
        (1, 1): "dependent_to_dependent",
        (1, 2): "dependent_to_dead",
    }
    transition_event_counts = {
        label: int(((panel.origin == origin) & (panel.destination == destination)).sum())
        for (origin, destination), label in transition_labels.items()
    }
    diagnostics = {
        "n": len(panel),
        "patient_n": int(panel.patient_id.nunique()),
        "parameters": len(params),
        "categories": {"stay": int(counts[0]), "switch": int(counts[1]), "death": int(counts[2])},
        "design_condition": condition,
        "max_score_per_patient": max_score,
        "converged": True,
        "covariance": "patient-cluster robust sandwich",
        "transition_event_counts": transition_event_counts,
        "zero_event_transitions": [label for label, count in transition_event_counts.items() if count == 0],
        "complete_separation_suspected": bool(np.abs(params).max() > 20),
        "warnings": [],
    }
    if diagnostics["zero_event_transitions"]:
        diagnostics["warnings"].append("zero-event transition observed")
    if diagnostics["complete_separation_suspected"]:
        diagnostics["warnings"].append("large coefficient; inspect for complete separation")
    return TransitionFit(terms, params, covariance, diagnostics)


def _probabilities(fit: TransitionFit, x: np.ndarray) -> np.ndarray:
    k = x.shape[1]
    beta = fit.params.reshape(k, 2, order="F")
    return softmax(np.column_stack([np.zeros(len(x)), x @ beta]), axis=1)


def _vectors(terms: list[str], exposure: str) -> dict[str, np.ndarray]:
    def vector(*names: str) -> np.ndarray:
        result = np.zeros(len(terms))
        for name in names:
            result[terms.index(name)] = 1
        return result

    other = f"other:{exposure}"
    other_delta = f"other:origin_1_x_{exposure}"
    death = f"dead:{exposure}"
    death_delta = f"dead:origin_1_x_{exposure}"
    return {
        "independent_to_dependent": vector(other),
        "dependent_to_independent": vector(other, other_delta),
        "independent_to_dead": vector(death),
        "dependent_to_dead": vector(death, death_delta),
    }


def _pool_contrast(fits: list[TransitionFit], vectors: list[np.ndarray]) -> dict:
    from ..imputation import pool_scalar

    q, u = [], []
    for fit, vector in zip(fits, vectors, strict=True):
        q.append(float(vector @ fit.params))
        u.append(float(vector @ fit.covariance @ vector))
    return pool_scalar(q, u)


def transition_effects(fits: list[TransitionFit]) -> tuple[pd.DataFrame, dict]:
    rows = []
    joint = {}
    for exposure in ("wmh_ml", "gm119_ml"):
        vectors = _vectors(fits[0].terms, exposure)
        labels = list(vectors)
        q_vectors = [vectors[label] for label in labels]
        estimates, covariances = [], []
        for label, vector in zip(labels, q_vectors, strict=True):
            pooled = _pool_contrast(fits, [vector] * len(fits))
            rows.append({
                "exposure": exposure,
                "transition": label,
                "estimate": pooled["estimate"],
                "lower": pooled["lower"],
                "upper": pooled["upper"],
                "ratio": float(np.exp(pooled["estimate"])),
                "ratio_lower": float(np.exp(pooled["lower"])),
                "ratio_upper": float(np.exp(pooled["upper"])),
                "p": pooled["p"],
                "m": pooled["m"],
                "scale": (
                    "relative transition probability ratio per 1 SD of log1p(WMH mL)"
                    if exposure == "wmh_ml" else
                    "relative transition probability ratio per 1 SD of raw GM119 mL"
                ),
            })
        for fit in fits:
            estimates.append(np.asarray([v @ fit.params for v in q_vectors]))
            covariances.append(np.asarray([[a @ fit.covariance @ b for b in q_vectors] for a in q_vectors]))
        joint[exposure] = joint_test(estimates, covariances)
    return pd.DataFrame(rows), joint


def _scenario_values(data: pd.DataFrame) -> dict[str, tuple[float, float]]:
    w = data.wmh_ml.quantile([.25, .5, .75]).to_numpy(float)
    g = data.gm119_ml.quantile([.25, .5, .75]).to_numpy(float)
    return {
        "wmh_q25": (float(w[0]), float(g[1])),
        "wmh_q50": (float(w[1]), float(g[1])),
        "wmh_q75": (float(w[2]), float(g[1])),
        "gm_q25": (float(w[1]), float(g[0])),
        "gm_q50": (float(w[1]), float(g[1])),
        "gm_q75": (float(w[1]), float(g[2])),
    }


def _scenario_design_matrices(
    data: pd.DataFrame,
    design: StudyDesign,
    wmh: float,
    gm: float,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Build each interval/origin design once for a standardized scenario."""
    working = data.copy()
    working["wmh_ml"] = wmh
    working["gm119_ml"] = gm
    n = len(working)
    result: list[tuple[np.ndarray, np.ndarray]] = []
    for interval_index, _ in enumerate(INTERVALS):
        interval_matrices = []
        for origin in (0, 1):
            panel = pd.DataFrame({
                "row_id": np.arange(n),
                "interval_index": interval_index,
                "origin": origin,
            })
            x, _ = _design_matrix(panel, working, design)
            interval_matrices.append(x)
        result.append((interval_matrices[0], interval_matrices[1]))
    return result


def _probabilities_from_params(params: np.ndarray, x: np.ndarray) -> np.ndarray:
    k = x.shape[1]
    beta = np.asarray(params).reshape(k, 2, order="F")
    return softmax(np.column_stack([np.zeros(len(x)), x @ beta]), axis=1)


def _scenario_occupancy(
    params: np.ndarray,
    matrices: list[tuple[np.ndarray, np.ndarray]],
) -> dict[int, np.ndarray]:
    """Propagate all three states through the five intervals once."""
    n = matrices[0][0].shape[0]
    occupied = np.zeros((n, 3), dtype=float)
    occupied[:, 0] = 1.0
    result: dict[int, np.ndarray] = {}
    for interval_index, (x_independent, x_dependent) in enumerate(matrices):
        independent = _probabilities_from_params(params, x_independent)
        dependent = _probabilities_from_params(params, x_dependent)[:, [1, 0, 2]]
        next_state = occupied[:, 0, None] * independent + occupied[:, 1, None] * dependent
        next_state[:, 2] += occupied[:, 2]
        occupied = next_state
        end = INTERVALS[interval_index][1]
        if end in REPORT_MONTHS:
            result[end] = occupied.mean(axis=0)
    return result


def _scenario_occupancy_draws(
    params: np.ndarray,
    matrices: list[tuple[np.ndarray, np.ndarray]],
) -> dict[int, np.ndarray]:
    """Propagate a batch of coefficient draws without Python-level draw loops."""
    if params.ndim != 2:
        raise ValueError("Parameter draws must be a two-dimensional array")
    n = matrices[0][0].shape[0]
    parameter_count = matrices[0][0].shape[1]
    if params.shape[1] != parameter_count * 2:
        raise ValueError("Parameter draw width does not match transition design")
    # Each row is flattened in Fortran order from the k-by-2 statsmodels table:
    # all switch coefficients followed by all death coefficients.
    beta = params.reshape(params.shape[0], 2, parameter_count).transpose(0, 2, 1)
    occupied = np.zeros((len(params), n, 3), dtype=float)
    occupied[:, :, 0] = 1.0
    result: dict[int, np.ndarray] = {}
    for interval_index, (x_independent, x_dependent) in enumerate(matrices):
        logits_independent = np.einsum("nk,dkc->dnc", x_independent, beta)
        logits_dependent = np.einsum("nk,dkc->dnc", x_dependent, beta)
        zeros = np.zeros((len(params), n, 1), dtype=float)
        independent = softmax(np.concatenate([zeros, logits_independent], axis=2), axis=2)
        dependent = softmax(np.concatenate([zeros, logits_dependent], axis=2), axis=2)[:, :, [1, 0, 2]]
        next_state = occupied[:, :, 0, None] * independent + occupied[:, :, 1, None] * dependent
        next_state[:, :, 2] += occupied[:, :, 2]
        occupied = next_state
        end = INTERVALS[interval_index][1]
        if end in REPORT_MONTHS:
            result[end] = occupied.mean(axis=1)
    return result


def state_probabilities(
    fits: list[TransitionFit],
    data_sets: list[pd.DataFrame],
    design: StudyDesign,
    draws: int = 16,
    seed: int = 20261012,
) -> pd.DataFrame:
    from ..imputation import pool_scalar

    scenarios = _scenario_values(data_sets[0])
    rows = []
    states = ("independent", "dependent", "dead")
    for scenario, (wmh, gm) in scenarios.items():
        point_by_fit: list[dict[int, np.ndarray]] = []
        within_by_fit: list[dict[int, np.ndarray]] = []
        for index, fit in enumerate(fits):
            matrices = _scenario_design_matrices(data_sets[index], design, wmh, gm)
            point = _scenario_occupancy(fit.params, matrices)
            point_by_fit.append(point)
            if draws <= 1:
                within_by_fit.append({month: np.zeros(3) for month in REPORT_MONTHS})
                continue
            rng = np.random.default_rng(seed + index * 1000)
            sampled = rng.multivariate_normal(fit.params, fit.covariance, size=draws)
            simulated_by_month = _scenario_occupancy_draws(sampled, matrices)
            within_by_fit.append({
                month: np.var(simulated_by_month[month], axis=0, ddof=1)
                for month in REPORT_MONTHS
            })
        for month in REPORT_MONTHS:
            for state_index, state in enumerate(states):
                q = [float(point[month][state_index]) for point in point_by_fit]
                u = [float(within[month][state_index]) for within in within_by_fit]
                pooled = pool_scalar(q, u)
                rows.append({
                    "scenario": scenario,
                    "month": month,
                    "state": state,
                    "wmh_ml": wmh,
                    "gm119_ml": gm,
                    "probability": pooled["estimate"],
                    "lower": pooled["lower"],
                    "upper": pooled["upper"],
                    "ci_method": "Rubin variance; parametric delta simulation within imputation",
                })
    return pd.DataFrame(rows)


def mrs_distribution(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for month in MONTHS:
        raw = _number(data, f"mrs{month}")
        state = _number(data, f"state{month}")
        score = raw.mask(raw.isna() & state.eq(2), 6)
        denominator = int(score.notna().sum())
        for mrs in range(7):
            n = int(score.eq(mrs).sum())
            rows.append({"month": month, "mrs": mrs, "n": n,
                         "denominator_known": denominator,
                         "proportion": n / denominator if denominator else np.nan})
        rows.append({"month": month, "mrs": "unknown", "n": int(score.isna().sum()),
                     "denominator_known": denominator, "proportion": np.nan})
    return pd.DataFrame(rows)


def _plot_outputs(folder: Path, probabilities: pd.DataFrame, effects: pd.DataFrame, mrs: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for state, ax in zip(("independent", "dependent", "dead"), axes, strict=True):
        for scenario in ("wmh_q25", "wmh_q75", "gm_q25", "gm_q75"):
            frame = probabilities[(probabilities.scenario == scenario) & (probabilities.state == state)]
            if frame.empty:
                continue
            ax.plot(frame.month, frame.probability, marker="o", label=scenario)
            ax.fill_between(frame.month, frame.lower, frame.upper, alpha=.12)
        ax.set_title(state)
        ax.set_xlabel("Months after stroke")
        ax.set_xticks(REPORT_MONTHS)
    axes[0].set_ylabel("Standardized state probability")
    axes[-1].legend(fontsize=8, loc="best")
    fig.suptitle("Longitudinal functional-state trajectories")
    fig.tight_layout()
    fig.savefig(folder / "functional_state_trajectory.png", dpi=180)
    plt.close(fig)

    if not effects.empty:
        fig, ax = plt.subplots(figsize=(9, 4.5))
        y = np.arange(len(effects))
        ratio = effects.ratio.to_numpy(float)
        lower = effects.ratio_lower.to_numpy(float)
        upper = effects.ratio_upper.to_numpy(float)
        ax.errorbar(ratio, y, xerr=[ratio-lower, upper-ratio], fmt="o", color="#1769a4")
        ax.axvline(1, color="gray", linestyle="--")
        ax.set_yticks(y)
        ax.set_yticklabels(effects.exposure + " | " + effects.transition)
        ax.set_xlabel("Relative transition probability ratio per 1 SD")
        fig.tight_layout()
        fig.savefig(folder / "transition_effects.png", dpi=180)
        plt.close(fig)

    known = mrs[mrs.mrs != "unknown"].copy()
    if not known.empty:
        pivot = known.pivot(index="month", columns="mrs", values="proportion").fillna(0)
        ax = pivot.plot(kind="bar", stacked=True, figsize=(9, 4), colormap="viridis")
        ax.set_xlabel("Months after stroke")
        ax.set_ylabel("Proportion among known mRS/death states")
        ax.set_title("Descriptive mRS distribution")
        ax.legend(title="mRS", ncol=4)
        fig = ax.get_figure()
        fig.tight_layout()
        fig.savefig(folder / "mrs_trajectory.png", dpi=180)
        plt.close(fig)


def prepare(cfg: dict) -> dict:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    output_root = root(cfg)
    folder = output_root / "runs" / stamp
    folder.mkdir(parents=True, exist_ok=False)
    local = copy.deepcopy(cfg)
    local["_out"] = str(folder)
    state = {"study": "recovery_trajectory", "run": stamp, "mode": cfg["mode"],
             "contract": CONTRACT, "path": str(folder), "status": "PREPARING",
             "analysis_scope": "trajectory_only", "output_root": str(output_root),
             "source_code_sha256": code_hashes()}
    dump_json(folder / "config_snapshot.json", local)
    dump_json(folder / "status.json", state)
    try:
        fields = source_fields()
        inv, owners = [], {}
        if cfg["inputs"].get("clinical_csv"):
            source = resolve(cfg, cfg["inputs"]["clinical_csv"])
            names = {str(c).casefold() for c in read_csv(source).columns}
            owners = {s: str(source) for s in fields if s.casefold() in names}
            dump_json(folder / "source_files.json", [{"path": str(source)}])
        else:
            issues = []
            inv = inventory(local, fields=fields)
            owners, issues = resolve_owners(inv, local, fields=fields)
            dump_json(folder / "source_inventory.json", inv)
            dump_json(folder / "source_issues.json", issues)
            if issues:
                raise DataError("Conflicting source owners; see source_issues.json")
            extract(local, fields=fields)
        invalid_masks: dict[str, pd.Series] = {}
        clinical, field_audit = read_clinical(local, invalid_masks=invalid_masks)
        field_audit = field_audit.loc[field_audit.source.isin(fields)].copy()
        field_audit["file"] = field_audit.source.map(owners).fillna("")
        field_audit.to_csv(folder / "field_audit.csv", index=False)
        images = read_imaging(local, include_gm=True)
        master = derive(clinical.merge(images, on="patient_id", how="left", validate="one_to_one"))
        base, exclusions, flow = build_cohort(master, "recovery")
        master.to_csv(folder / "master.csv", index=False)
        base.to_csv(folder / "base_eligible.csv", index=False)
        exclusions.to_csv(folder / "base_exclusions.csv", index=False)
        flow.to_csv(folder / "base_flow.csv", index=False)
        panel, interval_audit = transition_panel(base)
        panel.to_csv(folder / "transition_panel.csv", index=False)
        interval_audit.to_csv(folder / "transition_counts.csv", index=False)
        cohort = cohort_audit(base, "recovery")
        state_counts = {
            str(month): {str(k): int(base[f"state{month}"].eq(k).sum()) for k in (0, 1, 2)} |
            {"unknown": int(base[f"state{month}"].isna().sum())}
            for month in MONTHS
        }
        audit = {
            "clinical_n": len(clinical),
            "mapped_image_n": len(images),
            "base_eligible_n": len(base),
            "transition_rows_n": len(panel),
            "state_counts": state_counts,
            "intervals": interval_audit.to_dict("records"),
            "cohort": cohort,
            "source_absent": field_audit.loc[~field_audit.present, "source"].tolist(),
            "source_invalid": field_audit.loc[field_audit.invalid.gt(0), ["source", "invalid"]].to_dict("records"),
            "future_recurrence_fields_requested": sorted(
                s for s in fields if "STROKE" in s.upper() and s.upper() != "H_STROKE"
            ),
            "analysis_scope": "trajectory_only",
        }
        audit["status"] = "PREPARED" if len(base) and len(panel) and not interval_audit.death_to_alive_n.any() else "REVIEW_REQUIRED"
        dump_json(folder / "audit.json", audit)
        state.update(status=audit["status"], audit=audit)
    except (DataError, OSError, ValueError, KeyError) as exc:
        state.update(status="INPUTS_REQUIRED", error=str(exc))
        (folder / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
    dump_json(folder / "status.json", state)
    return state


def analyse(state: dict) -> dict:
    folder = Path(state["path"])
    import json

    cfg = json.loads((folder / "config_snapshot.json").read_text(encoding="utf-8"))
    options = settings(cfg)
    base = pd.read_csv(folder / "base_eligible.csv", dtype={"patient_id": str})
    panel = pd.read_csv(folder / "transition_panel.csv", dtype={"patient_id": str})
    design = StudyDesign.freeze(base, trajectory_spec())
    datasets, mi = impute(base, design, options)
    dump_json(folder / "imputation.json", mi)
    fits: list[TransitionFit] = []
    diagnostics = []
    for index, completed in enumerate(datasets, 1):
        fit = fit_transition(completed, panel, design)
        fits.append(fit)
        diagnostics.append(fit.diagnostics)
        if index % 10 == 0:
            print(f"  trajectory: fitted {index}/{len(datasets)}", flush=True)
    dump_json(folder / "fit_diagnostics.json", diagnostics)
    effects, joint = transition_effects(fits)
    effects.to_csv(folder / "transition_effects.csv", index=False)
    dump_json(folder / "joint_tests.json", joint)
    probability_draws = int(
        cfg.get("analysis", {}).get(
            "trajectory_probability_draws",
            cfg.get("studies", {}).get("trajectory_probability_draws", 16),
        )
    )
    probabilities = state_probabilities(fits, datasets, design, draws=probability_draws)
    probabilities.to_csv(folder / "state_probabilities_36_48_60.csv", index=False)
    mrs = mrs_distribution(base)
    mrs.to_csv(folder / "mrs_trajectory_distribution.csv", index=False)
    _plot_outputs(folder, probabilities, effects, mrs)
    state["analysis"] = {
        "status": "ESTIMATED",
        "imputations": len(fits),
        "transition_rows": len(panel),
        "trajectory_probability_draws": probability_draws,
        "joint_tests": joint,
        "outputs": [
            "transition_effects.csv", "state_probabilities_36_48_60.csv",
            "functional_state_trajectory.png", "transition_effects.png",
            "mrs_trajectory_distribution.csv", "mrs_trajectory.png",
        ],
    }
    state["status"] = "COMPLETED"
    state["finished_utc"] = datetime.now(UTC).isoformat()
    dump_json(folder / "status.json", state)
    dump_json(Path(state["output_root"]) / "latest_results.json",
              {"path": str(folder), "run": state["run"], "status": state["status"]})
    return state


def run(cfg: dict, through: str = "analyse") -> dict:
    if through not in {"prepare", "analyse"}:
        raise DataError("recovery-trajectory supports prepare or analyse")
    state = prepare(cfg)
    if through == "prepare" or state.get("status") != "PREPARED":
        return state
    state["status"] = "ANALYSING"
    dump_json(Path(state["path"]) / "status.json", state)
    try:
        return analyse(state)
    except (DataError, OSError, ValueError, KeyError, np.linalg.LinAlgError) as exc:
        state.update(status="PARTIAL", error=str(exc))
        (Path(state["path"]) / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
        dump_json(Path(state["path"]) / "status.json", state)
        return state
