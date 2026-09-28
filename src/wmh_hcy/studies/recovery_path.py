"""Exploratory three-month recovery landmark and subsequent all-stroke analysis.

This extension has its own runs and never changes the saved four-study contract.
Patient-level data remain inside the configured, ignored output directory.
"""
from __future__ import annotations

import copy
import json
import traceback
from datetime import UTC, datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..common import DataError, dump_json, outdir, read_csv, resolve
from ..imaging import read_imaging
from ..sas_extract import extract, inventory, resolve_owners
from .analyses import run_one
from .data import build_cohort, cohort_audit, read_clinical
from .design import StudyDesign
from .models import Fit
from .pooling import contrast
from .registry import SOURCES, ModelSpec, primary_spec
from .runner import code_hashes, settings, study_sources

HORIZON = 1825
MONTHS = (3, 12, 24, 36, 48, 60)
STROKE_YEARS = (2, 3, 4, 5)
CONTRACT = "recovery_path_all_stroke_20260928_v1"


def source_fields() -> dict:
    """The old recovery extraction stays unchanged; this is an explicit extension."""
    fields = dict(study_sources("recovery"))
    names = {"ONSET_D", "F3_DATE", "F12_DATE", "D_DEATH_D", "F3_DEATH_D",
             "F6_DEATH_D", "F12_DEATH_D"}
    for year in STROKE_YEARS:
        source = "Y5_STROKE" if year == 5 else f"y{year}_STROKE"
        names.update((source, source + ("_DD" if year == 5 else "_dd")))
    return fields | {s: SOURCES[s] for s in names}


def root(cfg: dict) -> Path:
    return outdir(cfg) / "studies" / "01_recovery_path"


def path_settings(cfg: dict) -> dict:
    values = settings(cfg)
    path = cfg.get("recovery_path", {})
    for key in ("prediction_bootstrap", "mediation_bootstrap"):
        values[key] = int(path.get(key, 100))
    return values


def _number(data: pd.DataFrame, name: str) -> pd.Series:
    return pd.to_numeric(data.get(name, pd.Series(np.nan, index=data.index)), errors="coerce")


def recurrence_spec() -> ModelSpec:
    return primary_spec("recovery").variant(
        "stroke_recurrence", family="cox", outcome="event_type",
        splines=("age", "wmh_ml"), primary=("wmh_ml", "wmh_ml_rcs"), tier="exploratory")


def first_stroke_conflicts(data: pd.DataFrame) -> pd.Series:
    """Compare only observed yearly records; absence is reported, never fabricated."""
    event = _number(data, "y5_stroke_event")
    day = _number(data, "y5_stroke_day")
    bad = pd.Series(False, index=data.index)
    for year in (2, 3, 4):
        earlier = _number(data, f"y{year}_stroke_event")
        earlier_day = _number(data, f"y{year}_stroke_day")
        observed = earlier.isin((0, 1)) & np.isfinite(earlier_day) & earlier_day.ge(0)
        bad |= observed & earlier.eq(1) & (~event.eq(1) | day.ne(earlier_day))
        bad |= observed & earlier.eq(0) & event.eq(1) & day.le(earlier_day.clip(upper=365*year))
        bad |= observed & earlier.eq(1) & earlier_day.gt(365*year)
    return bad


def observation_after_known_death(data: pd.DataFrame) -> pd.Series:
    """Reject clear follow-up beyond a visit confirming death, without dating death."""
    first = pd.Series(np.nan, index=data.index)
    for month in MONTHS[1:]:
        first.loc[first.isna() & _number(data, f"state{month}").eq(2)] = month
    stop = _number(data, "y5_stroke_day")
    return first.notna() & stop.gt(first*365/12)


def landmark(base: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Make the actual-visit risk set without extending loss to follow-up."""
    d = base.copy()
    d["entry"] = _number(d, "visit3_day")
    event, day = _number(d, "y5_stroke_event"), _number(d, "y5_stroke_day")
    d["exit"] = day.clip(upper=HORIZON)
    d["event_type"] = (event.eq(1) & day.le(HORIZON)).astype(int)
    checks = {
        "valid_actual_visit": np.isfinite(d.entry) & d.entry.gt(0) & d.entry.lt(HORIZON),
        "known_first_stroke_status_and_time": (event.isin((0, 1)) & np.isfinite(day) & day.ge(0)
                                               & day.mod(1).eq(0) & ~(event.eq(1) & day.gt(HORIZON))),
        "yearly_first_event_consistent": ~first_stroke_conflicts(d),
        "no_first_recurrence_before_or_on_visit": ~(event.eq(1) & day.le(d.entry)),
        "no_observation_after_dated_death": _number(d, "death_day").isna() | day.le(_number(d, "death_day")),
        "no_observation_beyond_confirmed_death_visit": ~observation_after_known_death(d),
        "observation_after_visit": d.exit.gt(d.entry),
    }
    keep = pd.Series(True, index=d.index)
    flow = []
    for label, check in checks.items():
        rejected = keep & ~check.fillna(False)
        keep &= check.fillna(False)
        flow.append({"step": label, "remaining": int(keep.sum()), "excluded_here": int(rejected.sum())})
    return d.loc[keep].reset_index(drop=True), pd.DataFrame(flow)


def chronology_audit(data: pd.DataFrame) -> dict:
    """An event inside an undated death interval has unknown temporal order."""
    death_month = pd.Series(np.nan, index=data.index)
    death_lower = pd.Series(np.nan, index=data.index)
    previous = pd.Series(0., index=data.index)
    ambiguous = pd.Series(False, index=data.index)
    contradictions = pd.Series(False, index=data.index)
    for month in MONTHS[1:]:
        state = _number(data, f"state{month}")
        death_flag = _number(data, f"death{month}")
        newly_dead = death_month.isna() & (state.eq(2) | death_flag.eq(1))
        event_day = _number(data, "y5_stroke_day")
        event = _number(data, "y5_stroke_event").eq(1)
        lower = previous.where(newly_dead)
        upper = float(month*365/12)
        contradictions |= newly_dead & event & event_day.gt(upper)
        if month == 12:
            exact = _number(data, "death_day")
            # Exact dates, when present, resolve this interval only.
            ambiguous |= newly_dead & exact.isna() & event & event_day.gt(lower) & event_day.le(upper)
            contradictions |= newly_dead & exact.notna() & event & event_day.gt(exact)
        else:
            ambiguous |= newly_dead & event & event_day.gt(lower) & event_day.le(upper)
        death_month.loc[newly_dead] = month
        death_lower.loc[newly_dead] = lower.loc[newly_dead]
        previous.loc[state.isin((0, 1)) | death_flag.eq(2)] = upper
    dated_death = _number(data, "death_day")
    contradictions |= dated_death.notna() & _number(data, "y5_stroke_day").gt(dated_death)
    early_stop = data.event_type.eq(0) & data.exit.lt(HORIZON)
    known_death_stop = (early_stop & death_month.notna() & data.exit.gt(death_lower)
                        & data.exit.le(death_month*365/12))
    return {"death_interval_ambiguous_event_n": int(ambiguous.sum()),
            "event_after_known_death_n": int(contradictions.sum()),
            "death_by_five_year_n": int(death_month.notna().sum()),
            "five_year_state_unknown_n": int(data.state60.isna().sum()),
            "early_censor_n": int(early_stop.sum()),
            "known_death_stop_n": int(known_death_stop.sum()),
            "unresolved_loss_n": int((early_stop & ~known_death_stop).sum())}


def prepare(cfg: dict) -> dict:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    folder = root(cfg) / "runs" / stamp
    folder.mkdir(parents=True, exist_ok=False)
    local = copy.deepcopy(cfg)
    local["_out"] = str(folder)
    state = {"study": "recovery_path", "run": stamp, "mode": cfg["mode"], "contract": CONTRACT,
             "path": str(folder), "status": "PREPARING", "statistics": path_settings(cfg),
             "source_code_sha256": code_hashes()}
    dump_json(folder / "config_snapshot.json", local)
    dump_json(folder / "status.json", state)
    try:
        fields = source_fields()
        if cfg["inputs"].get("clinical_csv"):
            source = resolve(cfg, cfg["inputs"]["clinical_csv"])
            raw_names = {n.casefold() for n in read_csv(source).columns}
            owners = {s: str(source) for s in fields if s.casefold() in raw_names}
        else:
            owners, issues = resolve_owners(inventory(local, fields=fields), local, fields=fields)
            dump_json(folder / "source_issues.json", issues)
            if issues:
                raise DataError("Conflicting source owners; see source_issues.json")
            extract(local, fields=fields)
        invalid_masks = {}
        clinical, field_audit = read_clinical(local, invalid_masks=invalid_masks)
        field_audit = field_audit.loc[field_audit.source.isin(fields)].copy()
        field_audit["file"] = field_audit.source.map(owners).fillna("")
        field_audit.to_csv(folder / "field_audit.csv", index=False)
        required = {"ONSET_D", "F3_DATE", "Y5_STROKE", "Y5_STROKE_DD"}
        missing = sorted(required - set(field_audit.loc[field_audit.present, "source"]))
        if missing:
            raise DataError("Missing recurrence-path required fields: " + ", ".join(missing))
        images = read_imaging(local, include_gm=True)
        from .data import derive
        master = derive(clinical.merge(images, on="patient_id", how="left", validate="one_to_one"))
        base, exclusions, base_flow = build_cohort(master, "recovery")
        base_member = clinical.patient_id.isin(base.patient_id)
        field_audit["invalid_base_eligible"] = field_audit.source.map(
            {name: int((mask & base_member).sum()) for name, mask in invalid_masks.items()}).fillna(0).astype(int)
        field_audit.to_csv(folder / "field_audit.csv", index=False)
        data, event_flow = landmark(base)
        master.to_csv(folder / "master.csv", index=False)
        base.to_csv(folder / "base_eligible.csv", index=False)
        data.to_csv(folder / "landmark.csv", index=False)
        exclusions.to_csv(folder / "base_exclusions.csv", index=False)
        base_flow.to_csv(folder / "base_flow.csv", index=False)
        event_flow.to_csv(folder / "event_flow.csv", index=False)
        audit = {"clinical_n": len(clinical), "exact_id_intersection": int(clinical.patient_id.isin(images.patient_id).sum()),
                 "base_eligible_n": len(base), "landmark_n": len(data), "events": int(data.event_type.sum()),
                 "base_state_known_n": int(base.state60.notna().sum()),
                 "landmark_state_known_n": int(data.state60.notna().sum()),
                 "landmark_state_counts": {str(k): int(data.state60.eq(k).sum()) for k in (0, 1, 2)},
                 "event_flow": event_flow.to_dict("records"),
                 "source_absent": field_audit.loc[~field_audit.present, "source"].tolist(),
                 "source_invalid": field_audit.loc[field_audit.invalid.gt(0), ["source", "invalid"]].to_dict("records"),
                 "source_invalid_base_eligible": field_audit.loc[field_audit.invalid_base_eligible.gt(0),
                                                                  ["source", "invalid_base_eligible"]].to_dict("records"),
                 "source_groups": {Path(path).name: ", ".join(group.source) for path, group in
                                   field_audit.loc[field_audit.present & field_audit.source.ne("code_n")].groupby("file")},
                 "flow_excluded": {row.step: int(row.excluded_here) for row in base_flow.itertuples()
                                   if row.excluded_here},
                 "yearly_fields_observed": {str(y): int(_number(base, f"y{y}_stroke_event").notna().sum()) for y in (2, 3, 4)},
                 "chronology": chronology_audit(data)}
        if len(data):
            covariates = cohort_audit(data, "recovery")
            audit["covariate_missing"] = covariates["covariate_missing"]
            audit["covariates_entirely_missing"] = covariates["covariates_entirely_missing"]
            try:
                frozen = StudyDesign.freeze(data, recurrence_spec())
                filled = data[list(frozen.spec.predictors)].copy()
                for col in filled:
                    code = frozen.coding[col]
                    fill = filled[col].mode().iloc[0] if code["kind"] == "category" else filled[col].median()
                    filled[col] = filled[col].fillna(fill)
                audit["model_parameters"] = frozen.transform(filled).shape[1]
                dump_json(folder / "frozen_recurrence_design.json", frozen.coding)
            except (DataError, IndexError) as exc:
                audit["design_blocker"] = str(exc)
        required_invalid = field_audit.loc[field_audit.source.isin(required) & field_audit.invalid_base_eligible.gt(0), "source"].tolist()
        audit["required_invalid_fields"] = required_invalid
        state["audit"] = audit
        if (not len(data) or required_invalid or audit["chronology"]["event_after_known_death_n"]
                or audit.get("covariates_entirely_missing") or audit.get("design_blocker")):
            state["status"] = "REVIEW_REQUIRED"
        else:
            state["status"] = "PREPARED"
        dump_json(folder / "audit.json", audit)
    except (DataError, OSError, ValueError, KeyError) as exc:
        state.update(status="INPUTS_REQUIRED", error=str(exc))
        (folder / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
    dump_json(folder / "status.json", state)
    return state


def _saved_fits(path: Path) -> list[Fit]:
    saved = np.load(path / "pooled_inputs.npz")
    return [Fit(saved["terms"].tolist(), p, u, {}) for p, u in zip(saved["params"], saved["covariance"], strict=True)]


def recurrence_curve(data: pd.DataFrame, design: StudyDesign, model_path: Path, output: Path) -> None:
    fits = _saved_fits(model_path)
    support = np.quantile(data.wmh_ml, [.05, .5, .95])
    grid = np.linspace(support[0], support[2], 81)
    rows = []
    for value in grid:
        vector = design.contrast(data, {"wmh_ml": value}, {"wmh_ml": support[1]})
        pooled = contrast(fits, vector)
        rows.append({"wmh_ml": value, "HR": np.exp(pooled["estimate"]),
                     "lower": np.exp(pooled["lower"]), "upper": np.exp(pooled["upper"]),
                     "reference_wmh_ml": support[1]})
    curve = pd.DataFrame(rows)
    curve.to_csv(output / "wmh_recurrence_curve.csv", index=False)
    fig, (ax, rug) = plt.subplots(2, 1, figsize=(8, 5), sharex=True, height_ratios=[4, 1])
    ax.plot(curve.wmh_ml, curve.HR, color="#1769a4")
    ax.fill_between(curve.wmh_ml, curve.lower, curve.upper, alpha=.18, color="#1769a4")
    ax.axhline(1, color="gray", linestyle="--")
    ax.set(ylabel="Adjusted first-stroke HR vs median WMH", title="Exploratory 3-month landmark association")
    rug.hist(data.wmh_ml, bins=np.linspace(support[0], support[2], 31), color="#768896")
    rug.set(xlabel="Baseline corrected WMH volume (mL); 5th–95th percentile support", ylabel="N")
    fig.tight_layout()
    fig.savefig(output / "wmh_recurrence_curve.png", dpi=180)
    plt.close(fig)


def analyse(state: dict) -> dict:
    """All new estimates are exploratory; failures remain visible and isolated."""
    folder = Path(state["path"])
    data = pd.read_csv(folder / "landmark.csv", dtype={"patient_id": str})
    base = pd.read_csv(folder / "base_eligible.csv", dtype={"patient_id": str})
    cfg = json.loads((folder / "config_snapshot.json").read_text(encoding="utf-8"))
    options = path_settings(cfg)
    state["status"] = "ANALYSING"
    dump_json(folder / "status.json", state)
    outcomes = []
    recurrence, design = run_one(data, recurrence_spec(), folder, options)
    outcomes.append(recurrence)
    if recurrence["status"] == "ESTIMATED":
        recurrence_curve(data, design, folder / "stroke_recurrence", folder)
    shape_spec = primary_spec("recovery").variant(
        "functional_shape", splines=("age", "wmh_ml", "gm119_ml"),
        primary=("dependent:wmh_ml_rcs", "dependent:gm119_ml_rcs"), tier="exploratory")
    shape, _ = run_one(base, shape_spec, folder, options)
    outcomes.append(shape)
    if shape["status"] == "ESTIMATED":
        from .pooling import terms_test
        fits = _saved_fits(folder / "functional_shape")
        dump_json(folder / "functional_shape_tests.json", {
            "wmh_log1p_nonlinearity": terms_test(fits, ("dependent:wmh_ml_rcs",)),
            "gm_raw_nonlinearity": terms_test(fits, ("dependent:gm119_ml_rcs",)),
            "note": "Exploratory shape checks do not replace the original main model"})
    logged = base.assign(gm119_log=np.log(base.gm119_ml))
    log_spec = primary_spec("recovery").variant(
        "functional_gm_log", exposures=("wmh_ml", "gm119_log"),
        primary=("dependent:gm119_log",), tier="exploratory")
    gm_log, _ = run_one(logged, log_spec, folder, options)
    outcomes.append(gm_log)
    from .recovery_mediation import estimate, feasibility
    gate = feasibility(data, state["audit"]["chronology"])
    dump_json(folder / "mediation_gate.json", gate)
    try:
        mediation = estimate(data, folder, options, gate)
    except (DataError, ValueError, np.linalg.LinAlgError) as exc:
        mediation = {"status": "NOT_ESTIMABLE", "reasons": [str(exc)], "gate": gate}
        dump_json(folder / "mediation.json", mediation)
    from .recovery_prediction import validate
    try:
        validate(base, folder / "prediction", options)
        prediction_status = "ESTIMATED_INTERNAL_VALIDATION"
    except (DataError, ValueError, np.linalg.LinAlgError) as exc:
        prediction_status = "NOT_ESTIMABLE"
        dump_json(folder / "prediction_failure.json", {"reason": str(exc)})
    pd.DataFrame(outcomes).to_csv(folder / "association_results.csv", index=False)
    state["analysis"] = {"recurrence": recurrence["status"], "functional_shape": shape["status"],
                         "functional_gm_log": gm_log["status"], "mediation": mediation["status"],
                         "prediction": prediction_status}
    if recurrence["status"] != "ESTIMATED":
        state["status"] = "PRIMARY_NOT_ESTIMABLE"
    elif shape["status"] != "ESTIMATED" or gm_log["status"] != "ESTIMATED" or prediction_status == "NOT_ESTIMABLE":
        state["status"] = "PARTIAL"
    else:
        state["status"] = "COMPLETED"
    state["finished_utc"] = datetime.now(UTC).isoformat()
    dump_json(folder / "status.json", state)
    pointer = {"path": str(folder), "run": state["run"], "status": state["status"]}
    dump_json(root(cfg) / "latest_path_attempt.json", pointer)
    if state["status"] == "COMPLETED":
        dump_json(root(cfg) / "latest_path_results.json", pointer)
    return state


def run(cfg: dict, through: str = "analyse") -> dict:
    if through not in {"prepare", "analyse"}:
        raise DataError("recovery-path supports prepare or analyse")
    state = prepare(cfg)
    if through == "prepare" or state["status"] != "PREPARED":
        return state
    try:
        return analyse(state)
    except (DataError, OSError, ValueError, KeyError, np.linalg.LinAlgError) as exc:
        state.update(status="FAILED", error=str(exc))
        folder = Path(state["path"])
        (folder / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
        dump_json(folder / "status.json", state)
        dump_json(root(cfg) / "latest_path_attempt.json",
                  {"path": str(folder), "run": state["run"], "status": state["status"]})
        return state
