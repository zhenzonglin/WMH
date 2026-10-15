"""Read and validate explicitly selected saved runs; never refit or mutate inputs.

Trajectory exponentiated coefficients are multinomial relative risk ratios:
the ratio P(destination)/P(staying in the origin state), compared per 1 SD of
exposure. They are neither a simple transition probability ratio nor a hazard
ratio. Patient-cluster sandwich inference is retained exactly as saved.
"""
from __future__ import annotations

import hashlib
import json
from itertools import pairwise, product
from pathlib import Path

import numpy as np
import pandas as pd

CONTINUOUS = ("age", "wmh_ml", "gm119_ml", "icv_ml", "lesion_ml", "nihss")
CATEGORICAL = ("sex", "smoking", "drinking", "hypertension", "diabetes", "prior_stroke",
               "education", "pre_mrs", "toast", "mrs3")
STATES = ("independent", "dependent", "dead")
MONTHS = (3, 12, 24, 36, 48, 60)
TRANSITIONS = ("independent_to_dependent", "dependent_to_independent",
               "independent_to_dead", "dependent_to_dead")
COUNT_TRANSITIONS = ("independent_to_independent", "independent_to_dependent", "independent_to_dead",
                     "dependent_to_independent", "dependent_to_dependent", "dependent_to_dead")
PRIMARY_TERMS = tuple(f"{state}:{exposure}" for state in ("dependent", "dead")
                      for exposure in ("wmh_ml", "gm119_ml"))
RECOVERY_ANALYSES = {"primary", "complete_case", "month12", "month24", "month36", "month48",
                     "actual_improvement", "raw_wmh", "dependent_or_dead", "observation_weighted"}
RATIO_COLUMNS = ("ratio", "ratio_lower", "ratio_upper", "p")
COUNT_COLUMNS = ("interval_index", "from_month", "to_month", "eligible_start_n",
                 "known_adjacent_state_n", "transition_rows_n", "unknown_start_n", "unknown_end_n",
                 "conflict_n", "death_start_n", "death_to_alive_n")


def _json(path: Path) -> dict | list:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read saved JSON artifact: {path.name}") from exc


def _object(path: Path) -> dict:
    value = _json(path)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path.name}")  # noqa: TRY004 - saved-artifact validation API
    return value


def _csv(path: Path, required: tuple[str, ...], optional: tuple[str, ...] = ()) -> pd.DataFrame:
    try:
        frame = pd.read_csv(path, usecols=lambda c: c in (*required, *optional),
                            dtype={"patient_id": str})
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read saved CSV artifact: {path.name}") from exc
    missing = set(required) - set(frame)
    if missing:
        raise ValueError(f"{path.name} lacks required columns: {', '.join(sorted(missing))}")
    return frame


def _numeric(frame: pd.DataFrame, columns: tuple[str, ...], name: str, *, missing: bool = False) -> None:
    for column in columns:
        try:
            values = pd.to_numeric(frame[column], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name}: nonnumeric {column}") from exc
        if np.isinf(values.to_numpy(float)).any() or (not missing and values.isna().any()):
            raise ValueError(f"{name}: missing or nonfinite {column}")
        frame[column] = values


def _count(value: object, label: str) -> int:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid saved count: {label}") from exc
    if not np.isfinite(number) or number < 0 or number != int(number):
        raise ValueError(f"Invalid saved count: {label}")
    return int(number)


def _counts(frame: pd.DataFrame, columns: tuple[str, ...], name: str) -> None:
    _numeric(frame, columns, name)
    for column in columns:
        values = frame[column].to_numpy(float)
        if (values < 0).any() or (values != np.floor(values)).any():
            raise ValueError(f"{name}: invalid count {column}")


def _p_value(value: object, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid saved P value: {label}") from exc
    if not np.isfinite(number) or not 0 <= number <= 1:
        raise ValueError(f"Invalid saved P value: {label}")
    return number


def _identity(run: Path, study: str, allow_demo: bool) -> tuple[Path, dict]:
    run = Path(run).resolve()
    status = _object(run / "status.json")
    if status.get("study") != study or status.get("run") != run.name:
        raise ValueError(f"Run identity mismatch for {study}")
    if not status.get("path") or Path(status["path"]).resolve() != run:
        raise ValueError(f"Saved path mismatch for {study}")
    mode = status.get("mode")
    if mode != "real":
        named_demo = any("demo" in part.casefold() or "synthetic" in part.casefold() for part in run.parts)
        if not (allow_demo and mode in {"demo", "synthetic"} and named_demo):
            raise ValueError("Publication input requires mode=real; demo requires explicitly named synthetic input")
    elif allow_demo:
        raise ValueError("Demo rendering cannot consume a real run")
    if status.get("status") not in {"COMPLETED", "PARTIAL"}:
        raise ValueError(f"Saved {study} run is not an analysis result: {status.get('status')}")
    return run, status


def _ratios(frame: pd.DataFrame, name: str) -> None:
    _numeric(frame, RATIO_COLUMNS, name)
    if ((frame[list(RATIO_COLUMNS[:3])] <= 0).any().any()
            or (frame.ratio_lower > frame.ratio).any() or (frame.ratio > frame.ratio_upper).any()
            or not frame.p.between(0, 1).all()):
        raise ValueError(f"{name}: malformed ratio, confidence interval, or P value")


def _probabilities(frame: pd.DataFrame, keys: list[str], name: str) -> None:
    _numeric(frame, ("wmh_ml", "gm119_ml", "probability", "lower", "upper"), name)
    if frame.duplicated([*keys, "state"]).any():
        raise ValueError(f"{name}: duplicate probability scenario/state")
    if set(frame.state) != set(STATES) or not frame.probability.between(0, 1).all():
        raise ValueError(f"{name}: invalid state or probability")
    if (frame.lower > frame.probability).any() or (frame.upper < frame.probability).any():
        raise ValueError(f"{name}: malformed confidence interval")
    for _, group in frame.groupby(keys, dropna=False):
        if len(group) != 3 or not np.isclose(group.probability.sum(), 1, atol=1e-4, rtol=0):
            raise ValueError(f"{name}: state probabilities do not form a complete unit-sum triplet")


def _recovery_cohort(run: Path, audit: dict) -> pd.DataFrame:
    cohort = _csv(run / "eligible.csv", ("patient_id", "state60", *CONTINUOUS, *CATEGORICAL))
    if cohort.patient_id.isna().any() or cohort.patient_id.str.strip().eq("").any():
        raise ValueError("Recovery cohort contains missing patient identifiers")
    if cohort.patient_id.duplicated().any():
        raise ValueError("Recovery cohort contains duplicate patient identifiers")
    cohort = cohort.drop(columns="patient_id")
    _numeric(cohort, ("state60", *CONTINUOUS, *CATEGORICAL), "eligible.csv", missing=True)
    if not cohort.state60.dropna().isin((0, 1, 2)).all():
        raise ValueError("Recovery cohort contains invalid state60 codes")
    actual = {"eligible_n": len(cohort), "outcome_observed_n": int(cohort.state60.notna().sum()),
              "independent": int(cohort.state60.eq(0).sum()), "dependent": int(cohort.state60.eq(1).sum()),
              "dead": int(cohort.state60.eq(2).sum()), "unknown": int(cohort.state60.isna().sum())}
    for key, value in actual.items():
        if key not in audit or _count(audit[key], key) != value:
            raise ValueError(f"Recovery cohort/audit mismatch: {key}")
    return cohort


def _recovery_flow(flow: pd.DataFrame, audit: dict) -> None:
    _counts(flow, ("remaining", "excluded_here"), "cohort_flow.csv")
    if flow.step.isna().any() or flow.step.duplicated().any():
        raise ValueError("Recovery flow contains missing or duplicated steps")
    if len(flow) < 2 or flow.iloc[-1].step != "known_five_year_functional_state":
        raise ValueError("Recovery flow lacks the final known-five-year-outcome step")
    first = flow.iloc[0]
    initial = int(first.remaining + first.excluded_here)
    if "clinical_n" in audit and initial != _count(audit["clinical_n"], "clinical_n"):
        raise ValueError("Recovery initial flow count differs from clinical audit")
    previous = initial
    for row in flow.itertuples(index=False):
        if previous - row.excluded_here != row.remaining:
            raise ValueError("Recovery flow exclusions do not reconcile sequentially")
        previous = row.remaining
    last = flow.iloc[-1]
    if (int(flow.iloc[-2].remaining) != audit["eligible_n"]
            or int(last.remaining) != audit["outcome_observed_n"]
            or int(last.excluded_here) != audit["unknown"]):
        raise ValueError("Recovery flow final eligible/known/unknown counts disagree with audit")


def load_recovery(run: Path, allow_demo: bool = False) -> dict:
    """Load pre-imputation TableOne inputs and saved functional-model estimates.

    The eligible cohort includes missing five-year outcomes. Patient identifiers
    are used only to validate uniqueness, then removed from the returned frame.
    """
    run, status = _identity(run, "recovery", allow_demo)
    audit = _object(run / "audit.json")
    cohort = _recovery_cohort(run, audit)
    flow = _csv(run / "cohort_flow.csv", ("step", "remaining", "excluded_here"))
    _recovery_flow(flow, audit)
    results = _csv(run / "results.csv", ("analysis", "status"),
                   ("n", "p", "reason", "family", "tier", "p_bh_secondary"))
    if results.analysis.isna().any() or results.analysis.duplicated().any():
        raise ValueError("Recovery results contain missing or duplicate analyses")
    if "primary" not in set(results.analysis):
        raise ValueError("Primary result is missing")
    coefficients, details = {}, {}
    for index, row in results.iterrows():
        analysis = str(row.analysis)
        if "/" in analysis or "\\" in analysis or analysis not in RECOVERY_ANALYSES:
            raise ValueError("Invalid saved analysis name")
        path = run / analysis / "result.json"
        detail = _object(path) if path.is_file() else row.dropna().to_dict()
        if detail.get("analysis", analysis) != analysis or detail.get("status") != row.status:
            raise ValueError(f"Saved result identity/status mismatch: {analysis}")
        details[analysis] = detail
        if row.status != "ESTIMATED":
            for field in ("p", "p_bh_secondary", "estimate", "lower", "upper", "ratio", "ratio_lower",
                          "ratio_upper", "se", "statistic", "df", "df1", "df2"):
                detail.pop(field, None)
            for field in ("p", "p_bh_secondary"):
                if field in results:
                    results.loc[index, field] = np.nan
            if analysis == "primary":
                raise ValueError("Primary recovery model is not ESTIMATED")
            continue
        n = _count(detail.get("n"), f"{analysis}.n")
        if "n" in row and pd.notna(row["n"]) and _count(row["n"], f"{analysis}.results.n") != n:
            raise ValueError(f"Result sample sizes differ: {analysis}")
        if analysis == "primary" and n != audit["outcome_observed_n"]:
            raise ValueError("Primary model N differs from observed five-year cohort")
        detail_p = _p_value(detail.get("p"), f"{analysis}.result.p")
        summary_p = _p_value(row.get("p"), f"{analysis}.results.p")
        if not np.isclose(detail_p, summary_p, rtol=1e-6, atol=1e-300):
            raise ValueError(f"Result P values differ: {analysis}")
        table = _csv(run / analysis / "coefficients.csv", ("term", *RATIO_COLUMNS),
                     ("estimate", "lower", "upper", "se", "df", "m", "scale"))
        if table.term.isna().any() or table.term.duplicated().any():
            raise ValueError(f"{analysis}: duplicate or missing coefficient terms")
        required = ("wmh_ml", "gm119_ml") if analysis == "dependent_or_dead" else PRIMARY_TERMS
        if not set(required).issubset(table.term):
            raise ValueError(f"{analysis}: incomplete imaging coefficient set")
        # Publication artifacts display only these predeclared imaging contrasts.
        # Exponentiating an extreme nuisance intercept can overflow despite finite
        # fitted log coefficients; it must not invalidate unrelated saved contrasts.
        table = table.loc[table.term.isin(required)].copy().reset_index(drop=True)
        _ratios(table, f"{analysis}/coefficients.csv")
        coefficients[analysis] = table
    design = _object(run / "primary" / "frozen_design.json")
    for exposure in ("wmh_ml", "gm119_ml"):
        code = design.get(exposure, {})
        try:
            scale = float(code["scale"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Missing frozen exposure scale: {exposure}") from exc
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"Invalid frozen exposure scale: {exposure}")
    probabilities = _csv(run / "primary" / "standardized_states.csv",
                         ("wmh_ml", "gm119_ml", "state", "probability", "lower", "upper"), ("ci",))
    _probabilities(probabilities, ["wmh_ml", "gm119_ml"], "standardized_states.csv")
    if (probabilities.lower < 0).any() or (probabilities.upper > 1).any():
        raise ValueError("Primary standardized probability confidence intervals must be within [0,1]")
    pairs = set(zip(probabilities.wmh_ml, probabilities.gm119_ml, strict=True))
    if (len(probabilities) != 27 or probabilities.wmh_ml.nunique() != 3
            or probabilities.gm119_ml.nunique() != 3
            or pairs != set(product(probabilities.wmh_ml.unique(), probabilities.gm119_ml.unique()))):
        raise ValueError("Primary standardized probabilities require the complete 9-by-3 grid")
    known = cohort.loc[cohort.state60.notna()]
    for exposure in ("wmh_ml", "gm119_ml"):
        expected_values = known[exposure].quantile([.25, .5, .75]).to_numpy(float)
        if not np.allclose(np.sort(probabilities[exposure].unique()), expected_values, atol=1e-4, rtol=0):
            raise ValueError(f"Primary standardized {exposure} scenarios differ from observed-cohort quartiles")
    return {"status": status, "audit": audit, "cohort": cohort, "flow": flow, "results": results,
            "coefficients": coefficients, "result_details": details, "design": design,
            "probabilities": probabilities}


def load_trajectory(run: Path, allow_demo: bool = False) -> dict:
    """Load saved transition RRRs (destination versus stay in origin) and occupancy."""
    run, status = _identity(run, "recovery_trajectory", allow_demo)
    if status.get("analysis_scope") != "trajectory_only":
        raise ValueError("Trajectory input has unexpected analysis scope")
    analysis = status.get("analysis", {})
    if analysis.get("status") != "ESTIMATED":
        raise ValueError("Trajectory model is not ESTIMATED")
    audit = _object(run / "audit.json")
    if audit.get("future_recurrence_fields_requested") != []:
        raise ValueError("Trajectory audit must confirm no future recurrence fields were requested")
    base_n = _count(audit.get("base_eligible_n"), "base_eligible_n")
    state_counts = audit.get("state_counts", {})
    for month in MONTHS:
        counts = state_counts.get(str(month), {})
        if set(counts) != {"0", "1", "2", "unknown"}:
            raise ValueError(f"Incomplete saved state counts at month {month}")
        if sum(_count(counts[k], f"month{month}.{k}") for k in counts) != base_n:
            raise ValueError(f"State counts differ from base cohort at month {month}")
    effects = _csv(run / "transition_effects.csv", ("exposure", "transition", *RATIO_COLUMNS, "m"),
                   ("estimate", "lower", "upper", "scale"))
    expected = set(product(("wmh_ml", "gm119_ml"), TRANSITIONS))
    if effects.duplicated(["exposure", "transition"]).any() or set(
            zip(effects.exposure, effects.transition, strict=True)) != expected:
        raise ValueError("Trajectory effects require eight unique imaging transitions")
    _ratios(effects, "transition_effects.csv")
    _counts(effects, ("m",), "transition_effects.csv")
    m = _count(analysis.get("imputations"), "trajectory.imputations")
    if m == 0 or not effects.m.eq(m).all():
        raise ValueError("Trajectory imputation counts disagree")
    probabilities = _csv(run / "state_probabilities_36_48_60.csv",
                         ("scenario", "month", "state", "wmh_ml", "gm119_ml", "probability", "lower", "upper"),
                         ("ci_method",))
    _counts(probabilities, ("month",), "state_probabilities_36_48_60.csv")
    _probabilities(probabilities, ["scenario", "month"], "state_probabilities_36_48_60.csv")
    scenarios = ("wmh_q25", "wmh_q50", "wmh_q75", "gm_q25", "gm_q50", "gm_q75")
    actual = set(zip(probabilities.scenario, probabilities.month, probabilities.state, strict=True))
    if len(probabilities) != 54 or actual != set(product(scenarios, (36, 48, 60), STATES)):
        raise ValueError("Trajectory probabilities require the complete 6-by-3-by-3 grid")
    for _, group in probabilities.groupby("scenario"):
        if group.wmh_ml.nunique() != 1 or group.gm119_ml.nunique() != 1:
            raise ValueError("Trajectory exposure scenario changes between states or months")
    shared_columns = ["wmh_ml", "gm119_ml", "probability", "lower", "upper"]
    medians = [probabilities.loc[probabilities.scenario.eq(name)].sort_values(["month", "state"])
               for name in ("wmh_q50", "gm_q50")]
    if not np.allclose(medians[0][shared_columns], medians[1][shared_columns], atol=1e-6, rtol=0):
        raise ValueError("Trajectory shared median scenarios disagree")
    counts = _csv(run / "transition_counts.csv", (*COUNT_COLUMNS, "death_absorbing_check",
                  *(name + "_n" for name in COUNT_TRANSITIONS)))
    count_columns = (*COUNT_COLUMNS, *(name + "_n" for name in COUNT_TRANSITIONS))
    _counts(counts, count_columns, "transition_counts.csv")
    expected_intervals = {(i, start, end) for i, (start, end) in enumerate(pairwise(MONTHS))}
    actual_intervals = set(zip(counts.interval_index, counts.from_month, counts.to_month, strict=True))
    if len(counts) != 5 or actual_intervals != expected_intervals:
        raise ValueError("Trajectory interval count grid is incomplete or duplicated")
    if counts.death_to_alive_n.any() or not counts.death_absorbing_check.astype(str).str.lower().eq("true").all():
        raise ValueError("Saved trajectory audit failed death absorption")
    if not counts[[name + "_n" for name in COUNT_TRANSITIONS]].sum(axis=1).eq(counts.transition_rows_n).all():
        raise ValueError("Transition event counts do not match interval row counts")
    total = int(counts.transition_rows_n.sum())
    if (total != _count(audit.get("transition_rows_n"), "transition_rows_n")
            or total != _count(analysis.get("transition_rows"), "analysis.transition_rows")):
        raise ValueError("Trajectory row counts disagree across saved artifacts")
    mrs = _csv(run / "mrs_trajectory_distribution.csv", ("month", "mrs", "n", "denominator_known", "proportion"))
    mrs["mrs"] = mrs.mrs.astype(str)
    _counts(mrs, ("month", "n", "denominator_known"), "mrs_trajectory_distribution.csv")
    _numeric(mrs, ("proportion",), "mrs_trajectory_distribution.csv", missing=True)
    if (mrs.duplicated(["month", "mrs"]).any() or set(zip(mrs.month, mrs.mrs, strict=True))
            != set(product(MONTHS, (*map(str, range(7)), "unknown")))):
        raise ValueError("Descriptive mRS grid is incomplete or duplicated")
    for month, group in mrs.groupby("month"):
        known = group.loc[group.mrs.ne("unknown")]
        unknown = group.loc[group.mrs.eq("unknown")]
        denominator = int(known.n.sum())
        if (denominator == 0 or not group.denominator_known.eq(denominator).all()
                or int(unknown.n.iloc[0]) + denominator != base_n
                or not unknown.proportion.isna().all()
                or not np.allclose(known.proportion, known.n / denominator, atol=1e-4, rtol=0)):
            raise ValueError(f"Invalid descriptive mRS denominator/proportions at month {month}")
        collapsed = {
            "0": int(group.loc[group.mrs.isin(("0", "1", "2")), "n"].sum()),
            "1": int(group.loc[group.mrs.isin(("3", "4", "5")), "n"].sum()),
            "2": int(group.loc[group.mrs.eq("6"), "n"].sum()),
            "unknown": int(unknown.n.iloc[0]),
        }
        if collapsed != state_counts[str(month)]:
            raise ValueError(f"Descriptive mRS counts differ from state audit at month {month}")
    diagnostics = _json(run / "fit_diagnostics.json")
    if not isinstance(diagnostics, list) or len(diagnostics) != m:
        raise ValueError("Trajectory fit diagnostics do not match imputation count")
    for item in diagnostics:
        if not isinstance(item, dict) or item.get("converged") is not True or _count(item.get("n"), "fit.n") != total:
            raise ValueError("Trajectory saved fit is incomplete, unconverged, or has a mismatched row count")
    warnings = []
    draws = analysis.get("trajectory_probability_draws")
    if draws is None:
        warnings.append("Saved trajectory probability simulation draw count is unknown; interval stability unverified.")
    elif _count(draws, "trajectory_probability_draws") < 200:
        warnings.append(f"Saved trajectory probability intervals use {draws} draws per imputation (<200); "
                        "Monte Carlo stability requires separate review. Saved intervals were not changed.")
    if (probabilities.lower < 0).any() or (probabilities.upper > 1).any():
        warnings.append("Saved normal-scale probability interval extends outside [0,1]; source values retained.")
    return {"status": status, "audit": audit, "effects": effects, "probabilities": probabilities,
            "counts": counts, "mrs": mrs, "diagnostics": diagnostics, "warnings": warnings}


def source_manifest(run: Path, filenames: list[str]) -> list[dict]:
    """Hash only named source artifacts; never collect rows or scan a run folder."""
    run = Path(run).resolve()
    manifest = []
    for filename in filenames:
        relative = Path(filename)
        path = (run / relative).resolve()
        if relative.is_absolute() or ".." in relative.parts or not path.is_relative_to(run):
            raise ValueError("Manifest artifact must be a relative path inside the selected run")
        if not path.is_file():
            raise ValueError(f"Missing manifest source artifact: {relative.as_posix()}")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        manifest.append({"artifact": relative.as_posix(), "sha256": digest.hexdigest(), "bytes": path.stat().st_size})
    return manifest
