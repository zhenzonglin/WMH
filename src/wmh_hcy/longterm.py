"""Cumulative 2-5 year extensions, isolated from the frozen one-year analysis."""

from __future__ import annotations

import copy
from datetime import UTC, datetime

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .common import DataError, dump_json, outdir, read_json, record_run
from .fields import LONGTERM_YEARS
from .imaging import available_volume

CONTRACT = "cumulative_first_is_event_or_censor_days_from_onset_20260917"
FAMILIES = ("H2", "H3", "H4", "H4_T1")


def values(data: pd.DataFrame, name: str) -> pd.Series:
    return pd.to_numeric(data.get(name, pd.Series(np.nan, index=data.index)), errors="coerce")


def endpoint_columns(year: int) -> tuple[str, str]:
    return f"y{year}_is_event", f"y{year}_is_day"


def consistency_issues(data: pd.DataFrame, year: int) -> pd.Series:
    """Flag conflicting first-event records; never merge or replace endpoints."""
    event, stop = (values(data, c) for c in endpoint_columns(year))
    bad = pd.Series(False, index=data.index)
    for earlier in range(1, year):
        names = ("is_event", "is_day") if earlier == 1 else endpoint_columns(earlier)
        prev, time = (values(data, c) for c in names)
        known = prev.eq(1) & np.isfinite(time) & time.between(0, 365 * earlier)
        bad |= known & (event.eq(0) | (event.eq(1) & stop.ne(time)))
        prior_stop = values(data, "last_contact_day") if earlier == 1 else time
        bad |= prev.eq(0) & event.eq(1) & stop.le(prior_stop.clip(upper=365 * earlier))
    return bad


def audit_longterm(data: pd.DataFrame, years: list[int]) -> pd.DataFrame:
    rows = []
    absent = set(data.attrs.get("absent_source_columns", []))
    for year in years:
        event, stop = (values(data, c) for c in endpoint_columns(year))
        mrs = values(data, f"mrs{12 * year}")
        state = event.isin([0, 1])
        time = np.isfinite(stop) & stop.ge(0)
        source = f"Y{year}_IS" if year == 5 else f"y{year}_IS"
        missing = [
            c for c in [source, source + ("_DD" if year == 5 else "_dd"), f"m{12 * year}_mrs"] if c in absent
        ]
        rows.append(
            {
                "year": year,
                "horizon_days": 365 * year,
                "patients": len(data),
                "missing_fields": ";".join(missing),
                "event_1": int(event.eq(1).sum()),
                "event_0": int(event.eq(0).sum()),
                "event_unknown_or_invalid": int((~state).sum()),
                "event_with_time": int((event.eq(1) & time).sum()),
                "non_event_with_time": int((event.eq(0) & time).sum()),
                "non_event_before_horizon": int((event.eq(0) & time & stop.lt(365 * year)).sum()),
                "event_after_horizon": int((event.eq(1) & time & stop.gt(365 * year)).sum()),
                "cross_year_conflicts": int(consistency_issues(data, year).sum()),
                "mrs_observed_0_6": int(mrs.isin(range(7)).sum()),
                "mrs_6": int(mrs.eq(6).sum()),
                "mrs_invalid": int((mrs.notna() & ~mrs.isin(range(7))).sum()),
            }
        )
    return pd.DataFrame(rows)


def apply_rules(
    data: pd.DataFrame, rules: list[tuple[str, pd.Series]]
) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    keep = pd.Series(True, index=data.index)
    exclusions = data[["patient_id"]].copy()
    exclusions["exclusion_reason"] = ""
    flow = [{"step": "clinical_image_union", "remaining": len(data), "excluded_here": 0}]
    for name, rule in rules:
        rule = rule.fillna(False)
        rejected = keep & ~rule
        exclusions.loc[rejected, "exclusion_reason"] = name
        keep &= rule
        flow.append({"step": name, "remaining": int(keep.sum()), "excluded_here": int(rejected.sum())})
    return data.loc[keep].copy(), exclusions, flow


def build_longterm_cohorts(data: pd.DataFrame, year: int) -> tuple[dict, dict, dict]:
    if year not in LONGTERM_YEARS:
        raise DataError("Long-term year must be 2, 3, 4 or 5")
    d = data.copy()
    horizon = 365 * year
    event, stop = (values(d, c) for c in endpoint_columns(year))
    d["entry"] = d.sample_day
    d["exit"] = stop.clip(upper=horizon)
    d["event_type"] = (event.eq(1) & stop.le(horizon)).astype(int)
    # No invented death dates or nominal censoring dates. All exits come from IS_DD.
    baseline = [
        ("adult_ischemic_stroke", d.age.ge(18) & d.diagnosis.eq(1)),
        ("available_wmh_icv", d.image_valid.fillna(False).astype(bool)),
        ("observed_positive_baseline_hcy", np.isfinite(d.hcy) & d.hcy.gt(0)),
        ("known_baseline_measurement_time", np.isfinite(d.entry) & d.entry.ge(0) & d.entry.lt(horizon)),
    ]
    death = values(d, "death_day")  # Only existing, explicitly dated deaths.
    conflict = consistency_issues(d, year)
    rules = baseline + [
        ("known_longterm_is_status", event.isin([0, 1])),
        ("observed_event_or_censor_time", np.isfinite(stop) & stop.ge(0)),
        ("consistent_cumulative_first_event", ~conflict),
        ("not_observed_after_known_death", death.isna() | stop.le(death)),
        ("observed_followup_after_entry", d.exit.gt(d.entry)),
    ]
    main, exclusions, flow = apply_rules(d, rules)
    m3_rules = rules + [
        ("observed_positive_month3_hcy", np.isfinite(d.hcy3) & d.hcy3.gt(0)),
        (
            "observed_month3_sample_before_exit",
            np.isfinite(d.sample3_day) & d.sample3_day.ge(d.entry) & d.sample3_day.lt(d.exit),
        ),
    ]
    m3, ex3, flow3 = apply_rules(d, m3_rules)
    m3["entry"] = m3.sample3_day

    # Functional samples do not require the recurrence endpoint or month-3 blood draw.
    outcome = f"mrs{year * 12}"
    supplied = values(d, outcome)
    known_dead = death.notna() & death.between(0, horizon)
    impossible_score = known_dead & supplied.isin(range(6))
    for earlier in range(2, year):
        impossible_score |= values(d, f"mrs{12 * earlier}").eq(6) & supplied.isin(range(6))
    d[outcome] = supplied.mask(known_dead & supplied.isna(), 6)
    frules = baseline + [
        ("alive_at_baseline_sampling", death.isna() | death.gt(d.entry)),
        ("mrs_consistent_with_known_death", ~impossible_score),
        ("observed_year_mrs_0_6", d[outcome].isin(range(7))),
        ("available_acute_lesion", available_volume(d.lesion_ml, allow_zero=True)),
    ]
    func, exf, flowf = apply_rules(d, frules)
    t1, ext, flowt = apply_rules(d, frules + [("available_t1_gm119", available_volume(d.gm119_ml))])
    return (
        {"H2": main, "H3": m3, "H4": func, "H4_T1": t1},
        {"H2": exclusions, "H3": ex3, "H4": exf, "H4_T1": ext},
        {"H2": flow, "H3": flow3, "H4": flowf, "H4_T1": flowt},
    )


def summarize_estimates(rows: list[dict]) -> pd.DataFrame:
    """Fixed four-year Holm families, including unestimated years."""
    table = pd.DataFrame(rows)
    for family in FAMILIES:
        mask = table.family.eq(family)
        indices = table.index[mask].tolist()
        if not indices:
            continue
        p = pd.to_numeric(table.loc[indices, "p"], errors="coerce")
        filled = p.fillna(1).tolist() + [1.0] * (4 - len(indices))
        adjusted = multipletests(filled, method="holm")[1][: len(indices)]
        table.loc[indices, "p_holm_4"] = np.where(p.notna(), adjusted, np.nan)
    for key in ["estimate", "lower", "upper"]:
        with np.errstate(over="ignore"):
            table["ratio" if key == "estimate" else f"ratio_{key}"] = np.exp(table[key].astype(float))
    table["interpretation"] = "supplementary_exploratory; pointwise_95pct_CI; four_year_Holm_within_family"
    return table


def run_longterm(
    cfg: dict, years: list[int] | None = None, through: str = "prepare", only: str | None = None
) -> dict:
    from .analysis import run_functional, run_survival
    from .design import Design
    from .harmonize import harmonize
    from .imaging import read_imaging
    from .longterm_reporting import report_longterm
    from .sas_extract import extract

    years = sorted(set(years or LONGTERM_YEARS))
    if any(year not in LONGTERM_YEARS for year in years) or only not in {None, "H2", "H3", "H4"}:
        raise DataError("longterm supports years 2-5 and hypotheses H2, H3 or H4")
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    target = outdir(cfg) / "longterm"
    root = target / "results" / stamp
    work = copy.deepcopy(cfg)
    work["_out"] = str(root / "inputs_derived")
    status = {
        "mode": cfg["mode"],
        "contract": CONTRACT,
        "years": years,
        "result_dir": str(root),
        "status": "RUNNING",
        "analyses": {},
        "absolute_risk_status": "NOT_REQUESTED_NO_LONGTERM_DEATH_TIMES",
    }
    dump_json(root / "status.json", status)
    dump_json(root / "config_snapshot.json", cfg)
    dump_json(target / "latest_results.json", {"path": str(root)})
    try:
        if not cfg["inputs"].get("clinical_csv"):
            print("Long-term: extracting supplied SAS fields", flush=True)
            extract(work)
        clinical = harmonize(work, require_one_year=False)
    except (DataError, OSError, ValueError) as exc:
        status.update(status="FAILED", reason=str(exc))
        dump_json(root / "status.json", status)
        raise
    audit = audit_longterm(clinical, years)
    audit.to_csv(root / "field_audit.csv", index=False)
    status["field_audit"] = audit.to_dict("records")
    print(
        audit[
            ["year", "event_1", "event_0", "non_event_with_time", "mrs_observed_0_6", "cross_year_conflicts"]
        ].to_string(index=False),
        flush=True,
    )
    if through == "audit":
        status["status"] = "AUDIT_COMPLETED"
        dump_json(root / "status.json", status)
        return status
    try:
        images = read_imaging(work)
        master = clinical.merge(images, on="patient_id", how="left", validate="one_to_one")
    except (DataError, OSError, ValueError) as exc:
        status.update(status="FAILED", reason=str(exc))
        dump_json(root / "status.json", status)
        raise
    counts, rows = [], []
    for year in years:
        cohorts, exclusions, flows = build_longterm_cohorts(master, year)
        for family in FAMILIES:
            frame = cohorts[family]
            folder = root / f"year{year}" / family
            folder.mkdir(parents=True, exist_ok=True)
            frame.to_csv(folder / "cohort.csv", index=False)
            exclusions[family].to_csv(folder / "exclusions.csv", index=False)
            pd.DataFrame(flows[family]).to_csv(folder / "flow.csv", index=False)
            count = {
                "year": year,
                "family": family,
                "n": len(frame),
                "events": int(frame.event_type.eq(1).sum()) if family in {"H2", "H3"} else None,
                "mrs_6": int(frame[f"mrs{12 * year}"].eq(6).sum()) if family.startswith("H4") else None,
            }
            counts.append(count)
            selected = (only is None or family.startswith(only)) and (
                family != "H4_T1" or cfg["analysis"].get("run_t1_extension", True)
            )
            result = {"status": "PREPARED" if through == "prepare" else "NOT_REQUESTED", **count}
            if through in {"analyse", "report"} and selected:
                print(f"Long-term analysis: year {year} {family} n={len(frame)}", flush=True)
                try:
                    if family in {"H2", "H3"}:
                        result.update(
                            run_survival(
                                frame,
                                cfg,
                                folder,
                                f"{family}_year{year}",
                                spec=Design(month3=family == "H3"),
                                kind="month3" if family == "H3" else "main",
                                risk=False,
                                death_auxiliaries=False,
                            )
                        )
                    else:
                        result.update(
                            run_functional(
                                frame,
                                cfg,
                                folder,
                                kind="functional_t1" if family == "H4_T1" else "functional",
                                outcome=f"mrs{12 * year}",
                                survival_auxiliaries=False,
                            )
                        )
                except (DataError, ValueError, np.linalg.LinAlgError) as exc:
                    result.update(status="NOT_ESTIMABLE", reason=str(exc))
            result.update(
                year=year,
                horizon_days=365 * year,
                family=family,
                absolute_risk_status="NOT_REQUESTED_NO_LONGTERM_DEATH_TIMES",
            )
            status["analyses"][f"year{year}/{family}"] = result
            estimate = read_json(folder / "hypothesis_estimate.json")
            rows.append(
                {
                    **count,
                    "status": result["status"],
                    "reason": result.get("reason", ""),
                    "estimate": estimate.get("estimate", np.nan),
                    "lower": estimate.get("lower", np.nan),
                    "upper": estimate.get("upper", np.nan),
                    "p": estimate.get("p", np.nan),
                    "estimand": estimate.get("estimand", "not_estimated"),
                }
            )
            dump_json(root / "status.json", status)
    pd.DataFrame(counts).to_csv(root / "cohort_counts.csv", index=False)
    summarize_estimates(rows).to_csv(root / "longterm_summary.csv", index=False)
    status["status"] = (
        "PREPARED"
        if through == "prepare"
        else (
            "COMPLETED_WITH_MODEL_FAILURES"
            if any(r["status"] == "NOT_ESTIMABLE" for r in rows)
            else "COMPLETED"
        )
    )
    dump_json(root / "status.json", status)
    if through == "report":
        status["report"] = str(report_longterm(root))
        dump_json(root / "status.json", status)
    record_run(work, "longterm", {"result_dir": str(root), "status": status["status"]})
    return status
