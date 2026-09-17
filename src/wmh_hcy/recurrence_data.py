"""One five-year cohort; every horizon truncates the same supplied endpoint."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .adjustment import adjustment_columns
from .common import DataError
from .fields import FIELDS
from .imputation import nelson_aalen_increment
from .longterm import apply_rules, consistency_issues, values

CONTRACT = "recurrence_v3_5y_shared_cohort_20260917"
HORIZONS = {3: 90, 6: 180, 12: 365, 24: 730, 36: 1095, 48: 1460, 60: 1825}
REQUIRED = (
    "code_n", "AGE", "GENDER", "D_DIAG", "BSL_HCY", "BSL_B12", "BSL_B9",
    "BSL_CYSC", "H_SMK", "H_DRINK", "H_HYPT", "H_DIAB", "H_STROKE",
    "ONSET_D", "I_BLDSAMP_DT", "Y5_IS", "Y5_IS_DD",
)
OPTIONAL = (
    "A_NIHSS", "IMG_C_TOAST", "y1_is", "y1_is_dd",
    "y2_IS", "y2_IS_dd", "y3_IS", "y3_IS_dd", "y4_IS", "y4_IS_dd",
    "D_DEATH", "D_DEATH_D", "F3_DATE", "F6_DATE", "F12_DATE",
    "F3_DEATH", "F3_DEATH_D", "F6_DEATH", "F6_DEATH_D", "F12_DEATH", "F12_DEATH_D",
)
SOURCES = {k: v for k, v in FIELDS.items() if k in set(REQUIRED + OPTIONAL)}


def field_audit(clinical: pd.DataFrame) -> pd.DataFrame:
    absent = set(clinical.attrs.get("absent_source_columns", []))
    rows = []
    for source, (name, _, _) in SOURCES.items():
        rows.append({"source": source, "canonical": name, "required": source in REQUIRED,
                     "present": source not in absent,
                     "observed": int(clinical[name].notna().sum()), "patients": len(clinical)})
    return pd.DataFrame(rows)


def build_cohort(data: pd.DataFrame):
    d = data.copy()
    event, stop = values(d, "y5_is_event"), values(d, "y5_is_day")
    d["entry"] = values(d, "sample_day")
    d["exit"] = stop.clip(upper=1825)
    d["event_type"] = (event.eq(1) & stop.le(1825)).astype(int)
    death = values(d, "death_day")
    check = d.copy()
    # When supplied, year-one DD also gives actual non-event observation time.
    # Fall back only to the existing year-one contact record for consistency checks.
    y1_time = values(check, "is_day")
    check["last_contact_day"] = y1_time.where(np.isfinite(y1_time) & y1_time.ge(0),
                                             values(check, "last_contact_day"))
    conflict = consistency_issues(check, 5)
    rules = [
        ("adult_ischemic_stroke", d.age.ge(18) & d.diagnosis.eq(1)),
        ("available_wmh_icv", d.image_valid.astype("boolean").fillna(False)),
        ("observed_positive_baseline_hcy", np.isfinite(d.hcy) & d.hcy.gt(0)),
        ("valid_baseline_sample_day", np.isfinite(d.entry) & d.entry.ge(0) & d.entry.lt(1825)),
        ("known_five_year_event_status", event.isin([0, 1])),
        ("valid_actual_event_or_censor_day", np.isfinite(stop) & stop.ge(0)),
        ("consistent_cumulative_first_event", ~conflict),
        ("not_observed_after_known_death", death.isna() | stop.le(death)),
        ("no_recurrence_before_or_on_entry", ~(event.eq(1) & stop.le(d.entry))),
        ("observed_followup_after_entry", d.exit.gt(d.entry)),
    ]
    cohort, exclusions, flow = apply_rules(d, rules)
    audit = {
        "clinical_image_records": len(d), "five_year_cohort": len(cohort),
        "cross_year_conflicts": int(conflict.sum()),
        "event_after_1825": int((event.eq(1) & stop.gt(1825)).sum()),
        "non_event_early_censor": int((event.eq(0) & stop.between(0, 1825, inclusive="left")).sum()),
        "missing_or_invalid_event_state": int((~event.isin([0, 1])).sum()),
        "invalid_or_missing_time": int((~np.isfinite(stop) | stop.lt(0)).sum()),
        "known_death_conflicts": int((death.notna() & stop.gt(death)).sum()),
        "endpoint_source": "Y5_IS / Y5_IS_DD; no inferred long-term death times",
        "censor_contract": "Supplied DD stops observation at first IS, death, last observation or horizon",
    }
    return cohort.reset_index(drop=True), exclusions, pd.DataFrame(flow), audit


def truncate(data: pd.DataFrame, month: int) -> pd.DataFrame:
    cutoff = HORIZONS[month]
    d = data.loc[data.entry.lt(cutoff)].copy()
    d["event_type"] = (d.event_type.eq(1) & d.exit.le(cutoff)).astype(int)
    d["exit"] = d.exit.clip(upper=cutoff)
    if d.exit.le(d.entry).any():
        raise DataError("Invalid (entry, exit] interval in the shared baseline cohort")
    return d


def counts(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for month, cutoff in HORIZONS.items():
        d = truncate(data, month)
        rows.append({"month": month, "cutoff_days": cutoff, "baseline_n": len(data),
                     "n": len(d), "entry_at_or_after_cutoff": len(data)-len(d),
                     "events": int(d.event_type.eq(1).sum()),
                     "censored_before_cutoff": int((d.event_type.eq(0) & d.exit.lt(cutoff)).sum())})
    return pd.DataFrame(rows)


def horizon_auxiliaries(data: pd.DataFrame) -> pd.DataFrame:
    """Entry-adjusted Nelson-Aalen auxiliaries, not fabricated event dates."""
    columns = {}
    for month in HORIZONS:
        d = truncate(data, month)
        eligible = data.index.isin(d.index)
        event = pd.Series(0.0, index=data.index)
        hazard = event.copy()
        event.loc[d.index] = d.event_type.eq(1).astype(float)
        hazard.loc[d.index] = nelson_aalen_increment(d, 1)
        columns[f"_eligible_m{month}"] = eligible.astype(float)
        columns[f"_event_m{month}"] = event
        columns[f"_hazard_m{month}"] = hazard
    return pd.DataFrame(columns, index=data.index)


def validate_covariates(data: pd.DataFrame) -> None:
    for name in adjustment_columns():
        x = values(data, name)
        if not np.isfinite(x).any():
            raise DataError(f"{name}: entirely unavailable in the five-year cohort; cannot impute")
        if (x.notna() & ~np.isfinite(x)).any():
            raise DataError(f"{name}: infinite covariate values")
        if name in {"b12", "folate", "cysc"} and x.le(0).any():
            raise DataError(f"{name}: nonpositive observed values need source review")


def split_at_90(data: pd.DataFrame) -> pd.DataFrame:
    """Nonoverlapping (entry, exit] rows; an event at day 90 belongs to early."""
    early = data.loc[data.entry.lt(90)].copy()
    early["event_type"] = (early.event_type.eq(1) & early.exit.le(90)).astype(int)
    early["exit"] = early.exit.clip(upper=90)
    early["late"] = 0.0
    late = data.loc[data.exit.gt(90)].copy()
    late["entry"] = late.entry.clip(lower=90)
    late["late"] = 1.0
    return pd.concat([early, late], ignore_index=True)
