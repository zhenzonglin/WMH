"""Private study-specific data preparation with explicit missingness and time rules."""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..common import DataError, outdir, read_csv, read_json, resolve, unique_ids
from ..harmonize import sas_date
from ..longterm import apply_rules, consistency_issues
from .registry import COVARIATES, IMAGE_COLUMNS, SOURCES, UNITS, UNKNOWN


def number(data, name):
    return pd.to_numeric(data.get(name, pd.Series(np.nan, index=data.index)), errors="coerce")


def read_clinical(cfg):
    supplied = cfg["inputs"].get("clinical_csv")
    path = resolve(cfg, supplied) if supplied else outdir(cfg) / "extracted/clinical_raw.csv"
    if not path.is_file():
        raise DataError("Clinical input missing; supply the SAS directory or an explicit whitelist CSV")
    raw = read_csv(path)
    spellings = {c.casefold(): c for c in raw}
    if len(spellings) != len(raw.columns):
        raise DataError("Case-insensitive duplicate source columns")
    raw = raw.rename(columns={spellings[s.casefold()]: s for s in SOURCES if s.casefold() in spellings})
    unique_ids(raw, "code_n", str(path))
    fmt_path = cfg["inputs"].get("clinical_formats_json")
    formats = read_json(resolve(cfg, fmt_path)) if fmt_path else read_json(outdir(cfg) / "extracted/formats.json")
    overrides = cfg.get("sas", {}).get("date_formats", {})
    columns, audits = {}, []
    for source, (name, kind, codes) in SOURCES.items():
        present = source in raw
        bad = 0
        issue = ""
        if not present:
            value = pd.Series(pd.NaT if kind in {"date", "datetime"} else np.nan, index=raw.index)
        elif kind == "id":
            value = raw[source].astype(str)
        elif kind in {"date", "datetime"}:
            try:
                value = sas_date(raw[source], formats.get(source, ""), overrides.get(source))
            except DataError as exc:
                # Isolate a date error to studies needing this field; never guess its unit.
                value = pd.Series(pd.NaT, index=raw.index)
                bad = int(raw[source].ne("").sum())
                issue = str(exc)
        else:
            text = raw[source].astype("string").str.strip()
            missing = text.str.match(r"^\.?[A-Z_]$|^\.$|^$", na=False)
            value = pd.to_numeric(text.mask(missing), errors="coerce").astype(float)
            known_unknown = value.isin(UNKNOWN.get(source, []))
            invalid = (~missing & value.isna()) | (value.notna() & ~np.isfinite(value))
            if codes:
                invalid |= value.notna() & ~value.isin(codes) & ~known_unknown
            if name in {"nihss", "pre_mrs", "mrs3", "discharge_mrs"}:
                upper = 42 if name == "nihss" else 5
                invalid |= value.notna() & (~value.between(0, upper) | value.mod(1).ne(0))
            if name.startswith("mrs"):
                invalid |= value.notna() & (~value.between(0, 6) | value.mod(1).ne(0))
            # Positivity is a measurement contract, not percentile trimming.
            if name in {"cysc", "cysc3", "tg", "hdl", "ldl", "apo_ai", "bmi"}:
                invalid |= value.notna() & value.le(0)
            if name in {"cer16", "cer20", "cer24", "cer241", "uacr0", "uacr3", "cec"}:
                invalid |= value.notna() & value.lt(0)
            if "sbp" in name or "dbp" in name:
                invalid |= value.notna() & value.le(0)
            bad = int(invalid.sum())
            value = value.mask(invalid | known_unknown)
        columns[name] = value
        audits.append({"source": source, "canonical": name, "present": present,
                       "observed": int(value.notna().sum()), "invalid": bad,
                       "unit": UNITS.get(name, "dictionary"), "patients": len(raw), "issue": issue})
    data = pd.DataFrame(columns)
    # Parent "no medication" answers are structural negatives, not missing-child imputation.
    contradiction = data.prior_lipid_med.eq(1) & data.prior_statin.eq(1)
    data["medication_conflict"] = contradiction
    data["prior_statin"] = data.prior_statin.mask(contradiction)
    data.loc[data.prior_lipid_med.eq(1) & ~contradiction, "prior_statin"] = 0
    for child in ("discharge_acei", "discharge_arb"):
        contradiction = data.discharge_bp_med.eq(1) & data[child].eq(1)
        data["medication_conflict"] |= contradiction
        data[child] = data[child].mask(contradiction)
        data.loc[data.discharge_bp_med.eq(1) & ~contradiction, child] = 0
    data["raas"] = np.where(data.discharge_acei.eq(1) | data.discharge_arb.eq(1), 1,
                            np.where(data.discharge_acei.eq(0) & data.discharge_arb.eq(0), 0, np.nan))
    onset = data.onset_date.dt.normalize()
    for date, day in [("sample_date", "sample_day"), ("visit3_date", "visit3_day"),
                      ("visit12_date", "visit12_day")]:
        data[day] = (data[date].dt.normalize() - onset).dt.days.astype(float)
    death_days = []
    for flag, date in [("hospital_death", "hospital_death_date")] + [
        (f"death{m}", f"death{m}_date") for m in (3, 6, 12)
    ]:
        death_days.append(((data[date].dt.normalize() - onset).dt.days).where(data[flag].eq(1)))
    data["death_day"] = pd.concat(death_days, axis=1).min(axis=1)
    for month in (3, 12):
        data = add_bp(data, month)
    return derive(data), pd.DataFrame(audits)


def add_bp(data, month):
    d = data.copy()
    left, right = number(d, f"lsbp{month}"), number(d, f"rsbp{month}")
    ldbp, rdbp = number(d, f"ldbp{month}"), number(d, f"rdbp{month}")
    # Logical errors remain visible and cannot drive choice of the other arm silently.
    invalid = (left.notna() & ldbp.notna() & left.le(ldbp)) | (
        right.notna() & rdbp.notna() & right.le(rdbp))
    d[f"bp{month}_conflict"] = invalid
    choose_left = left.notna() & (right.isna() | left.ge(right))
    d[f"sbp{month}"] = left.where(choose_left, right).mask(invalid)
    d[f"dbp{month}"] = ldbp.where(choose_left, rdbp).mask(invalid)
    d[f"sbp{month}_mean"] = pd.concat([left, right], axis=1).mean(axis=1).mask(invalid)
    d[f"bp{month}_arm"] = np.where(d[f"sbp{month}"].isna(), "missing", np.where(choose_left, "left", "right"))
    return d


def functional_state(data, month):
    """Known death is absorbing forward only. Never carry an old living score forward."""
    score = number(data, f"mrs{month}")
    dead = number(data, "hospital_death").eq(1)
    bad = pd.Series(False, index=data.index)
    for prior in (3, 6, 12, 24, 36, 48, 60):
        if prior > month:
            continue
        flag, mrs = number(data, f"death{prior}"), number(data, f"mrs{prior}")
        alive = flag.eq(2) | mrs.between(0, 5)
        dead_now = flag.eq(1) | mrs.eq(6)
        bad |= (dead & alive) | (flag.eq(1) & mrs.between(0, 5)) | (flag.eq(2) & mrs.eq(6))
        dead |= dead_now
    state = pd.Series(np.nan, index=data.index)
    state.loc[score.between(0, 2)] = 0
    state.loc[score.between(3, 5)] = 1
    state.loc[dead] = 2
    state.loc[bad] = np.nan
    return state, bad


def derive(data):
    d = data.copy()
    c16, c24, c241 = (number(d, n) for n in ("cer16", "cer24", "cer241"))
    d["cer_ratio"] = np.log2(c16.where(c16.gt(0)) / c24.where(c24.gt(0)))
    d["cer_alt"] = np.log2(c241.where(c241.gt(0)) / c24.where(c24.gt(0)))
    d["log_c16"] = np.log2(c16.where(c16.gt(0)))
    d["log_c24"] = np.log2(c24.where(c24.gt(0)))
    a, b = number(d, "uacr0"), number(d, "uacr3")
    d["albuminuria"] = (a.ge(3).astype(int) + 2*b.ge(3).astype(int)).where(a.notna() & b.notna())
    d["uacr0_log"], d["uacr3_log"] = np.log1p(a), np.log1p(b)
    for month in (3, 12, 24, 36, 48, 60):
        d[f"state{month}"], d[f"state{month}_conflict"] = functional_state(d, month)
    pv, deep = number(d, "fazekas_pv"), number(d, "fazekas_deep")
    d["severe_wmh"] = (pv.eq(3) | deep.eq(3)).astype(float).where(
        pv.eq(3) | deep.eq(3) | (pv.notna() & deep.notna()))
    if "wmh_ml" in d:
        d["log_wmh"] = np.log1p(d.wmh_ml.where(d.wmh_ml.ge(0)))
        d["log_lesion"] = np.log1p(d.lesion_ml.where(d.lesion_ml.ge(0)))
    return d


def volume_ok(data, name):
    v = number(data, name)
    valid = np.isfinite(v) & (v.gt(0) if name in {"icv_ml", "gm119_ml"} else v.ge(0))
    if name != "icv_ml":
        valid &= v.le(number(data, "icv_ml"))
    return valid


def build_cohort(master, study, month=3):
    d = derive(master)
    rules = [("adult_ischemic_stroke", d.age.ge(18) & d.diagnosis.eq(1)),
             ("available_true_icv", volume_ok(d, "icv_ml"))]
    images = ("gm119_ml", "lesion_ml") if study == "cec" else ("wmh_ml",)
    if study == "recovery":
        images += ("gm119_ml", "lesion_ml")
    if study == "ceramide":
        images += ("lesion_ml",)
    rules += [(f"available_{c}", volume_ok(d, c)) for c in images]
    if study in {"recovery", "bp", "kidney"}:
        alive = (number(d, f"death{month}").eq(2) | number(d, f"mrs{month}").between(0, 5))
        rules.append((f"alive_at_month{month}", alive & ~d[f"state{month}_conflict"]
                      & ~number(d, f"death{month}").eq(1)))
    if study == "recovery":
        rules.append(("independent_at_month3", d.mrs3.between(0, 2)))
    elif study == "bp":
        d["entry"] = number(d, f"visit{month}_day")
        # Explicit exposure substitution for the independent month12 landmark.
        d["sbp3"] = number(d, f"sbp{month}")
        if month == 12:
            d["mrs3"] = number(d, "mrs12")
        event, stop = number(d, "y5_is_event"), number(d, "y5_is_day")
        d["exit"] = stop.clip(upper=1825)
        d["event_type"] = (event.eq(1) & stop.le(1825)).astype(int)
        check = d.copy()
        check["last_contact_day"] = number(d, "is_day")
        conflict = consistency_issues(check, 5)
        rules += [
            ("observed_recovery_sbp", np.isfinite(d.sbp3) & d.sbp3.gt(0)),
            ("known_actual_visit_day", np.isfinite(d.entry) & d.entry.gt(0) & d.entry.lt(1825)),
            ("known_five_year_status_and_time", event.isin([0, 1]) & np.isfinite(stop) & stop.ge(0)),
            ("consistent_first_recurrence", ~conflict),
            ("no_recurrence_before_or_on_visit", ~(event.eq(1) & stop.le(d.entry))),
            ("no_observation_after_dated_death", d.death_day.isna() | stop.le(d.death_day)),
            ("followup_after_visit", d.exit.gt(d.entry)),
        ]
    elif study == "ceramide":
        rules += [("observed_positive_cer16_and_cer24", np.isfinite(d.cer_ratio)),
                  ("known_baseline_sample_time", np.isfinite(d.sample_day) & d.sample_day.ge(0))]
    elif study == "cec":
        rules += [("observed_baseline_cec", np.isfinite(d.cec) & d.cec.ge(0)),
                  ("known_baseline_sample_time", np.isfinite(d.sample_day) & d.sample_day.ge(0))]
    elif study == "kidney":
        rules.append(("both_uacr_observed", np.isfinite(d.albuminuria)))
    if study in {"recovery", "kidney"}:
        rules.append(("no_five_year_state_contradiction", ~d.state60_conflict))
    eligible, exclusions, flow = apply_rules(d, rules)
    outcome = "event_type" if study == "bp" else "log_wmh" if study == "ceramide" else (
        "gm119_ml" if study == "cec" else "state60")
    if study in {"recovery", "kidney"}:
        flow.append({"step": "known_five_year_functional_state", "remaining": int(eligible[outcome].notna().sum()),
                     "excluded_here": int(eligible[outcome].isna().sum())})
    return eligible.reset_index(drop=True), exclusions, pd.DataFrame(flow)


def cohort_audit(data, study):
    result = {"eligible_n": len(data), "outcome_observed_n": len(data), "covariates_assessed": bool(len(data)),
              "covariate_missing": {c: int(data[c].isna().sum()) for c in COVARIATES[study]},
              "covariates_entirely_missing": [c for c in COVARIATES[study] if len(data) and data[c].notna().sum() == 0]}
    if study == "bp":
        result.update(events=int(data.event_type.sum()),
                      censored=int(data.event_type.eq(0).sum()),
                      early_censor=int((data.event_type.eq(0) & data.exit.lt(1825)).sum()))
    elif study in {"recovery", "kidney"}:
        result.update(outcome_observed_n=int(data.state60.notna().sum()),
                      independent=int(data.state60.eq(0).sum()), dependent=int(data.state60.eq(1).sum()),
                      dead=int(data.state60.eq(2).sum()), unknown=int(data.state60.isna().sum()))
    for c in IMAGE_COLUMNS:
        result[f"available_{c}"] = int(volume_ok(data, c).sum())
    return result
