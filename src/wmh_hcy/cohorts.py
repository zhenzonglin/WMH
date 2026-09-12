"""Calendar-day delayed entry and one frozen ischemic endpoint."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .common import DataError, dump_json, outdir, record_run
from .harmonize import harmonize
from .imaging import read_imaging


def build_cohorts(data: pd.DataFrame, horizon: float = 365) -> tuple[dict, pd.DataFrame, list]:
    d = data.copy()
    d["entry"] = d[["sample_day", "mri_day"]].max(axis=1, skipna=False)
    d["event_day"] = d.is_day.where(d.is_event.eq(1))
    d["censor_day"] = d.last_contact_day.clip(upper=horizon)
    # A confirmed event/death is itself an observed contact, even if visit date is absent.
    d["censor_day"] = pd.concat([d.censor_day, d.event_day.clip(upper=horizon),
                                 d.death_day.clip(upper=horizon)], axis=1).max(axis=1)
    d["exit"] = pd.concat([d.censor_day, d.event_day, d.death_day,
                            pd.Series(horizon, index=d.index)], axis=1).min(axis=1)
    d["event_type"] = np.select([
        d.is_event.eq(1) & d.event_day.eq(d.exit),
        d.death_day.eq(d.exit),
    ], [1, 2], default=0).astype(int)
    # Documented ischemic event on the death date has priority. This does not infer fatal recurrence.
    conditions = [
        ("adult_ischemic_stroke", d.age.ge(18) & d.diagnosis.eq(1)),
        ("valid_wmh_icv_qc", d.image_valid.fillna(False).astype(bool)),
        ("observed_positive_baseline_hcy", np.isfinite(d.hcy) & d.hcy.gt(0)),
        ("known_baseline_measurement_time", d.entry.notna() & d.sample_day.ge(0) & d.mri_day.ge(0)),
        ("known_endpoint_status", d.is_event.isin([0, 1])),
        ("known_event_date_if_event", ~d.is_event.eq(1) | d.event_day.notna()),
        ("event_within_one_year", ~d.is_event.eq(1) | d.event_day.between(0, horizon)),
        ("known_death_date_if_death", ~d.death_date_missing.fillna(False).astype(bool)),
        ("valid_event_chronology", ~(d.event_day.notna() & d.death_day.notna()
                                     & d.event_day.gt(d.death_day))),
        ("no_event_before_or_on_entry", d.event_day.isna() | d.event_day.gt(d.entry)),
        ("alive_at_entry", d.death_day.isna() | d.death_day.gt(d.entry)),
        ("observed_followup_after_entry", d.censor_day.notna() & d.exit.gt(d.entry)
                                        & d.entry.lt(horizon)),
    ]
    keep = pd.Series(True, index=d.index)
    audit = d[["patient_id"]].copy()
    audit["exclusion_reason"] = ""
    flow = [{"step": "clinical_records", "remaining": len(d), "excluded_here": 0}]
    for name, rule in conditions:
        rule = rule.fillna(False)
        excluded = keep & ~rule
        audit.loc[excluded, "exclusion_reason"] = name
        keep &= rule
        flow.append({"step": name, "remaining": int(keep.sum()), "excluded_here": int(excluded.sum())})
    main = d.loc[keep].copy()
    xsec = d.loc[conditions[0][1] & conditions[1][1] & conditions[2][1]].copy()
    m3 = main.loc[main.sample3_day.notna() & np.isfinite(main.hcy3) & main.hcy3.gt(0)
                   & main.sample3_day.ge(main.entry) & main.sample3_day.lt(horizon)
                   & main.sample3_day.lt(main.exit)].copy()
    m3["entry"] = m3.sample3_day
    # Actual sample dates only; no nominal day-90 substitution or carry-forward.
    func = main.copy()
    func.loc[func.death_day.notna() & func.death_day.le(horizon), "mrs12"] = 6
    func = func.loc[func.mrs12.isin(range(7))].copy()
    if "t1_qc" in func:
        t1 = func.loc[func.gm119_ml.notna() & func.gm119_ml.gt(0) & func.t1_qc.eq("pass")].copy()
    else:
        t1 = func.iloc[:0].copy()
    return {"main": main, "cross_sectional": xsec, "month3": m3,
            "functional": func, "functional_t1": t1}, audit, flow


def prepare(cfg: dict) -> dict:
    clinical = harmonize(cfg)
    imaging = read_imaging(cfg)
    d = clinical.merge(imaging, on="patient_id", how="left", validate="one_to_one")
    out = outdir(cfg) / "prepared"
    d.to_csv(out / "analysis_master.csv", index=False)
    cohorts, audit, flow = build_cohorts(d, cfg["analysis"]["horizon"])
    audit.to_csv(out / "exclusions.csv", index=False)
    pd.DataFrame(flow).to_csv(out / "flow.csv", index=False)
    summary = {}
    for name, frame in cohorts.items():
        frame.to_csv(out / f"cohort_{name}.csv", index=False)
        summary[name] = {"n": len(frame), "ischemic_events": int(frame.event_type.eq(1).sum()),
                         "deaths_first": int(frame.event_type.eq(2).sum())}
    dump_json(out / "cohort_summary.json", summary)
    record_run(cfg, "prepare", summary)
    if cohorts["main"].empty:
        raise DataError("Main cohort empty. Inspect prepared/flow.csv and exclusions.csv")
    return summary
