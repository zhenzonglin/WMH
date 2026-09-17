"""Read-only time summaries from an existing run; no identifiers or row output."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .common import DataError, read_json


def timing_summary(root: Path) -> list[str]:
    wanted = {"entry", "exit", "event_type", "sample_day", "onset_date", "sample_date"}
    path = root / "cohort.csv"
    if not path.is_file():
        raise DataError("Prepared recurrence cohort absent in the selected run")
    data = pd.read_csv(path, usecols=lambda name: name in wanted)
    if wanted - set(data):
        raise DataError(f"Prepared timing fields absent: {sorted(wanted-set(data))}")
    for col in ["entry", "exit", "event_type", "sample_day"]:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    entry, stop, event = data.entry, data.exit, data.event_type
    early_censor = event.eq(0) & stop.lt(1825)
    lines = ["RECURRENCE V3 | TIMING CHECK | run=" + root.name,
             "READ ONLY: local date/time columns -> aggregate counts; no IDs, row values, edits or fits.",
             (f"N={len(data)}; IS={int(event.eq(1).sum())}; early_censor={int(early_censor.sum())}; "
              f"event_free_at_day1825={int((event.eq(0) & stop.eq(1825)).sum())}")]
    q = entry.quantile([0, .25, .5, .75, 1]).tolist()
    lines.append("ENTRY DAYS FROM ONSET: min/Q1/median/Q3/max=" + "/".join(f"{x:g}" for x in q))
    lines.append("Entry interval (days)      N    IS_by_5y")
    for low, high, label in [(0, 8, "0-7"), (8, 15, "8-14"), (15, 31, "15-30"),
                             (31, 90, "31-89"), (90, 180, "90-179"),
                             (180, 365, "180-364"), (365, np.inf, ">=365")]:
        mask = entry.ge(low) & entry.lt(high)
        lines.append(f"{label:<23}{int(mask.sum()):>5}{int((mask & event.eq(1)).sum()):>12}")
    invalid_entry = ~np.isfinite(entry) | entry.lt(0) | entry.mod(1).ne(0)
    onset = pd.to_datetime(data.onset_date, errors="coerce", format="mixed")
    sample = pd.to_datetime(data.sample_date, errors="coerce", format="mixed")
    calculated = (sample.dt.normalize()-onset.dt.normalize()).dt.days
    comparable = calculated.notna() & entry.notna()
    mismatch = comparable & calculated.ne(entry)
    lines.append(f"Invalid/noninteger entry days={int(invalid_entry.sum())}; "
                 f"unreadable canonical dates={int(calculated.isna().sum())}")
    lines.append(f"Canonical date subtraction vs entry: checked={int(comparable.sum())}; "
                 f"mismatches={int(mismatch.sum())}; sample_day vs entry mismatches="
                 f"{int((data.sample_day.isna() | data.sample_day.ne(entry)).sum())}")
    lines.append(f"Invalid follow-up intervals={int((~np.isfinite(stop) | stop.le(entry)).sum())}")
    manifest = read_json(root / "inputs_derived/extracted/manifest.json")
    formats = read_json(root / "inputs_derived/extracted/formats.json")
    cfg = read_json(root / "config_snapshot.json")
    overrides = cfg.get("sas", {}).get("date_formats", {})
    for source in ["ONSET_D", "I_BLDSAMP_DT"]:
        owner = manifest.get("variable_owner", {}).get(source)
        # Only metadata and source filenames; never print raw dates or patient IDs.
        basename = str(owner).replace("\\", "/").rsplit("/", 1)[-1] if owner else "CSV/metadata_unavailable"
        lines.append(f"{source}: source={basename}; format={formats.get(source) or 'not_recorded'}; "
                     f"override={overrides.get(source) or 'none'}")
    audit = read_json(root / "endpoint_audit.json")
    if audit:
        lines.append("Pre-exclusion endpoint audit: cross_year_conflicts="
                     f"{audit.get('cross_year_conflicts', 'NA')}; known_death_conflicts="
                     f"{audit.get('known_death_conflicts', 'NA')}")
    lines.append("Late entries are flags only: unchanged eligibility; zero mismatches does not validate source dates.")
    lines.append("Early censoring does not distinguish death, loss to follow-up or shorter observed follow-up.")
    return lines
