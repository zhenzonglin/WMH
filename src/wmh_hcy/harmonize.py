"""Explicit units/codes/dates. No inferred medication exposure or endpoint reconstruction."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .common import DataError, dump_json, outdir, read_csv, read_json, resolve, unique_ids
from .fields import FIELDS, LONGTERM_FIELDS, UNKNOWN_98


def sas_date(series: pd.Series, fmt: str, override: str | None = None) -> pd.Series:
    clean = series.astype("string").str.strip().replace({"": pd.NA, ".": pd.NA})
    clean = clean.mask(clean.str.match(r"^\.?[A-Z_]$", na=False), pd.NA)
    nums = pd.to_numeric(clean, errors="coerce")
    numeric = nums.notna()
    result = pd.Series(pd.NaT, index=series.index, dtype="datetime64[ns]")
    if numeric.any():
        kind = override
        f = fmt.upper()
        if kind is None:
            if any(x in f for x in ["DATETIME", "E8601DT", "B8601DT"]):
                kind = "datetime"
            elif any(x in f for x in ["DATE", "YYMMDD", "MMDDYY", "DDMMYY", "E8601DA"]):
                kind = "date"
        if kind not in {"date", "datetime"}:
            raise DataError(f"Unformatted numeric date in {series.name}; set sas.date_formats explicitly")
        result.loc[numeric] = pd.to_datetime(nums.loc[numeric], unit="D" if kind == "date" else "s",
                                             origin="1960-01-01", errors="coerce")
    strings = ~numeric & clean.notna()
    if strings.any():
        result.loc[strings] = pd.to_datetime(clean.loc[strings], errors="coerce", format="mixed")
        if result.loc[strings].isna().any():
            raise DataError(f"Unparseable date text in {series.name}")
    if result.loc[numeric].isna().any():
        raise DataError(f"Out-of-range SAS dates in {series.name}")
    return result


def harmonize(cfg: dict, require_one_year: bool = True, *, source_fields=None,
              required_sources=None) -> pd.DataFrame:
    supplied = cfg["inputs"].get("clinical_csv")
    raw_path = resolve(cfg, supplied) if supplied else outdir(cfg) / "extracted/clinical_raw.csv"
    if not raw_path.is_file():
        raise DataError("Clinical input missing. Run extract or explicitly set inputs.clinical_csv")
    raw = read_csv(raw_path)
    spelling = {c.lower(): c for c in raw}
    raw = raw.rename(columns={spelling[c.lower()]: c for c in FIELDS if c.lower() in spelling})
    unique_ids(raw, "code_n", str(raw_path))
    fp = cfg["inputs"].get("clinical_formats_json")
    formats = read_json(resolve(cfg, fp)) if fp else read_json(outdir(cfg) / "extracted/formats.json")
    overrides = cfg["sas"].get("date_formats", {})
    result = pd.DataFrame(index=raw.index)
    issues = []
    absent = []
    for source, (name, kind, codes) in FIELDS.items():
        if source not in raw or (source_fields is not None and source not in source_fields):
            result[name] = pd.NaT if kind in {"date", "datetime"} else np.nan
            absent.append(source)
            continue
        if kind == "id":
            result[name] = raw[source].astype(str)
        elif kind in {"date", "datetime"}:
            result[name] = sas_date(raw[source], formats.get(source, ""), overrides.get(source))
        else:
            x = pd.to_numeric(raw[source], errors="coerce")
            if source in UNKNOWN_98:
                x = x.mask(x == 98)
            if codes:
                invalid = x.notna() & ~x.isin(codes)
                if invalid.any():
                    if source in LONGTERM_FIELDS:
                        issues.append({"source": source, "invalid_optional_values": int(invalid.sum())})
                    else:
                        raise DataError(f"{source}: invalid codes {sorted(x[invalid].unique())}; no guessing")
            bad_text = raw[source].ne("") & x.isna()
            permitted_missing = raw[source].str.match(r"^\.?[A-Z_]$|^\.$|^$", na=False)
            if source in UNKNOWN_98:
                permitted_missing |= raw[source].isin(["98", "98.0"])
            if (bad_text & ~permitted_missing).any():
                raise DataError(f"{source}: unparsed nonnumeric values; explicit correction required")
            if bad_text.any():
                issues.append({"source": source, "missing_or_unparsed_count": int(bad_text.sum())})
            result[name] = x.astype(float)
    required = ["code_n", "AGE", "D_DIAG", "BSL_HCY", "y1_is", "y1_is_dd",
                "ONSET_D", "I_BLDSAMP_DT"]
    if not require_one_year:
        required = [c for c in required if c not in {"y1_is", "y1_is_dd"}]
    if required_sources is not None:
        required = list(required_sources)
    if set(required) & set(absent):
        raise DataError(f"Missing required columns: {sorted(set(required) & set(absent))}")
    for c, lower, upper in [("nihss", 0, 42), ("mrs12", 0, 5), ("pre_mrs", 0, 5)]:
        invalid = result[c].notna() & (~result[c].between(lower, upper) | (result[c] % 1 != 0))
        if invalid.any():
            raise DataError(f"{c}: score outside supplied dictionary coding")
    # Long-term scores are validated per year, so a new optional field cannot
    # remove cases from the existing one-year analysis.
    # Day precision is the common denominator; same-day order remains unknown.
    onset = result.onset_date.dt.normalize()
    result["sample_day"] = (result.sample_date.dt.normalize() - onset).dt.days.astype(float)
    result["sample3_day"] = (result.sample3_date.dt.normalize() - onset).dt.days.astype(float)
    contacts = pd.concat([(result[f"visit{v}_date"].dt.normalize() - onset).dt.days.where(
                          result[f"death{v}"].isin([1, 2]))
                          for v in [3, 6, 12]], axis=1)
    result["last_contact_day"] = contacts.max(axis=1)
    death_days = []
    death_flags = []
    for flag, date in [("hospital_death", "hospital_death_date")] + [
        (f"death{v}", f"death{v}_date") for v in [3, 6, 12]
    ]:
        confirmed = result[flag].eq(1)
        death_flags.append(confirmed)
        dd = (result[date].dt.normalize() - onset).dt.days.astype(float).where(confirmed)
        death_days.append(dd)
    result["death_confirmed"] = pd.concat(death_flags, axis=1).any(axis=1)
    result["death_day"] = pd.concat(death_days, axis=1).min(axis=1)
    result["death_date_missing"] = result.death_confirmed & result.death_day.isna()
    result.attrs["absent_source_columns"] = absent
    out = outdir(cfg)
    dump_json(out / "prepared/harmonization.json", {"absent_source_columns": absent, "issues": issues,
              "endpoint": ("Y5_IS + Y5_IS_DD only; annual records checked for consistency"
                           if required_sources is not None and "Y5_IS" in required_sources
                           else "y1_is + y1_is_dd only" if require_one_year
                           else "supplied year2-year5 IS + IS_DD"),
              "date_resolution": "calendar days"})
    (out / "prepared").mkdir(parents=True, exist_ok=True)
    result.to_csv(out / "prepared/clinical.csv", index=False)
    return result
