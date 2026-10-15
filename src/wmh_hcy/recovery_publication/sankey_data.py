"""Aggregate saved, observed mRS paths without refitting or exporting patients.

Missing visits remain a distinct node. The only derived mRS value is 6 when
the saved, audited functional state establishes death and the raw score is
missing. This retains the source workflow's absorbing death state, including
death evidence at visits that are not displayed in the Sankey figure.
"""
from __future__ import annotations

from itertools import pairwise, product
from pathlib import Path

import numpy as np
import pandas as pd

from .inputs import _count, _csv, _identity, _numeric, _object

MONTHS = (3, 12, 24, 36, 48, 60)
MRS_STATES = (*map(str, range(7)), "missing")
SOURCE_FILES = ("status.json", "audit.json", "base_eligible.csv")
REQUIRED_COLUMNS = ("patient_id", *(f"{prefix}{month}{suffix}" for month in MONTHS
                    for prefix, suffix in (("mrs", ""), ("state", ""), ("state", "_conflict"))))


def _observed_scores(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Validate score/state agreement and return a temporary anonymous matrix."""
    missing_columns = set(REQUIRED_COLUMNS) - set(frame)
    if missing_columns:
        raise ValueError(f"base_eligible.csv lacks required columns: {', '.join(sorted(missing_columns))}")
    ids = frame.patient_id.astype("string")
    if ids.isna().any() or ids.str.strip().eq("").any():
        raise ValueError("Sankey cohort contains missing patient identifiers")
    if ids.duplicated().any():
        raise ValueError("Sankey cohort contains duplicate patient identifiers")
    values = frame.loc[:, [f"{prefix}{month}" for month in MONTHS for prefix in ("mrs", "state")]].copy()
    _numeric(values, tuple(values.columns), "base_eligible.csv", missing=True)
    matrix = pd.DataFrame(index=frame.index)
    derived_death = {}
    previously_dead = pd.Series(False, index=frame.index)
    for month in MONTHS:
        score, state = values[f"mrs{month}"], values[f"state{month}"]
        if not score.dropna().isin(range(7)).all():
            raise ValueError(f"Invalid mRS score at month {month}; expected integers 0 through 6")
        if not state.dropna().isin((0, 1, 2)).all():
            raise ValueError(f"Invalid saved functional state at month {month}")
        conflict = frame[f"state{month}_conflict"].astype("string").str.strip().str.lower()
        if conflict.isna().any() or not conflict.isin(("true", "false", "1", "0")).all():
            raise ValueError(f"Invalid saved conflict flag at month {month}")
        if conflict.isin(("true", "1")).any():
            raise ValueError(f"Saved functional-state contradiction at month {month}")
        compatible = (
            (state.eq(0) & score.isin((0, 1, 2)))
            | (state.eq(1) & score.isin((3, 4, 5)))
            | (state.eq(2) & (score.eq(6) | score.isna()))
            | (state.isna() & score.isna())
        )
        if not compatible.all():
            raise ValueError(f"Saved mRS/functional-state mismatch at month {month}")
        if (previously_dead & ~state.eq(2)).any():
            raise ValueError(f"Death absorption violated at month {month}")
        previously_dead |= state.eq(2)
        derived = state.eq(2) & score.isna()
        derived_death[str(month)] = int(derived.sum())
        display = score.mask(derived, 6)
        matrix[month] = display.astype("Int64").astype("string").fillna("missing")
    if not matrix[3].isin(("0", "1", "2")).all():
        raise ValueError("Sankey cohort must include only patients independent at month 3 (mRS 0-2)")
    return matrix, derived_death


def aggregate_mrs(frame: pd.DataFrame, audit: dict) -> dict:
    """Return full-cohort counts and adjacent-visit links; never patient records.

    All eight states, including missing, have explicit rows at every visit.
    The link table retains all 64 state pairs per interval for transparent
    zero-count auditing. Every denominator is the complete starting cohort.
    """
    n = len(frame)
    if not n or _count(audit.get("base_eligible_n"), "base_eligible_n") != n:
        raise ValueError("Sankey cohort size differs from base_eligible_n, or cohort is empty")
    matrix, derived_death = _observed_scores(frame)
    saved_counts = audit.get("state_counts")
    if not isinstance(saved_counts, dict):
        raise ValueError("Sankey audit lacks state_counts")  # noqa: TRY004 - saved-artifact validation API
    node_rows, link_rows = [], []
    for month in MONTHS:
        counts = matrix[month].value_counts().reindex(MRS_STATES, fill_value=0)
        collapsed = {"0": int(counts.loc[["0", "1", "2"]].sum()),
                     "1": int(counts.loc[["3", "4", "5"]].sum()),
                     "2": int(counts["6"]), "unknown": int(counts["missing"])}
        saved = saved_counts.get(str(month))
        if not isinstance(saved, dict) or set(saved) != set(collapsed):
            raise ValueError(f"Sankey audit has incomplete state counts at month {month}")
        if any(_count(saved[key], f"month{month}.{key}") != count for key, count in collapsed.items()):
            raise ValueError(f"Sankey cohort/state-count audit mismatch at month {month}")
        node_rows.extend({"month": month, "state": state, "n": int(counts[state]),
                          "proportion": float(counts[state] / n)} for state in MRS_STATES)
    for start, end in pairwise(MONTHS):
        counts = matrix.groupby([start, end], dropna=False).size()
        link_rows.extend({"from_month": start, "to_month": end, "source": source, "target": target,
                          "n": int(counts.get((source, target), 0)),
                          "proportion": float(counts.get((source, target), 0) / n)}
                         for source, target in product(MRS_STATES, repeat=2))
    nodes, links = pd.DataFrame(node_rows), pd.DataFrame(link_rows)
    # Verify both interval totals and each node's incoming/outgoing margins.
    for start, end in pairwise(MONTHS):
        interval = links.loc[links.from_month.eq(start) & links.to_month.eq(end)]
        if int(interval.n.sum()) != n:
            raise ValueError("Sankey interval does not preserve the complete cohort")
        for column, month in (("source", start), ("target", end)):
            totals = interval.groupby(column).n.sum().reindex(MRS_STATES)
            expected = nodes.loc[nodes.month.eq(month)].set_index("state").n.reindex(MRS_STATES)
            if not np.array_equal(totals.to_numpy(), expected.to_numpy()):
                raise ValueError("Sankey node/link margin mismatch")
    return {"n": n, "months": MONTHS, "nodes": nodes, "links": links,
            "validation": {"death_absorbing": True, "all_visits_include_full_cohort": True,
                           "all_intervals_include_full_cohort": True,
                           "living_scores_carried_forward": 0,
                           "death_state_to_mrs6_n": derived_death,
                           "missing_n_by_month": {str(month): int(matrix[month].eq("missing").sum())
                                                  for month in MONTHS}}}


def load_sankey(run: Path, allow_demo: bool = False) -> dict:
    """Read only three named artifacts from one exact saved trajectory run."""
    run, status = _identity(run, "recovery_trajectory", allow_demo)
    if status.get("analysis_scope") != "trajectory_only":
        raise ValueError("Sankey input must be an isolated trajectory run")
    audit = _object(run / "audit.json")
    if audit.get("future_recurrence_fields_requested") != []:
        raise ValueError("Sankey audit must confirm no future recurrence fields were requested")
    frame = _csv(run / "base_eligible.csv", REQUIRED_COLUMNS)
    return {"status": status, "audit": audit, **aggregate_mrs(frame, audit)}
