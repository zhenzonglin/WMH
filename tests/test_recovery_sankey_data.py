"""Synthetic patient paths verify descriptive aggregation and strict source audits."""
from __future__ import annotations

import json
from itertools import pairwise

import numpy as np
import pandas as pd
import pytest

from wmh_hcy.recovery_publication.inputs import source_manifest
from wmh_hcy.recovery_publication.sankey_data import (
    MONTHS,
    MRS_STATES,
    SOURCE_FILES,
    aggregate_mrs,
    load_sankey,
)


def fixture_cohort() -> tuple[pd.DataFrame, dict]:
    paths = np.array([
        [0, 1, 3, 0, 0, 2],
        [1, 6, np.nan, np.nan, np.nan, np.nan],
        [2, np.nan, 2, np.nan, 2, np.nan],
        [0, 2, 4, 5, 6, np.nan],
        [1, 0, 2, 3, 4, 5],
        [2, np.nan, np.nan, 6, np.nan, np.nan],
    ], dtype=float)
    frame = pd.DataFrame({"patient_id": [f"{i:04d}" for i in range(6)]})
    death = np.zeros(6, dtype=bool)
    counts = {}
    for index, month in enumerate(MONTHS):
        score = paths[:, index]
        death |= score == 6
        state = np.where(np.isnan(score), np.nan, np.where(score <= 2, 0, 1))
        state[death] = 2
        frame[f"mrs{month}"] = score
        frame[f"state{month}"] = state
        frame[f"state{month}_conflict"] = False
        counts[str(month)] = {str(k): int((state == k).sum()) for k in (0, 1, 2)} | {
            "unknown": int(np.isnan(state).sum())}
    return frame, {"base_eligible_n": 6, "state_counts": counts,
                   "future_recurrence_fields_requested": []}


def fixture_run(tmp_path):
    run = tmp_path / "synthetic_trajectory" / "runs" / "20261015T120001000000Z"
    run.mkdir(parents=True)
    frame, audit = fixture_cohort()
    frame.to_csv(run / "base_eligible.csv", index=False)
    (run / "audit.json").write_text(json.dumps(audit), encoding="utf-8")
    (run / "status.json").write_text(json.dumps({
        "study": "recovery_trajectory", "run": run.name, "path": str(run),
        "mode": "synthetic", "status": "COMPLETED", "analysis_scope": "trajectory_only",
    }), encoding="utf-8")
    return run


def test_full_cohort_and_every_node_and_link_margin_are_preserved():
    frame, audit = fixture_cohort()
    original = frame.copy(deep=True)
    graph = aggregate_mrs(frame, audit)
    pd.testing.assert_frame_equal(frame, original)
    assert graph["n"] == 6
    assert graph["months"] == (3, 12, 24, 36, 48, 60)
    nodes, links = graph["nodes"], graph["links"]
    assert len(nodes) == 48 and len(links) == 320
    assert nodes.groupby("month").n.sum().eq(6).all()
    assert np.allclose(nodes.groupby("month").proportion.sum(), 1)
    for start, end in pairwise(MONTHS):
        interval = links.loc[links.from_month.eq(start)]
        assert interval.to_month.eq(end).all()
        assert interval.n.sum() == 6
        assert interval.proportion.sum() == pytest.approx(1)
        for column, month in (("source", start), ("target", end)):
            observed = interval.groupby(column).n.sum().reindex(MRS_STATES)
            expected = nodes.loc[nodes.month.eq(month)].set_index("state").n.reindex(MRS_STATES)
            assert np.array_equal(observed.to_numpy(), expected.to_numpy())
    assert "patient_id" not in nodes and "patient_id" not in links
    assert not any(isinstance(value, pd.DataFrame) and len(value) == 6 for value in graph.values())


def test_missing_visit_and_return_are_preserved_without_living_forward_fill():
    frame, audit = fixture_cohort()
    graph = aggregate_mrs(frame, audit)
    links = graph["links"].set_index(["from_month", "to_month", "source", "target"])
    assert links.loc[(3, 12, "2", "missing"), "n"] == 2
    assert links.loc[(12, 24, "missing", "2"), "n"] == 1
    assert links.loc[(12, 24, "missing", "missing"), "n"] == 1
    assert links.loc[(24, 36, "3", "0"), "n"] == 1  # Observed recovery is retained.
    assert graph["validation"]["living_scores_carried_forward"] == 0
    assert graph["validation"]["missing_n_by_month"] == {
        "3": 0, "12": 2, "24": 1, "36": 1, "48": 0, "60": 1}


def test_saved_death_with_missing_raw_score_is_mrs6_not_missing():
    frame, audit = fixture_cohort()
    graph = aggregate_mrs(frame, audit)
    assert graph["validation"]["death_state_to_mrs6_n"] == {
        "3": 0, "12": 0, "24": 1, "36": 1, "48": 2, "60": 3}
    links = graph["links"]
    assert links.loc[links.source.eq("6") & links.target.ne("6"), "n"].eq(0).all()
    assert graph["validation"]["death_absorbing"] is True


def test_saved_death_evidence_before_display_visit_does_not_require_raw_mrs6():
    frame, audit = fixture_cohort()
    frame.loc[1, "mrs12"] = np.nan
    graph = aggregate_mrs(frame, audit)
    assert graph["validation"]["death_state_to_mrs6_n"]["12"] == 1
    assert graph["nodes"].query("month == 12 and state == '6'").n.item() == 1


@pytest.mark.parametrize("score,state", [(0., 0.), (np.nan, np.nan)])
def test_alive_or_unknown_after_saved_death_is_rejected(score, state):
    frame, audit = fixture_cohort()
    frame.loc[1, ["mrs24", "state24"]] = [score, state]
    with pytest.raises(ValueError, match="Death absorption"):
        aggregate_mrs(frame, audit)


@pytest.mark.parametrize("score", [-1, 7, 2.5, np.inf, "invalid"])
def test_invalid_mrs_values_are_not_replaced_with_missing(score):
    frame, audit = fixture_cohort()
    frame["mrs12"] = frame.mrs12.astype(object)
    frame.loc[0, "mrs12"] = score
    with pytest.raises(ValueError, match="mRS|nonnumeric|nonfinite"):
        aggregate_mrs(frame, audit)


@pytest.mark.parametrize("score,state", [(6, 0), (3, 0), (1, 1), (np.nan, 0), (2, np.nan), (2, 2)])
def test_raw_score_and_saved_state_must_agree(score, state):
    frame, audit = fixture_cohort()
    frame.loc[0, ["mrs12", "state12"]] = [score, state]
    with pytest.raises(ValueError, match="mRS/functional-state mismatch"):
        aggregate_mrs(frame, audit)


@pytest.mark.parametrize("flag", [True, "true", "1", "unexpected", None])
def test_saved_conflicts_or_invalid_flags_are_never_silently_dropped(flag):
    frame, audit = fixture_cohort()
    frame["state12_conflict"] = frame.state12_conflict.astype(object)
    frame.loc[0, "state12_conflict"] = flag
    with pytest.raises(ValueError, match="contradiction|conflict flag"):
        aggregate_mrs(frame, audit)


def test_initial_cohort_requires_independence():
    frame, audit = fixture_cohort()
    frame.loc[0, ["mrs3", "state3"]] = [3, 1]
    with pytest.raises(ValueError, match="independent at month 3"):
        aggregate_mrs(frame, audit)


@pytest.mark.parametrize("kind", ["total", "state", "missing_month", "missing_state"])
def test_exact_saved_audit_counts_are_required(kind):
    frame, audit = fixture_cohort()
    if kind == "total":
        audit["base_eligible_n"] += 1
    elif kind == "state":
        audit["state_counts"]["60"]["unknown"] += 1
    elif kind == "missing_month":
        del audit["state_counts"]["60"]
    else:
        del audit["state_counts"]["60"]["2"]
    with pytest.raises(ValueError, match="cohort size|audit"):
        aggregate_mrs(frame, audit)


@pytest.mark.parametrize("identifier", ["0000", "", None])
def test_duplicate_or_missing_identifiers_fail(identifier):
    frame, audit = fixture_cohort()
    frame.loc[1, "patient_id"] = identifier
    with pytest.raises(ValueError, match="patient identifiers"):
        aggregate_mrs(frame, audit)


def test_reader_preserves_leading_zero_ids_and_never_changes_source(tmp_path):
    run = fixture_run(tmp_path)
    frame = pd.read_csv(run / "base_eligible.csv", dtype={"patient_id": str})
    frame.loc[0:1, "patient_id"] = ["001", "1"]
    frame["unrelated_patient_secret"] = "not read or exported"
    frame.to_csv(run / "base_eligible.csv", index=False)
    before = source_manifest(run, list(SOURCE_FILES))
    graph = load_sankey(run, allow_demo=True)
    assert graph["n"] == 6
    assert source_manifest(run, list(SOURCE_FILES)) == before
    assert all("patient_id" not in value for value in graph.values() if isinstance(value, pd.DataFrame))
    assert "not read or exported" not in str(graph)


@pytest.mark.parametrize("key,value", [("study", "recovery"), ("run", "wrong"),
                                     ("path", "wrong"), ("status", "FAILED"),
                                     ("analysis_scope", "recurrence")])
def test_reader_requires_exact_trajectory_run_identity(tmp_path, key, value):
    run = fixture_run(tmp_path)
    status = json.loads((run / "status.json").read_text(encoding="utf-8"))
    status[key] = value
    (run / "status.json").write_text(json.dumps(status), encoding="utf-8")
    with pytest.raises(ValueError):
        load_sankey(run, allow_demo=True)


def test_demo_needs_explicit_opt_in_and_cannot_consume_real_data(tmp_path):
    run = fixture_run(tmp_path)
    with pytest.raises(ValueError, match="mode=real"):
        load_sankey(run)
    status = json.loads((run / "status.json").read_text(encoding="utf-8"))
    status["mode"] = "real"
    (run / "status.json").write_text(json.dumps(status), encoding="utf-8")
    with pytest.raises(ValueError, match="cannot consume a real run"):
        load_sankey(run, allow_demo=True)
    assert load_sankey(run)["n"] == 6


def test_missing_required_columns_fail_before_plotting(tmp_path):
    run = fixture_run(tmp_path)
    frame = pd.read_csv(run / "base_eligible.csv").drop(columns="mrs60")
    frame.to_csv(run / "base_eligible.csv", index=False)
    with pytest.raises(ValueError, match="mrs60"):
        load_sankey(run, allow_demo=True)
