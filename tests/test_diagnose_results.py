"""Standalone reporter tests: no analysis dependencies or patient data."""

import csv
import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "diagnose_results", Path(__file__).resolve().parents[1] / "scripts/diagnose_results.py"
)
REPORTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORTER)


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def test_page_one_groups_failures_and_does_not_expose_details(tmp_path):
    folder = tmp_path / "results/run"
    write_json(
        folder / "02_recurrence/risk_diagnostics.json",
        {
            "status": "BOOTSTRAP_UNSTABLE",
            "valid": 150,
            "requested": 200,
            "failures": [
                {"imputation": 0, "reason": "Singular matrix"},
                {"imputation": 0, "reason": "secret patient path /private/123"},
            ],
        },
    )
    write_json(
        tmp_path / "prepared/cohort_summary.json",
        {"main": {"n": 100, "ischemic_events": 10, "deaths_first": 3}},
    )
    status = {"analyses": {"02_recurrence": {"n": 100, "ischemic_events": 10}}}
    text = "\n".join(REPORTER.page_one(folder, tmp_path, status))
    assert "150/200" in text and "singular_matrix:1" in text and "other:1" in text
    assert "100 / 10 / 3" in text
    assert "secret" not in text and "/private" not in text


def test_missing_and_malformed_json_are_explicit(tmp_path):
    assert REPORTER.read_json(tmp_path / "absent.json")["_unavailable"] == "MISSING"
    (tmp_path / "bad.json").write_text("not json", encoding="utf-8")
    assert REPORTER.read_json(tmp_path / "bad.json")["_unavailable"] == "UNREADABLE"
    text = "\n".join(REPORTER.page_two(tmp_path))
    assert "MI: MISSING/UNREADABLE" in text
    assert "MAR or proportional odds" in text


def test_imputation_trace_and_fit_ranges(tmp_path):
    write_json(
        tmp_path / "02_recurrence/imputation.json",
        {
            "m": 2,
            "missing": {"b12": 12},
            "imputed_mean_trace_transformed_scale": [
                {"iteration": 1, "imputation": 0, "b12": 7},
                {"iteration": 2, "imputation": 0, "b12": 8},
                {"iteration": 2, "imputation": 1, "b12": 9},
            ],
        },
    )
    write_json(
        tmp_path / "02_recurrence/diagnostics.json",
        [
            {
                "events": 10,
                "parameters": 5,
                "rank": 5,
                "max_abs_gradient": 0.1,
                "scaled_condition_number": 20,
                "warnings": [],
            },
            {
                "events": 10,
                "parameters": 5,
                "rank": 5,
                "max_abs_gradient": 0.2,
                "scaled_condition_number": 25,
                "warnings": ["warning"],
            },
        ],
    )
    text = "\n".join(REPORTER.page_two(tmp_path))
    assert "first/last iteration mean=7/8.5" in text
    assert "last-iteration MI range=8..9" in text
    assert "condition=20..25" in text and "warning fits=1" in text


def test_current_prepared_count_mismatch_is_visible(tmp_path):
    write_json(tmp_path / "prepared/cohort_summary.json", {"main": {"n": 50}})
    status = {"analyses": {"02_recurrence": {"n": 100}}}
    text = "\n".join(REPORTER.page_one(tmp_path, tmp_path, status))
    assert "prepared count differs" in text


def test_command_reads_only_aggregate_files_and_writes_nothing(tmp_path, monkeypatch, capsys):
    root = tmp_path / "results/test-run"
    write_json(root / "status.json", {"mode": "synthetic", "analyses": {}})
    write_json(tmp_path / "latest_results.json", {"path": str(root)})
    patient = root / "analysis_participants.csv"
    patient.write_text("patient_id\nPRIVATE_ID_123\n", encoding="utf-8")
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    monkeypatch.setattr("sys.argv", ["diagnose_results.py", "--output-dir", str(tmp_path)])
    assert REPORTER.main() == 0
    text = capsys.readouterr().out
    assert "PRIVATE_ID_123" not in text and "mode=synthetic" in text
    assert before == {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}


def test_nonfinite_numbers_do_not_look_valid():
    assert REPORTER.span([None, float("nan"), float("inf")]) == "NA"
    assert REPORTER.span([1, 2, None]) == "1..2"


def test_threshold_diagnostics_and_no_raw_warning_text(tmp_path):
    folder = tmp_path / "04_function"
    write_json(folder / "diagnostics.json", [{"converged": True, "warnings": ["PRIVATE_DETAIL"]}])
    (folder / "threshold_diagnostics.csv").write_text(
        "threshold,term,estimate,failure\n1,H_x_W,0.1,\n2,H_x_W,0.3,\n3,,,PRIVATE_DETAIL\n", encoding="utf-8"
    )
    text = "\n".join(REPORTER.page_two(tmp_path))
    assert "converged=1; warning fits=1" in text
    assert "HW thresholds=1,2; coefficient range=0.1..0.3; failed rows=1" in text
    assert "PRIVATE_DETAIL" not in text


def test_moved_latest_directory_and_page_two(tmp_path, monkeypatch, capsys):
    root = tmp_path / "results/moved-run"
    write_json(root / "status.json", {"mode": "synthetic", "analyses": {}})
    write_json(tmp_path / "latest_results.json", {"path": str(tmp_path / "old/moved-run")})
    monkeypatch.setattr("sys.argv", ["diagnose_results.py", "--output-dir", str(tmp_path), "--page", "2"])
    assert REPORTER.main() == 0
    text = capsys.readouterr().out
    assert "run=moved-run" in text and "[2/2]" in text


def test_longterm_screenshots_use_separate_pointer(tmp_path, monkeypatch, capsys):
    root = tmp_path / "longterm/results/check"
    write_json(root / "status.json", {"mode": "synthetic", "status": "PREPARED", "analyses": {
        "year2/H2": {"n": 120, "events": 20, "status": "PREPARED"}}})
    write_json(tmp_path / "longterm/latest_results.json", {"path": str(root)})
    monkeypatch.setattr("sys.argv", ["diagnose_results.py", "--output-dir", str(tmp_path), "--longterm"])
    assert REPORTER.main() == 0
    text = capsys.readouterr().out
    assert "year2/H2" in text and "120" in text and "SUPPLEMENTARY" in text


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def page_three_fixture(tmp_path):
    root = tmp_path / "results/test-run"
    rows = [
        {"patient_id": "PRIVATE_001", "entry": "3", "exit": "100", "event_type": "1.0"},
        {"patient_id": "PRIVATE_002", "entry": "3", "exit": "120", "event_type": "2"},
        {"patient_id": "PRIVATE_003", "entry": "3", "exit": "365", "event_type": "0"},
    ]
    for row in rows:
        row.update({name: "1.0" for name in REPORTER.CORE_LEVELS})
    rows[1]["sex"] = "2"
    rows[1]["hypertension"] = ""
    rows[2]["smoking"] = "PRIVATE_BAD_CODE"
    for label, name in [("H2", "main"), ("H3", "month3")]:
        write_csv(root / REPORTER.MODELS[label] / "analysis_participants.csv", rows,
                  REPORTER.PARTICIPANT_FIELDS)
        write_csv(tmp_path / f"prepared/cohort_{name}.csv", rows,
                  REPORTER.PARTICIPANT_FIELDS + tuple(REPORTER.CORE_LEVELS))
    status = {"mode": "synthetic", "analyses": {
        REPORTER.MODELS[label]: {"n": 3} for label in ("H2", "H3")}}
    write_json(root / "status.json", status)
    write_json(tmp_path / "latest_results.json", {"path": str(root)})
    return root, rows, status


def test_page_three_counts_privacy_and_readonly_command(tmp_path, monkeypatch, capsys):
    root, _, _ = page_three_fixture(tmp_path)
    write_csv(root / "02_recurrence/death_coefficients.csv",
              [{"estimate": "-2", "se": "0.5"}, {"estimate": "0.2", "se": "1.5"}],
              ("estimate", "se"))
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    monkeypatch.setattr("sys.argv", ["diagnose_results.py", "--output-dir", str(tmp_path), "--page", "3"])
    assert REPORTER.main() == 0
    text = capsys.readouterr().out
    assert "H2: MATCH; N=3; IS=1; death=1; censored=1" in text
    assert "sex: 1:2/1/0  2:1/0/1" in text
    assert "MISS:1/0/1" in text and "INVALID:1/0/0" in text
    assert "max|beta|=2; max SE=1.5" in text
    assert "local patient rows summarized" in text
    assert "PRIVATE" not in text
    assert before == {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}


@pytest.mark.parametrize("change", ["id", "entry", "exit", "event_type", "duplicate", "status_n"])
def test_page_three_refuses_stale_or_duplicate_cohort(tmp_path, change):
    root, rows, status = page_three_fixture(tmp_path)
    if change == "id":
        rows[0]["patient_id"] = "PRIVATE_NEW_ID"
    elif change in ("entry", "exit", "event_type"):
        rows[0][change] = {"entry": "4", "exit": "101", "event_type": "2"}[change]
    elif change == "duplicate":
        rows[0]["patient_id"] = rows[1]["patient_id"]
    else:
        status["analyses"]["02_recurrence"]["n"] = 4
    write_csv(tmp_path / "prepared/cohort_main.csv", rows,
              REPORTER.PARTICIPANT_FIELDS + tuple(REPORTER.CORE_LEVELS))
    text = "\n".join(REPORTER.page_three(root, tmp_path, status))
    h2 = text.split("H2:", 1)[1].split("H3:", 1)[0]
    assert "MATCH;" not in h2 and "sex:" not in h2
    assert "PRIVATE" not in text
    assert "H3: MATCH;" in text


def test_page_three_missing_inputs_are_explicit(tmp_path):
    text = "\n".join(REPORTER.page_three(tmp_path, tmp_path, {}))
    assert text.count("UNAVAILABLE/INVALID inputs") == 2


def test_page_three_preserves_leading_zero_ids():
    rows = [{"patient_id": value, "entry": "1", "exit": "365", "event_type": "0"}
            for value in ("001", "1")]
    assert len(REPORTER.participant_index(rows)) == 2


def test_page_three_rejects_longterm(monkeypatch):
    monkeypatch.setattr("sys.argv", ["diagnose_results.py", "--longterm", "--page", "3"])
    with pytest.raises(SystemExit) as exc:
        REPORTER.main()
    assert exc.value.code == 2


def test_risk_failure_stage_is_visible_without_raw_error(tmp_path):
    write_json(tmp_path / "03_month3_update/risk_failure.json",
               {"stage": "death_fit", "imputation": 2, "reason": "PRIVATE_DETAIL"})
    write_json(tmp_path / "02_recurrence/risk_diagnostics.json",
               {"failures": [{"stage": "death_fit"}, {"stage": "PRIVATE_DETAIL"}]})
    text = "\n".join(REPORTER.page_one(tmp_path, tmp_path, {}))
    assert "failure stage=death_fit; imputation=2" in text
    assert "failure stages: death_fit:1, unrecorded:1" in text
    assert "PRIVATE_DETAIL" not in text
