"""Read-only aggregate review for one recovery-trajectory run.

The script never refits models and never prints patient identifiers.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def read_json(path: Path) -> dict:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def require_run(path: Path) -> tuple[dict, dict]:
    if not path.is_dir():
        raise ValueError(f"Run directory absent: {path}")
    status = read_json(path / "status.json")
    if status.get("study") != "recovery_trajectory":
        raise ValueError("Expected a recovery_trajectory status.json")
    return status, read_json(path / "audit.json")


def page_status(run: Path, status: dict, audit: dict) -> list[str]:
    lines = [
        "[1/4] RECOVERY TRAJECTORY | STATUS",
        f"run={status.get('run', '')}",
        f"directory={run}",
        f"status={status.get('status', '')}; scope={status.get('analysis_scope', '')}",
        f"error={status.get('error', '') or 'none'}",
        f"base_eligible_n={audit.get('base_eligible_n', 'NA')}; transition_rows_n={audit.get('transition_rows_n', 'NA')}",
        f"future_recurrence_fields_requested={audit.get('future_recurrence_fields_requested', [])}",
    ]
    state_counts = audit.get("state_counts", {})
    if state_counts:
        lines.append("state counts by visit (independent/dependent/dead/unknown):")
        for month, counts in state_counts.items():
            lines.append(f"  month{month}: {counts.get('0', 0)}/{counts.get('1', 0)}/{counts.get('2', 0)}/{counts.get('unknown', 0)}")
    return lines


def page_transitions(run: Path) -> list[str]:
    lines = ["[2/4] TRANSITION AUDIT AND EFFECTS"]
    counts = run / "transition_counts.csv"
    if counts.is_file():
        lines.extend(pd.read_csv(counts).to_string(index=False).splitlines())
    effects = run / "transition_effects.csv"
    if effects.is_file():
        lines.append("\ntransition effects:")
        lines.extend(pd.read_csv(effects).to_string(index=False).splitlines())
    joint = read_json(run / "joint_tests.json")
    if joint:
        lines.append(f"\njoint tests: {json.dumps(joint, ensure_ascii=False)}")
    return lines


def page_probabilities(run: Path) -> list[str]:
    lines = ["[3/4] STANDARDIZED STATE PROBABILITIES"]
    path = run / "state_probabilities_36_48_60.csv"
    if path.is_file():
        lines.extend(pd.read_csv(path).to_string(index=False).splitlines())
    else:
        lines.append("state probability file absent")
    return lines


def page_descriptive(run: Path) -> list[str]:
    lines = ["[4/4] DESCRIPTIVE mRS AND ARTIFACTS"]
    path = run / "mrs_trajectory_distribution.csv"
    if path.is_file():
        frame = pd.read_csv(path)
        lines.extend(frame.to_string(index=False).splitlines())
    lines.append("\nartifacts:")
    for name in (
        "status.json", "audit.json", "transition_counts.csv", "transition_effects.csv",
        "state_probabilities_36_48_60.csv", "functional_state_trajectory.png",
        "transition_effects.png", "mrs_trajectory_distribution.csv", "mrs_trajectory.png",
        "imputation.json", "fit_diagnostics.json",
    ):
        lines.append(f"  {'EXISTS' if (run / name).exists() else 'MISSING'}: {run / name}")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path, help="Exact saved recovery-trajectory directory")
    parser.add_argument("--no-pause", action="store_true", help="Print all pages without waiting")
    args = parser.parse_args()
    run = args.run.resolve()
    status, audit = require_run(run)
    pages = [page_status(run, status, audit), page_transitions(run), page_probabilities(run), page_descriptive(run)]
    for index, lines in enumerate(pages):
        if index and not args.no_pause:
            input("\n截图后按回车显示下一页（Ctrl+C退出）：")
        print("\n".join(lines))


if __name__ == "__main__":
    main()
