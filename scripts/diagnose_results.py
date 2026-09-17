"""Read-only screenshot diagnostics; Python standard library only."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import Counter
from pathlib import Path

MODELS = {"H1": "01_structure", "H2": "02_recurrence", "H3": "03_month3_update", "H4": "04_function"}
CORE_LEVELS = {
    "sex": (1, 2), "smoking": (1, 2, 3, 4), "drinking": (1, 2, 3, 4),
    "hypertension": (1, 2), "diabetes": (1, 2), "prior_stroke": (1, 2),
}
PARTICIPANT_FIELDS = ("patient_id", "entry", "exit", "event_type")


def read_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except FileNotFoundError:
        return {"_unavailable": "MISSING"}
    except (OSError, ValueError):
        return {"_unavailable": "UNREADABLE"}


def numbers(values) -> list[float]:
    result = []
    for value in values:
        try:
            number = float(value)
            if math.isfinite(number):
                result.append(number)
        except (TypeError, ValueError):
            continue
    return result


def span(values) -> str:
    usable = numbers(values)
    if not usable:
        return "NA"
    lo, hi = min(usable), max(usable)
    return f"{lo:.4g}" if lo == hi else f"{lo:.4g}..{hi:.4g}"


def reason_group(reason: object) -> str:
    """Fixed labels only: error text may contain local paths or identifiers."""
    text = str(reason).lower()
    for needle, label in [
        ("singular", "singular_matrix"),
        ("rank-deficient", "rank_deficient"),
        ("convergence", "nonconvergence"),
        ("no events for cause 2", "no_death_events"),
        ("no events for cause 1", "no_ischemic_events"),
        ("invalid cox", "invalid_cox_estimate"),
        ("non-finite", "nonfinite_matrix"),
        ("invalid cumulative", "invalid_cumulative_risk"),
    ]:
        if needle in text:
            return label
    return "other" if text not in {"", "none", "nan"} else "NA"


def records(value: object) -> list[dict]:
    return [row for row in value if isinstance(row, dict)] if isinstance(value, list) else []


def page_one(root: Path, output: Path, status: dict) -> list[str]:
    analyses = status.get("analyses", {})
    lines = [
        "[1/2] COUNTS AND ABSOLUTE-RISK FAILURES",
        "Model    n       ischemic events    model status / absolute-risk status",
    ]
    for label, folder in MODELS.items():
        row = analyses.get(folder, {})
        lines.append(
            f"{label:<7} {row.get('n', 'NA')!s:<7} {row.get('ischemic_events', 'NA')!s:<18} "
            f"{row.get('status', 'MISSING')} / {row.get('absolute_risk_status', 'not requested')}"
        )
    cohort = read_json(output / "prepared/cohort_summary.json")
    lines.append("Current prepared cohorts: n / ischemic-first / death-first")
    if isinstance(cohort, dict) and "_unavailable" not in cohort:
        for key, label in [("main", "H2"), ("month3", "H3")]:
            row = cohort.get(key, {})
            lines.append(
                f"  {label}: {row.get('n', 'NA')} / {row.get('ischemic_events', 'NA')} / "
                f"{row.get('deaths_first', 'NA')}"
            )
            saved = analyses.get(MODELS[label], {}).get("n")
            if saved is not None and row.get("n") != saved:
                lines.append(f"  WARNING: {label} prepared count differs from this model run.")
    else:
        lines.append("  MISSING/UNREADABLE; death counts cannot be inferred.")
    lines.append("Prepared counts are current files, not a frozen model-run snapshot.")
    for label in ["H2", "H3"]:
        folder = root / MODELS[label]
        row = analyses.get(MODELS[label], {})
        diag = read_json(folder / "risk_diagnostics.json")
        lines.append(f"{label} risk error: {reason_group(row.get('absolute_risk_reason'))}")
        failure = read_json(folder / "risk_failure.json")
        if isinstance(failure, dict) and "_unavailable" not in failure:
            stage = failure.get("stage")
            known = {"death_fit", "death_pooling", "risk_bootstrap", "month3_M0_M1_risk_update"}
            lines.append(f"  failure stage={stage if stage in known else 'other'}; "
                         f"imputation={span([failure.get('imputation')])}")
        if not isinstance(diag, dict) or "_unavailable" in diag:
            lines.append(
                f"  risk_diagnostics.json: {diag.get('_unavailable', 'INVALID') if isinstance(diag, dict) else 'INVALID'}"
            )
        else:
            failures = records(diag.get("failures"))
            counts = Counter(reason_group(r.get("reason")) for r in failures)
            lines.append(
                f"  bootstrap={diag.get('status', 'NA')}; "
                f"valid/requested at stop={diag.get('valid', 'NA')}/{diag.get('requested', 'NA')}"
            )
            lines.append(
                f"  recorded failures={len(failures)}; "
                + (", ".join(f"{k}:{v}" for k, v in counts.most_common(5)) or "none")
            )
            by_m = Counter(r.get("imputation", "NA") for r in failures)
            lines.append(
                "  failures by imputation (index starts at 0): "
                + (", ".join(f"{k}:{v}" for k, v in by_m.most_common(5)) or "none")
            )
            stages = Counter(r.get("stage") if r.get("stage") in
                             {"ischemic_fit", "death_fit", "risk_prediction"} else "unrecorded"
                             for r in failures)
            if failures:
                lines.append("  failure stages: " + ", ".join(f"{k}:{v}" for k, v in stages.items()))
        lines.append(
            "  death_coefficients.csv: "
            + ("PRESENT" if (folder / "death_coefficients.csv").is_file() else "ABSENT")
        )
    lines.append("MISSING diagnostic files do not by themselves identify the failing model.")
    return lines


def page_two(root: Path) -> list[str]:
    lines = [
        "[2/2] IMPUTATION AND MODEL CHECKS",
        "MI: missing counts; imputed-mean traces use transformed scales.",
    ]
    for label, folder in MODELS.items():
        mi = read_json(root / folder / "imputation.json")
        if not isinstance(mi, dict) or "_unavailable" in mi:
            lines.append(f"{label} MI: MISSING/UNREADABLE")
            continue
        missing = mi.get("missing", {})
        ranked = sorted(
            ((k, int(v)) for k, v in missing.items() if numbers([v]) and float(v) > 0),
            key=lambda item: item[1],
            reverse=True,
        )[:3]
        trace = records(mi.get("imputed_mean_trace_transformed_scale"))
        iterations = numbers(r.get("iteration") for r in trace)
        lines.append(
            f"{label} MI: m={mi.get('m', 'NA')}; iterations={max(iterations, default=0):g}; "
            "missing=" + (", ".join(f"{k}:{v}" for k, v in ranked) or "none recorded")
        )
        if ranked and iterations:
            variable = ranked[0][0]
            first = numbers(r.get(variable) for r in trace if r.get("iteration") == min(iterations))
            last = numbers(r.get(variable) for r in trace if r.get("iteration") == max(iterations))
            a = f"{statistics.mean(first):.5g}" if first else "NA"
            b = f"{statistics.mean(last):.5g}" if last else "NA"
            lines.append(
                f"  {variable}: first/last iteration mean={a}/{b}; last-iteration MI range={span(last)}"
            )
    lines.append("Cox diagnostics: ranges across imputed fits; PH p-values are descriptive only.")
    for label in ["H2", "H3"]:
        fits = records(read_json(root / MODELS[label] / "diagnostics.json"))
        if not fits:
            lines.append(f"{label}: diagnostics MISSING/UNREADABLE")
            continue
        lines.append(
            f"{label}: fits={len(fits)}; events={span(r.get('events') for r in fits)}; "
            f"parameters={span(r.get('parameters') for r in fits)}; rank={span(r.get('rank') for r in fits)}"
        )
        lines.append(
            f"  condition={span(r.get('scaled_condition_number') for r in fits)}; "
            f"max|gradient|={span(r.get('max_abs_gradient') for r in fits)}; "
            f"warning fits={sum(bool(r.get('warnings')) for r in fits)}"
        )
        ph = []
        for term in ["H", "W", "H_x_W"]:
            values = [
                r.get("schoenfeld_time_correlations", {}).get(term, {}).get("p_descriptive") for r in fits
            ]
            ph.append(f"{term}={span(values)}")
        lines.append("  PH descriptive p range: " + "; ".join(ph))
    for name, folder in [("H4", "04_function"), ("H4+T1", "04_function/t1_extension")]:
        fits = records(read_json(root / folder / "diagnostics.json"))
        if not fits:
            lines.append(f"{name}: diagnostics MISSING/UNREADABLE")
            continue
        lines.append(
            f"{name}: fits={len(fits)}; converged={sum(r.get('converged') is True for r in fits)}; "
            f"warning fits={sum(bool(r.get('warnings')) for r in fits)}"
        )
        path = root / folder / "threshold_diagnostics.csv"
        if path.is_file():
            with path.open(encoding="utf-8-sig", newline="") as handle:
                thresholds = list(csv.DictReader(handle))
            hw = [r for r in thresholds if r.get("term") == "H_x_W"]
            cuts = sorted({r.get("threshold", "?") for r in hw})
            lines.append(
                f"  HW thresholds={','.join(cuts) or 'NA'}; coefficient range={span(r.get('estimate') for r in hw)}; "
                f"failed rows={sum(bool(r.get('failure')) for r in thresholds)}"
            )
        else:
            lines.append("  threshold_diagnostics.csv: MISSING")
    lines.append(
        "These summaries screen for problems; they do not certify convergence, MAR or proportional odds."
    )
    return lines


def read_selected_csv(path: Path, fields: tuple[str, ...]) -> list[dict]:
    """Read only named columns; callers must never print raw rows or errors."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not set(fields).issubset(reader.fieldnames or []):
            raise ValueError("Required columns absent")
        return [{name: row[name] for name in fields} for row in reader]


def participant_index(rows: list[dict]) -> dict[str, tuple[float, float, float]]:
    index = {}
    for row in rows:
        patient = row["patient_id"]
        times = numbers(row[name] for name in PARTICIPANT_FIELDS[1:])
        if (not patient or patient in index or len(times) != 3
                or times[0] < 0 or times[1] <= times[0] or times[2] not in (0, 1, 2)):
            raise ValueError("Invalid participant record")
        index[patient] = tuple(times)
    return index


def page_three(root: Path, output: Path, status: dict) -> list[str]:
    lines = [
        "[3] DEATH-EVENT DISTRIBUTION (LOCAL ROWS -> AGGREGATE COUNTS ONLY)",
        "Check: exact IDs + entry/exit + event type against the saved model-run list.",
        "Cells are observed category: N/IS-first/death-first; MISS includes empty/NaN.",
    ]
    for label, name in [("H2", "main"), ("H3", "month3")]:
        folder = root / MODELS[label]
        try:
            frozen = read_selected_csv(folder / "analysis_participants.csv", PARTICIPANT_FIELDS)
            current = read_selected_csv(
                output / f"prepared/cohort_{name}.csv", PARTICIPANT_FIELDS + tuple(CORE_LEVELS)
            )
            frozen_index = participant_index(frozen)
            current_index = participant_index(current)
        except (OSError, ValueError, KeyError, TypeError, csv.Error):
            lines.append(f"{label}: UNAVAILABLE/INVALID inputs; no category counts printed.")
            continue
        saved_n = status.get("analyses", {}).get(MODELS[label], {}).get("n")
        if not frozen_index or frozen_index != current_index or saved_n != len(frozen_index):
            lines.append(f"{label}: MISMATCH/EMPTY run list, current cohort or saved n; counts withheld.")
            continue
        event_counts = Counter(int(item[2]) for item in frozen_index.values())
        lines.append(
            f"{label}: MATCH; N={len(current)}; IS={event_counts[1]}; death={event_counts[2]}; "
            f"censored={event_counts[0]}"
        )
        for variable, levels in CORE_LEVELS.items():
            cells = {str(level): Counter() for level in levels}
            cells.update(MISS=Counter(), INVALID=Counter())
            for row in current:
                raw = row[variable]
                value = numbers([raw])
                if value and value[0] in levels:
                    key = str(int(value[0]))
                elif raw is None or raw.strip().lower() in ("", "nan", "na", "none"):
                    key = "MISS"
                else:
                    key = "INVALID"
                cells[key][int(float(row["event_type"]))] += 1
            text = []
            for key, counts in cells.items():
                if key in ("MISS", "INVALID") and not counts:
                    continue
                text.append(f"{key}:{sum(counts.values())}/{counts[1]}/{counts[2]}")
            lines.append(f"  {variable}: " + "  ".join(text))
        path = folder / "death_coefficients.csv"
        try:
            coefficients = read_selected_csv(path, ("estimate", "se"))
            max_beta = max((abs(v) for v in numbers(r["estimate"] for r in coefficients)), default=None)
            max_se = max(numbers(r["se"] for r in coefficients), default=None)
            lines.append(
                f"  Saved death model: {len(coefficients)} terms; "
                f"max|beta|={span([max_beta])}; max SE={span([max_se])}"
            )
        except (OSError, ValueError, KeyError, TypeError, csv.Error):
            lines.append("  Saved death model: MISSING/UNREADABLE")
    lines.extend([
        "Matching verifies membership/times/events only; covariates are CURRENT, pre-imputation.",
        "Few/zero deaths in a cell flag sparse data, not proof of Cox separation.",
        "No fitting, imputation, variable selection, file writes or patient IDs in this output.",
    ])
    return lines


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--page", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "outputs/real"
    )
    parser.add_argument("--results", type=Path, help="Optional exact result batch directory")
    parser.add_argument("--longterm", action="store_true", help="Summarize the separate 2-5 year extension")
    args = parser.parse_args()
    if args.longterm and args.page == 3:
        parser.error("--page 3 checks the one-year H2/H3 cohorts only; omit --longterm.")
    if args.longterm:
        args.output_dir = args.output_dir / "longterm"
    root = args.results
    if root is None:
        latest = read_json(args.output_dir / "latest_results.json")
        if not isinstance(latest, dict) or not latest.get("path"):
            print("No readable latest_results.json. Supply --output-dir or --results.")
            return 2
        root = Path(latest["path"])
        if not root.is_dir():
            root = args.output_dir / "results" / root.name
    status = read_json(root / "status.json")
    if not isinstance(status, dict) or "_unavailable" in status:
        print("No readable status.json in the selected result directory.")
        return 2
    print(f"WMH DIAGNOSTICS | run={root.name} | mode={status.get('mode', 'NA')}")
    print("READ ONLY: local patient rows summarized; no IDs printed; no fitting or file writes."
          if args.page == 3 else "READ ONLY: aggregate files; no patient rows; no analysis is rerun.")
    if args.longterm:
        print("SUPPLEMENTARY 2-5 YEAR ANALYSES | " + str(status.get("status", "NA")))
        if args.page == 1:
            print("Year / model      n     IS events   status")
            for name, row in status.get("analyses", {}).items():
                print(
                    f"{name:<16} {row.get('n', 'NA')!s:<6} {row.get('events', 'NA')!s:<11} {row.get('status', 'NA')}"
                )
            for row in status.get("field_audit", []):
                print(
                    f"Y{row['year']}: IS_dd present for non-events={row['non_event_with_time']}; "
                    f"cross-year conflicts={row['cross_year_conflicts']}; missing={row.get('missing_fields') or 'none'}"
                )
        else:
            print("Year / model      ratio [pointwise 95% CI]       P / Holm P (4 years)")
            path = root / "longterm_summary.csv"
            if path.is_file():
                with path.open(encoding="utf-8-sig", newline="") as handle:
                    for row in csv.DictReader(handle):
                        print(
                            f"Y{row['year']} {row['family']:<7} "
                            f"{span([row['ratio']])} [{span([row['ratio_lower']])}, {span([row['ratio_upper']])}] "
                            f"P={span([row['p']])} / {span([row['p_holm_4']])}"
                        )
            else:
                print("longterm_summary.csv: MISSING")
            for name, row in status.get("analyses", {}).items():
                if row.get("status") == "NOT_ESTIMABLE":
                    print(f"{name} failure: {reason_group(row.get('reason'))}")
        print("No long-term competing-death absolute risk is calculated.")
        return 0
    if args.page == 1:
        lines = page_one(root, args.output_dir, status)
    elif args.page == 2:
        lines = page_two(root)
    else:
        lines = page_three(root, args.output_dir, status)
    for line in lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
