"""Read-only, two-page aggregate diagnostics; Python standard library only."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import Counter
from pathlib import Path

MODELS = {"H1": "01_structure", "H2": "02_recurrence", "H3": "03_month3_update", "H4": "04_function"}


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
            f"{label:<7} {str(row.get('n', 'NA')):<7} {str(row.get('ischemic_events', 'NA')):<18} "
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--page", type=int, choices=[1, 2], default=1)
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "outputs/real"
    )
    parser.add_argument("--results", type=Path, help="Optional exact result batch directory")
    args = parser.parse_args()
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
    print("READ ONLY: aggregate files; no patient rows; no analysis is rerun.")
    for line in page_one(root, args.output_dir, status) if args.page == 1 else page_two(root):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
