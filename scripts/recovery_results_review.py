"""Review saved recovery primary and recovery-path extension separately, read only."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from recovery_path_review import (
    build_pages,
    fit_lines,
    fmt,
    local_columns,
    read_json,
    resolve_run,
    run_directories,
    states,
)

from wmh_hcy.common import DataError

DEFAULT_ROOT = Path(__file__).resolve().parents[1] / "outputs/real/studies"
CONTINUOUS = ("age", "wmh_ml", "gm119_ml", "icv_ml", "lesion_ml", "nihss")
CATEGORICAL = ("sex", "smoking", "drinking", "hypertension", "diabetes", "prior_stroke",
               "education", "pre_mrs", "toast", "mrs3")
READ_ERRORS = (DataError, OSError, ValueError, KeyError, AttributeError, TypeError, np.linalg.LinAlgError)


def resolve_recovery(root: Path, explicit: Path | None = None) -> tuple[Path, dict]:
    if explicit is not None:
        run = explicit.resolve()
    else:
        run = None
        for candidate in run_directories(root):
            state = read_json(candidate / "status.json") or {}
            has_model = any((candidate / file).is_file() for file in (
                "primary/result.json", "primary/failure.txt", "results.csv"))
            if has_model or state.get("status") in {
                "ANALYSING", "COMPLETED", "PRIMARY_NOT_ESTIMABLE", "FAILED", "PARTIAL",
            }:
                run = candidate
                break
        if run is None:
            raise ValueError("No saved recovery analysis attempt; pass --recovery-run with its exact directory")
    state = read_json(run / "status.json")
    if not isinstance(state, dict) or state.get("study") != "recovery":
        raise ValueError("Expected a saved recovery status.json")
    if state.get("run") != run.name or Path(state.get("path", "")).resolve() != run:
        raise ValueError("Recovery run identity/path differs from status.json")
    if state.get("mode") != "real":
        raise ValueError("Recovery review requires mode=real")
    return run, state


def recovery_cohort(run: Path, state: dict) -> pd.DataFrame:
    data = local_columns(run, "eligible.csv", ("state60", *CONTINUOUS, *CATEGORICAL))
    if "state60" not in data:
        raise ValueError("Same-run eligible.csv/state60 is missing")
    y = pd.to_numeric(data.state60, errors="raise")
    if not y.dropna().isin([0, 1, 2]).all():
        raise ValueError("Unexpected saved state60 codes")
    audit = read_json(run / "audit.json") or state.get("audit", {})
    actual = {"eligible_n": len(data), "outcome_observed_n": int(y.notna().sum()),
              "independent": int(y.eq(0).sum()), "dependent": int(y.eq(1).sum()),
              "dead": int(y.eq(2).sum()), "unknown": int(y.isna().sum())}
    for key, count in actual.items():
        if key in audit and int(audit[key]) != count:
            raise ValueError(f"Same-run {key} differs between cohort and audit")
    primary = read_json(run / "primary/result.json") or {}
    if primary.get("status") == "ESTIMATED" and int(primary.get("n", -1)) != actual["outcome_observed_n"]:
        raise ValueError("Primary model N differs from same-run observed five-year cohort")
    return data


def continuous_cell(values: pd.Series) -> str:
    observed = pd.to_numeric(values, errors="coerce").dropna()
    if observed.empty:
        return "NA"
    q25, q50, q75 = observed.quantile([.25, .5, .75])
    return f"{fmt(q50)}[{fmt(q25)},{fmt(q75)}] n={len(observed)}"


def cohort_page(run: Path, state: dict, data: pd.DataFrame) -> list[str]:
    y = pd.to_numeric(data.state60)
    known = data.loc[y.notna()]
    lines = [f"Overall status={state.get('status')}; mode={state.get('mode')}; error={state.get('error', 'none')}",
             f"Eligible={len(data)}; five-year known={len(known)}; unknown={int(y.isna().sum())}",
             "state60 independent/dependent/dead/unknown=" + states(data),
             "Continuous: median[Q25,Q75], observed N before MI.",
             "  Columns: ALL KNOWN | INDEPENDENT | DEPENDENT | DEAD"]
    for name in CONTINUOUS:
        if name not in known:
            lines.append(name + ": COLUMN_MISSING")
            continue
        cells = [continuous_cell(known[name])]
        cells += [continuous_cell(known.loc[known.state60.eq(k), name]) for k in (0, 1, 2)]
        lines.append(f"  {name}: " + " | ".join(cells))
    flow = run / "cohort_flow.csv"
    if flow.is_file():
        lines.append("Sequential exclusions (nonzero steps):")
        for row in pd.read_csv(flow).itertuples():
            if row.excluded_here:
                lines.append(f"  {row.step}: excluded={row.excluded_here}; remaining={row.remaining}")
    coding = read_json(run / "primary/frozen_design.json") or {}
    lines += [f"WMH effect: 1 SD log(1+mL)={fmt(coding.get('wmh_ml', {}).get('scale'))}",
              f"GM119 effect: 1 SD raw volume={fmt(coding.get('gm119_ml', {}).get('scale'))} mL",
              "ICV adjusted separately; original function cohort is distinct from recurrence landmark cohort."]
    return lines


def categorical_page(run: Path, state: dict, data: pd.DataFrame) -> list[str]:
    known = data.loc[data.state60.notna()]
    lines = ["Categorical baseline among five-year known outcomes, before MI.",
             "Codes retain source meanings; percent denominators exclude missing values.",
             "  Each level: ALL KNOWN | INDEPENDENT | DEPENDENT | DEAD n(%)"]
    for name in CATEGORICAL:
        if name not in known:
            lines.append(name + ": COLUMN_MISSING")
            continue
        values = pd.to_numeric(known[name], errors="coerce")
        groups = [values, *(values.loc[known.state60.eq(k)] for k in (0, 1, 2))]
        missing = "/".join(str(int(g.isna().sum())) for g in groups)
        lines.append(f"{name}: missing ALL/I/D/dead={missing}")
        for level in sorted(values.dropna().unique()):
            cells = [f"{int(g.eq(level).sum())}({100*g.eq(level).sum()/g.notna().sum():.1f}%)"
                     if g.notna().any() else "NA" for g in groups]
            lines.append(f"  code={level:g}: " + " | ".join(cells))
    return lines


def model_page(run: Path, names: tuple[str, ...]) -> list[str]:
    lines = ["Same-run coefficients; multinomial RRR vs independent, binary OR labeled separately."]
    table = local_columns(run, "results.csv", ("analysis", "status", "n", "p", "reason"))
    for name in names:
        result = read_json(run / name / "result.json")
        if not isinstance(result, dict):
            rows = table.loc[table.analysis.eq(name)] if "analysis" in table else pd.DataFrame()
            if len(rows) > 1:
                raise ValueError(f"Duplicate results rows: {name}")
            result = rows.iloc[0].to_dict() if len(rows) else {}
        lines.append(f"{name}: {result.get('status', 'MISSING')}; n={result.get('n', 'NA')}; joint P={fmt(result.get('p'))}")
        reason = result.get("reason")
        if isinstance(reason, str) and reason:
            lines.append("  Reason: " + reason)
        if result.get("status") != "ESTIMATED":
            continue
        coefficients = local_columns(run / name, "coefficients.csv",
                                     ("term", "ratio", "ratio_lower", "ratio_upper", "p"))
        if "term" not in coefficients:
            lines.append("  Coefficients: MISSING")
            continue
        prefixes = ("",) if name == "dependent_or_dead" else ("dependent:", "dead:")
        for prefix in prefixes:
            for exposure in ("wmh_ml", "wmh_raw_ml", "gm119_ml"):
                rows = coefficients.loc[coefficients.term.eq(prefix + exposure)]
                for row in rows.itertuples():
                    scale = "OR" if name == "dependent_or_dead" else "RRR"
                    lines.append(f"  {row.term}: {scale}={fmt(row.ratio)} [{fmt(row.ratio_lower)},{fmt(row.ratio_upper)}]; P={fmt(row.p)}")
    lines.append("Joint significance does not imply both exposures significant; RRR is not HR.")
    return lines


def diagnostics_page(run: Path, state: dict, data: pd.DataFrame) -> list[str]:
    lines = []
    for name in ("primary", "month12", "month24", "observation_weighted"):
        mi = read_json(run / name / "imputation.json") or {}
        lines.append(f"{name} MI: M={mi.get('m', 'NA')}; iterations={mi.get('iterations', 'NA')}; seed={mi.get('seed', 'NA')}")
        missing = {k: v for k, v in mi.get("missing", {}).items() if v}
        lines.append("  Missing: " + (", ".join(f"{k}={v}" for k, v in missing.items()) or "none/unsaved"))
        traces = pd.DataFrame(mi.get("transformed_mean_traces", []))
        if "iteration" in traces:
            for name in missing:
                if name not in traces:
                    continue
                first = traces.loc[traces.iteration.eq(traces.iteration.min()), name].dropna()
                last = traces.loc[traces.iteration.eq(traces.iteration.max()), name].dropna()
                lines.append(f"  {name} transformed imputed-value means: {fmt(first.mean())}[{fmt(first.min())},{fmt(first.max())}] -> {fmt(last.mean())}[{fmt(last.min())},{fmt(last.max())}]")
    weights = read_json(run / "observation_weighted/observation_weights.json")
    if isinstance(weights, list) and weights:
        lines.append("Observation weights across imputations: median [min,max]")
        for key in ("minimum_probability", "weight_max", "weight_p99", "ess"):
            values = pd.to_numeric(pd.Series([w.get(key) for w in weights]), errors="coerce").dropna()
            lines.append(f"  {key}={fmt(values.median())} [{fmt(values.min())},{fmt(values.max())}]")
    else:
        lines.append("Observation-weight diagnostics: NOT SAVED")
    lines.append("MI means alone do not prove MAR or adequate imputation; weighting does not adjust recurrence.")
    return lines


def states_page(run: Path, state: dict, data: pd.DataFrame) -> list[str]:
    lines = ["Standardized five-year states (%); same-run fixed WMH/GM scenarios.",
             "WMH mL / GM mL: independent [95%CI] | dependent [95%CI] | dead [95%CI]"]
    table = local_columns(run / "primary", "standardized_states.csv",
                          ("wmh_ml", "gm119_ml", "state", "probability", "lower", "upper"))
    if not table.empty:
        for (wmh, gm), group in table.groupby(["wmh_ml", "gm119_ml"], sort=True):
            cells = []
            for name in ("independent", "dependent", "dead"):
                selected = group.loc[group.state.eq(name)]
                if len(selected) != 1:
                    raise ValueError("Missing or duplicated standardized-state scenario")
                row = selected.iloc[0]
                cells.append(f"{fmt(100*row.probability)}[{fmt(100*row.lower)},{fmt(100*row.upper)}]")
            lines.append(f"  {fmt(wmh)}/{fmt(gm)}: " + " | ".join(cells))
    else:
        lines.append("Standardized-state estimates: NOT SAVED")
    for name in ("primary", "month36", "month48", "observation_weighted"):
        lines.append(name + " fit diagnostics:")
        lines += fit_lines(run / name)
    lines.append("Same-run primary figure and historical report:")
    for name in ("primary_result.png", "report.html"):
        path = run / name
        lines.append(f"  {'EXISTS' if path.is_file() else 'NOT GENERATED'}: {path}")
    lines.append("State probabilities are model standardized associations, not a validated individual risk threshold.")
    return lines


def recovery_pages(run: Path, state: dict) -> list[tuple[str, list[str]]]:
    try:
        data = recovery_cohort(run, state)
    except READ_ERRORS as exc:
        # Do not continue displaying coefficients when same-run provenance checks failed.
        return [("主分析同次运行核验", [f"READ_ERROR: {exc}; recovery results withheld"])]
    pages = []
    functions = (
        ("主分析病例与连续基线", lambda: cohort_page(run, state, data)),
        ("主分析分类基线", lambda: categorical_page(run, state, data)),
        ("主分析及各月份功能状态", lambda: model_page(run, ("primary", "month12", "month24", "month36", "month48"))),
        ("主分析敏感性估计", lambda: model_page(run, ("complete_case", "actual_improvement", "raw_wmh", "dependent_or_dead", "observation_weighted"))),
        ("主分析插补与权重诊断", lambda: diagnostics_page(run, state, data)),
        ("主分析状态概率与模型诊断", lambda: states_page(run, state, data)),
    )
    for title, function in functions:
        try:
            lines = function()
        except READ_ERRORS as exc:
            lines = [f"READ_ERROR: {exc}; other pages remain available"]
        pages.append((title, lines))
    return pages


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--path-run", required=True, type=Path, help="Exact saved recovery-path directory")
    parser.add_argument("--recovery-run", type=Path, help="Exact primary recovery directory; default newest analysis attempt")
    parser.add_argument("--page", choices=["all", *(str(i) for i in range(1, 13))], default="all")
    parser.add_argument("--no-pause", action="store_true")
    args = parser.parse_args()
    try:
        path_run, path_state = resolve_run(args.root / "01_recovery_path", args.path_run)
        pages = []
        try:
            recovery_run, recovery_state = resolve_recovery(args.root / "01_recovery", args.recovery_run)
            pages += [(title, lines, "RECOVERY", recovery_run)
                      for title, lines in recovery_pages(recovery_run, recovery_state)]
        except READ_ERRORS as exc:
            pages.append(("主分析运行核验", [f"READ_ERROR: {exc}; recovery results withheld"], "RECOVERY", None))
        pages += [(title, lines, "RECOVERY_PATH", path_run) for title, lines in build_pages(path_run, path_state)]
    except READ_ERRORS as exc:
        parser.exit(2, f"Stopped: {exc}\n")
    if args.page != "all" and int(args.page) > len(pages):
        parser.exit(2, f"Only {len(pages)} pages are available for these saved runs\n")
    selected = range(len(pages)) if args.page == "all" else [int(args.page)-1]
    for position, index in enumerate(selected):
        if position and not args.no_pause and sys.stdin.isatty() and sys.stdout.isatty():
            try:
                input("\n截图后按回车显示下一页（Ctrl+C退出）：")
            except EOFError:
                pass
            except KeyboardInterrupt:
                print("\n查看结束；结果文件未修改。")
                return
        title, lines, section, run = pages[index]
        print(f"\n[{index+1}/{len(pages)}] {title} | {section} | RUN={run.name if run else 'UNAVAILABLE'}")
        print(f"Directory: {run or 'UNAVAILABLE'}")
        print("READ ONLY: same-run saved results; aggregates only; no fitting, report generation or edits.")
        print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
