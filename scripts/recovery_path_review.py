"""Print six screenshot-sized pages from one saved recovery-path run, read only."""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from wmh_hcy.common import DataError
from wmh_hcy.imputation import pool_scalar
from wmh_hcy.studies.design import StudyDesign
from wmh_hcy.studies.pooling import joint_test
from wmh_hcy.studies.recovery_path import _wmh_curve_contrast
from wmh_hcy.studies.registry import ModelSpec

MODELS = ("stroke_recurrence", "functional_shape", "functional_gm_log")
DEFAULT_ROOT = Path(__file__).resolve().parents[1] / "outputs/real/studies/01_recovery_path"


def read_json(path: Path) -> dict | list | None:
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def fmt(value: object) -> str:
    try:
        number = float(value)
        return f"{number:.4g}" if math.isfinite(number) else "NA"
    except (TypeError, ValueError):
        return "NA"


def resolve_run(root: Path, explicit: Path | None = None) -> tuple[Path, dict]:
    root = root.resolve()
    if explicit is None:
        pointer = read_json(root / "latest_path_attempt.json")
        if not isinstance(pointer, dict):
            raise ValueError("No latest_path_attempt.json; pass --run with the exact saved run")
        run = Path(pointer["path"]).resolve()
        if run.parent != (root / "runs").resolve() or pointer.get("run") != run.name:
            raise ValueError("Latest-attempt pointer does not match this recovery-path run root")
    else:
        run = explicit.resolve()
    state = read_json(run / "status.json")
    if not isinstance(state, dict) or state.get("study") != "recovery_path":
        raise ValueError("Expected a saved recovery_path status.json")
    if state.get("run") != run.name or Path(state.get("path", "")).resolve() != run:
        raise ValueError("Run identity/path differs from status.json")
    if state.get("mode") != "real":
        raise ValueError("This workstation review requires mode=real")
    return run, state


def local_columns(run: Path, filename: str, names: tuple[str, ...]) -> pd.DataFrame:
    path = run / filename
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_csv(path, usecols=lambda name: name in names)


def states(data: pd.DataFrame) -> str:
    if "state60" not in data:
        return "UNAVAILABLE"
    y = pd.to_numeric(data.state60, errors="coerce")
    counts = [int(y.eq(k).sum()) for k in (0, 1, 2)]
    return "/".join(map(str, counts)) + f"/unknown={int(y.isna().sum())}"


def overview(run: Path, state: dict) -> list[str]:
    audit = read_json(run / "audit.json") or state.get("audit", {})
    lines = [f"Overall status={state.get('status')}; mode={state.get('mode')}",
             f"Components: {state.get('analysis', {})}",
             f"Error: {state.get('error', 'none')}"]
    for key in ("clinical_n", "exact_id_intersection", "base_eligible_n", "landmark_n", "events"):
        lines.append(f"  {key}={audit.get(key, 'UNAVAILABLE')}")
    lines.append("state60: independent/dependent/dead/unknown")
    for name in ("base_eligible", "landmark"):
        lines.append(f"  {name}: {states(local_columns(run, name+'.csv', ('state60',)))}")
    lines.append("Sequential exclusions (nonzero steps only):")
    for name in ("base_flow", "event_flow"):
        path = run / (name + ".csv")
        if path.is_file():
            for row in pd.read_csv(path).itertuples():
                if row.excluded_here:
                    lines.append(f"  {row.step}: excluded={row.excluded_here}; remaining={row.remaining}")
    lines.append("Chronology counts:")
    for key, value in audit.get("chronology", {}).items():
        lines.append(f"  {key}={value}")
    return lines


def saved_inputs(path: Path) -> tuple[list[str], np.ndarray, np.ndarray]:
    with np.load(path / "pooled_inputs.npz", allow_pickle=False) as saved:
        terms, params, covariance = saved["terms"].tolist(), saved["params"], saved["covariance"]
    if (params.ndim != 2 or not len(params) or params.shape[1] != len(terms)
            or covariance.shape != (len(params), len(terms), len(terms))
            or not np.isfinite(params).all() or not np.isfinite(covariance).all()):
        raise ValueError("Invalid saved coefficient/covariance dimensions or values")
    return terms, params, covariance


def saved_contrast(path: Path, value: float, reference: float) -> dict:
    terms, params, covariance = saved_inputs(path)
    coding = read_json(path / "frozen_design.json")
    if not isinstance(coding, dict):
        raise TypeError("Saved recurrence frozen_design.json is missing")
    definition = read_json(path / "model_definition.json")
    if not isinstance(definition, dict) or definition.get("family") != "cox":
        raise TypeError("Saved recurrence Cox model_definition.json is missing")
    design = StudyDesign(ModelSpec(**definition), coding)
    vector = _wmh_curve_contrast(design, terms, value, reference)
    variances = np.einsum("i,mij,j->m", vector, covariance, vector)
    pooled = pool_scalar((params @ vector).tolist(), variances.tolist())
    return {key: float(np.exp(pooled[key])) for key in ("estimate", "lower", "upper")}


def fit_lines(path: Path) -> list[str]:
    diagnostics = read_json(path / "fit_diagnostics.json")
    if not isinstance(diagnostics, list) or not diagnostics:
        return ["  Fit diagnostics: UNAVAILABLE"]
    lines = [f"  Fits saved={len(diagnostics)}; converged={sum(d.get('converged') is True for d in diagnostics)}"]
    for key in ("design_condition", "max_score_per_patient", "max_standardized_case_influence"):
        values = [float(d[key]) for d in diagnostics if key in d and np.isfinite(d[key])]
        if values:
            lines.append(f"  {key}: median={fmt(np.median(values))}; max={fmt(max(values))}")
    for key in ("max_scaled_score_per_event", "scaled_information_condition"):
        values = [d["numerical"][key] for d in diagnostics if key in d.get("numerical", {})]
        if values:
            lines.append(f"  {key}: max={fmt(max(values))}")
    warnings = Counter(w for d in diagnostics for w in d.get("warnings", []))
    lines.append(f"  Saved warnings={sum(warnings.values())}; distinct={len(warnings)}")
    return lines


def recurrence(run: Path, state: dict) -> list[str]:
    path = run / "stroke_recurrence"
    result = read_json(path / "result.json") or {}
    lines = ["First ANY-type recurrent stroke after actual F3_DATE; cause-specific Cox.",
             f"Status={result.get('status', 'MISSING')}; n={result.get('n', 'NA')}; M={result.get('imputations', 'NA')}",
             f"WMH overall 2-df P={fmt(result.get('p'))}; terms={result.get('primary_terms', 'NA')}"]
    if result.get("reason"):
        lines.append("Reason: " + result["reason"])
    if (path / "pooled_inputs.npz").is_file():
        terms, params, covariance = saved_inputs(path)
        j = terms.index("wmh_ml_rcs")
        nonlinear = joint_test(params[:, [j]], covariance[:, [j]][:, :, [j]])
        lines.append(f"WMH nonlinearity P={fmt(nonlinear['p'])}")
        volumes = local_columns(run, "landmark.csv", ("wmh_ml",)).wmh_ml.dropna()
        q25, q50, q75 = volumes.quantile([.25, .5, .75]).tolist()
        lines.append(f"WMH Q25/Q50/Q75 (mL)={fmt(q25)}/{fmt(q50)}/{fmt(q75)}")
        for label, value in (("Q25", q25), ("Q75", q75)):
            effect = saved_contrast(path, value, q50)
            lines.append(f"  {label} vs Q50: HR={fmt(effect['estimate'])} [{fmt(effect['lower'])}, {fmt(effect['upper'])}]")
        code = read_json(path / "frozen_design.json")
        lines.append(f"GM119 effect scale: 1 SD={fmt(code.get('gm119_ml', {}).get('scale'))} mL")
    coeff = path / "coefficients.csv"
    if coeff.is_file():
        table = pd.read_csv(coeff)
        for row in table.loc[table.term.eq("gm119_ml")].itertuples():
            lines.append(f"GM119 per 1 SD: HR={fmt(row.ratio)} [{fmt(row.ratio_lower)}, {fmt(row.ratio_upper)}]; P={fmt(row.p)}")
    lines += fit_lines(path)
    diagnostics = read_json(path / "fit_diagnostics.json") or []
    for term in ("wmh_ml", "wmh_ml_rcs"):
        values = [d["PH_descriptive"][term] for d in diagnostics if term in d.get("PH_descriptive", {})]
        values = [v for v in values if np.isfinite(v)]
        if values:
            lines.append(f"PH descriptive {term}: P median={fmt(np.median(values))}; range={fmt(min(values))}..{fmt(max(values))}")
    lines.append("PH P ranges are per-imputation descriptive checks, not a pooled PH test.")
    return lines


def functional(run: Path, state: dict) -> list[str]:
    lines = ["Exploratory functional-shape/log checks; historical primary model is separate."]
    shape = read_json(run / "functional_shape_tests.json") or {}
    for key in ("wmh_log1p_nonlinearity", "gm_raw_nonlinearity"):
        lines.append(f"Dependent vs independent {key}: P={fmt(shape.get(key, {}).get('p'))}")
    selected = {f"{s}:{term}" for s in ("dependent", "dead")
                for term in ("wmh_ml", "wmh_ml_rcs", "gm119_ml", "gm119_ml_rcs", "gm119_log")}
    for model in ("functional_shape", "functional_gm_log"):
        path = run / model
        result = read_json(path / "result.json") or {}
        lines.append(f"{model}: {result.get('status', 'MISSING')}; n={result.get('n', 'NA')}; saved-test P={fmt(result.get('p'))}")
        lines.append("  Tested terms: " + result.get("primary_terms", "UNAVAILABLE"))
        if result.get("reason"):
            lines.append("  Reason: " + result["reason"])
        if (path / "coefficients.csv").is_file():
            for row in pd.read_csv(path / "coefficients.csv").itertuples():
                if row.term in selected:
                    lines.append(f"  {row.term}: beta={fmt(row.estimate)} [{fmt(row.lower)}, {fmt(row.upper)}]; P={fmt(row.p)}")
        lines += fit_lines(path)[:3]
    lines.append("Spline basis betas are not standalone clinical HRs; gm119_log is natural log(mL).")
    return lines


def imputation_and_baseline(run: Path, state: dict) -> list[str]:
    names = ("age", "wmh_ml", "gm119_ml", "icv_ml", "lesion_ml", "nihss")
    lines = ["Continuous baseline: median [Q25,Q75]; observed values, before covariate MI.",
             "variable                 base_eligible                 landmark"]
    frames = [local_columns(run, filename, names) for filename in ("base_eligible.csv", "landmark.csv")]
    for name in names:
        cells = []
        for frame in frames:
            values = pd.to_numeric(frame[name], errors="coerce").dropna() if name in frame else pd.Series(dtype=float)
            q = values.quantile([.25, .5, .75])
            cells.append(f"{fmt(q.loc[.5])} [{fmt(q.loc[.25])},{fmt(q.loc[.75])}] n={len(values)}")
        lines.append(f"{name:<12} {cells[0]:<30} {cells[1]}")
    for model in MODELS:
        mi = read_json(run / model / "imputation.json") or {}
        lines.append(f"{model}: M={mi.get('m', 'NA')}; iterations={mi.get('iterations', 'NA')}; seed={mi.get('seed', 'NA')}")
        missing = {k: v for k, v in mi.get("missing", {}).items() if v}
        lines.append("  Missing: " + (", ".join(f"{k}={v}" for k, v in missing.items()) or "none/unsaved"))
        traces = mi.get("transformed_mean_traces", [])
        if traces:
            first, last = min(t["iteration"] for t in traces), max(t["iteration"] for t in traces)
            for name in missing:
                means = [np.mean([t[name] for t in traces if t["iteration"] == iteration and name in t])
                         for iteration in (first, last)]
                lines.append(f"  {name} transformed imputed-value mean: first={fmt(means[0])}; last={fmt(means[1])}")
    lines.append("Mean traces alone do not establish adequate imputation or the MAR assumption.")
    return lines


def mediation(run: Path, state: dict) -> list[str]:
    gate = read_json(run / "mediation_gate.json") or {}
    result = read_json(run / "mediation.json") or {}
    lines = [f"Gate={gate.get('status', 'MISSING')}; mediation={result.get('status', 'MISSING')}"]
    for key, value in gate.get("checks", {}).items():
        lines.append(f"  {key}={value}")
    reasons = list(dict.fromkeys(gate.get("reasons", []) + result.get("reasons", [])))
    for reason in reasons:
        lines.append("Reason: " + reason)
    if "contrasts" in result:
        lines.append(f"Exposure contrasts (high-risk, low-risk): {result['contrasts']}")
    if "bootstrap_attempts" in result:
        lines.append(f"Patient bootstrap: success={result.get('bootstrap_success', 'NA')}/{result['bootstrap_attempts']}")
    failures = result.get("bootstrap_failures", [])
    if failures:
        lines.append(f"Bootstrap failures={len(failures)}; reasons={dict(Counter(f['reason'] for f in failures))}")
    table = run / "mediation_effects.csv"
    if table.is_file():
        lines.append("Indirect probability difference (percentage points; bootstrap 95% interval):")
        for row in pd.read_csv(table).itertuples():
            lines.append(f"  {row.exposure}/{row.outcome}: {fmt(100*row.estimate)} [{fmt(100*row.lower)}, {fmt(100*row.upper)}]")
    else:
        lines.append("mediation_effects.csv: NOT GENERATED; this is not an estimate of zero.")
    return lines


def prediction(run: Path, state: dict) -> list[str]:
    result = read_json(run / "prediction/validation.json") or {}
    failure = read_json(run / "prediction_failure.json") or {}
    lines = ["Internal validation: apparent -> optimism-corrected; no clinical cutoff."]
    if failure:
        lines.append("Prediction failure: " + failure.get("reason", str(failure)))
    for name in ("full", "wmh_only"):
        if name not in result:
            lines.append(f"{name}: NO SAVED VALIDATION")
            continue
        model = result.get(name, {})
        lines.append(f"{name}: M={model.get('m', 'NA')}; bootstrap success={model.get('bootstrap_success', 'NA')}/{model.get('bootstrap_attempts', 'NA')}")
        for outcome in ("dependent", "dependent_or_dead"):
            target = model.get(outcome, {})
            apparent = target.get("apparent", {})
            lines.append(f"  {outcome}: n={apparent.get('n', 'NA')}; events={apparent.get('events', 'NA')}; observed/mean risk={fmt(apparent.get('observed_risk'))}/{fmt(apparent.get('mean_predicted_risk'))}")
            for metric in ("auc", "brier", "calibration_intercept", "calibration_slope"):
                lines.append(f"    {metric}: {fmt(apparent.get(metric))} -> {fmt(target.get('optimism_corrected_'+metric))}")
        failures = model.get("bootstrap_failures", [])
        if failures:
            lines.append(f"  Bootstrap failures={len(failures)}; example={failures[0].get('reason')}")
    lines.append("Bootstrap optimism is conditional on the first imputed dataset; apparent plot is uncorrected.")
    lines.append("Please screenshot these existing image files separately:")
    for name in ("wmh_recurrence_curve.png", "prediction/apparent_calibration.png"):
        path = run / name
        lines.append(f"  {'EXISTS' if path.is_file() else 'NOT GENERATED'}: {path}")
    return lines


def build_pages(run: Path, state: dict) -> list[tuple[str, list[str]]]:
    pages = []
    for title, function in (("运行与病例数量", overview), ("全部卒中复发关联", recurrence),
                            ("五年功能模型形状", functional), ("基线与插补诊断", imputation_and_baseline),
                            ("复发路径可行性与估计", mediation), ("预测验证与图片路径", prediction)):
        try:
            lines = function(run, state)
        except (DataError, OSError, ValueError, KeyError, AttributeError, TypeError, np.linalg.LinAlgError) as exc:
            lines = [f"READ_ERROR: {exc}; other pages remain available"]
        pages.append((title, lines))
    return pages


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--run", type=Path, help="Read this exact saved run instead of latest attempt")
    parser.add_argument("--page", choices=["all", "1", "2", "3", "4", "5", "6"], default="all")
    parser.add_argument("--no-pause", action="store_true", help="Print all selected pages without terminal pauses")
    args = parser.parse_args()
    try:
        run, state = resolve_run(args.root, args.run)
        pages = build_pages(run, state)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"Stopped: {exc}\n")
    selected = range(6) if args.page == "all" else [int(args.page)-1]
    for position, index in enumerate(selected):
        if position and not args.no_pause and sys.stdin.isatty() and sys.stdout.isatty():
            try:
                input("\n截图后按回车显示下一页（Ctrl+C退出）：")
            except EOFError:
                pass
            except KeyboardInterrupt:
                print("\n查看结束；结果文件未修改。")
                return
        title, lines = pages[index]
        print(f"\n[{index+1}/6] {title} | RUN={state['run']}")
        print("READ ONLY: saved same-run artifacts; aggregates only; no fitting or edits.")
        print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
