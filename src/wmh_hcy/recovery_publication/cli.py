"""Build manuscript assets from saved runs, without importing analysis runners."""
from __future__ import annotations

import argparse
import hashlib
import html
import importlib.metadata
import json
import platform
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .figures import (
    STYLES,
    dependence_points,
    figure1,
    figure2,
    flow_rows,
    save_figure,
    trajectory_figure,
    transition_figure,
)
from .inputs import load_recovery, load_trajectory, source_manifest
from .tables import build_baseline_tables, export_table

ANALYSIS_LABELS = {
    "primary": "Primary five-year analysis", "month12": "12 months", "month24": "24 months",
    "month36": "36 months", "month48": "48 months", "complete_case": "Complete cases",
    "actual_improvement": "Improvement from discharge to three months",
    "raw_wmh": "Original WMH measurement", "dependent_or_dead": "Dependence or death combined",
    "observation_weighted": "Outcome-observation probability weighting",
}
ADJUSTMENT = (
    "Both imaging measures entered jointly; adjusted for age, sex, education, smoking, alcohol use, "
    "hypertension, diabetes, previous stroke, pre-stroke mRS, admission NIHSS, TOAST subtype, "
    "acute infarct volume, three-month mRS, and intracranial volume."
)
BMJ_LINKS = ["https://authors.bmj.com/writing-and-formatting/formatting-your-paper/",
             "https://www.bmj.com/about-bmj/resources-authors/article-submission?change_country=1",
             "https://www.bmj.com/about-bmj/resources-authors/house-style"]


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def fmt_p(value: float) -> str:
    if pd.isna(value):
        return "Not available"
    return "<0.001" if value < .001 else f"{value:.3f}"


def estimate_text(row) -> str:
    return f"{row.ratio:.3f} ({row.ratio_lower:.3f}–{row.ratio_upper:.3f})"


def result_rows(main: dict, analyses: list[str]) -> list[list[str]]:
    rows = []
    for name in analyses:
        detail = main["result_details"].get(name)
        if detail is None:
            raise ValueError(f"Required saved analysis absent: {name}; no substitution is permitted")
        n = str(int(detail["n"])) if pd.notna(detail.get("n")) else "Not available"
        if detail["status"] != "ESTIMATED":
            rows.append([ANALYSIS_LABELS[name], n, "Not estimable", "—", "—", "—"])
            continue
        coefficients = main["coefficients"][name].set_index("term")
        for state in ([""] if name == "dependent_or_dead" else ["dependent:", "dead:"]):
            for exposure, measure in (("wmh_ml", "WMH"), ("gm119_ml", "GMV")):
                row = coefficients.loc[state + exposure]
                comparison = {"": "Dependence/death vs independence (OR)",
                              "dependent:": "Dependence vs independence (RRR)",
                              "dead:": "Death vs independence (RRR)"}[state]
                rows.append([ANALYSIS_LABELS[name], n, comparison, measure,
                             estimate_text(row), fmt_p(row.p)])
    return rows


def result_tables(main: dict, out: Path, *, demo: bool = False) -> list[Path]:
    wmh_sd = main["design"]["wmh_ml"]["scale"]
    gm_sd = main["design"]["gm119_ml"]["scale"]
    notes = [ADJUSTMENT,
             f"WMH: per 1 SD of log(1 + mL) ({wmh_sd:.4f}); GMV: per 1 SD of raw mL ({gm_sd:.2f} mL).",
             "RRR is the ratio of P(outcome)/P(independence) per SD; OR is used only for the binary outcome.",
             "Values and pointwise 95% confidence intervals are copied from saved, pooled model outputs."]
    rows = result_rows(main, ["primary"])
    files = export_table(out, "table2", "Table 2  Imaging measures and five-year functional status",
                         ["Outcome contrast", "Measure per SD", "RRR (95% CI)", "P value"],
                         [row[2:] for row in rows],
                         [f'Primary analysis n = {int(main["audit"]["outcome_observed_n"]):,}.', *notes,
                          ("The prespecified primary joint test assesses the two imaging coefficients in the "
                          f'dependence equation: P {fmt_p(main["result_details"]["primary"]["p"])}. '
                          "A significant joint test does not imply that both measures are significant.")], demo=demo)
    for stem, title, names in [
        ("supplement_table2", "Supplementary Table 2  Functional status at earlier visits",
         ["month12", "month24", "month36", "month48"]),
        ("supplement_table3", "Supplementary Table 3  Prespecified sensitivity analyses",
         ["complete_case", "actual_improvement", "raw_wmh", "dependent_or_dead", "observation_weighted"]),
    ]:
        reasons = [f'{ANALYSIS_LABELS[n]}: {main["result_details"][n].get("reason", "Reason not saved")}'
                   for n in names if main["result_details"].get(n, {}).get("status") != "ESTIMATED"]
        extra = []
        if stem == "supplement_table2":
            for n in names:
                saved = main["results"].loc[main["results"].analysis == n]
                if not saved.empty and saved.iloc[0]["status"] == "ESTIMATED":
                    row = saved.iloc[0]
                    adjusted = fmt_p(row.get("p_bh_secondary", float("nan")))
                    extra.append(f'{ANALYSIS_LABELS[n]}: joint P {fmt_p(row.p)}; '
                                 f'BH-adjusted joint P {adjusted} (saved four-visit family).')
        else:
            extra.append("The original WMH sensitivity analysis retains the primary model scaling constant.")
        files.extend(export_table(out, stem, title,
                                  ["Analysis", "n", "Outcome contrast", "Measure", "Effect (95% CI)", "P value"],
                                  result_rows(main, names), notes + [
                                      "Coefficient P values are unadjusted. Not estimable is not a null finding.",
                                      *extra, *reasons], demo=demo))
    return files


def trajectory_tables(data: dict, out: Path, *, demo: bool = False) -> list[Path]:
    frame = data["effects"]
    rows = [[r.transition.replace("_to_", " → ").replace("_", " ").capitalize(),
             "WMH" if r.exposure == "wmh_ml" else "GMV", estimate_text(r), fmt_p(r.p)]
            for r in frame.itertuples(index=False)]
    files = export_table(out, "optional_table_s4", "Optional Supplementary Table 4  Functional transitions",
                         ["Transition", "Measure per SD", "RRR (95% CI)", "P value"], rows,
                         [ADJUSTMENT, ("Reference outcome is staying in the starting state; RRR compares "
                          "P(destination)/P(stay). These are not hazard ratios or simple transition risk ratios."),
                          ("Exploratory secondary analysis; unadjusted P values; patient-cluster covariance and "
                          "multiple-imputation uncertainty retained as saved.")], demo=demo)
    frame = data["probabilities"]
    rows = [[r.scenario, str(r.month), r.state,
             f"{100*r.probability:.2f} ({100*r.lower:.2f}–{100*r.upper:.2f})"]
            for r in frame.itertuples(index=False)]
    files.extend(export_table(out, "optional_table_s5", "Optional Supplementary Table 5  State probabilities",
                              ["Scenario", "Month", "State", "Probability % (95% CI)"], rows,
                              ["Saved pointwise intervals: Rubin variance with within-imputation parameter simulation.",
                               ("The other imaging measure is fixed at its median; these are standardised "
                               "associations, not validated individual risks."), *data["warnings"]], demo=demo))
    return files


def safe_output(path: Path, inputs: list[Path]) -> Path:
    path = path.resolve()
    if any(part in {"01_recovery", "01_recovery_path", "01_recovery_trajectory"} for part in path.parts):
        raise ValueError("Publication output must be separate from historical study directories")
    for source in inputs:
        source = source.resolve()
        if path == source or source in path.parents or path in source.parents:
            raise ValueError("Publication output overlaps an input run")
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError("Output already exists and is not empty; select a new directory")
    return path


def provenance(main_run: Path, trajectory_run: Path | None, main: dict) -> dict:
    names = ["status.json", "audit.json", "eligible.csv", "cohort_flow.csv", "results.csv",
             "primary/standardized_states.csv", "primary/frozen_design.json"]
    if not (main_run / "audit.json").exists():
        names.remove("audit.json")
    for name in main["result_details"]:
        for artifact in ("result.json", "coefficients.csv"):
            if (main_run / name / artifact).is_file():
                names.append(f"{name}/{artifact}")
    sources = {"recovery": {"run": main_run.name, "files": source_manifest(main_run, names)}}
    if trajectory_run:
        names = ["status.json", "audit.json", "transition_effects.csv", "transition_counts.csv",
                 "state_probabilities_36_48_60.csv", "mrs_trajectory_distribution.csv", "fit_diagnostics.json"]
        sources["trajectory"] = {"run": trajectory_run.name,
                                 "files": source_manifest(trajectory_run, names)}
    return sources


def write_gallery(out: Path, styles: list[str], demo: bool, trajectory: bool) -> None:
    intro = ("SYNTHETIC STYLE PREVIEW — NOT RESEARCH RESULTS" if demo else
             "Saved-run figures. All styles use identical estimates and confidence intervals.")
    parts = ["<!doctype html><html lang='en'><meta charset='utf-8'><title>Recovery figure styles</title>",
             ("<style>body{font:16px Arial,sans-serif;margin:32px;color:#17212a}section{margin:32px 0}"
             "img{width:100%;max-width:780px;border:1px solid #ddd}a{color:#0066a0}"
             ".choices{display:flex;flex-wrap:wrap;gap:32px}.choice{max-width:800px}</style>"),
             "<h1>Recovery manuscript figure styles</h1><p>" + intro + "</p>", "<div class='choices'>"]
    figures = [("figure2", "Five-year dependence probability"), ("figure1", "Participant flow")]
    if trajectory:
        figures += [("optional_figure_s1", "Longitudinal state probabilities"),
                    ("optional_figure_s2", "Transition associations")]
    for style in styles:
        parts.append(f"<div class='choice'><h2>{style.capitalize()}</h2>"
                     f"<p>{html.escape(STYLES[style]['description'])}</p>")
        for stem, title in figures:
            parts.append(f"<section><h3>{title}</h3><a href='{style}/{stem}.pdf'>Vector PDF</a> · "
                         f"<a href='{style}/{stem}.svg'>Editable SVG</a><br>"
                         f"<img src='{style}/{stem}.png' alt='{html.escape(title)}'></section>")
        parts.append("</div>")
    parts.append("</div><h2>Editable tables</h2>")
    for p in sorted((out / "tables").glob("*.html")):
        parts.append(f"<p><a href='tables/{p.name}'>{p.stem}</a> · "
                     f"<a href='tables/{p.stem}.docx'>Word table</a></p>")
    (out / "style_gallery.html").write_text("\n".join(parts), encoding="utf-8")


def build(main_run: Path, trajectory_run: Path | None, out: Path, *, styles: list[str],
          preview_only: bool = False, demo: bool = False, education_labels: dict | None = None) -> Path:
    main = load_recovery(main_run, allow_demo=demo)
    trajectory = load_trajectory(trajectory_run, allow_demo=demo) if trajectory_run else None
    if not styles or any(style not in STYLES for style in styles):
        raise ValueError("Select at least one recognised publication style")
    if trajectory and main["audit"]["eligible_n"] != trajectory["audit"]["base_eligible_n"]:
        raise ValueError("Recovery and trajectory base cohort sizes differ; select compatible saved runs")
    if trajectory:
        for key, state in (("independent", "0"), ("dependent", "1"), ("dead", "2"), ("unknown", "unknown")):
            if main["audit"][key] != trajectory["audit"]["state_counts"]["60"][state]:
                raise ValueError("Recovery and trajectory five-year state counts differ; review saved runs")
    # Validate every requested result table before creating any deliverables.
    result_rows(main, list(ANALYSIS_LABELS))
    out = safe_output(out, [main_run] + ([trajectory_run] if trajectory_run else []))
    before = provenance(main_run, trajectory_run, main)
    out.mkdir(parents=True, exist_ok=True)
    dump(out / "status.json", {"status": "BUILDING", "mode": "synthetic" if demo else "real"})
    try:
        table_dir = out / "tables"
        build_baseline_tables(main["cohort"], table_dir, education_labels, demo=demo)
        result_tables(main, table_dir, demo=demo)
        source_dir = out / "source_data"
        source_dir.mkdir()
        pd.DataFrame(flow_rows(main)).to_csv(source_dir / "figure1_flow.csv", index=False)
        dependence_points(main).to_csv(source_dir / "figure2_probabilities.csv", index=False)
        for name, frame in main["coefficients"].items():
            frame.to_csv(source_dir / f"{name}_imaging_effects.csv", index=False)
        main["results"].to_csv(source_dir / "analysis_statuses.csv", index=False)
        if trajectory:
            trajectory_tables(trajectory, table_dir, demo=demo)
            trajectory["effects"].to_csv(source_dir / "trajectory_effects.csv", index=False)
            trajectory["probabilities"].to_csv(source_dir / "trajectory_probabilities.csv", index=False)
        for style in styles:
            for stem, plot, data in [("figure1", figure1, main), ("figure2", figure2, main)]:
                save_figure(plot(data, style), out / style, stem, preview_only=preview_only, demo=demo)
            if trajectory:
                for stem, plot in [("optional_figure_s1", trajectory_figure),
                                   ("optional_figure_s2", transition_figure)]:
                    save_figure(plot(trajectory, style), out / style, stem, preview_only=preview_only, demo=demo)
        write_gallery(out, styles, demo, bool(trajectory))
        warnings = ["Education labels require confirmation against the source codebook."] if not education_labels else []
        warnings += trajectory["warnings"] if trajectory else []
        warnings += ["No model was refitted. Rendering does not validate model assumptions or causal interpretation."]
        after = provenance(main_run, trajectory_run, main)
        if before != after:
            raise ValueError("Input artifacts changed during export; publication package cannot be certified")
        versions = {package: importlib.metadata.version(package)
                    for package in ("pandas", "numpy", "matplotlib", "pubfig", "tableone", "python-docx")}
        caption = (
            "Figure 1. Participant selection. Exclusion counts are sequential and come from the saved "
            "recovery cohort audit. Three-month survival and independence exclusions are eligibility "
            "failures, not assumed numbers of deaths or dependent participants.\n\n"
            "Figure 2. Standardised five-year probability of functional dependence at the 25th, 50th and "
            "75th percentiles of baseline WMH volume, with GMV fixed at its median. Other covariates "
            "retain the primary cohort's empirical distribution. Bars indicate saved pointwise 95% "
            "confidence intervals from multiple-imputation pooling and logit delta transformation. "
            "The three display points are not clinical categories or cut-offs.\n\n"
            "Abbreviations: WMH, white matter hyperintensity; GMV, grey matter volume (GM119); "
            "mRS, modified Rankin scale; RRR, relative risk ratio; OR, odds ratio; CI, confidence interval.\n"
        )
        if trajectory:
            caption += (
                "\nOptional Figure S1. Saved standardised state probabilities at 36, 48 and 60 months. "
                "Circles/solid lines and squares/dashed lines show the 25th and 75th percentile scenarios. "
                "The other imaging measure is fixed at its median. Bars are saved pointwise 95% CIs "
                "using Rubin variance and within-imputation parameter simulation; check the recorded "
                "simulation draw count in manifest.json. Connecting lines do not estimate unsaved visits.\n\n"
                "Optional Figure S2. Imaging associations with adjacent-visit transitions. RRRs compare "
                "P(destination)/P(staying in the starting state), per SD of log(1 + WMH) or raw GMV. "
                "Intervals retain saved patient-cluster and multiple-imputation uncertainty.\n")
        if demo:
            caption = "SYNTHETIC STYLE PREVIEW — NOT RESEARCH RESULTS\n\n" + caption
        (out / "figure_legends.txt").write_text(caption, encoding="utf-8")
        files = [{"artifact": str(p.relative_to(out)).replace("\\", "/"), "bytes": p.stat().st_size,
                  "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                 for p in sorted(out.rglob("*")) if p.is_file() and p.name != "status.json"]
        dump(out / "manifest.json", {"mode": "synthetic" if demo else "real", "sources": before,
              "styles": styles, "python": platform.python_version(), "packages": versions,
              "exporter_source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                         for p in sorted(Path(__file__).parent.glob("*.py"))},
              "no_refitting": True, "patient_data_exported": False, "preview_only": preview_only,
              "bmj_sources": BMJ_LINKS, "requirements_checked": "2026-10-15", "warnings": warnings,
              "outputs": files})
        dump(out / "status.json", {"status": "COMPLETED", "mode": "synthetic" if demo else "real",
                                   "scope": "publication_only", "warnings": warnings})
    except Exception as exc:
        dump(out / "status.json", {"status": "FAILED", "error_type": type(exc).__name__})
        raise
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Recovery manuscript figures and TableOne tables; no refits")
    parser.add_argument("--recovery-run", type=Path, help="Exact saved recovery run, never a latest pointer")
    parser.add_argument("--trajectory-run", type=Path, help="Optional exact saved trajectory run")
    parser.add_argument("--out", type=Path, help="New empty publication output directory")
    parser.add_argument("--style", choices=["all", *STYLES], default="all")
    parser.add_argument("--preview-only", action="store_true", help="Omit large 1200 dpi TIFF; retain PDF/SVG/PNG")
    parser.add_argument("--education-labels", type=Path, help='Verified local JSON map: {"1": "label", ...}')
    parser.add_argument("--demo", action="store_true", help="Synthetic style preview; no real run allowed")
    args = parser.parse_args(argv)
    if args.demo and (args.recovery_run or args.trajectory_run or args.education_labels):
        parser.error("--demo cannot be combined with real inputs or a clinical label dictionary")
    if not args.demo and not args.recovery_run:
        parser.error("--recovery-run is required (or use --demo for synthetic style previews)")
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    out = args.out or Path("outputs/publication/recovery") / ("demo_" + stamp if args.demo else stamp)
    styles = list(STYLES) if args.style == "all" else [args.style]
    try:
        if args.demo:
            from .demo import create_demo
            with tempfile.TemporaryDirectory(prefix="recovery_publication_synthetic_") as temporary:
                main_run, trajectory_run = create_demo(Path(temporary))
                build(main_run, trajectory_run, out, styles=styles, preview_only=args.preview_only, demo=True)
        else:
            labels = json.loads(args.education_labels.read_text(encoding="utf-8")) if args.education_labels else None
            build(args.recovery_run, args.trajectory_run, out, styles=styles,
                  preview_only=args.preview_only, education_labels=labels)
    except (ValueError, OSError, ImportError) as exc:
        print(f"Stopped: {exc}", file=sys.stderr)
        return 2
    print(f"Completed publication export: {out.resolve()}")
    print("Open style_gallery.html to compare styles; figures/tables use saved estimates without refitting.")
    return 0
