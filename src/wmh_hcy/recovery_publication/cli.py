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

from .effect_figures import (
    EFFECT_VIEWS,
    combined_effect_figure,
    combined_effect_rows,
    trajectory_effect_figure,
)
from .figures import (
    STYLES,
    dependence_points,
    figure2,
    primary_probability_rows,
    save_figure,
    trajectory_figure,
)
from .inputs import load_recovery, load_trajectory, source_manifest
from .tables import build_baseline_tables, export_table
from .theme import FONT_SIZE, resolved_font

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


def probability_interval_table(main: dict, out: Path, *, demo: bool = False) -> list[Path]:
    """Retain each component's saved CI alongside compositional graphics."""
    frame = primary_probability_rows(main)
    rows = [[f'{r.wmh_ml:.2f}', f'{r.gm119_ml:.2f}', r.state.capitalize(),
             f'{100*r.probability:.2f} ({100*r.lower:.2f}–{100*r.upper:.2f})']
            for r in frame.itertuples(index=False)]
    return export_table(out, "display_probability_intervals", "Display companion: five-year state probabilities",
                        ["WMH (mL)", "GMV (mL)", "State", "Probability % (95% CI)"], rows,
                        ["Saved model-standardised probabilities at the two imaging measures' quartile scenarios.",
                         "Other covariates retain the primary cohort's empirical distribution.",
                         "Intervals are pointwise for individual components, not cumulative stack boundaries.",
                         "Scenarios are not patient subgroups, interaction tests, or individual risk thresholds."],
                        demo=demo)


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


def provenance(main_run: Path, trajectory_run: Path | None) -> dict:
    names = ["status.json", "audit.json", "eligible.csv", "cohort_flow.csv", "results.csv",
             "primary/standardized_states.csv", "primary/frozen_design.json"]
    if not (main_run / "audit.json").exists():
        names.remove("audit.json")
    for name in ANALYSIS_LABELS:
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


def write_gallery(out: Path, styles: list[str], effect_views: list[str], demo: bool, trajectory: bool) -> None:
    from PIL import Image
    intro = ("SYNTHETIC STYLE PREVIEW — NOT RESEARCH RESULTS" if demo else
             "Saved-run results only — no statistical refitting.")
    parts = ["<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>Recovery统计图型比较</title>",
             ("<style>body{font:16px Arial,sans-serif;margin:32px;color:#17212a;line-height:1.5}"
              "main{max-width:1800px;margin:auto}section{margin:32px 0}"
              "img{height:auto;max-width:100%;display:block;margin:20px auto}a{color:#0066a0}"
              ".notice{padding:16px;background:#fff3df}"
              ".choices{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,690px),1fr));gap:24px}"
              ".choice{padding:16px;border:1px solid #ddd;border-radius:6px}h3{margin-top:0}"
              "nav a{margin-right:18px}.choice p{font-size:14px}</style><main>"),
             "<h1>Recovery：不同统计展示方式</h1><p class='notice'>" + intro + "</p>",
             ("<p>三套图型用于比较并选择。RRR各版本表达相同估计；概率图版本使用不同的已保存场景，"
              "详见图注。可分别选择主图与补图，无需将所有版本放入论文。</p>"),
             "<p>Seaborn + SciencePlots · 按投稿物理宽度预览（96 CSS dpi）；点击图片查看完整分辨率。</p>",
             ("<nav><a href='#effects'>主要关联</a><a href='#probabilities'>五年概率</a>"
              "<a href='#tables'>三线表</a></nav>")]

    def group(anchor: str, title: str, description: str, cards: list[tuple[str, str, str]]) -> None:
        parts.append(f"<section id='{anchor}'><h2>{title}</h2><p>{description}</p><div class='choices'>")
        for label, stem, detail in cards:
            with Image.open(out / (stem + ".png")) as raster:
                css_width = raster.width / raster.info.get("dpi", (300, 300))[0] * 96
            parts.append(f"<article class='choice'><h3>{html.escape(label)}</h3>"
                         f"<p>{html.escape(detail)}</p><a href='{stem}.pdf'>PDF</a> · "
                         f"<a href='{stem}.svg'>SVG</a> · <a href='{stem}.png'>查看完整图片</a>"
                         f"<a href='{stem}.png'><img loading='lazy' src='{stem}.png' style='width:{css_width:.1f}px' "
                         f"alt='{html.escape(title + ': ' + label)}'></a></article>")
        parts.append("</div></section>")

    effect_names = {"forest": "A · 森林图", "heatmap": "B · 效应热图", "diverging": "C · 发散效应条形图"}
    probability_names = {"estimates": "A · 概率及置信区间", "heatmap": "B · 概率热图",
                         "composition": "C · 三状态组成图"}
    group("effects", "主要与补充关联：同一张RRR图",
          "五年主分析置于首位；36、48个月为随访节点分析；五年结局可观测加权为敏感性分析。"
          "所有版本包含这四项分析，区分依赖/死亡结局与WMH/GMV。",
          [(effect_names[v], f"effects/{v}/combined_rrr", EFFECT_VIEWS[v]["description"]) for v in effect_views])
    group("probabilities", "五年功能状态概率",
          "A：GM中位数下的三个WMH场景及95%CI。B：WMH×GM九种组合的依赖概率。"
          "C：九种组合的三状态组成；分量置信区间见配套表。",
          [(probability_names[s], f"{s}/figure2", STYLES[s]["description"]) for s in styles])
    if trajectory:
        group("trajectories", "纵向标准化状态概率",
              "A：第25/75百分位场景的访视折线与点式95%CI带。B：三个影像水平×三个访视的概率热图。"
              "C：四个影像场景的三状态堆积图；各场景均将另一影像指标固定在中位数。",
              [(probability_names[s], f"{s}/optional_figure_s1", STYLES[s]["description"]) for s in styles])
        group("transitions", "纵向转移关联",
              "四种转移分别呈现；每项RRR的参照为保持其起始状态，不能与五年结局的RRR混作同一结局。",
              [(effect_names[v], f"effects/{v}/trajectory_rrr", EFFECT_VIEWS[v]["description"])
               for v in effect_views])
    parts.append("<section id='tables'><h2>可编辑三线表与区间配套表</h2>"
                 "<p>原有主表与补表保留分析记录。图型选择不会删除表中的未估计状态或其他已保存分析。</p>")
    for p in sorted((out / "tables").glob("*.html")):
        parts.append(f"<p><a href='tables/{p.name}'>{p.stem}</a> · "
                     f"<a href='tables/{p.stem}.docx'>Word table</a></p>")
    parts.append("</section><p><a href='figure_legends.txt'>完整英文图注</a> · "
                 "<a href='manifest.json'>来源与导出清单</a></p></main></html>")
    (out / "style_gallery.html").write_text("\n".join(parts), encoding="utf-8")


def figure_legends(styles: list[str], effect_views: list[str], trajectory: bool, *, demo: bool = False) -> str:
    entries = ["SYNTHETIC STYLE PREVIEW — NOT RESEARCH RESULTS"] if demo else []
    entries += [
        ("Combined associations. Five-year primary analysis (emphasised), secondary 36- and 48-month "
        "analyses, and five-year outcome-observation-weighted sensitivity analysis. Each panel shows "
        "an imaging measure and a dependent-versus-independent or dead-versus-independent contrast. "
        "RRR compares P(outcome)/P(independence) per 1 SD of log(1 + WMH volume) or raw GMV. "
        "The primary model is the comparison reference, not a meta-analysis summary. All confidence "
        "intervals are saved pointwise 95% intervals. No between-analysis differences are tested. "
        "Exact RRRs and intervals are supplied in Tables 2, S2 and S3; the forest and diverging views "
        "display the intervals, while the heatmap gives an overview of point estimates."),
        ADJUSTMENT,
    ]
    effect_notes = {
        "forest": "Forest view: horizontal axes are logarithmic; the null value is RRR = 1. "
                  "Primary estimates are visually emphasised. Intervals retain the saved pooled bounds.",
        "heatmap": "Effect heatmap: colours encode log(RRR), centred at zero (RRR = 1); cells retain "
                   "the original RRR. Corresponding 95% CIs are given in the companion tables. "
                   "Colour does not indicate statistical significance.",
        "diverging": "Diverging bars: the axis shows 100 × (RRR − 1), with the saved interval endpoints "
                     "transformed identically. This is the percentage change in P(outcome)/P(reference) "
                     "per SD, not a percentage-point change in outcome probability or an individual risk. "
                     "Original RRRs and intervals remain in the companion tables. The null value is zero.",
    }
    entries.extend(effect_notes[v] for v in effect_views)
    entries.append("Any not-estimable result has no plotted estimate and is explicitly labelled. "
                   "It is not a null estimate. Detailed statuses remain in the companion tables.")
    primary_notes = {
        "estimates": "Figure 2, interval view. Standardised five-year dependence probabilities at the "
                     "25th, 50th and 75th percentiles of WMH, with GMV at its median. The x axis is raw "
                     "WMH mL; points are discrete saved scenarios, not a continuous fitted curve.",
        "heatmap": "Figure 2, probability heatmap. Dependence probabilities for the nine combinations "
                   "of WMH and GMV at their 25th, 50th and 75th percentiles. Cell labels show probability "
                   "alone; saved pointwise 95% CIs are provided in display_probability_intervals. "
                   "These model-standardised scenarios are not observed subgroups "
                   "and do not constitute an interaction test.",
        "composition": "Figure 2, composition view. Independent, dependent and dead probabilities for "
                       "the nine joint imaging scenarios, faceted by GMV. Each bar sums to 100%. "
                       "The display_probability_intervals companion table provides every component's "
                       "saved 95% CI; component intervals are not cumulative stack-boundary intervals.",
    }
    entries.extend(primary_notes[s] for s in styles)
    entries.append("Five-year probability scenarios retain other covariates' empirical distribution in "
                   "the primary cohort. Saved intervals retain multiple-imputation pooling and logit delta "
                   "transformation. Quantile scenarios are not clinical cut-offs or validated individual risks.")
    if trajectory:
        entries += [
            ("Longitudinal probabilities, optional Figure S1. Saved standardised state probabilities at "
            "36, 48 and 60 months. The interval view connects the 25th and 75th percentile scenarios at "
            "saved visits only, with shaded pointwise 95% intervals (not simultaneous confidence bands). "
            "The heatmap shows all three quantiles and visits separately for each "
            "imaging measure and state. Composition panels show the 25th/75th percentile scenarios of "
            "each measure; each bar sums to 100%. The other imaging measure is fixed at its median. "
            "Heatmap colour scales are shared between imaging measures within each state, but differ "
            "between states. Line-plot state panels also use different y-axis ranges. "
            "Stack components are labelled through the legend; numerical values are in the companion table. "
            "Optional Table S5 contains all saved component intervals. Intervals are pointwise and use "
            "Rubin variance with within-imputation parameter simulation; inspect the draw-count warning "
            "in manifest.json. Lines do not estimate unsaved visits."),
            ("Trajectory transition associations. Saved RRRs compare P(destination)/P(staying in the "
            "starting state), per SD of log(1 + WMH) or raw GMV. These are neither hazard ratios nor "
            "simple transition risk ratios. Intervals retain saved patient-cluster and multiple-imputation "
            "uncertainty. The forest, heatmap and diverging view definitions above apply with this "
            "transition-specific reference. These exploratory secondary estimates are separate from "
            "the five-year outcome contrasts."),
        ]
    entries.append("Abbreviations: WMH, white matter hyperintensity; GMV, grey matter volume (GM119); "
                   "RRR, relative risk ratio; SD, standard deviation; CI, confidence interval.")
    return "\n\n".join(entries) + "\n"


def build(main_run: Path, trajectory_run: Path | None, out: Path, *, styles: list[str],
          effect_views: list[str] | None = None, preview_only: bool = False,
          demo: bool = False, education_labels: dict | None = None) -> Path:
    before = provenance(main_run, trajectory_run)
    main = load_recovery(main_run, allow_demo=demo)
    trajectory = load_trajectory(trajectory_run, allow_demo=demo) if trajectory_run else None
    if not styles or any(style not in STYLES for style in styles):
        raise ValueError("Select at least one recognised publication style")
    effect_views = list(EFFECT_VIEWS) if effect_views is None else effect_views
    if not effect_views or any(view not in EFFECT_VIEWS for view in effect_views):
        raise ValueError("Select at least one recognised effect view")
    if trajectory and main["audit"]["eligible_n"] != trajectory["audit"]["base_eligible_n"]:
        raise ValueError("Recovery and trajectory base cohort sizes differ; select compatible saved runs")
    if trajectory:
        for key, state in (("independent", "0"), ("dependent", "1"), ("dead", "2"), ("unknown", "unknown")):
            if main["audit"][key] != trajectory["audit"]["state_counts"]["60"][state]:
                raise ValueError("Recovery and trajectory five-year state counts differ; review saved runs")
    # Validate every requested result table before creating any deliverables.
    result_rows(main, list(ANALYSIS_LABELS))
    out = safe_output(out, [main_run] + ([trajectory_run] if trajectory_run else []))
    if before != provenance(main_run, trajectory_run):
        raise ValueError("Input artifacts changed while loading; retry with stable saved runs")
    out.mkdir(parents=True, exist_ok=True)
    dump(out / "status.json", {"status": "BUILDING", "mode": "synthetic" if demo else "real"})
    try:
        table_dir = out / "tables"
        build_baseline_tables(main["cohort"], table_dir, education_labels, demo=demo)
        result_tables(main, table_dir, demo=demo)
        probability_interval_table(main, table_dir, demo=demo)
        source_dir = out / "source_data"
        source_dir.mkdir()
        dependence_points(main).to_csv(source_dir / "figure2_probabilities.csv", index=False)
        primary_probability_rows(main).to_csv(source_dir / "primary_probability_grid.csv", index=False)
        combined_effect_rows(main).to_csv(source_dir / "combined_rrr.csv", index=False)
        for name, frame in main["coefficients"].items():
            frame.to_csv(source_dir / f"{name}_imaging_effects.csv", index=False)
        main["results"].to_csv(source_dir / "analysis_statuses.csv", index=False)
        if trajectory:
            trajectory_tables(trajectory, table_dir, demo=demo)
            trajectory["effects"].to_csv(source_dir / "trajectory_effects.csv", index=False)
            trajectory["probabilities"].to_csv(source_dir / "trajectory_probabilities.csv", index=False)
        for style in styles:
            save_figure(figure2(main, style), out / style, "figure2", preview_only=preview_only, demo=demo)
            if trajectory:
                save_figure(trajectory_figure(trajectory, style), out / style, "optional_figure_s1",
                            preview_only=preview_only, demo=demo)
        for view in effect_views:
            save_figure(combined_effect_figure(main, view), out / "effects" / view, "combined_rrr",
                        preview_only=preview_only, demo=demo)
            if trajectory:
                save_figure(trajectory_effect_figure(trajectory, view), out / "effects" / view,
                            "trajectory_rrr", preview_only=preview_only, demo=demo)
        write_gallery(out, styles, effect_views, demo, bool(trajectory))
        warnings = ["Education labels require confirmation against the source codebook."] if not education_labels else []
        warnings += trajectory["warnings"] if trajectory else []
        warnings += ["No model was refitted. Rendering does not validate model assumptions or causal interpretation."]
        after = provenance(main_run, trajectory_run)
        if before != after:
            raise ValueError("Input artifacts changed during export; publication package cannot be certified")
        versions = {package: importlib.metadata.version(package)
                    for package in ("pandas", "numpy", "matplotlib", "seaborn", "SciencePlots",
                                    "tableone", "python-docx", "openpyxl", "tabulate", "lxml", "et_xmlfile")}
        caption = figure_legends(styles, effect_views, bool(trajectory), demo=demo)
        (out / "figure_legends.txt").write_text(caption, encoding="utf-8")
        files = [{"artifact": str(p.relative_to(out)).replace("\\", "/"), "bytes": p.stat().st_size,
                  "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                 for p in sorted(out.rglob("*")) if p.is_file() and p.name != "status.json"]
        dump(out / "manifest.json", {"mode": "synthetic" if demo else "real", "sources": before,
              "styles": styles, "effect_views": effect_views,
              "combined_rrr_analyses": ["primary", "month36", "month48", "observation_weighted"],
              "python": platform.python_version(), "packages": versions,
              "figure_engine": "Seaborn + SciencePlots with Matplotlib export",
              "font": {"family": resolved_font()["family"], "size_pt": FONT_SIZE,
                       "file": Path(resolved_font()["path"]).name,
                       "sha256": hashlib.sha256(Path(resolved_font()["path"]).read_bytes()).hexdigest()},
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
    parser.add_argument("--style", choices=["all", *STYLES], default="all",
                        help="Probability representation; replaces the former classic/colour/compact palettes")
    parser.add_argument("--effect-view", choices=["all", *EFFECT_VIEWS], default="all",
                        help="Representation of combined RRRs and optional trajectory RRRs")
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
    effect_views = list(EFFECT_VIEWS) if args.effect_view == "all" else [args.effect_view]
    try:
        if args.demo:
            from .demo import create_demo
            with tempfile.TemporaryDirectory(prefix="recovery_publication_synthetic_") as temporary:
                main_run, trajectory_run = create_demo(Path(temporary))
                build(main_run, trajectory_run, out, styles=styles, effect_views=effect_views,
                      preview_only=args.preview_only, demo=True)
        else:
            labels = json.loads(args.education_labels.read_text(encoding="utf-8")) if args.education_labels else None
            build(args.recovery_run, args.trajectory_run, out, styles=styles, effect_views=effect_views,
                  preview_only=args.preview_only, education_labels=labels)
    except (ValueError, OSError, ImportError) as exc:
        print(f"Stopped: {exc}", file=sys.stderr)
        return 2
    print(f"Completed publication export: {out.resolve()}")
    print("Open style_gallery.html to compare styles; figures/tables use saved estimates without refitting.")
    return 0
