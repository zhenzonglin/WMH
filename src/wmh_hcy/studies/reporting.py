"""Patient-free HTML reports, plots and compact screen-readable diagnostic summaries."""
from __future__ import annotations

import html
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ..common import read_json
from .registry import STUDIES, TITLES


def table(frame):
    return frame.to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")


def report(folder, state, results):
    primary = folder / "primary"
    figure = None
    figure_label = "SYNTHETIC DEMONSTRATION" if state["mode"] == "synthetic" else "OBSERVATIONAL ANALYSIS"
    if (primary / "clinical_contrasts.csv").is_file():
        d = pd.read_csv(primary / "clinical_contrasts.csv")
        d = d.loc[d.sbp.eq(130) & d.status.eq("ESTIMATED")]
        if len(d):
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.errorbar(d.HR, d.wmh_percentile, xerr=[d.HR-d.HR_lower, d.HR_upper-d.HR], fmt="o", capsize=3)
            ax.axvline(1, color="gray", linestyle="--")
            ax.set(xlabel="HR: SBP 130 versus 140 mmHg (95% CI)", ylabel="WMH percentile")
            ax.set_title(figure_label)
            figure = "primary_result.png"
            fig.tight_layout()
            fig.savefig(folder / figure, dpi=160)
            plt.close(fig)
    elif (primary / "standardized_states.csv").is_file():
        d = pd.read_csv(primary / "standardized_states.csv")
        labels = [c for c in d if c not in {"state", "probability", "lower", "upper", "ci"}]
        d["scenario"] = d[labels].apply(lambda r: "; ".join(f"{k}={v:.2f}" for k, v in r.items()), axis=1)
        fig, ax = plt.subplots(figsize=(9, 5))
        for label, part in d.groupby("state", sort=False):
            ax.errorbar(part.probability*100, part.scenario,
                        xerr=[(part.probability-part.lower)*100, (part.upper-part.probability)*100],
                        fmt="o", capsize=2, label=label, alpha=.8)
        ax.set(xlabel="Standardized five-year state probability (%)", ylabel="Prespecified display scenario")
        ax.legend()
        ax.set_title(figure_label)
        figure = "primary_result.png"
        fig.tight_layout()
        fig.savefig(folder / figure, dpi=160)
        plt.close(fig)
    elif (primary / "coefficients.csv").is_file():
        d = pd.read_csv(primary / "coefficients.csv")
        target = "cer_ratio" if state["study"] == "ceramide" else "cec"
        d = d.loc[d.term.eq(target)]
        if len(d):
            fig, ax = plt.subplots(figsize=(7, 3))
            ax.errorbar(d.estimate.to_numpy(), [0], xerr=[(d.estimate-d.lower).to_numpy(), (d.upper-d.estimate).to_numpy()], fmt="o", capsize=4)
            ax.axvline(0, color="gray", linestyle="--")
            ax.set(yticks=[0], yticklabels=[target], xlabel="Adjusted difference and 95% CI\n" +
                   ("log(1 + WMH mL), per doubling of C16/C24" if target == "cer_ratio" else "GM119 mL, per 1 SD CEC"))
            ax.set_title(figure_label)
            figure = "primary_result.png"
            fig.tight_layout()
            fig.savefig(folder / figure, dpi=160)
            plt.close(fig)
    mode = "SYNTHETIC DEMONSTRATION — NOT PATIENT RESULTS" if state["mode"] == "synthetic" else "OBSERVATIONAL ANALYSIS"
    sections = [f"<h1>{TITLES[state['study']]}</h1><p><strong>{mode}</strong></p>",
                "<p>五项独立队列；不要求Hcy。患者行和ID仅保存在本地分析目录，不进入本报告。</p>",
                "<h2>病例流程</h2>", table(pd.read_csv(folder / "cohort_flow.csv")),
                "<h2>预设检验</h2>", table(results),
                "<p>多分类系数取指数表示 P(该状态)/P(独立) 的比值之比，不能标为HR。联合检验显著不代表每个系数都显著。</p>"]
    if figure:
        sections.append(f'<img src="{figure}" style="max-width:100%" alt="主要结果图">')
    for name in ("coefficients", "clinical_contrasts", "standardized_states"):
        if (primary / f"{name}.csv").is_file():
            sections.extend([f"<h2>{name}</h2>", table(pd.read_csv(primary / f"{name}.csv"))])
    diag = read_json(primary / "fit_diagnostics.json")
    if diag:
        sections.extend(["<h2>诊断摘要</h2>", (f"<p>完成模型 {len(diag)}；参数数 {diag[0]['parameters']}；"
                         f"最大设计条件数 {max(d['design_condition'] for d in diag):.3g}。完整诊断见 fit_diagnostics.json。</p>")])
    sections.append("<p>缺失处理依赖MAR假设。功能状态图为固定访视结局的标准化概率，不是首次失能的累积发生率。"
                    "血压分析是实测血压的条件关联。次要分析报告研究内BH-FDR；总汇报另列固定五项Holm校正。</p>")
    write_html(folder / "report.html", "".join(sections))


def write_html(path, body):
    path.write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>CNSR-III studies</title>'
                    '<style>body{font-family:Arial,"Microsoft YaHei",sans-serif;margin:35px;color:#213547;line-height:1.5}'
                    'table{border-collapse:collapse;font-size:13px;display:block;overflow:auto}'
                    'td,th{border:1px solid #ddd;padding:7px}th{background:#f2f4f5}h1,h2{color:#152b3c}</style>'
                    + body + '</html>', encoding="utf-8")


def summary(cfg, page=1):
    from ..common import outdir
    from .runner import read_results, study_root
    rows = read_results(cfg)
    path = outdir(cfg) / "studies"
    path.mkdir(parents=True, exist_ok=True)
    rows.to_csv(path / "summary.csv", index=False)
    links = "".join(f'<li>{html.escape(s)}: <a href="{Path(p).as_posix()}/report.html">报告</a></li>'
                    for s, p in zip(rows.study, rows.get("path", pd.Series([""]*5)), strict=True) if isinstance(p, str) and p)
    write_html(path / "summary.html", f"<h1>五项独立研究汇总</h1><p>mode={cfg['mode']}；未估计项不当作阴性结果。</p>"
               +table(rows)+"<ul>"+links+"</ul>")
    print(f"CNSR-III STUDIES | mode={cfg['mode']} | page={page}")
    if page == 1:
        keep = [c for c in ["study", "status", "n", "parameters", "df1", "p", "p_holm_five"] if c in rows]
        print(rows[keep].to_string(index=False))
        print("Holm family fixed at five primary tests; missing tests not reported as p=1.")
    else:
        for study in STUDIES:
            pointer = read_json(study_root(cfg, study) / "latest_prepared.json")
            state = read_json(Path(pointer["path"]) / "status.json") if pointer else {}
            print_audit(state)
    print(f"Report: {path / 'summary.html'}")
    return path / "summary.html"


def print_audit(state):
    print(f"\n{state.get('study', '?')} | {state.get('mode', '?')} | {state.get('status', 'NOT_RUN')}")
    if "error" in state:
        print(state["error"])
        return
    a = state.get("audit", {})
    keys = ["clinical_n", "exact_id_intersection", "eligible_n", "outcome_observed_n", "independent",
            "dependent", "dead", "unknown", "events", "early_censor", "model_parameters"]
    print(" | ".join(f"{k}={a[k]}" for k in keys if k in a))
    print("Entirely missing covariates:", ", ".join(a.get("covariates_entirely_missing", [])) or "none")
    print("Missing counts:", "; ".join(f"{k}:{v}" for k, v in a.get("covariate_missing", {}).items() if v) or "none")
    print("Invalid field counts:", a.get("invalid_fields", []))
    if a.get("unit_conflicts"):
        print("Unit conflicts:", a["unit_conflicts"])
    print("Exclusions:", "; ".join(f"{k}:{v}" for k, v in a.get("flow_excluded", {}).items()) or "none")
    print("Source files / fields:")
    import textwrap
    for source, variables in a.get("source_groups", {}).items():
        print(textwrap.fill(f"  {source}: {variables}", width=115, subsequent_indent="    "))
    if a.get("primary_design_blocker"):
        print("Design blocker:", a["primary_design_blocker"])
    print("Field provenance and all absent optional/required sources: field_audit.csv")
    print("Run:", state.get("path", ""))
