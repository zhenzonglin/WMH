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
    if (primary / "continuous_sbp.csv").is_file():
        d = pd.read_csv(primary / "continuous_sbp.csv")
        d = d.loc[d.status.eq("ESTIMATED")]
        if len(d):
            fig, (ax, hist) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, height_ratios=[3, 1])
            distribution = pd.read_csv(primary / "sbp_distribution.csv")
            for percentile, part in d.groupby("wmh_percentile", sort=True):
                line, = ax.plot(part.sbp, part.HR, label=f"WMH percentile {percentile}")
                ax.fill_between(part.sbp, part.HR_lower, part.HR_upper, color=line.get_color(), alpha=.15)
                values = distribution.loc[distribution.wmh_percentile.eq(percentile)]
                hist.stairs(values.n, [*values.sbp_left, values.sbp_right.iloc[-1]], color=line.get_color())
            ax.axhline(1, color="gray", linestyle="--")
            ax.set(yscale="log", ylabel="HR versus SBP 140 mmHg (log scale)", title=figure_label)
            ax.legend()
            hist.set(xlabel="Recovery SBP (mmHg), continuous", ylabel="Local N")
            fig.text(.5, .01, "Pointwise 95% CIs; WMH display neighborhoods; no extrapolation beyond supported SBP", ha="center", fontsize=8)
            figure = "primary_result.png"
            fig.tight_layout(rect=(0, .03, 1, 1))
            fig.savefig(folder / figure, dpi=160)
            plt.close(fig)
    elif (primary / "kidney_interaction_curves.csv").is_file():
        d = pd.read_csv(primary / "kidney_interaction_curves.csv")
        d = d.loc[d.status.eq("ESTIMATED")]
        if len(d):
            fig, (ax, diff) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
            for name, label in (("both_low", "Both UACR <3 mg/mmol"), ("persistent", "Both UACR >=3 mg/mmol")):
                line, = ax.plot(d.wmh_ml, 100*d[name+"_probability"], label=label)
                ax.fill_between(d.wmh_ml, 100*d[name+"_lower"], 100*d[name+"_upper"], color=line.get_color(), alpha=.15)
            ax.set(ylabel="Five-year dependence probability (%)", title=figure_label)
            ax.legend()
            diff.plot(d.wmh_ml, 100*d.difference_estimate)
            diff.fill_between(d.wmh_ml, 100*d.difference_lower, 100*d.difference_upper, alpha=.15)
            diff.axhline(0, color="gray", linestyle="--")
            diff.set(xlabel="Baseline WMH (mL), continuous", ylabel="Persistent minus both-low\n(percentage points)")
            fig.text(.5, .01, "Pointwise 95% CIs; standardized visit-state probabilities; primary test is on the relative scale", ha="center", fontsize=8)
            figure = "primary_result.png"
            fig.tight_layout(rect=(0, .03, 1, 1))
            fig.savefig(folder / figure, dpi=160)
            plt.close(fig)
    elif (folder / "cec_hdl_comparison.csv").is_file():
        d = pd.read_csv(folder / "cec_hdl_comparison.csv")
        d = d.loc[d.status.eq("ESTIMATED")].sort_values("analysis", ascending=False)
        if len(d):
            fig, ax = plt.subplots(figsize=(8, 3))
            labels = d.analysis.map({"primary": "M1: with HDL-C (primary)", "without_hdl_same_sample": "M0: without HDL-C"})
            ax.errorbar(d.estimate, range(len(d)), xerr=[d.estimate-d.lower, d.upper-d.estimate], fmt="o", capsize=4)
            ax.axvline(0, color="gray", linestyle="--")
            ax.set(yticks=range(len(d)), yticklabels=labels, xlabel="GM119 difference (mL) per 1 SD CEC; 95% CI", title=figure_label)
            fig.text(.5, .01, "Same patients, CEC scale and completed covariates; coefficient changes are not mediation estimates", ha="center", fontsize=8)
            figure = "primary_result.png"
            fig.tight_layout(rect=(0, .05, 1, 1))
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
        target = "cec"
        d = d.loc[d.term.eq(target)]
        if len(d):
            fig, ax = plt.subplots(figsize=(7, 3))
            ax.errorbar(d.estimate.to_numpy(), [0], xerr=[(d.estimate-d.lower).to_numpy(), (d.upper-d.estimate).to_numpy()], fmt="o", capsize=4)
            ax.axvline(0, color="gray", linestyle="--")
            ax.set(yticks=[0], yticklabels=[target], xlabel="Adjusted difference and 95% CI\nGM119 mL, per 1 SD CEC")
            ax.set_title(figure_label)
            figure = "primary_result.png"
            fig.tight_layout()
            fig.savefig(folder / figure, dpi=160)
            plt.close(fig)
    mode = "SYNTHETIC DEMONSTRATION — NOT PATIENT RESULTS" if state["mode"] == "synthetic" else "OBSERVATIONAL ANALYSIS"
    sections = [f"<h1>{TITLES[state['study']]}</h1><p><strong>{mode}</strong></p>",
                "<p>四项独立队列；不要求Hcy。患者行和ID仅保存在本地分析目录，不进入本报告。</p>",
                "<h2>病例流程</h2>", table(pd.read_csv(folder / "cohort_flow.csv")),
                "<h2>预设检验</h2>", table(results),
                "<p>多分类系数取指数表示 P(该状态)/P(独立) 的比值之比，不能标为HR。联合检验显著不代表每个系数都显著。</p>"]
    if figure:
        sections.append(f'<img src="{figure}" style="max-width:100%" alt="主要结果图">')
    if state["study"] == "bp":
        sections.append("<p>SBP始终连续建模；主图展示样条关联。140 mmHg仅为参照，固定点对比是补充表。阴影为逐点区间。</p>")
    if state["study"] == "kidney":
        sections.append("<p>主要检验：依赖方程中持续白蛋白尿×标准化WMH，1自由度。三个类别交互的整体检验为次要。"
                        "模型含四类UACR、WMH及三个交互，在依赖和死亡方程中均估计。主要比值表示两组WMH斜率的相对概率比之比。"
                        "绝对概率差是补充展示；两种尺度不要求得到相同交互结论。</p>")
    if (folder / "cec_hdl_comparison.csv").is_file():
        sections.extend(["<h2>CEC与HDL-C的配对模型</h2>", table(pd.read_csv(folder / "cec_hdl_comparison.csv")),
                         ("<p>主要模型调整HDL-C；M0使用同一患者、同一套完成数据与固定CEC尺度。仅删除模型中的HDL-C项。"
                          "比较系数和区间，不以P值一显著一不显著判断模型差异。</p>")])
    for name in ("coefficients", "clinical_contrasts", "standardized_states", "kidney_interaction_curves"):
        if (primary / f"{name}.csv").is_file():
            sections.extend([f"<h2>{name}</h2>", table(pd.read_csv(primary / f"{name}.csv"))])
    diag = read_json(primary / "fit_diagnostics.json")
    if diag:
        sections.extend(["<h2>诊断摘要</h2>", (f"<p>完成模型 {len(diag)}；参数数 {diag[0]['parameters']}；"
                         f"最大设计条件数 {max(d['design_condition'] for d in diag):.3g}。完整诊断见 fit_diagnostics.json。</p>")])
    sections.append("<p>缺失处理依赖MAR假设。功能状态图为固定访视结局的标准化概率，不是首次失能的累积发生率。"
                    "血压分析是实测血压的条件关联。次要分析报告研究内BH-FDR；总汇报另列本版固定四项Holm校正。</p>")
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
                    for s, p in zip(rows.study, rows.get("path", pd.Series([""]*len(STUDIES))), strict=True) if isinstance(p, str) and p)
    write_html(path / "summary.html", f"<h1>四项独立研究汇总</h1><p>mode={cfg['mode']}；未估计项不当作阴性结果。</p>"
               +table(rows)+"<ul>"+links+"</ul>")
    print(f"CNSR-III STUDIES | mode={cfg['mode']} | page={page}")
    from . import CONTRACT
    print(f"Contract={CONTRACT}; PREVIOUS_VERSION results are retained but excluded from current Holm pooling.")
    if page == 1:
        keep = [c for c in ["study", "status", "n", "parameters", "df1", "p", "p_holm_four"] if c in rows]
        print(rows[keep].to_string(index=False))
        print("Current Holm family fixed at four primary tests; missing tests not reported as p=1.")
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
    if a.get("eligible_n") == 0:
        print("Covariate missingness: NOT ASSESSED (empty eligible cohort; not evidence of absent source fields)")
    else:
        print("Entirely missing covariates:", ", ".join(a.get("covariates_entirely_missing", [])) or "none")
        print("Missing counts:", "; ".join(f"{k}:{v}" for k, v in a.get("covariate_missing", {}).items() if v) or "none")
    print("Invalid field counts (all extracted):", a.get("invalid_fields", []))
    if "invalid_fields_eligible" in a:
        print("Invalid field counts (eligible):", a["invalid_fields_eligible"])
        print("Invalid primary fields requiring review:", a.get("invalid_primary_covariates_require_review", []))
    if a.get("chd_rules"):
        print("CHD skip-rule audit (eligible):", a["chd_rules"]["eligible"])
        if a["chd_rules"]["eligible"].get("conflicts"):
            print("CHD rule conflict blocks BP fitting; local patient details: chd_conflicts.csv")
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
