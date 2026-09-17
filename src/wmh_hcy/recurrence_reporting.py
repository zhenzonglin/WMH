"""Compact, aggregate recurrence report and relative-association figures."""
from __future__ import annotations

import html

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator, NullFormatter, ScalarFormatter

from .common import read_json


def report_recurrence(root):
    status = read_json(root / "status.json")
    synthetic = status["mode"] == "synthetic"
    title = "SYNTHETIC VALIDATION ONLY" if synthetic else "OBSERVATIONAL ASSOCIATION"
    tables = {}
    names = ["cohort_counts", "cohort_flow", "horizon_results", "clinical_contrasts",
             "robustness_results", "diagnostic_summary"]
    for name in names:
        path = root / f"{name}.csv"
        try:
            tables[name] = pd.read_csv(path) if path.is_file() else pd.DataFrame()
        except pd.errors.EmptyDataError:
            tables[name] = pd.DataFrame()
    figures = root / "figures"
    figures.mkdir(exist_ok=True)

    def save(fig, name):
        fig.suptitle(title, fontsize=11, color="#a32c24" if synthetic else "#324f66")
        fig.tight_layout()
        for suffix in ["png", "svg", "pdf"]:
            fig.savefig(figures / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
        plt.close(fig)

    h = tables["horizon_results"]
    fig, ax = plt.subplots(figsize=(8, 4.7))
    for y, row in h.iterrows():
        if np.isfinite(row.HR) and np.isfinite(row.HR_lower) and np.isfinite(row.HR_upper):
            color = "#a64a25" if row.month == 60 else "#26718b"
            ax.errorbar(row.HR, y, xerr=[[row.HR-row.HR_lower], [row.HR_upper-row.HR]],
                        fmt="o", color=color, capsize=3)
        else:
            ax.text(.02, y, "Not estimable", transform=ax.get_yaxis_transform())
    ax.set_yticks(range(len(h)), [f"{int(m)} months" + (" (primary)" if m == 60 else "") for m in h.month])
    ax.set_xscale("log")
    ax.axvline(1, color="gray", linestyle="--")
    ax.set_xlabel("Hcy x WMH: ratio of hazard ratios (pointwise 95% CI)")
    ax.invert_yaxis()
    save(fig, "horizon_interactions")

    curve_path = root / "clinical_contrast_curve.csv"
    if curve_path.is_file():
        c = pd.read_csv(curve_path)
        fig, ax = plt.subplots(figsize=(8, 4.7))
        ax.plot(c.wmh_ml, c.HR, color="#26718b")
        ax.fill_between(c.wmh_ml, c.HR_lower, c.HR_upper, color="#26718b", alpha=.18)
        ax.axhline(1, color="gray", linestyle="--")
        ax.set_xlabel("Baseline WMH volume (mL)")
        ax.set_ylabel("HR: Hcy 15 vs 10 micromol/L (pointwise 95% CI)")
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        if not c.status.eq("ESTIMATED").any():
            ax.text(.5, .5, "10 and 15 outside local observed support", ha="center", transform=ax.transAxes)
        save(fig, "conditional_hcy_contrast")

    def table(name, columns=None):
        d = tables[name]
        if columns:
            d = d[[c for c in columns if c in d]]
        return d.to_html(index=False, float_format=lambda x: f"{x:.4g}", na_rep="—", escape=True)

    primary = h.loc[h.month.eq(60)]
    if len(primary) and primary.iloc[0].status == "ESTIMATED":
        r = primary.iloc[0]
        interpretation = (f"五年主要交互HR之比为 {r.HR:.3f}（95%CI {r.HR_lower:.3f}–{r.HR_upper:.3f}），"
                          f"双侧P={r.p:.4g}。数值>1表示WMH更高时Hcy的相对风险关联增强，<1表示减弱。"
                          "须结合区间及分段检验解读，不能解释为降低Hcy的治疗效果。")
    else:
        interpretation = "五年主模型未成功估计；请查看状态和失败原因，不以其他月份替代主要结论。"
    d = tables["diagnostic_summary"]
    diag_html = "无成功拟合的诊断记录。"
    if not d.empty:
        summary = d.groupby("analysis", sort=False).agg(
            fits=("imputation", "count"), events=("events", "first"), parameters=("parameters", "first"),
            max_information_condition=("information_condition", "max"),
            max_scaled_score=("max_scaled_score_per_event", "max"),
        ).reset_index()
        diag_html = summary.to_html(index=False, float_format=lambda x: f"{x:.4g}")
    time_path = root / "robustness/time_varying/time_interactions.csv"
    time_html = pd.read_csv(time_path).to_html(index=False, float_format=lambda x: f"{x:.4g}") \
        if time_path.is_file() else "时间分段模型未成功估计。"
    cols = ["month", "n", "events", "status", "HR", "HR_lower", "HR_upper", "p", "p_holm_6", "reason"]
    body = f"""<h1>Hcy–WMH五年复发分析 · v3</h1><strong>{title}</strong>
<p>方案修订于2026-09-17，已查看既往一年结果。五年是修订后的主要分析，不称为原先预注册的确认性检验。</p>
<p>状态：{html.escape(status['status'])}。进入时间为基线采血日；全部时点使用同一个五年事件/实际删失记录。</p>
<h2>1. 病例与事件数量</h2>{table('cohort_counts')}{table('cohort_flow')}
<h2>2. 五年主要交互与各月份敏感性分析</h2><p>{interpretation}</p>{table('horizon_results', cols)}
<p>主要检验为60个月双侧P值；其余六个月份使用固定六项Holm校正。置信区间为逐项95%区间。</p>
<img src="figures/horizon_interactions.png" width="820">
<h2>3. 临床对比：同一WMH水平下Hcy 15对10的HR</h2>
{table('clinical_contrasts', ['wmh_quantile','wmh_ml','local_n','hcy_local_p05','hcy_local_p95','status','HR','HR_lower','HR_upper'])}
<p>WMH百分位是展示点。每个点在相邻±10个百分点的WMH观测窗核验Hcy的5%–95%范围；不支持的点留空。曲线不是绝对复发概率。</p>
"""
    if curve_path.is_file():
        body += '<img src="figures/conditional_hcy_contrast.png" width="820">'
    body += f"""<h2>4. 五年稳健性与时间变化</h2>
{table('robustness_results', ['analysis','n','events','status','HR','HR_lower','HR_upper','p','reason'])}
<p>clinical_same_subset_core与clinical_extended使用同一子样本及完成数据。time_varying行检验交互系数的晚期–早期差异。</p>
{time_html}<h2>5. 数值诊断</h2>{diag_html}
<p>残差时间相关只作描述，不能以收敛通过替代比例风险判断。完整诊断在diagnostic_summary.csv及各模型diagnostics.json。</p>
<p>本流程不拟合死亡回归，不计算竞争风险绝对概率，不调用bootstrap。各时点共用一次核心插补；原始WMH与临床扩展各作独立兼容插补。失败时保留原因，不删变量重试。</p>
"""
    page = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>五年复发分析</title><style>' \
        'body{font-family:Arial,"Microsoft YaHei",sans-serif;margin:36px;color:#213747;line-height:1.65}' \
        'table{border-collapse:collapse;font-size:13px;margin:18px 0}td,th{border:1px solid #ccd4dc;padding:6px}' \
        'th{background:#eef3f7}img{max-width:100%}h2{margin-top:30px}</style><body>' + body + '</body></html>'
    target = root / "report.html"
    target.write_text(page, encoding="utf-8")
    (root / "SUMMARY.md").write_text(
        f"# 五年复发 v3\n\n{title}\n\n状态：{status['status']}\n\n{interpretation}\n\n"
        "打开同目录 report.html 查看汇总；cohort_counts.csv 为实际各时点例数和事件数。\n"
        "本报告与旧版一年结果独立；不产生绝对风险或治疗建议。\n", encoding="utf-8")
    return target
