"""Readable HTML/Markdown reports and exportable statistical figures."""
from __future__ import annotations

import html
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import DataError, outdir, read_json, record_run
from .fields import CATEGORIES, CORE


def baseline_table(data: pd.DataFrame) -> pd.DataFrame:
    """Describe observed baseline data; no P-value screening and no identifiers."""
    rows = []
    for col in dict.fromkeys(["hcy", "wmh_ml", *CORE, "wmh_raw_ml", "lesion_ml", "gm119_ml"]):
        if col not in data:
            continue
        observed = data[col].dropna()
        common = {"variable": col, "n_observed": len(observed), "n_missing": int(data[col].isna().sum())}
        if col in CATEGORIES:
            for category in CATEGORIES[col]:
                n = int(observed.eq(category).sum())
                value = f"{n} ({100*n/len(observed):.1f}%)" if len(observed) else "not measured"
                rows.append({**common, "category": str(category), "summary": value})
        else:
            q = observed.quantile([.25, .5, .75])
            value = f"{q.loc[.5]:.2f} [{q.loc[.25]:.2f}, {q.loc[.75]:.2f}]" if len(observed) else "not measured"
            rows.append({**common, "category": "median [Q1, Q3]", "summary": value})
    return pd.DataFrame(rows)


def save_figure(fig, folder: Path, name: str) -> None:
    for ext in ["png", "svg", "pdf"]:
        fig.savefig(folder / f"{name}.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def report(cfg: dict) -> Path:
    latest = read_json(outdir(cfg) / "latest_results.json")
    if not latest:
        raise DataError("No analysis result manifest. Run analyse first")
    root = Path(latest["path"])
    status = read_json(root / "status.json")
    folder = root / "report"
    folder.mkdir(exist_ok=True)
    synthetic = cfg["mode"] == "synthetic"
    banner = "SYNTHETIC DATA — SOFTWARE VALIDATION ONLY" if synthetic else "OBSERVATIONAL ANALYSIS"
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    sections = [f"<h1>Hcy–WMH分析报告</h1><h2>{banner}</h2>",
                "<p>H1结构关联 → H2复发交互（主要）→ H3三个月风险更新 → H4功能交互。</p>",
                "<p>主要终点：y1_is及y1_is_dd。曲线在观测数据支持范围内展示。</p>"]
    if synthetic:
        sections.append("<p><strong>本报告中的人数、事件率、效应及P值均来自合成数据，不能用于论文结果。</strong></p>")
    hypotheses = pd.read_csv(root / "hypothesis_summary.csv")
    sections.append("<h2>四条假设的预设检验</h2><p>H2：双侧α=0.05；H1、H3、H4：固定三项Holm校正，结合效应方向判定。区间为逐项95%区间。</p>"
                    + hypotheses.to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}"))
    curve_path = root / "01_structure/structural_curve.csv"
    if curve_path.is_file():
        d = pd.read_csv(curve_path)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(d.hcy, d.estimate, color="#246A80", linewidth=2)
        ax.fill_between(d.hcy, d.lower, d.upper, color="#246A80", alpha=.18)
        ax.set(xlabel="Baseline Hcy (µmol/L)", ylabel="Adjusted mean log(1 + WMH mL)")
        ax.set_title(banner, fontsize=10)
        save_figure(fig, folder, "H1_structure")
        sections.append('<h2>H1：调整后结构曲线</h2><img src="H1_structure.png" width="850">')
    rows = [{"analysis": name, **values} for name, values in status.get("analyses", {}).items()]
    statuses = pd.DataFrame(rows)
    sections.append("<h2>分析状态</h2>" + statuses.to_html(index=False, escape=True))
    flow = pd.read_csv(outdir(cfg) / "prepared/flow.csv")
    sections.append("<h2>病例流程</h2>" + flow.to_html(index=False, escape=True))
    prepared = pd.read_csv(outdir(cfg) / "prepared/cohort_main.csv", dtype={"patient_id": str})
    baseline = baseline_table(prepared)
    baseline.to_csv(root / "baseline_table.csv", index=False)
    sections.append("<h2>观察到的基线特征</h2><p>分类百分比分母为非缺失者；编码含义见字段字典。</p>"
                    + baseline.to_html(index=False, escape=True))
    sections.append("<p>H和W系数是样条的线性分量，其单独HR不代表整个暴露范围的固定效应。"
                    "主要解释H×W以及连续风险曲线。</p>")
    missing = pd.read_csv(root / "missingness.csv")
    show = missing.loc[missing.fraction.gt(0)].sort_values("fraction").tail(15)
    if len(show):
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh(show.variable, 100*show.fraction, color="#477B94")
        ax.set_xlabel("Missing (%)")
        ax.set_title(banner, fontsize=10)
        save_figure(fig, folder, "missingness")
    curves = root / "02_recurrence/risk_curves.csv"
    if curves.is_file():
        d = pd.read_csv(curves)
        fig, ax = plt.subplots(figsize=(8, 5))
        for q, sub in d.groupby("wmh_quantile"):
            ax.plot(sub.hcy, 100*sub.risk, label=f"WMH percentile {q*100:.0f}", linewidth=2)
        interval_path = root / "02_recurrence/risk_intervals.csv"
        if interval_path.is_file() and interval_path.stat().st_size > 5:
            ci = pd.read_csv(interval_path)
            for row in ci.itertuples():
                ax.errorbar(row.hcy, 100*row.risk, yerr=[[100*(row.risk-row.lower)],
                            [100*(row.upper-row.risk)]], fmt="o", color="#263A4B", capsize=4)
        ax.set(xlabel="Hcy (µmol/L)", ylabel="1-year ischemic stroke cumulative incidence (%)")
        ax.set_title(banner, fontsize=10)
        ax.legend(frameon=False)
        save_figure(fig, folder, "absolute_risk")
        sections.append("<h2>H2：绝对风险</h2><p>区间标注于固定比较点10和15 μmol/L。</p>"
                        '<img src="absolute_risk.png" width="850">')
    forest = []
    for name in status.get("analyses", {}):
        path = root / name / "coefficients.csv"
        if path.is_file():
            d = pd.read_csv(path)
            selected = d.loc[d.term.isin(["H_x_W", "gm119_100ml"])]
            sections.append(f"<h2>{html.escape(name)}</h2>" + selected.to_html(index=False, escape=True,
                                                                                 float_format=lambda x: f"{x:.4g}"))
            if "HR" in selected:
                for row in selected.loc[selected.term.eq("H_x_W")].to_dict("records"):
                    forest.append({"analysis": name, **row})
    if forest:
        f = pd.DataFrame(forest)
        fig, ax = plt.subplots(figsize=(8, max(3, 0.42*len(f))))
        y = np.arange(len(f))
        ax.errorbar(f.HR, y, xerr=[f.HR-f.HR_lower, f.HR_upper-f.HR], fmt="o", color="#246A80")
        ax.axvline(1, color="#999999", linestyle="--")
        ax.set_yticks(y, f.analysis)
        ax.set_xscale("log")
        ax.set_xlabel("Hcy × WMH interaction: ratio of hazard ratios (95% CI)")
        ax.set_title(banner, fontsize=10)
        save_figure(fig, folder, "interaction_forest")
        sections.append('<h2>交互估计</h2><img src="interaction_forest.png" width="850">')
    update = root / "03_month3_update/risk_update_M0_M1.csv"
    if update.is_file():
        u = pd.read_csv(update)
        sections.append("<h2>H3：加入三个月Hcy前后的模型风险</h2><p>M0与M1使用相同患者和其他协变量；差异为三个月Hcy的两个样条项。下表为拟合风险分布，模型校准与外部预测性能须另行验证。</p>"
                        + u[["M0_fitted_risk", "M1_fitted_risk"]].describe().to_html())
    update_curve = root / "03_month3_update/risk_curves.csv"
    if update_curve.is_file():
        d = pd.read_csv(update_curve)
        fig, ax = plt.subplots(figsize=(8, 5))
        for q, sub in d.groupby("wmh_quantile"):
            ax.plot(sub.hcy, 100*sub.risk, label=f"WMH percentile {q*100:.0f}", linewidth=2)
        ax.set(xlabel="Month-3 Hcy (µmol/L)", ylabel="Post-measurement ischemic stroke incidence (%)")
        ax.set_title(banner, fontsize=10)
        ax.legend(frameon=False)
        save_figure(fig, folder, "H3_updated_risk")
        sections.append('<h2>H3：复测后的风险曲线</h2><img src="H3_updated_risk.png" width="850">')
    sections.append("<h2>模型审核</h2><p>按各假设的diagnostics.json核验收敛、比例风险及比例优势假设；"
                    "按adjustment_certificate.json追溯变量选择。插补假设为MAR，并用完整病例结果评估稳定性。</p>")
    target = folder / "report.html"
    target.write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>Hcy WMH</title>'
                      '<style>body{font:16px sans-serif;max-width:1200px;margin:40px auto;color:#20303a}'
                      'table{border-collapse:collapse;font-size:13px}td,th{padding:6px;border:1px solid #ddd}'
                      'h2{margin-top:36px}img{max-width:100%}</style>' + "\n".join(sections) + '</html>',
                      encoding="utf-8")
    (folder / "SUMMARY.md").write_text(f"# {banner}\n\n主终点：y1_is + y1_is_dd。\n\n"
        + "\n".join(f"- {r['analysis']}: {r['status']}" for r in rows)
        + "\n\n详见hypothesis_summary.csv、report.html与各假设hypothesis_estimate.json。\n", encoding="utf-8")
    record_run(cfg, "report", {"report": str(target)})
    return target
