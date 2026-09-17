"""Reports for the separate, supplementary multi-year analysis."""

from __future__ import annotations

import html
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator, NullLocator

from .common import read_json


def report_longterm(root: Path) -> Path:
    summary = pd.read_csv(root / "longterm_summary.csv")
    counts = pd.read_csv(root / "cohort_counts.csv")
    audit = pd.read_csv(root / "field_audit.csv")
    status = read_json(root / "status.json")
    figures = root / "figures"
    figures.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2), sharey=True)
    names = {
        "H2": "Baseline Hcy x WMH (HR ratio)",
        "H3": "Month-3 Hcy 15 vs 10 (HR)",
        "H4": "Hcy x WMH on mRS (OR ratio)",
        "H4_T1": "mRS with T1 (OR ratio)",
    }
    for ax, (family, title) in zip(axes, names.items(), strict=True):
        rows = summary.loc[summary.family.eq(family)]
        for r in rows.itertuples():
            usable = np.isfinite([r.ratio, r.ratio_lower, r.ratio_upper]).all()
            if usable and 0 < r.ratio_lower <= r.ratio <= r.ratio_upper:
                ax.errorbar(
                    r.ratio,
                    r.year,
                    xerr=[[r.ratio - r.ratio_lower], [r.ratio_upper - r.ratio]],
                    fmt="o",
                    color="#17657a",
                    capsize=3,
                )
            else:
                ax.text(0.02, (r.year - 1.5) / 4, "not estimated", transform=ax.transAxes, fontsize=8)
        ax.axvline(1, color="gray", linestyle="--")
        ax.set_xscale("log")
        low, high = ax.get_xlim()
        ticks = MaxNLocator(nbins=4).tick_values(low, high)
        ax.xaxis.set_major_locator(FixedLocator([v for v in ticks if low <= v <= high and v > 0]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:g}"))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_title(title, fontsize=10)
        ax.set_ylim(1.5, 5.5)
        ax.set_yticks([2, 3, 4, 5])
        ax.set_xlabel("Estimate and pointwise 95% CI")
    axes[0].set_ylabel("Cumulative follow-up year")
    fig.suptitle("SUPPLEMENTARY OBSERVATIONAL ANALYSES | " + str(status.get("mode", "unknown")))
    fig.tight_layout()
    fig.savefig(figures / "longterm_forest.png", dpi=160)
    plt.close(fig)
    sections = [
        '<meta charset="utf-8"><title>Hcy–WMH长期随访补充分析</title>',
        (
            "<style>body{font:16px sans-serif;margin:32px;color:#173245}table{border-collapse:collapse;"
            "font-size:13px}td,th{padding:6px;border:1px solid #ccc}img{max-width:100%}</style>"
        ),
        "<h1>Hcy–WMH：2—5年补充分析</h1>",
        (
            f"<p>运行模式：{html.escape(str(status.get('mode')))}；"
            f"状态：{html.escape(str(status.get('status')))}</p>"
        ),
        (
            "<p>本扩展在一年结果已知后增加，属于探索性补充分析。各年使用累计首次缺血性事件，"
            "病例与事件有重叠，不是四次独立验证。每项假设固定四年Holm校正；置信区间为逐项95%区间。</p>"
        ),
        (
            "<p>Cox按IS_DD记录的事件或删失时间、实际采血日及各年行政截止分析；没有把未复发者一律记为完整随访。"
            "当前没有完整长期死亡日期，未计算死亡竞争风险绝对发生概率；mRS的6分不能用来推算死亡日期。</p>"
        ),
        "<h2>字段与随访时间审计</h2>",
        audit.to_html(index=False, escape=True),
        "<h2>各分析病例及事件数</h2>",
        counts.to_html(index=False, escape=True),
        "<h2>全部年份的假设估计</h2>",
        summary.to_html(index=False, escape=True),
        (
            "<p>H2、H4及H4_T1的ratio为交互比值，不能解释为高Hcy组相对低Hcy组的整体HR/OR。"
            "H3为三个月Hcy 15对10的条件HR。不同年份使用各自队列的WMH标准差，转换参数保存在design.json。</p>"
        ),
        '<img src="figures/longterm_forest.png" alt="长期随访森林图">',
        "<h2>模型状态与诊断入口</h2>",
        (
            "<p>每个year*/H*目录中包含flow.csv、coefficients.csv、imputation.json、diagnostics.json。"
            "收敛、比例风险和mRS阈值检查需与估计精度一起审阅。缺失结局与失访不按无事件处理。</p>"
        ),
    ]
    models = [{"analysis": name, **row} for name, row in status.get("analyses", {}).items()]
    sections.append(pd.DataFrame(models).to_html(index=False, escape=True))
    target = root / "report.html"
    target.write_text("\n".join(sections), encoding="utf-8")
    return target
