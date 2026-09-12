# H3：三个月Hcy提供后续风险更新

**研究假设：** 在基线Hcy、WMH及C3相同条件下，三个月Hcy为15 μmol/L者的后续缺血性复发风险高于10 μmol/L者。

**零假设：** Δ3 = 0

**预期方向：** Δ3 > 0

**模型：** 实际三个月采血日起的原因别Cox

**检验：** 三个月Hcy非线性15对10的条件logHR；完整协方差Wald检验

**协变量：** C3 + 基线Hcy RCS + 基线Hcy×WMH + 两次采血时点 + ICV

**检验家族：** H1/H3/H4固定三项Holm

运行：`uv run wmh-hcy analyse --hypothesis H3`。完整研究按H1–H4依次运行。单独运行仍保留原多重比较家族。

输出：`hypothesis_estimate.json`、`coefficients.csv`、病例清单及模型/插补诊断；研究级判定见`hypothesis_summary.csv`。
