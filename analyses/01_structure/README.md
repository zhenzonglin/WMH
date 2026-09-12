# H1：Hcy与背景白质损伤

**研究假设：** 在C0及测量因素相同条件下，基线Hcy为15 μmol/L者的WMH负担高于10 μmol/L者。

**零假设：** Δ1 = 0

**预期方向：** Δ1 > 0

**模型：** HC3稳健标准误多变量线性回归

**检验：** fH(log2(15)) - fH(log2(10))

**协变量：** C0 + ICV + 采血距发病天数

**检验家族：** H1/H3/H4固定三项Holm

运行：`uv run wmh-hcy analyse --hypothesis H1`。完整研究按H1–H4依次运行。单独运行仍保留原多重比较家族。

输出：`hypothesis_estimate.json`、`coefficients.csv`、病例清单及模型/插补诊断；研究级判定见`hypothesis_summary.csv`。
