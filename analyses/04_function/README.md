# H4：Hcy–WMH交互与一年功能结局

**研究假设：** 在背景因素、本次卒中负担和卒中前功能相同条件下，WMH越重，Hcy升高所对应的12个月较差mRS共同OR越大。

**零假设：** δHW = 0

**预期方向：** δHW > 0

**模型：** 12个月mRS 0–6有序logistic回归

**检验：** H×W的1自由度Wald检验；效应为共同OR之比

**协变量：** C4；另加ICV及采血时点，T1子样本补充GM119

**检验家族：** H1/H3/H4固定三项Holm

运行：`uv run wmh-hcy analyse --hypothesis H4`。完整研究按H1–H4依次运行。单独运行仍保留原多重比较家族。

输出：`hypothesis_estimate.json`、`coefficients.csv`、病例清单及模型/插补诊断；研究级判定见`hypothesis_summary.csv`。
