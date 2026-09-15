# H2：WMH修饰Hcy与复发的关联

**研究假设：** 基线WMH越重，Hcy升高所对应的一年缺血性复发HR越大，即存在正向交互。

**零假设：** βHW = 0

**预期方向：** βHW > 0

**模型：** 延迟入组原因别Cox；分别建模缺血性卒中和死亡

**时间：** 以发病日为时间尺度，从实际基线采血日进入风险集；不使用发病到MRI的天数。

**检验：** H×W的1自由度Wald检验；补充WMH分位点的15对10绝对风险差

**协变量：** C0 + ICV + 采血距发病天数

**检验家族：** 唯一主要假设；双侧α=0.05

运行：`uv run wmh-hcy analyse --hypothesis H2`。完整研究按H1–H4依次运行。单独运行仍保留原多重比较家族。

输出：`hypothesis_estimate.json`、`coefficients.csv`、病例清单及模型/插补诊断；研究级判定见`hypothesis_summary.csv`。
