# 独立纵向功能轨迹补丁

`recovery-trajectory` 是 recovery 的独立次要分析。它从3个月时功能独立者开始，读取基线影像、recovery协变量、重复mRS和死亡字段，拟合3–12、12–24、24–36、36–48、48–60个月相邻访视的离散时间多状态模型；正式汇总节点为36、48和60个月。它不读取卒中复发字段，也不会调用或重写 `01_recovery`、`01_recovery_path`、血压、CEC或脑肾分析。

在工作站的仓库根目录中执行：

```bash
wmh-study recovery-trajectory --through analyse
```

也可以先只准备并审计输入：

```bash
wmh-study recovery-trajectory --through prepare
```

结果写入：

```text
outputs/real/studies/01_recovery_trajectory/runs/<timestamp>/
```

查看已保存运行时使用精确目录：

```bash
python scripts/recovery_trajectory_review.py \
  --run outputs/real/studies/01_recovery_trajectory/runs/<timestamp> \
  --no-pause
```

每次运行会保存状态、输入审计、相邻区间转移计数、转移效应、36/48/60个月标准化状态概率、mRS描述性分布、插补记录、模型诊断及图形。缺失访视不做前向填充；只有相邻两端状态明确的区间进入转移模型。死亡是吸收状态，任何死亡后出现的存活状态会使运行进入 `REVIEW_REQUIRED`，并保留审计原因。

状态为 `PARTIAL` 或 `REVIEW_REQUIRED` 时，应先查看该运行的 `status.json`、`audit.json`、`failure.txt` 和 `fit_diagnostics.json`。补丁不会自动删协变量、替换结局或用既有运行结果补齐不可估计的转移。
