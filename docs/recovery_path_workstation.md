# 早期功能独立后的远期失能：全部卒中复发路径扩展

这是查看原五年功能结果之后提出的**探索性修订**。原recovery运行、结果指针和报告均保留。血压、CEC和脑肾代码不删除，但默认不执行。本扩展不生成Word或HTML。

## 一条命令

本地代码修订同步到工作站可拉取的仓库之后，在有原始SAS、影像和`config/workstation.local.yml`的工作站运行：

```bash
cd /data/usersdir/linzhenzong/WMH
conda activate wmh-hcy
git pull --ff-only
python -m pip install --no-deps -e .
wmh-study recovery-path --through analyse
```

只准备与审计时用`wmh-study recovery-path --through prepare`。两条命令都创建新的时间戳目录；分析命令会自己重新审计和准备，无需先运行prepare。实际患者人数和估计仅由工作站生成。

## 固定定义

- 基础队列：成人缺血性卒中、3个月存活且mRS 0–2，校正后全脑WMH、GM119、真实ICV和急性梗死体积可用。无需Hcy、CEC、UACR。
- 复发路径队列：实际`F3_DATE`进入风险集；入组前或当天首次复发者不进入。`Y5_STROKE/Y5_STROKE_DD`提供**任何类型卒中首次复发**及首次事件或有效观察终止日；右删失保留实际结束日，行政上限1825天。
- 2–4年累计记录只核查首次事件一致性，不要求每人所有年份字段齐全。死亡是复发模型的退出风险集事件；Cox的HR不是累计绝对复发概率。
- Cox同时调整WMH和GM119及原recovery的固定协变量。WMH为`log(1+mL)`三节点样条，GM119保留原mL每1 SD，年龄三节点样条。输出WMH原始mL轴上相对中位数的HR曲线、完整协方差95%区间及观测5–95百分位支持范围。
- 原五年三状态模型不覆盖。本扩展对WMH、GM119增加样条形状检验及GM119-log形式作为探索性敏感性分析；不按P值替换原主模型。

## 中介与预测边界

`mediation_gate.json`先核查年度生死状态、无明确死亡解释的提前失访、实际入组日和死亡—复发次序。中途mRS缺失但死亡字段明确存活时，仍能用于生死区间；五年依赖状态缺失者不补为独立，存活者结局可观测概率作为权重。若时间信息不足，`mediation.json`标记`NOT_ESTIMABLE`并列出原因，仍保留影像—复发和影像—五年功能关联。通过后才计算按随访区间的探索性复发过程间接概率差：WMH第75相对25百分位、GM119第25相对75百分位，分别给出独立、依赖、死亡及依赖/死亡合并结果。协变量50份插补；默认100次患者重抽样，每次重新插补后给出百分位区间。影像不是随机分配，不能把此分解写为已证明的因果中介，也不报告“中介比例”。长期死亡日期不造假。

`prediction/`只用3个月时已知变量，分别校验五年依赖概率（死亡单列）与依赖或死亡概率。比较仅WMH与完整临床模型，保存校准、AUC、Brier与重抽样乐观偏差。预测主拟合使用全部插补集；自助验证目前条件于第一份完成数据，因此其乐观偏差不包含插补不确定性。没有指定额外随访行动、误判代价或外部验证队列，本轮**不提供可执行WMH cutoff或个体风险阈值**。

## 输出与状态

新运行位于`outputs/real/studies/01_recovery_path/runs/<时间戳>/`。重点看`status.json`、`audit.json`、`base_flow.csv`、`event_flow.csv`、`stroke_recurrence/`、`wmh_recurrence_curve.csv/.png`、`functional_shape/`、`functional_shape_tests.json`、`functional_gm_log/`、`mediation_gate.json`、`mediation.json`及`prediction/validation.json`。`PARTIAL`表示某项必需补充分析未估计；每个模型目录的`failure.txt`和`result.json`给出原因。`NOT_ESTIMABLE`的中介是预设可行性门控，不等同于零效应。

患者级CSV留在工作站的忽略目录，不进入Git。不得把合成演示人数、效应或P值写为真实研究结果。
