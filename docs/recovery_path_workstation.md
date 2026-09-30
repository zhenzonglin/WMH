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

## 按页查看并截图

需要同时汇报原recovery主分析和recovery-path补充分析时，使用组合查看脚本，并明确指定补充分析目录。例如本次研究者提供的新运行：

```bash
python scripts/recovery_results_review.py \
  --path-run outputs/real/studies/01_recovery_path/runs/20260929T081847414216Z
```

组合脚本通常显示12页：前6页为recovery主分析的病例与基线、主检验与各月份估计、敏感性估计、插补/权重诊断及标准化状态概率；后6页为明确指定的recovery-path运行。每页打印独立的研究名称、运行编号和绝对目录。主分析默认扫描`01_recovery/runs/`，选取最新有模型保存记录或已进入分析阶段的运行；新的仅审计/准备目录不会替代已分析目录。可以用`--recovery-run <确切主分析目录>`固定主分析。若同次运行的队列、审计与主模型人数不匹配，则保留错误并停止展示该主分析数值，补充分析查看仍可继续。没有远程通道时，本机不能据此核验工作站上的实际数值，须将工作站输出截图反馈。

只查看扩展时可运行`python scripts/recovery_path_review.py --run <确切运行目录>`。省略`--run`时直接选取磁盘上最新时间戳目录，不依赖可能滞后的结果指针。它分6页显示状态与人数、复发关联、功能模型形状、基线与插补诊断、中介、预测验证。两种脚本均支持交互终端中截图后按回车继续、`--page 2`单独查看某页、`--no-pause`一次输出全部页。

这是只读查看，不重新拟合、插补或修改结果，不打印患者ID或逐人数据，不加载血压、CEC或脑肾结果。WMH四分位对比及非线性检验由已保存的系数、完整协方差和固定尺度计算。脚本不会用旧成功结果替换最新失败运行，也不会把缺失文件当作已生成。组合查看列出3张图片的实际存在状态和绝对路径：原recovery的`primary_result.png`，以及指定recovery-path的`wmh_recurrence_curve.png`和`prediction/apparent_calibration.png`。图片可另外截图。校准图为表观校准，乐观偏差修正指标见文字页。原`report.html`只打印路径，本次不会重新生成报告。
