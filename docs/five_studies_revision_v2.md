# 五项研究第二版修订与验证

日期：2026-09-17。版本标识：`imaging_five_studies_20260917_v2`。本次修订落实血压连续建模、CEC与HDL-C的条件关联比较，以及脑肾损伤交互三个问题。未接入患者数据。

## 统计方案的具体变化

| 研究 | 第二版主要模型或输出 | 主要检验 |
| --- | --- | --- |
| 恢复期血压 | SBP连续样条及WMH交互；主图为连续SBP的HR曲线、逐点95%区间和血压分布 | 两个SBP样条与WMH交互系数的2自由度联合检验，保持不变 |
| CEC | M1调整HDL-C；新增M0仅去掉HDL-C，使用完全相同患者、完成数据和CEC标准差 | M1中CEC系数，保持不变；M0为敏感性分析 |
| 脑肾 | 四类白蛋白尿、连续WMH及三个类别交互，在依赖和死亡方程中均估计 | 依赖方程中持续升高×WMH的1自由度检验，替代旧版主效应检验 |

血压140 mmHg是计算参照，不是推荐目标。120、130、150相对140的固定点对比保留在补充表，不用于血压分类。连续曲线只在预设观测支持范围内显示，限制绘图不删除入模病例。

CEC的主要问题仍是在HDL-C相近时，CEC是否提供结构关联信息；HDL-C不因单因素P值进入或退出模型。M0和M1的系数变化不解释为中介比例，也不以两模型显著性不同判断效应差异。

脑肾研究保留三类别交互的3自由度整体检验作为次要分析，连续恢复期UACR×WMH并调整基线UACR作为次要分析；旧版不含交互的条件关联模型也保留为次要。依赖概率和概率差由同一多分类模型标准化，使用完整跨方程协方差及共享系数的协方差，不用bootstrap。相对概率尺度的交互与绝对概率差属于不同统计量。

## 运行与文件

工作站更新后先重新审计：

```bash
conda activate wmh-hcy
git pull --ff-only
python -m pip install --no-deps -e .
wmh-study audit --study all
```

确认对应审计后运行：

```bash
wmh-study run --study bp --through report
wmh-study run --study cec --through report
wmh-study run --study kidney --through report
wmh-study summary --page 1
```

需要完整的五项新版汇总时，研究recovery和ceramide也按现有入口重新运行。旧结果物理保留；旧合同显示`PREVIOUS_VERSION`，其P值不混入第二版固定五项Holm校正。未完成的新版研究仍占据预设检验家族位置，并显示未估计。

每个研究的新运行目录下重点查看：

- `02_bp/.../primary_result.png`、`primary/continuous_sbp.csv`、`primary/sbp_distribution.csv`、`primary/clinical_contrasts.csv`。
- `05_cec/.../cec_hdl_comparison.csv`、`primary_result.png`；M0插补诊断记录`reused_from: primary`。
- `06_kidney/.../primary/kidney_interaction_curves.csv`、`primary_result.png`；主要项为`dependent:albuminuria_3_x_wmh_ml`。

## 软件验证

全部验证使用合成数据。Python 3.11下143项测试通过；ruff和依赖一致性检查通过。新增测试覆盖预设白蛋白尿参照及三个交互、主要1自由度检验、概率差的数值梯度与完整协方差、血压连续曲线参照与观测支持、CEC共用插补和患者顺序，以及旧合同结果隔离。

五模块合成演示使用2400条生成记录、2份插补和2次迭代，完成45项估计：恢复10项、血压6项、神经酰胺9项、CEC9项、脑肾11项。演示参数不改变真实配置的50份插补和10次迭代。汇总与主要表图见[第二版合成示例](../examples/studies_demo_v2/summary.html)。示例只发布汇总，没有患者或模拟个体行。

新增脑肾交互Monte Carlo模拟采用每场景60次、每次1200条生成记录，零交互和已知交互共120次，均完成估计。完整数值见[第二版模拟表](five_studies_v2_simulation_validation.csv)。

| 依赖方程持续升高×WMH真系数 | 平均偏差 | 偏差MC标准误 | 95%区间覆盖率及Monte Carlo区间 | 拒绝率及Monte Carlo区间 |
| --- | ---: | ---: | --- | --- |
| 0 | -0.0206 | 0.0228 | 98.3%（91.1%–99.7%） | 1.7%（0.3%–8.9%） |
| 0.4 | 0.0121 | 0.0281 | 95.0%（86.3%–98.3%） | 46.7%（34.6%–59.1%） |

上述已知效应场景的检出率只有46.7%，不是实际队列的功效计算；模拟不保证真实数据模型正确或交互有足够精度。第一版的600次模拟及近似插补偏差记录继续保留在[第一版验证记录](five_studies_validation.md)，不重新标记为本轮执行。

复现本轮模拟：

```bash
python scripts/validate_studies.py --family kidney_interaction --repetitions 60 --n 1200 --output outputs/synthetic/kidney_validation_new
wmh-study demo --output outputs/synthetic/studies_revision_new --n 2400 --through report
```

修订Word共15页，经渲染逐页检查；五项合成主要图均已检查。此次验证环境为本机Windows Python 3.11，未在数据所在工作站运行真实SAS或重建Conda环境。真实例数、事件数、拟合稳定性及研究结果由工作站新运行产生。
