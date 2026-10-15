# Recovery 手稿图表的独立导出

此补丁对应 2026-10-15 提供的 Manuscript.docx 中黄色标记的图表。只读取保存的 recovery 结果和插补前队列，生成图表；不调用 audit、prepare、analyse 或既有研究入口，不重新插补、不重新拟合，也不改写历史运行或结果指针。原手稿不修改。

## 手稿编号与输入

| 手稿位置 | 输出 | 来源及口径 |
|---|---|---|
| Figure 1 | `figure1.*` | `cohort_flow.csv` 和 `audit.json`，顺序排除人数，不把未满足生存条件者全称为死亡 |
| Figure 2 | `figure2.*` | `primary/standardized_states.csv`，GM119固定中位数时，WMH第25/50/75百分位对应的五年依赖概率及95%CI |
| Table 1 | `tables/table1.*` | `eligible.csv`中`state60`已知者，插补前的主要基线特征 |
| Supplementary Table 1 | `tables/supplement_table1.*` | 同一人群，完整基线特征 |
| Table 2 | `tables/table2.*` | 同次运行`primary/coefficients.csv`的四项影像系数，RRR和95%CI |
| Supplementary Table 2 | `tables/supplement_table2.*` | 12/24/36/48个月已保存估计；未估计者保留状态；仅使用已保存的BH联合检验，不重新计算P值 |
| Supplementary Table 3 | `tables/supplement_table3.*` | 完整病例、实际改善、原始WMH、合并结局、可观测概率加权；RRR和二分类OR分别标注 |

主基线表采用年龄、性别、血管病史、卒中前和3个月mRS、NIHSS、WMH、GM119、ICV、急性病灶体积；补充表增加吸烟、饮酒、教育和TOAST。选择按手稿用途固定，不按P值筛选。两张基线表调用真正的 Python `tableone.TableOne`；其余效应表直接排版保存的合并估计，不能用TableOne重新计算回归结果。

基线统计使用中位数［Q1，Q3］或n（%），按总体及五年独立、依赖、死亡分层。分类百分比分母为该列该变量非缺失人数，另列各组缺失人数；不自动增加基线P值。编码依据`docs/variable_dictionary.csv`；教育程度未在该字典中解释，默认显示中性类别，投稿前须补充核实标签。BMI、房颤等未保存于该recovery队列的变量不会从其他课题拼接。

## 工作站安装

在既有仓库中操作；如Git提示存在未提交改动，先保留并检查，不运行reset/clean。

```bash
cd /data/usersdir/linzhenzong/WMH
conda activate wmh-hcy
git fetch origin
git switch feat/recovery-publication-figures
git pull --ff-only origin feat/recovery-publication-figures
python -m pip install -c requirements-conda.txt -r requirements-publication.txt
python scripts/recovery_publication.py --help
```

如果首次切换提示找不到本地分支：

```bash
git switch --track origin/feat/recovery-publication-figures
```

项目统计环境为Python 3.11。增量依赖约束保留`requirements-conda.txt`中的版本；如工作站实际安装版本已偏离约束，pip可能调整为约束版本，应先核对环境。无需重新执行`wmh-study`或安装整个项目，脚本直接使用当前checkout的src。

## 先输出三套真实结果预览

```bash
RUN=outputs/real/studies/01_recovery/runs/20260918T115107666042Z
OUT=outputs/publication/recovery/review_$(date -u +%Y%m%dT%H%M%SZ)
python scripts/recovery_publication.py \
  --recovery-run "$RUN" \
  --style all --preview-only \
  --out "$OUT"
printf '预览目录：%s\n' "$OUT"
```

查看`style_gallery.html`，或分别查看`classic/figure2.png`、`colour/figure2.png`和`compact/figure2.png`。

| 参数 | 样式 |
|---|---|
| `--style classic` | 黑白竖向点估计和置信区间 |
| `--style colour` | 蓝色竖向点估计；可选多组图用蓝橙及点型/线型共同区分 |
| `--style compact` | 横向置信区间和独立数值栏；适合精确阅读 |
| `--style all` | 同时输出三套，统计数值完全相同 |

Figure 2仅展示三个已保存取值，不插值为连续风险曲线，不增设cutoff。`--preview-only`省略大尺寸TIFF，但保留PDF、SVG和PNG。脚本拒绝覆盖非空输出目录；每次使用新目录。

## 选定样式后导出投稿素材

以`colour`为例，删除`--preview-only`即可生成1200 dpi TIFF：

```bash
python scripts/recovery_publication.py \
  --recovery-run outputs/real/studies/01_recovery/runs/20260918T115107666042Z \
  --style colour \
  --out outputs/publication/recovery/colour_$(date -u +%Y%m%dT%H%M%SZ)
```

每幅图输出：矢量PDF、保留文字的SVG、300 dpi PNG和单层1200 dpi TIFF。表格输出可编辑Word三线表、CSV、booktabs LaTeX和HTML预览。三线表只有顶线、表头底线、底线，无竖线或内部网格；长表重复表头，数据行不跨页拆分。

## 可选轨迹补充图表

当前手稿的轨迹段落尚无黄色图表编号，因此默认不加入。若需要，显式追加`--trajectory-run`：

```bash
python scripts/recovery_publication.py \
  --recovery-run outputs/real/studies/01_recovery/runs/20260918T115107666042Z \
  --trajectory-run outputs/real/studies/01_recovery_trajectory/runs/20261011T171844379322Z \
  --style all --preview-only \
  --out outputs/publication/recovery/with_trajectory_$(date -u +%Y%m%dT%H%M%SZ)
```

额外生成可选Figure S1（WMH/GMV × 三状态概率轨迹）、Figure S2（四类转移的森林图）、Table S4（转移效应）和Table S5（36/48/60个月概率及CI）。编号标为optional，需在手稿中确认后引用。

轨迹图只连接保存的36/48/60个月；线条不代表估计了中间月份。转移模型的参照是“维持起始状态”，指数系数为P(目标状态)/P(维持原状态)之比的每SD变化，即多分类RRR，不能写成HR或简单概率比。独立起点的参照为保持独立，依赖起点的参照为保持依赖。

已保存轨迹概率CI使用每份插补内参数模拟及Rubin合并。原实现默认每份仅16次参数抽样；导出时记录实际抽样次数并提示稳定性复核，不通过绘图改算CI。工作站原CSV提供完整精度，优先于历史报告或截图中的转录值。

## 教育类别标签

核对原始数据字典后，在Git忽略的本地文件中保存5个编码的英文标签，例如`config/education_labels.local.yml.json`，内容为JSON对象、键为字符串`1`至`5`。不要按编码猜测学历含义。

```bash
python scripts/recovery_publication.py \
  --recovery-run "$RUN" --style colour \
  --education-labels config/education_labels.local.yml.json \
  --out outputs/publication/recovery/labelled_$(date -u +%Y%m%dT%H%M%SZ)
```

## 本地样式演示及验证

```bash
python scripts/recovery_publication.py --demo --style all --preview-only \
  --out outputs/publication/recovery/demo_$(date -u +%Y%m%dT%H%M%SZ)
ruff check src tests
python -m pytest -q
git diff --check
```

演示中的人数、估计、CI和P值均为合成样式数据，所有图和表均作标记，不能用于论文。真实模式不接受合成输入；演示模式不接受真实运行。测试涵盖TableOne缺失分母、数值与区间跨样式一致、未知结局保留、运行身份和概率网格、输入只读、输出不覆盖、死亡吸收、三线表XML、TIFF分辨率和数据不泄漏。缺少可选依赖时，专门的制图/制表测试会跳过；安装本文件的增量依赖后应完整执行。

输出`manifest.json`保存输入/输出SHA256、运行编号、依赖版本和校验说明；`source_data/`仅保存支持图表的汇总数值，绝不导出患者ID或逐人记录。`status.json=COMPLETED`表示图表导出完成，不表示统计模型假设、医学因果解释或投稿审查已通过。原始SAS、临床配置、患者CSV、图像和输出不进入Git。

## BMJ依据与需作者核对的内容

2026-10-15查阅的[BMJ Author Hub](https://authors.bmj.com/writing-and-formatting/formatting-your-paper/)要求图独立提交；非矢量图至少300 dpi、线图1200 dpi；TIFF为单层，多面板合成一个文件；表格为可编辑Word。另参阅[The BMJ submission](https://www.bmj.com/about-bmj/resources-authors/article-submission?change_country=1)和[house style](https://www.bmj.com/about-bmj/resources-authors/house-style)。173 mm宽、8–10 pt无衬线字体和配色是本项目设计选择，不冒称为BMJ唯一强制模板。最终以目标期刊及文章类型要求为准。

TableOne方法参考[官方API](https://tableone.readthedocs.io/en/latest/tableone.html)及Pollard等（JAMIA Open 2018，[doi:10.1093/jamiaopen/ooy012](https://doi.org/10.1093/jamiaopen/ooy012)）。

手稿内容仍需核对：Python 3.10的表述与本项目3.11环境不同；12/24个月不收敛的具体原因不能仅凭总事件数断定；轨迹效应应明确上述参照分母；教育编码标签尚待字典确认。本补丁不修改这些正文，也不增补未保存的统计结果。
