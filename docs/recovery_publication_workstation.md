# Recovery 手稿图表的独立导出

此补丁读取明确指定的 recovery 保存结果和插补前队列，生成手稿图表。按本次确认的图形范围，展示状态概率、五年主分析与相关补充分析的RRR，并可选展示纵向轨迹；不再生成病例流程图或OR图。表格的既有范围和统计口径保留。

导出不调用audit、prepare、analyse或既有研究入口，不重新插补、不重新拟合，也不改写历史运行、结果指针或手稿。真实结果以工作站保存文件为准，不从截图重建数值。

## 手稿编号与输入

| 手稿位置 | 输出 | 来源及口径 |
|---|---|---|
| 五年状态概率 | `<style>/figure2.*` | `primary/standardized_states.csv`；按图形分别展示GM中位数场景、3×3影像网格或九个场景的三状态组成，不重新预测 |
| 合并RRR图 | `effects/<view>/combined_rrr.*` | 五年主分析、36个月、48个月及五年可观测概率加权；分主分析、次要分析和敏感性分析三组 |
| Table 1 | `tables/table1.*` | `eligible.csv`中`state60`已知者，插补前的主要基线特征 |
| Supplementary Table 1 | `tables/supplement_table1.*` | 同一人群，完整基线特征 |
| Table 2 | `tables/table2.*` | 同次运行`primary/coefficients.csv`的四项影像系数，RRR和95%CI |
| Supplementary Table 2 | `tables/supplement_table2.*` | 12/24/36/48个月已保存估计；未估计者保留状态；仅使用已保存的BH联合检验，不重新计算P值 |
| Supplementary Table 3 | `tables/supplement_table3.*` | 完整病例、实际改善、原始WMH、合并结局、可观测概率加权；RRR和二分类OR分别标注 |
| 概率图的区间伴随表 | `tables/display_probability_intervals.*` | WMH与GM各三个取值组成九个场景，每个场景三种状态，共27项保存概率及95%CI |

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

绘图使用Seaborn 0.13.2和SciencePlots 2.2.2，采用`science`与`no-latex`样式，直接由Matplotlib导出文件；无需安装LaTeX。当前增量依赖已替换旧的绘图导出依赖，请使用本分支的`requirements-publication.txt`更新工作站环境。

## 先输出全部真实结果预览

```bash
RUN=outputs/real/studies/01_recovery/runs/20260918T115107666042Z
OUT=outputs/publication/recovery/review_$(date -u +%Y%m%dT%H%M%SZ)
python scripts/recovery_publication.py \
  --recovery-run "$RUN" \
  --style all --effect-view all --preview-only \
  --out "$OUT"
printf '预览目录：%s\n' "$OUT"
```

查看`style_gallery.html`。页面按“状态概率”“主分析与补充分析的RRR”“可选纵向轨迹”组织；每个问题下比较不同图形表达。预览保留各图的设计宽度，不把单栏小图拉伸到双栏图的宽度。不同版本不重新估计数据，也不依据P值筛选指标。

概率图由`--style`控制，默认`all`：

| 参数 | 图形表达 |
|---|---|
| `--style estimates` | GM固定中位数时，WMH第25/50/75百分位的五年依赖概率及95%CI，共三个点估计 |
| `--style heatmap` | WMH与GM第25/50/75百分位的3×3网格，单元格只标注九个场景的五年依赖概率点估计；95%CI见伴随表 |
| `--style composition` | 按GM三个取值分面，每个面板展示WMH三个取值的独立、依赖和死亡组成；采用低饱和度状态配色，无纹理填充，减少重复数值；九个场景的27项概率及95%CI见伴随表 |
| `--style all` | 同时输出以上三种图形 |

RRR图由`--effect-view`控制，默认`all`：

| 参数 | 图形表达 |
|---|---|
| `--effect-view forest` | 对数刻度森林图；参照线为RRR=1，点和区间展示原始RRR及95%CI，不另设重复数值栏 |
| `--effect-view heatmap` | 颜色编码log(RRR)，中性点对应RRR=1；单元格只标注原始RRR点估计，完整95%CI保留在伴随效应表 |
| `--effect-view diverging` | 发散条形长度为100×(RRR−1)，95%CI同步按100×(CI−1)转换并保留；原始RRR及区间见伴随效应表 |
| `--effect-view all` | 同时输出以上三种图形 |

发散条形表示相对比值的百分比变化，不是结局概率的百分点差。例如RRR=1.40对应比值增加40%，不能写成依赖概率增加40个百分点。热图和发散条形仅改变展示尺度，不改变合并估计或区间。

`combined_rrr`的森林图和发散条形图用两个面板分别展示依赖和死亡结局，每个面板按WMH、GM分组。各影像组内先列五年主分析，再列36/48个月次要分析，最后列五年可观测概率加权敏感性分析。热图则以四项分析为行、结局与影像指标的四种组合为列。主分析标签以同字号加粗，并在森林图中使用方形点。三种图均包括上述四项分析的全部16项RRR，不以显著性取舍；精确数值和完整95%CI保留在原Table 2、Supplementary Table 2和3中。完整病例、实际改善、原始WMH及12/24个月不进入这张图；不另外生成加权分析图，也不生成二分类OR图。表格内容保持不变。

概率图只展示保存的取值，不插值为连续风险曲线，不增设cutoff。组成图各状态概率总和为1；各状态的点态区间单独给出，不把各区间端点相加当作联合置信区间。`--preview-only`省略大尺寸TIFF，但保留PDF、SVG和PNG。脚本拒绝覆盖非空输出目录；每次使用新目录。

三种概率图的展示范围不同：点区间图聚焦GM中位数下的三个依赖概率；热图扩展至九个影像组合；组成图显示这些组合的全部三种状态。重叠场景必须使用相同保存数值，不能把新增面板解释为新增统计分析。原有五张表加区间伴随表，共六张；加入轨迹时增加两张表，共八张。

旧版`--style classic`、`colour`和`compact`已停用，调用会报错并列出可用图型。它们原本只是视觉样式；应按本次希望表达的统计问题，明确改用`estimates`、`heatmap`、`composition`或`all`，不要在旧输出目录中混放新版本。

## 选定样式后导出投稿素材

以概率点区间图和RRR森林图为例，删除`--preview-only`即可生成1200 dpi TIFF：

```bash
python scripts/recovery_publication.py \
  --recovery-run outputs/real/studies/01_recovery/runs/20260918T115107666042Z \
  --style estimates --effect-view forest \
  --out outputs/publication/recovery/estimates_forest_$(date -u +%Y%m%dT%H%M%SZ)
```

每幅图输出：矢量PDF、保留文字的SVG、300 dpi PNG和单层1200 dpi TIFF。表格输出可编辑Word三线表、CSV、booktabs LaTeX和HTML预览。三线表只有顶线、表头底线、底线，无竖线或内部网格；长表重复表头，数据行不跨页拆分。

正式图内文字统一使用9 pt Arial，仅保留辨认面板、坐标和图例所需文字；不加图内总标题或解释性页脚。没有Arial的工作站依次使用Liberation Sans、Helvetica或DejaVu Sans，实际选择的字体记入`manifest.json`。完整定义、参照组、区间方法和解释限制另存于`figure_legends.txt`，用于手稿图注。主分析以同字号加粗突出，不通过放大字体改变层级。演示输出仍保留合成数据标记。

不同内容使用不同的物理尺寸：独立的五年概率点图为85×68 mm，五年概率热图为110×90 mm，复杂多面板图使用178 mm宽度。热图简化为点估计标注时，95%CI仍完整保留在伴随表中，不删除、不替换，也不重新计算区间。

主要输出路径：

```text
<out>/style_gallery.html
<out>/estimates/figure2.*
<out>/heatmap/figure2.*
<out>/composition/figure2.*
<out>/effects/forest/combined_rrr.*
<out>/effects/heatmap/combined_rrr.*
<out>/effects/diverging/combined_rrr.*
<out>/tables/table1.*
<out>/tables/table2.*
<out>/tables/supplement_table1.*
<out>/tables/supplement_table2.*
<out>/tables/supplement_table3.*
<out>/tables/display_probability_intervals.*
```

只指定一个`--style`或`--effect-view`时，只生成相应图形子目录；两组选项相互独立。

## 可选轨迹补充图表

当前手稿的轨迹段落尚无黄色图表编号，因此默认不加入。若需要，显式追加`--trajectory-run`：

```bash
python scripts/recovery_publication.py \
  --recovery-run outputs/real/studies/01_recovery/runs/20260918T115107666042Z \
  --trajectory-run outputs/real/studies/01_recovery_trajectory/runs/20261011T171844379322Z \
  --style all --effect-view all --preview-only \
  --out outputs/publication/recovery/with_trajectory_$(date -u +%Y%m%dT%H%M%SZ)
```

额外生成`<style>/optional_figure_s1.*`（WMH/GMV与三状态概率轨迹）、`effects/<view>/trajectory_rrr.*`（四类转移效应）、Table S4（转移效应）和Table S5（36/48/60个月概率及CI）。可选图表的手稿编号需在正文中确认后引用。

轨迹概率同样提供三类表达。`estimates`使用Seaborn访视折线，并叠加保存的点态95%CI带；`heatmap`按两个影像指标×三种状态形成六个面板，每个面板为三个影像分位值×三个随访月份，单元格只标注概率点估计；`composition`用四个面板展示低/高WMH和低/高GM场景，每个场景展示36/48/60个月的三状态组成，另一影像指标固定中位数。完整区间保留在Table S5中。轨迹RRR另提供`forest`、`heatmap`、`diverging`三类表达。

轨迹折线调用`seaborn.lineplot`时指定`estimator=None`和`errorbar=None`，禁止Seaborn重新汇总或生成区间。CI带仅由已保存的lower/upper连接，表示各访视点态区间的可视化；不代表估计了连续时间轨迹，也不是同时置信带。

轨迹图只展示保存的36/48/60个月；若使用连接线，线条不代表估计了中间月份。转移模型的参照是“维持起始状态”，指数系数为P(目标状态)/P(维持原状态)之比的每SD变化，即多分类RRR，不能写成HR或简单概率比。独立起点的参照为保持独立，依赖起点的参照为保持依赖。因此轨迹RRR单独成图，不与参照为五年独立状态的`combined_rrr`混排。

已保存轨迹概率CI使用每份插补内参数模拟及Rubin合并。原实现默认每份仅16次参数抽样；导出时记录实际抽样次数并提示稳定性复核，不通过绘图改算CI。工作站原CSV提供完整精度，优先于历史报告或截图中的转录值。

## 教育类别标签

核对原始数据字典后，在Git忽略的本地文件中保存5个编码的英文标签，例如`config/education_labels.local.yml.json`，内容为JSON对象、键为字符串`1`至`5`。不要按编码猜测学历含义。

```bash
python scripts/recovery_publication.py \
  --recovery-run "$RUN" --style estimates --effect-view forest \
  --education-labels config/education_labels.local.yml.json \
  --out outputs/publication/recovery/labelled_$(date -u +%Y%m%dT%H%M%SZ)
```

## 本地样式演示及验证

```bash
python scripts/recovery_publication.py --demo --style all --effect-view all --preview-only \
  --out outputs/publication/recovery/demo_$(date -u +%Y%m%dT%H%M%SZ)
ruff check src tests
python -m pytest -q
git diff --check
```

演示中的人数、估计、CI和P值均为合成样式数据，所有图和表均作标记，不能用于论文。真实模式不接受合成输入；演示模式不接受真实运行。演示也必须写入新目录，以免与旧版图形混淆。

验证范围包括：TableOne缺失分母；不同图形中重叠场景的原始估计和区间一致；发散条形的尺度及区间转换；概率组成及27项伴随CI；合并RRR仅含四项指定分析、主分析在首位；无病例流程图、OR图或单独加权图；未知结局保留；运行身份和概率网格；输入只读；输出不覆盖；死亡吸收；三线表XML；TIFF分辨率及数据不泄漏。测试验证软件契约，不代替研究结果审阅。缺少可选依赖时，专门的制图/制表测试会跳过；安装本文件的增量依赖后应完整执行。

输出`manifest.json`保存输入/输出SHA256、运行编号、依赖版本和校验说明；`source_data/`仅保存支持图表的汇总数值，绝不导出患者ID或逐人记录。`status.json=COMPLETED`表示图表导出完成，不表示统计模型假设、医学因果解释或投稿审查已通过。原始SAS、临床配置、患者CSV、图像和输出不进入Git。

## BMJ依据与需作者核对的内容

2026-10-15查阅的[BMJ Author Hub](https://authors.bmj.com/writing-and-formatting/formatting-your-paper/)要求图独立提交；非矢量图至少300 dpi、线图1200 dpi；TIFF为单层，多面板合成一个文件；表格为可编辑Word。另参阅[The BMJ submission](https://www.bmj.com/about-bmj/resources-authors/article-submission?change_country=1)和[house style](https://www.bmj.com/about-bmj/resources-authors/house-style)。按内容选择85、110或178 mm图宽、统一9 pt字体和配色是本项目设计选择，不冒称为BMJ唯一强制模板。最终以目标期刊及文章类型要求为准。

TableOne方法参考[官方API](https://tableone.readthedocs.io/en/latest/tableone.html)及Pollard等（JAMIA Open 2018，[doi:10.1093/jamiaopen/ooy012](https://doi.org/10.1093/jamiaopen/ooy012)）。

手稿内容仍需核对：Python 3.10的表述与本项目3.11环境不同；12/24个月不收敛的具体原因不能仅凭总事件数断定；轨迹效应应明确上述参照分母；教育编码标签尚待字典确认。本补丁不修改这些正文，也不增补未保存的统计结果。
