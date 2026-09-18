# CNSR-III：影像与长期预后四项独立研究

当前研究为：早期功能独立后的远期失能、恢复期血压与WMH、CEC与灰质结构、脑肾微血管损伤。四个队列分别建立，不要求Hcy，不以其他研究的检测完整性限制入组。原编号01、02、05、06保留。

当前为第四版：按研究者决定移除IMG_ICAS调整项和Apo_AI扩展分析，两字段不再进入提取、插补、模型或缺失诊断。此前已取消神经酰胺研究和视觉评分定义的重度WMH子组。血压保留连续WMH交互及冠心病跳答核查；CEC保留同样本同插补的HDL-C调整前后比较；脑肾保留持续白蛋白尿×WMH主要检验。汇总为四项Holm校正。旧版结果保留并标为 `PREVIOUS_VERSION`。详见[第四版修订与验证](docs/four_studies_revision_v4.md)。

补充SAS字段后重新执行审计。随后可用 `wmh-study diagnose --page 1`、`--page 2`、`--page 3` 查看来源、候选名称、冠心病规则补齐和异常计数，方便截图反馈。参见[工作站操作](docs/four_studies_workstation.md)。

```bash
conda activate wmh-hcy
git pull --ff-only
python -m pip install --no-deps -e .
wmh-study audit --study all
# 核对对应队列审计后，例如先运行恢复研究：
wmh-study run --study recovery --through report
wmh-study summary --page 1
```

[四项统计方案](docs/four_studies_plan.md) · [Word方案](docs/CNSRIII_四项独立研究统计分析方案_20260918_v4.docx) · [工作站操作](docs/four_studies_workstation.md) · [本版验证](docs/four_studies_revision_v4.md) · [第三版历史合成演示](examples/studies_demo_v3/README.md)。复用现有多SAS、SuStaIn与精确ID配置；Python 3.11、Conda兼容。新结果写入 `outputs/real/studies/`。旧五项研究文档、示例及Hcy历史结果保留，以下内容是历史方案说明。

---

## 历史研究 Hcy–WMH五年首次缺血性复发分析 v3

这是基于既有WMH分割产物和已提供临床字段的纯Python分析项目。支持Ubuntu或WSL工作站，不依赖R，也不使用中心变量。

**仓库包含代码、方法学材料和合成数据验证。真实患者分析在数据所在工作站执行，患者数据及结果不随仓库发布。** 合成数据事件率、样本量、效应和P值不能用于论文结果。

## 历史研究入口：五年主分析，其他月份敏感性

2026-09-17修订后的核心问题：**基线Hcy与五年首次缺血性卒中复发的关联，是否因WMH负担而不同？** 已查看一年结果，这次修订不称为原先预注册的确认性分析。

```bash
conda activate wmh-hcy
git pull --ff-only
wmh-hcy recurrence --through prepare
python scripts/diagnose_recurrence.py
# 检查本次病例/事件数后：
wmh-hcy recurrence --through report
```

同一份五年终点记录截断为3、6、12、24、36、48、60个月。各时点共享核心插补和指标尺度；不拟合死亡模型、不计算绝对风险、不调用bootstrap。结果写入`outputs/real/recurrence_v3/`，旧结果保留。

[新版统计分析计划](docs/recurrence_v3_plan.md) · [新版工作站步骤](docs/recurrence_v3_workstation.md)。新分析只需本版字段，三个月血液、mRS及T1灰质资料不阻碍运行。以下保留旧版H1–H4的操作与背景，`run`和`longterm`仍属于旧方案。

## 旧版H1–H4：在另一台工作站开始

在Linux或WSL终端执行。已有Conda时可直接使用：

```bash
git clone https://github.com/zhenzonglin/WMH.git
cd WMH
conda env create -f environment.yml
conda activate wmh-hcy
python -m pip check
python -m pytest -q

# 先检查SAS字段和患者数；此时不需要影像目录。
wmh-hcy configure --sas-dir "/data/CNSRIII/SAS"
wmh-hcy audit

# 提供既有SuStaIn路径后，连接影像并准备队列。
wmh-hcy configure --sustain-dir "/data/SuStaIn"
wmh-hcy run --through prepare
```

Conda版不需要安装uv。`environment.yml`通过Conda安装Python和OpenMP运行库，再由环境内的pip安装`requirements-conda.txt`中的固定版本及本项目。审阅病例后，用`wmh-hcy run --through report`启动完整分析。**下文所有`uv run wmh-hcy ...`，在已激活的Conda环境中都可直接写成`wmh-hcy ...`。**

也可以继续使用[uv](https://docs.astral.sh/uv/getting-started/installation/)：

```bash
git clone https://github.com/zhenzonglin/WMH.git
cd WMH
uv python install 3.11
uv sync --frozen
uv run pytest -q

# 第一阶段：只指定SAS目录，检查字段与人数；此时不需要影像目录。
uv run wmh-hcy configure --sas-dir "/data/CNSRIII/SAS"
uv run wmh-hcy audit
```

查看`outputs/real/audit/SUMMARY.md`、`summary.json`与`report.html`。审计区分字段存在、有效测量人数、独立患者数、重复ID和各假设临床交集。

```bash
# 第二阶段：接入既有SuStaIn项目/derivatives目录，检查影像交集并准备队列。
uv run wmh-hcy configure --sustain-dir "/data/SuStaIn"
uv run wmh-hcy run --through prepare

# 第三阶段：审阅病例流程后，明确启动完整统计及报告。
uv run wmh-hcy run --through report
```

`configure`生成`config/workstation.local.yml`，其他命令自动读取它。该配置、原始数据、派生CSV和分析输出均不纳入Git。默认`run`仅运行至`prepare`。如果检查发现问题，会记录原因并停止后续阶段。**完整配置说明和问题处理见[工作站操作手册](docs/workstation_setup.md)。**

## 四个研究假设

先读`docs/statistical_analysis_plan.md`。研究和项目按以下顺序排列：

|目录|假设|主要检验|
|---|---|---|
|`analyses/01_structure`|H1：Hcy 15对10对应更高WMH|非线性对比Δ1|
|`analyses/02_recurrence`|H2：WMH增强Hcy与复发的关联|H×W；唯一主要检验|
|`analyses/03_month3_update`|H3：M3 Hcy在基线信息之外更新风险|M3 Hcy15对10条件logHR|
|`analyses/04_function`|H4：WMH增强Hcy与较差mRS的关联|有序模型H×W|

H2使用双侧α=0.05；H1、H3、H4使用固定三项Holm。B12作为调整因素；协变量按`adjustment.py`的角色规则选择，详见`docs/covariate_rationale.csv`。

### 2—5年补充分析

新增`y2_IS`至`Y5_IS`及对应`_dd`、`m24_mrs`至`m60_mrs`接口。在原Conda环境更新项目后运行：

```bash
wmh-hcy longterm --through prepare
wmh-hcy longterm --through report
```

各年份分别执行H2复发、H3三个月Hcy更新、H4有序mRS及T1补充；结果单独写入`outputs/real/longterm/`，原一年结果保留。长期Cox使用实际事件/删失天数，没有完整长期死亡日期时不计算竞争风险绝对发生概率。四年结果按补充探索性分析报告，每项假设固定四年Holm校正。[统计规则](docs/longterm_analysis_plan.md)及[操作说明](docs/workstation_setup.md)包含字段和截图命令。

不接入患者数据也可以运行演示：

```bash
uv run wmh-hcy demo
```

虚拟环境由项目独立管理，不修改系统Python。

`demo`生成明确标记的合成CSV，运行全部分析模块，输出`outputs/demo/latest_results.json`指向的HTML报告。演示使用2份插补和每份3次bootstrap，仅验证软件流程。真实配置使用50份插补、每份200次bootstrap，计算量明显更大。

## 接入工作站真实数据

先使用`configure`创建工作站配置；需要指定多个来源或高级设置时，编辑`config/workstation.local.yml`：

```yaml
inputs:
  sas_globs:
    - /actual/data/clinical.sas7bdat
    - /actual/data/laboratory.sas7bdat
  derivatives_root: /actual/Substain/derivatives
  id_map_csv: /actual/data/image_clinical_ids.csv
```

ID一致时，将`id_map_csv`留空。需要映射时，文件列为`participant_id,code_n`，双方均按字符串精确匹配。不存在明确映射时，程序不会猜测。

```bash
uv run wmh-hcy audit
uv run wmh-hcy extract
uv run wmh-hcy image-audit
uv run wmh-hcy prepare
uv run wmh-hcy analyse
uv run wmh-hcy report
```

单独运行一个假设：`uv run wmh-hcy analyse --hypothesis H3`。完整运行按H1–H4依次执行；单项运行仍保持固定多重比较家族。`hypothesis_summary.csv`是四条假设的统一结果表。

死亡模型或绝对风险计算失败时，可先运行`wmh-hcy check-fit --hypothesis H3`（也支持H2）。它使用现有prepared队列和配置中的插补份数，独立检查两类Cox模型，不做bootstrap、不覆盖正式结果；摘要保存在`outputs/real/diagnostics/`。备用求解保持模型公式与非惩罚估计目标，不能替代模型假设检查。详见[工作站操作手册](docs/workstation_setup.md)。

分析命令接受`--config config/其他.local.yml`。`audit`发现数据问题时返回退出码2，并保存报告。一个临床变量在多份SAS中出现时，用`variable_sources`指定一个来源文件名；同名文件用绝对路径区分。主终点始终只使用`y1_is`和`y1_is_dd`。

源SAS文件只读。CSV输出保留原编码，标准化步骤另存分析变量。SAS日期按照元数据区分从1960年起的天数或秒数。字段没有日期格式时需要在`sas.date_formats`明确指定，不根据数值大小猜测。

## 输出位置

- `outputs/real/extracted/`：白名单原始CSV、元数据和来源清单。
- `outputs/real/audit/`：SAS字段/人数审计、影像连接与体积可用性计数。影像人工质控标记不参与入组。
- `outputs/real/pipeline_status.json`：本次分阶段运行状态与停止原因。
- `outputs/real/prepared/`：主表、每个分析队列、病例流程及排除原因。
- `outputs/real/results/<UTC时间>/`：本次分析的系数、诊断和风险曲线。
- 每次分析创建独立结果目录，`latest_results.json`指向本次运行。失败记录不会被旧结果掩盖。
- 各分析都有单独状态。模型无法估计时给出原因，不自动删变量、不改终点。

## 必读方法文件

- `docs/statistical_analysis_plan.md`：研究问题、模型、时间和解释规则。
- `docs/data_contract.md`：SAS、日期、ID和影像接口。
- `docs/variable_dictionary.csv`：当前分析使用的38个源字段。
- `docs/covariate_rationale.csv`：协变量角色和文献对应。
- `docs/references.md`、`docs/evidence_matrix.csv`：文献及证据边界。
- `presentations/CNSRIII_Hcy_WMH_假设驱动统计方案_v2.pptx`：28页假设驱动方法学汇报。
- `docs/validation_report.md`：软件验证情况及未验证边界。

## 软件验证

```bash
uv run pytest -q
uv run ruff check src tests
uv run python scripts/validate_simulations.py
```

公共SAS测试样例来自pyreadstat上游，来源、哈希和许可证位于`tests/fixtures/sas/`。这些样例只验证SAS格式读取，不是本研究患者。合成数据生成器和数值验证均可在不访问真实数据的情况下运行。

## 研究解释

主要问题由H2检验，H1/H3/H4补充结构、复测和功能证据。分析报告同时列出效应方向、区间和固定阈值判定。不同假设使用各自的病例清单，H3不会从H2中删除早期复发者。模型假设由诊断文件供研究者审阅。

仓库发布当前v2分析主体与方法学PPT。历史本机归档和生成过程文件不属于运行依赖。

本项目不替代原研究伦理、二次分析授权和研究者对模型及结果的审阅。代码就绪与患者数据就绪分别记录。
