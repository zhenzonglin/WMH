# CNSR-III：Hcy–WMH四假设分析项目 v2

这是基于既有WMH分割产物和已提供临床字段的纯Python分析项目。支持Ubuntu或WSL工作站，不依赖R，也不使用中心变量。

**当前交付包含代码、方法学材料和合成数据验证。尚未读取或分析患者数据。** 合成数据事件率、样本量、效应和P值不能用于论文结果。

## 在另一台工作站开始

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
