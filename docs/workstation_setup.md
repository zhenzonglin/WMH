# 工作站接入：先核查SAS，再连接WMH，最后启动分析

适用环境：Linux或Windows中的WSL。数据使用`.sas7bdat`文件；代码由Python 3.11运行。Conda与uv两种安装方式使用同一组分析依赖版本，分别通过`requirements-conda.txt`与`uv.lock`安装。工作站无需R、MATLAB、GPU或SuStaIn模型运行环境；本项目读取已经产出的WMH结果。

## 1. 下载项目与安装环境

### 方式A：已有Conda的工作站

在克隆后的仓库根目录执行环境创建，保证环境文件中的相对路径能够找到依赖文件和本项目：

```bash
git clone https://github.com/zhenzonglin/WMH.git
cd WMH
conda env create -f environment.yml
conda activate wmh-hcy
python -m pip check
python -m pytest -q
wmh-hcy --help
```

如果此前已克隆项目，进入该目录执行`git pull --ff-only`后再创建环境。Conda环境名为`wmh-hcy`。安装过程使用Conda中的pip，不需要uv，也不会向base环境安装分析包。`libgomp`提供Linux下LightGBM需要的OpenMP运行库。

后续操作直接执行：

```bash
wmh-hcy configure --sas-dir "/data/CNSRIII/SAS"
wmh-hcy audit
wmh-hcy configure --sustain-dir "/data/SuStaIn"
wmh-hcy run --through prepare
# 审阅病例流程后执行：
wmh-hcy run --through report
```

下文的`uv run wmh-hcy ...`命令，在激活Conda后省略`uv run`即可。无法在批处理脚本中激活环境时，可以使用`conda run --no-capture-output -n wmh-hcy wmh-hcy audit`等对应命令。

环境文件使用Conda官方支持的pip依赖段，分析依赖从既有锁文件导出，避免两种安装方式分别选择不同版本。[Conda环境管理文档](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html)。当前完整安装验证针对Linux/WSL，未验证Windows原生Conda。

### 方式B：使用uv

```bash
# 如系统尚无这些程序，可先安装。
sudo apt-get update
sudo apt-get install -y git curl libgomp1

# 如尚未安装uv，按官方安装入口安装，然后重新打开终端。
curl -LsSf https://astral.sh/uv/install.sh | sh
source "$HOME/.local/bin/env"

git clone https://github.com/zhenzonglin/WMH.git
cd WMH
uv python install 3.11
uv sync --frozen
uv run pytest -q
```

已有uv时跳过安装步骤。[uv官方安装文档](https://docs.astral.sh/uv/getting-started/installation/)。后面的命令都在克隆后的`WMH`目录执行。WSL中的数据路径使用`/data/...`或`/mnt/d/...`等Linux路径。

## 2. 第一阶段：只检查临床SAS

```bash
uv run wmh-hcy configure --sas-dir "/data/CNSRIII/SAS"
uv run wmh-hcy audit
```

脚本递归发现该目录下的SAS文件，读取元数据并扫描38个既定字段的有效值。此阶段不读取影像，不拟合模型，也不插补数据。源文件只读。大数据目录首次扫描可能耗时。

`config/workstation.local.yml`保存本机配置。后续命令自动选用它；更新该配置时会在本地保留备份。不要把真实路径、患者CSV或影像文件提交到仓库。

|审计文件（位于`outputs/real/audit/`）|用途|
|---|---|
|`SUMMARY.md` / `summary.json`|总体状态、来源并集人数、字段缺口、阻断原因|
|`report.html`|可在浏览器中查看的汇总表|
|`files.csv`|各文件元数据行数、扫描行数、独立ID、重复行及空ID数|
|`field_coverage.csv`|38个字段是否存在、所选来源、有效值人数、预期单位/类型|
|`variables_by_file.csv`|每个字段逐文件的有效、缺失、非法编码计数|
|`hypothesis_counts.csv`|H1–H4影像连接前的临床完整观察值交集|
|`file_intersections.csv`|临床文件间共享独立ID数|
|`sas_inventory.json`|原变量名、标签、SAS日期格式、文件编码与来源|

**这些计数不是最终分析样本量。** 重复ID不会自动保留第一行；协变量散在缺失也不会仅因“完整病例计数较少”就被永久剔除。最终病例资格由影像体积可用性、时间规则和各假设的变量要求共同确定。

`READY_FOR_EXTRACTION`表示基础临床字段/编码检查可继续；`REVIEW_REQUIRED`表示需要按报告修正输入配置或数据问题。缺少可选子分析字段会在对应假设报告中列出。没有发现患者缺血事件时，事件时间全部为空本身不阻断提取，但不意味着能拟合复发模型。

2026-09-16起，按研究者决定移除发病到MRI的天数，白名单由39项改为38项。风险集进入日改为实际基线采血日；该MRI间隔不再提取、调整或用于病例筛选。升级后重新运行`wmh-hcy audit`与`wmh-hcy run --through prepare`以更新审计和病例清单；不要沿用升级前的派生队列。

### 多份SAS出现相同变量

在本机配置中为该变量指定一个来源，不会根据文件顺序自动覆盖：

```yaml
variable_sources:
  BSL_HCY: laboratory_baseline.sas7bdat
  y1_is: outcomes_one_year.sas7bdat
  y1_is_dd: outcomes_one_year.sas7bdat
```

不同目录中存在同名文件时使用绝对路径。SAS变量大小写按其语言规则匹配，并保留实际拼写；不把近似名称视为同一变量。每个被选来源必须具有唯一、非空的`code_n`。ID为字符串并保留前导零，不根据临床数据推测影像ID。

### 日期、编码与单位

默认采用文件自带编码和日期元数据。确认编码错误时才设置`sas.encoding`。未声明格式的数值日期，需要在`sas.date_formats`为对应变量明确指定`date`或`datetime`。不根据数值大小猜测日期类型。

只有原SAS把ID存成数字、且已知固定宽度时才设置`sas.id_width`。字符串ID不会被自动补零。特殊缺失符号在提取CSV及元数据中保留；分析时按既定规则处理。血液单位以`docs/variable_dictionary.csv`为准，程序不推断未记录的单位换算。

## 3. 第二阶段：提供SuStaIn路径并准备病例

```bash
uv run wmh-hcy configure --sustain-dir "/data/SuStaIn"
uv run wmh-hcy run --through prepare
```

接受原项目根目录，或直接包含`sub-*/wmh/`与`sub-*/t1/`的衍生目录。自动探查常见的`derivatives`、`outputs/derivatives`及`output/derivatives`位置；发现多个候选时，传入唯一的具体目录。

每个受试者使用既有全掩膜WMH记录，必要时从WMH掩膜计算体积；不以20个区域的简单加和替代全脑WMH。ICV使用DLMUSE 702定义。无需健康锚定，也不要求SuStaIn全部40个特征完整。T1灰质体积仅影响相应子分析。

临床与影像ID不完全相同，需要另行提供**已经确认的精确映射**：

```bash
uv run wmh-hcy configure --id-map "/data/mapping/image_clinical_ids.csv"
```

映射CSV列为`participant_id,code_n`，前者对应原影像记录（通常包含`sub-`前缀）。不能确认映射的病例不会模糊匹配。源路径变化可在本机配置`imaging.path_prefix_map`中指定旧前缀到新前缀的替换。

按研究者修订，不需要影像人工质控，也无需提供QC表。旧配置中即使仍有`require_qc: true`或`qc_csv`，也不会据此阻断。既有CSV中的审核标签不参与筛选，不改写成通过。WMH和ICV按数值可用性纳入，T1灰质或急性病灶体积仅影响对应子分析。

已完成SAS和影像接入的工作站，更新后可直接运行到分析报告：

```bash
git pull --ff-only
python -m pip install --no-deps -e .
wmh-hcy run --through report
```

该命令会重新提取、连接并生成病例队列，再按H1–H4分析，避免沿用曾因审核状态排除病例的旧队列。无需重建Conda环境。

`run --through prepare`依次执行：临床审计 → 白名单CSV提取 → 影像连接/体积可用性计数 → 各假设病例清单。若某步发现数据问题，后续步骤停止，状态保存在`outputs/real/pipeline_status.json`。

重点审阅：

- `audit/imaging_summary.json`：影像ID交集、全脑WMH与ICV可用性，`wmh_icv_eligible_matched`为两种体积均可用且ID匹配的人数；`manual_image_review_required`固定为false。
- `prepared/flow.csv`：各病例规则的纳入/排除数量。
- `prepared/`各队列和排除记录：实际用于H1–H4的病例及原因；包含患者级信息，留在工作站。

## 4. 第三阶段：统计分析和报告

审阅病例准备后执行：

```bash
# 全流程重新审计和准备后，按H1-H4顺序运行并生成报告。
uv run wmh-hcy run --through report

# 或使用已准备好的队列，单独执行H2。
uv run wmh-hcy analyse --hypothesis H2
uv run wmh-hcy report
```

主分析使用50份多重插补、每份200次bootstrap，完整运行可能较长。演示的低重复次数不作为正式研究配置。模型输出保留每项`ESTIMATED`或`NOT_ESTIMABLE`状态，并给出失败原因；总流程结束并不代表每个模型均可解释，需审阅`status.json`与诊断文件。

报告路径由`outputs/real/latest_results.json`指向。各次模型结果保存在`outputs/real/results/<UTC时间>/`，汇总表为`hypothesis_summary.csv`。`run`不带`--through`时默认停在病例准备，不会自动启动长时间统计。

修改源文件或输入配置后，重新运行`run --through prepare`，避免让旧CSV对应新配置。需要不同输出批次时，在本机配置中修改`output_dir`。

## 5. 无患者数据验证与后续协作

```bash
uv run wmh-hcy demo --n 700
```

演示数据以`SYN`标记，仅用于验证软件，保存在`examples/synthetic/`和`outputs/demo/`。真实入口遇到缺少数据会停止，不会自动改用演示数据。

第一轮完成后，可以提供`audit/SUMMARY.md`、`summary.json`和`hypothesis_counts.csv`的汇总内容，随后根据实际字段与病例交集修订配置。不要上传患者CSV或ID映射。要继续在工作站执行，可在该工作站打开克隆后的项目；当前电脑无法直接操作一台尚未连接的工作站。

更新代码使用`git pull --ff-only`。Conda用户激活`wmh-hcy`后执行`python -m pip install -r requirements-conda.txt`、`python -m pip install --no-deps -e .`和`python -m pip check`；uv用户执行`uv sync --frozen`。本机配置和输出已被Git忽略；不要使用强制添加将其纳入提交。

维护者修改`uv.lock`后，应同时运行`uv export --frozen --no-emit-project --no-hashes --no-annotate --output-file requirements-conda.txt`更新Conda依赖清单。该导出步骤仅用于维护仓库，工作站安装和分析不需要uv。

## 6. 只能截图时：诊断摘要

在已完成分析的工作站执行以下命令。摘要脚本只使用Python标准库，无需重新安装项目或重新运行分析。

```bash
cd /data/usersdir/linzhenzong/WMH
git pull --ff-only
python scripts/diagnose_results.py --page 1
```

第一张截图保留整个摘要，包含H1–H4模型状态、H2/H3病例及事件数量、绝对风险计算状态、bootstrap失败类型和死亡模型结果文件的存在情况。然后执行并截取第二页：

```bash
python scripts/diagnose_results.py --page 2
```

第二页包含插补份数及轨迹摘要、Cox模型矩阵及梯度诊断、描述性比例风险检查和mRS模型诊断。每页正常约30行；可将终端最大化后截图。摘要不读取患者CSV、不显示患者ID、不写入或修改结果。

默认读取本项目`outputs/real/latest_results.json`指向的结果。若输出位置自定义，加`--output-dir /实际输出目录`；若需检查指定批次，再加`--results /实际结果批次目录`。两页应使用相同路径参数。缺少诊断文件会显示MISSING，不需要为了生成摘要重新分析。

摘要中prepared病例计数来自当前准备目录，并非本批次模型保存的快照；与模型人数不一致时会提示。插补均值和描述性PH检验仅用于定位问题，不能单独判定插补收敛、MAR或比例优势假设成立。摘要不足以定位时，再根据具体行定向截取一个文件。

若H2绝对风险bootstrap不稳定或H3死亡模型报错，再执行：

```bash
python scripts/diagnose_results.py --page 3
```

第三页在工作站本地读取`prepared/cohort_main.csv`、`cohort_month3.csv`与该批次的`analysis_participants.csv`，先逐一核对患者ID、进入时间、结束时间、事件类型和模型人数。一致后，显示六个预设分类变量各类别的总人数、缺血事件数和死亡事件数，以及已保存死亡模型的系数和标准误范围摘要。输出仅为汇总，不显示患者ID、不联网、不写文件、不重新拟合。与前两页不同，此页需要读取本地患者行；只需截图终端输出。

MATCH仅表示成员、时间、事件与模型批次一致；分类变量取自当前prepare文件、尚未插补，不代表已核验全部历史协变量。少数或零死亡的类别是定位稀疏数据的线索，不能单独证明Cox分离，也不据此自动删除协变量。若显示MISMATCH或UNAVAILABLE，保留截图，不为获得诊断而重建队列或重跑模型。第三页目前只适用于一年H2/H3，不与`--longterm`组合。

## 7. 新增2、3、4、5年分析

将包含新字段的SAS文件放在已经配置的SAS目录（允许子目录）后更新项目，沿用当前Conda环境：

```bash
cd /data/usersdir/linzhenzong/WMH
git pull --ff-only
python -m pip install --no-deps -e .
wmh-hcy longterm --through prepare
```

此命令自动从多份SAS中重新提取字段，连接原WMH影像，并建立各年H2/H3/H4队列。不会运行统计模型，也不覆盖一年结果。若SAS路径已经改变，先使用`wmh-hcy configure --sas-dir "/实际SAS目录"`。新SAS同样必须有`code_n`；若仅有另一种脱敏编号，需要已确认的精确连接键，不能按行号拼接。

检查准备摘要后执行四年全部模型及报告：

```bash
wmh-hcy longterm --through report
```

该命令重新提取并准备本批次，随后按年份执行H2基线Hcy×WMH复发模型、H3三个月Hcy更新、H4对应年份mRS及T1补充模型。它不重跑H1或一年模型。若需分开运行，可用`--years 2 3`或`--hypothesis H2`；未运行年份不会缩小Holm四次检验家族。

结果位于`outputs/real/longterm/results/<UTC时间>/report.html`，新指针为`outputs/real/longterm/latest_results.json`。每次运行保留独立批次，原`outputs/real/latest_results.json`保持不变。每年目录含病例流程、模型系数、插补和诊断；汇总为`longterm_summary.csv`及`figures/longterm_forest.png`。失败模型会标记NOT_ESTIMABLE并保留原因，其余模型继续执行。

只能截图时，按页执行：

```bash
python scripts/diagnose_results.py --longterm --page 1
python scripts/diagnose_results.py --longterm --page 2
```

第一页面向病例和事件数量，第二页面向四年估计与Holm校正。缺失或失败状态也请保留在截图中。

白名单增加28个可选字段，总数66；旧的一年数据只有38项仍可运行。长期分析仅使用各年IS/IS_DD与mRS，STROKE/HS目前只提取。`_dd`是事件或删失时间，不能将0状态均视作随访完成。当前没有完整的长期死亡日期，因此不生成2—5年竞争风险绝对发生率。详细规则及统计层级见[长期分析方案](longterm_analysis_plan.md)。

软件演示（只生成合成数据，不读取真实配置中的患者文件）：

```bash
python examples/run_longterm_demo.py --n 600
```

合成结果单独保存在`outputs/demo_longterm/longterm/`，采用2份插补，不能用于研究推断。
