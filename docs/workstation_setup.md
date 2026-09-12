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

脚本递归发现该目录下的SAS文件，读取元数据并扫描39个既定字段的有效值。此阶段不读取影像，不拟合模型，也不插补数据。源文件只读。大数据目录首次扫描可能耗时。

`config/workstation.local.yml`保存本机配置。后续命令自动选用它；更新该配置时会在本地保留备份。不要把真实路径、患者CSV或影像文件提交到仓库。

|审计文件（位于`outputs/real/audit/`）|用途|
|---|---|
|`SUMMARY.md` / `summary.json`|总体状态、来源并集人数、字段缺口、阻断原因|
|`report.html`|可在浏览器中查看的汇总表|
|`files.csv`|各文件元数据行数、扫描行数、独立ID、重复行及空ID数|
|`field_coverage.csv`|39个字段是否存在、所选来源、有效值人数、预期单位/类型|
|`variables_by_file.csv`|每个字段逐文件的有效、缺失、非法编码计数|
|`hypothesis_counts.csv`|H1–H4影像连接前的临床完整观察值交集|
|`file_intersections.csv`|临床文件间共享独立ID数|
|`sas_inventory.json`|原变量名、标签、SAS日期格式、文件编码与来源|

**这些计数不是最终分析样本量。** 重复ID不会自动保留第一行；协变量散在缺失也不会仅因“完整病例计数较少”就被永久剔除。最终病例资格由影像QC、时间规则和各假设的变量要求共同确定。

`READY_FOR_EXTRACTION`表示基础临床字段/编码检查可继续；`REVIEW_REQUIRED`表示需要按报告修正输入配置或数据问题。缺少可选子分析字段会在对应假设报告中列出。没有发现患者缺血事件时，事件时间全部为空本身不阻断提取，但不意味着能拟合复发模型。

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

原有QC记录可从`derivatives/tables/qc_reviews.tsv`读取；也可显式传入模态QC CSV：

```bash
uv run wmh-hcy configure --qc-csv "/data/qc/modality_reviews.csv"
```

具体列和状态见`docs/data_contract.md`。未审核不会自动记为通过。影像未通过和未匹配的原因保存在本地输出。主队列不要求T1灰质或急性病灶子分析均完整。

`run --through prepare`依次执行：临床审计 → 白名单CSV提取 → 影像连接/QC计数 → 各假设病例清单。若某步发现阻断问题，后续步骤停止，状态保存在`outputs/real/pipeline_status.json`。

重点审阅：

- `audit/imaging_summary.json`：影像ID交集、全脑WMH与ICV可用性、QC通过人数。
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
