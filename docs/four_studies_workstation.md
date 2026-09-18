# 四项独立研究的工作站操作

新命令是 `wmh-study`。已有 `wmh-hcy`、一年分析和 `recurrence_v3` 的结果及指针均保留。无需重新整理原始SAS或SuStaIn目录。

当前为第四版 `imaging_four_studies_20260918_v4`：保留recovery、bp、cec、kidney，取消神经酰胺和视觉重度WMH子组。IMG_ICAS和Apo_AI不再进入任何本版研究的字段提取、插补或模型；不再要求补齐，也不再安排Apo_AI扩展。冠心病跳答核查、CEC主模型HDL-C调整保持不变。汇总使用四项Holm校正；旧合同结果标为PREVIOUS_VERSION。

## 第一步 更新代码和入口

在已有工作站仓库中执行：

```bash
cd /data/usersdir/linzhenzong/WMH
conda activate wmh-hcy
git pull --ff-only
python -m pip install --no-deps -e .
wmh-study --help
```

最后一步安装新增的命令入口，不更新已经锁定的分析依赖。新工作站可使用：

```bash
git clone https://github.com/zhenzonglin/WMH.git
cd WMH
conda env create -f environment.yml
conda activate wmh-hcy
```

若使用精简代码目录 `WMH-code-v3`，在该目录更新代码，再回到原 `WMH` 目录读取本地配置：

```bash
git -C /data/usersdir/linzhenzong/WMH-code-v3 pull --ff-only
cd /data/usersdir/linzhenzong/WMH
```

适用于已从 `WMH-code-v3` 执行过 `pip install --no-deps -e .` 的环境。纯Python代码修订无需重新安装依赖。不要在精简代码目录运行真实审计，否则不会读取原目录的 `config/workstation.local.yml`。首次精简克隆的旧Git兼容方式为 `clone --no-checkout --depth 1 --filter=blob:none` 后，在新目录依次执行 `sparse-checkout init --cone`、`sparse-checkout set src config scripts`、`checkout main`；避免在旧Git中直接使用 `clone --sparse`。

默认复用 `config/workstation.local.yml`。没有本地配置时，先用既有配置入口指定实际目录：

```bash
wmh-hcy configure --sas-dir "/实际SAS目录" --sustain-dir "/实际SuStaIn产物目录"
```

如患者编号需要映射，在本地配置 `inputs.id_map_csv` 指定含 `participant_id,code_n` 两列的CSV。字符串精确连接，保留前导零。重复ID或多文件同名变量有冲突时，脚本会停止对应研究；在 `variable_sources` 中明确来源，不通过文件顺序选择。

## 第二步 先审计 不拟合

```bash
wmh-study audit --study all
```

它递归扫描多份SAS，提取研究白名单CSV，连接既有影像，分别建立四个队列。每项显示来源、例数、结局数、缺失、排除和参数数。可以逐项运行，方便截图：

```bash
wmh-study audit --study recovery
wmh-study audit --study bp
wmh-study audit --study cec
wmh-study audit --study kidney
```

`PREPARED` 表示该队列已准备好，不代表已有统计结果。`REVIEW_REQUIRED` 会列出整列缺失、非法协变量编码或设计阻断；`INPUTS_REQUIRED` 表示文件、日期、ID或来源接口需要修正。一项失败不会停止其余项审计。详细字段来源及实际观测数在该运行的 `field_audit.csv`，病例流程在 `cohort_flow.csv`。

`wmh-study run --study all --through prepare` 与上述审计执行相同的数据准备步骤。每次准备生成新时间戳目录，不覆盖旧运行。

## 第三步 分项执行模型

```bash
wmh-study run --study recovery --through report
wmh-study run --study bp --through report
wmh-study run --study cec --through report
wmh-study run --study kidney --through report
wmh-study summary --page 1
```

如果已有同配置、未分析的准备快照，报告命令使用该快照；它校验派生文件哈希。原始数据若已更新，必须重新执行 `audit` 或 `--through prepare`。已经完成的分析再次执行时会建立新运行，不覆盖原结果。

默认每个研究模型插补50份、迭代10次。不同次要结局、扩展协变量或研究人群可能需要重新插补，以包含该模型实际结局；不会沿用Hcy插补公式。完整病例分析不插补。过程会打印插补和模型进度。真实研究不要采用演示的2份插补配置。

## 第四步 查看或截图反馈

如果审计显示 `REVIEW_REQUIRED`，先使用以下三页只读诊断，无需重复提取SAS或重新拟合：

```bash
wmh-study diagnose --page 1
wmh-study diagnose --page 2
wmh-study diagnose --page 3
```

第一页区分字段不存在与字段存在但无有效值；第二页在已保存的SAS列名和标签中查找冠心病和CEC候选名称；第三页只输出冠心病编码计数、CEC及恢复期胱抑素C异常类型及其在队列中的人数。候选名称不会自动替换正式变量。三页均不写文件、不改编码、不重新扫描SAS、不拟合模型、不打印ID或个体化验值；CSV输入若在审计后发生改变，第三页拒绝检查该来源。

空的主要队列显示协变量缺失“未评估”，不能据此认为原始年龄等所有字段都不存在。`M03_CYSC`非法值计数来自临床提取表，第三页另列其中进入肾脏分析队列的人数。`H_CHD`仅依据本次已确认规则补齐：H_HD=1补0，H_CHD_TP非空补1。冲突单列，其他缺失不当作0。原始值和推导值均保留；第3页列出规则补齐及冲突计数。

异常审计同时保存全部提取记录与本研究入组队列的计数，后者为 `field_audit.csv` 的 `invalid_eligible`。本次修复限定为肾脏研究的 `M03_CYSC`：只有入组队列中的异常值阻断，已按既定规则排除者的异常不阻断；该队列包含五年状态尚未知的入组者。这个恢复期指标不用于更大的基线结构子研究。其他字段保留原有核查范围，避免遗漏血压研究独立12个月队列等人群中的异常。CEC负值、非有限值或不可解析值仍要求确认检测/编码，零值沿用现有规则允许。不会自动改变阈值、将异常值插补或丢弃全部异常记录。旧审计没有队列内计数时，需重跑对应研究的 `audit`，只读诊断本身不改变状态。

```bash
wmh-study summary --page 1
wmh-study summary --page 2
```

第一页显示四项主要检验的样本量、参数数、原始P值和固定四项Holm校正。第二页显示独立队列审计。每个报告有明确的 `real` 或 `synthetic` 标记。

结果位置：

```text
outputs/real/studies/
  summary.html
  summary.csv
  01_recovery/runs/运行时间/
  02_bp/runs/运行时间/
  05_cec/runs/运行时间/
  06_kidney/runs/运行时间/
```

每次运行包括 `config_snapshot.json`、`status.json`、`field_audit.csv`、`cohort_flow.csv`、`audit.json`、`results.csv`、`report.html`、`primary_result.png`。每个模型目录包括定义、固定尺度、插补诊断、系数、完整诊断和失败原因。

重点查看：血压的 `primary/continuous_sbp.csv`、`primary/sbp_distribution.csv` 和连续主图，固定点对比保留在 `primary/clinical_contrasts.csv`；CEC的 `cec_hdl_comparison.csv`，两模型必须同人数、同插补；脑肾的 `primary/kidney_interaction_curves.csv`，包含依赖概率、概率差和条件相对概率比。脑肾主要P值对应交互。血压另有 `chd_rule_audit.json` 聚合规则计数和本地 `chd_conflicts.csv` 冲突清单。

`eligible.csv`、`master.csv`、`exclusions.csv` 和 `model_membership.csv` 含患者信息，仅供工作站本地追溯；不要上传GitHub。截图请优先使用上面的摘要命令。`NOT_ESTIMABLE` 是未估计，不能解释为阴性。

## 本机代码验证

```bash
python -m pytest -q
python -m ruff check src tests
wmh-study demo --output outputs/synthetic/my_new_demo --n 1600
python scripts/validate_studies.py --output outputs/synthetic/my_new_validation
```

演示和验证目录必须是新路径，避免覆写。演示数据由随机程序生成，不能用来推断真实研究人数、效应或检验效能。真实入口缺输入时不会自动使用演示数据。
