# 五年复发v3：工作站操作

此流程使用已有Conda环境和`config/workstation.local.yml`，无需重装依赖。所有真实患者数据留在工作站。

## 已有工作站的三步操作

```bash
cd /data/usersdir/linzhenzong/WMH
conda activate wmh-hcy
git pull --ff-only
python -m pip install --no-deps -e .

# 第一步：核验五年字段、连接WMH、输出病例与各月份事件数，不拟合模型。
wmh-hcy recurrence --through prepare

# 第二步：显示一屏汇总，适合截图。
python scripts/diagnose_recurrence.py

# 第三步：审计无问题后，运行五年主分析、各月份和限定稳健性分析。
wmh-hcy recurrence --through report
python scripts/diagnose_recurrence.py
```

两个入口都从只读输入重新准备一个独立运行目录，不覆盖旧运行。`prepare`不会触发插补或模型；`report`包含准备、共享插补、模型和报告。生产插补份数与迭代仍读取现有配置，默认50和10。`--hypothesis`不适用于v3。

## 新工作站

```bash
git clone https://github.com/zhenzonglin/WMH.git
cd WMH
conda env create -f environment.yml
conda activate wmh-hcy
wmh-hcy configure --sas-dir "/实际/SAS目录" --sustain-dir "/实际/SuStaIn目录"
# 如临床与影像ID不同，使用已知的一对一映射：
# wmh-hcy configure --id-map "/实际/映射.csv"
wmh-hcy recurrence --through prepare
```

既有多个SAS仍递归扫描，源变量大小写不敏感。遇到多个文件含同一使用变量，沿用`variable_sources`明确选择。不把不同患者的同名字段按行序拼接。

## 输出在哪里

路径由配置中的`output_dir`决定，默认：

```text
outputs/real/recurrence_v3/
  latest_run.json
  latest_results.json
  runs/<时间>/
    status.json
    field_audit.csv
    cohort_flow.csv
    cohort_counts.csv
    endpoint_audit.json
    horizon_results.csv
    clinical_contrasts.csv
    robustness_results.csv
    diagnostic_summary.csv
    report.html
    figures/
```

终端会打印本次运行路径。`PREPARED`表示队列准备完成；`COMPLETED`表示本轮模型均成功估计；`COMPLETED_WITH_MODEL_FAILURES`表示至少一个预设分析无法估计，查看汇总中的具体原因。通过数值拟合不代表假设成立。

截图脚本只读本轮聚合输出，不打印患者ID，不重新插补或拟合。`--config`可指定非默认配置；`--results`可指定某次运行目录。仅有截图时，先发这一屏即可。

现有`run`、`longterm`和`check-fit`继续对应旧方案。新版请明确使用`recurrence`。不要把旧版5857例和444事件写作本轮五年结果。

方法与解释：[统计分析计划](recurrence_v3_plan.md)。
