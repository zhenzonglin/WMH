# 独立新增mRS桑基图

本脚本只生成一张描述性桑基图。已确认的v3图形、其源代码及原分析结果不修改，也不重新拟合或插补。

## 人群和口径

- 使用指定recovery trajectory运行的`base_eligible.csv`：全部3个月功能独立者，当前研究为7,529人；包括五年结局未知者。
- 访视固定为3、12、24、36、48、60个月，不加入6个月。
- 节点为mRS 0–6和Missing。每个访视及相邻区间均使用完整入组人数作为分母。
- 相邻访视流量来自同一患者的记录，不能由各访视的边际人数推算。流带不表示潜在轨迹分型，也不能沿多段流带追踪某个具体患者。
- 已确认死亡为吸收状态；原始mRS缺失但保存状态已确认死亡时显示为6。其余缺失不填补，后续再次观察到mRS时允许从Missing返回已知状态。
- 若发现死亡后转活、死亡后变为未知、mRS与保存状态冲突、人数不符或重复ID，停止绘图并保存失败原因，不删人或修改源数据。

## 工作站执行

```bash
cd /data/usersdir/linzhenzong/WMH
conda activate wmh-hcy
git fetch origin
git switch feat/recovery-publication-figures
git pull --ff-only origin feat/recovery-publication-figures
python -m pip install -c requirements-conda.txt -r requirements-publication.txt

RUN=outputs/real/studies/01_recovery_trajectory/runs/20261011T171844379322Z
OUT=outputs/publication/recovery_sankey/$(date -u +%Y%m%dT%H%M%SZ)
python scripts/recovery_mrs_sankey.py --trajectory-run "$RUN" --out "$OUT"
printf '桑基图目录：%s\n' "$OUT"
```

本命令只读取该运行的`status.json`、`audit.json`和`base_eligible.csv`。不运行`wmh-study`，不读取SAS、复发字段、影像协变量或其他课题数据。源文件在读取前后和导出后核对SHA256。

## 输出

- `mrs_sankey.pdf`、`.svg`、`.png`和单层1200 dpi `.tiff`。
- `mrs_sankey_nodes.csv`：六次访视×八种状态，共48行，含人数及完整队列占比。
- `mrs_sankey_transitions.csv`：五个相邻区间×八乘八种转移，共320行，保留零计数用于核对。
- `audit.json`、`status.json`、`manifest.json`及`figure_legend.txt`。
- `preview.html`按实际投稿宽度展示；仅新的桑基图，不重新生成v3预览。

所有导出数据均为汇总数据，不含患者ID或逐人记录。列等距表示连续访视，不能将3–12个月的间距解释为与其余区间相同的实际历时。图形沿用v3的Seaborn/SciencePlots主题、Arial类9 pt字体和低饱和配色；新增颜色区分mRS等级，灰色区分死亡及缺失。

使用`--preview-only`可省略TIFF；省略`--out`时自动创建新的时间戳目录。脚本拒绝覆盖非空目录及历史研究目录。

## 合成预览及软件验证

```bash
python scripts/recovery_mrs_sankey.py --demo --preview-only \
  --out outputs/publication/recovery_sankey/demo_$(date -u +%Y%m%dT%H%M%SZ)
ruff check src tests
python -m pytest -q
git diff --check
```

合成预览使用虚构患者随访记录，图上明确标记，不是研究结果。真实数据只在工作站执行。软件测试覆盖流量守恒、缺失后再次观测、死亡吸收、稀有转移、源文件只读及既有v3结果不被覆盖。
