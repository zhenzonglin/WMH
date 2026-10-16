# 使用已导出的图表完成局部修订

`scripts/recovery_publication_revise.py`只读取下载后含有`v3/`的图表目录。它将原包复制到新的输出目录，然后更新两张基线表中的卒中前mRS展示、Figure 2概率热图及其说明。所有模型估计、P值、置信区间、其他图表和桑基图均保留。

```bash
python scripts/recovery_publication_revise.py \
  --export-root /path/to/all_20261015T123833Z \
  --out /path/to/all_20261015T123833Z_revised
```

- 卒中前mRS：从两张原表保存的整数频数计算中位数及第25、75百分位数，使用与原TableOne相同的线性顺序统计量定义。原始频数须完整且与每组人数一致；不根据已四舍五入的百分比推算人数，也不构造患者记录。
- 原TableOne生成的其他行和注释保留；新增脚注说明这一个汇总行的来源。将来直接从原队列导出时，TableOne也按中位数及四分位数汇总卒中前mRS。3个月mRS仍显示分类分布；回归变量编码不变。
- Figure 2：从`v3/source_data/primary_probability_grid.csv`读取九组影像场景的依赖、死亡概率；A为依赖，B为死亡，两侧使用同一绝对概率色阶。原95%CI不重新估计，继续保留在伴随表中。
- `trajectory_rrr_legend.txt`提供森林图的英文正式图注，A为WMH，B为GMV。

脚本核对原导出manifest中记录的输出校验值，拒绝覆盖原目录和非空输出目录。如果仅两张待修订基线Word表的校验值不同，则进一步逐格核对其内容与已校验CSV是否一致；内容一致时保留当前Word格式并记录校验差异，内容不一致则停止。输出`revision_manifest.json`记录全部原文件校验值、允许的修改、频数与分位数依据和代码校验值；`v3/manifest_before_revision.json`保留原工作站manifest。不同运行环境可能采用字体回退，修订图的实际字体单独记录。

此命令不访问患者级文件，不运行协变量插补、回归、联合检验或任何历史分析命令。数据不足或校验不符时停止并说明原因。
