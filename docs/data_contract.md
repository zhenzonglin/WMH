# 数据接口

## SAS和CSV

仅使用`variable_dictionary.csv`的源字段。SAS变量按大小写不敏感规则匹配，并记录实际拼写；提取CSV使用字典的标准拼写。直接提供CSV时需使用标准列名。主要ID为`code_n`。不读取全基因组、无关临床变量或任何假设存在的用药变量。

每个SAS来源内部一名患者一行。重复ID会停止提取，避免笛卡尔积。数值ID只能恢复其存储值，不能凭空恢复前导零；确有统一位数时在`sas.id_width`明确配置。输入文件中的字符ID保留原样。

输出UTF-8 CSV和JSON元数据。SAS特殊缺失值在原始提取中保留，在标准化中转为缺失，并记录来源。仅字典明确规定的字段把98解释为不详。Hcy=98等合法连续测量值不因数值相同被删除。

日期字符串需可解析。数值日期先用原始SAS格式识别DATE/DATETIME。格式缺失时通过`sas.date_formats: {I_BLDSAMP_DT: datetime}`显式配置。实际单位不符时修订映射，禁止混合单位直接分析。

## 影像输入

默认读取`inputs.derivatives_root`的`sub-*`目录。

- WMH：`wmh/wmh_features.json`中的`contralateral_correction.wmh_volume_after_correction_ml`、`wmh_volume_before_correction_ml`，或阶段状态的相应details。
- WMH记录缺失时：读取已存在的校正与原始二值掩膜。要求三维、可识别空间单位、有效affine，不自动阈值化概率图。
- ICV：`t1/t1_features.json`的`dlicv_icv.icv_label702_ml`；或`t1/nichart_tool_output/DLMUSE_Volumes.csv`的702列，单位mm³，除1000转换mL。
- GM119：119个原生标签体积之和；JSON未生成时使用项目携带的官方标签集合读取原始体积CSV。该集合提取自原Substain映射的`in_official_gm119`列。
- 急性梗死：既有`lesion/lesion_space-FLAIR.nii.gz`，或原participants.tsv中的`lesion_mask`。只提取体积，不重新配准。

不要求常模转换成功，不读取40特征入组标志。20区域WMH体积之和不作为全脑WMH。

也可提供`inputs.imaging_csv`，列如下：

```text
participant_id,wmh_ml,wmh_raw_ml,icv_ml,lesion_ml,gm119_ml,wmh_qc,icv_qc,t1_qc,lesion_qc
```

体积统一为mL。缺失T1或梗死体积只影响对应扩展分析。原始WMH缺失只影响该敏感性分析。

## 质控和映射

`wmh_qc/icv_qc/t1_qc/lesion_qc`取`pass/fail/unreviewed/stale`。原项目`tables/qc_reviews.tsv`的全局pass可用；全局fail或stale不会自动变为局部pass。可以用明确的模态复核CSV覆盖对应字段，列为`participant_id`及上述质控列。

主分析要求WMH和ICV质控通过。人工未审阅不会被解释为通过。程序的体积合理性检查不能替代图像阅片。因原T1常模失败而未生成完整JSON的病例，可从已有原始702体积恢复ICV，但仍需对应质控合格。

临床ID和影像ID不一致时，提供`participant_id,code_n`映射CSV。需要选择重复扫描时先明确一例一条映射，不取“第一条”。移动过的影像绝对路径可用`imaging.path_prefix_map`修改目录前缀。

## 事件与时间

主终点仅`y1_is`（0/1）及`y1_is_dd`（发病后天数）。正事件缺失日期、日期超过365天或事件晚于已知死亡均不能进入主时间结局模型。死亡原因不用于补造缺血性复发。

采用日精度。`entry=基线采血日−发病日`，风险区间为(entry, exit]。采血前或同日发生的事件单列排除；采血后3个月内的事件仍进入主分析。PHReg传入nextafter(entry,+∞)仅实现开放左端点，不创造临床事件时间。

不读取或推算发病到MRI的天数。既有合格MRI作为基线结构表型，入组时间仅由基线采血日期确定。

有效随访日期须有同次明确生存状态。已确认的缺血性事件或死亡本身也证明观察到相应日期。同日明确缺血性事件和死亡优先记录缺血性事件，不据此诊断致死性复发。

## 不支持的输入假设

没有中心信息、实验室结果回报日期或检测批次。`eGFR`歧义字段不使用。实际SAS列存在不代表患者确有测量，完全缺失的调整变量不能通过插补创造出来。

## v2按假设的数据要求

H1使用基线影像/代谢合格病例，不要求结局完整。H2从实际基线采血日进入风险集。H3使用实际M3采血日及同期HCY/B12/B9/CYSC，保留基线Hcy和WMH信息。H4要求mRS及合格急性损伤体积，T1另建子样本。H_STROKE进入v2核心调整集合。当前源字段白名单为38项。
