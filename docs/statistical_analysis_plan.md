# CNSR-III：Hcy–WMH假设驱动统计分析方案 v2

版本日期：2026-09-16。项目：`/home/zhenzong2/analysis/WMH`。修订：从数据要求中移除发病到MRI的天数，采用基线采血日入组。工作站已启动真实SAS审计，尚未开展患者级统计推断。

## 1. 研究主题和临床意义

**研究主题：背景白质损伤是否改变Hcy所代表的卒中后风险，以及三个月Hcy能否进一步更新风险判断。**

同一个Hcy值可以出现在不同B12、叶酸和肾功能背景下；同样的Hcy异常也可以发生在背景白质损伤程度不同的人。仅分析“Hcy高者结局更差”不能回答影像是否改变Hcy的临床含义。本研究将Hcy与定量WMH共同建模，先验证结构关联，再检验复发风险是否存在交互，最后用复测值和功能结局补充证据。

### 1.1 为什么提出H1

Northern Manhattan Study在无卒中人群中发现Hcy与定量WMH有关，为代谢指标与白质负担的关联提供直接依据。[NOMAS] SMART研究把Hcy与小血管病进展联系起来，提供纵向背景。[SMART] 两项研究支持在卒中患者中检验该关联，并不能替代本研究结果。本研究沿用定量影像问题，协变量选择独立采用预设因果准则；不沿用部分既有研究按P值保留变量的做法。

B12参与Hcy再甲基化及神经髓鞘相关生物过程。B12与白质病灶也存在人体观察性证据。[B12WMH; NIH] 因此，B12在本研究的主要角色是调整Hcy解读的背景因素。H1检验的是在B12等背景相近时Hcy与WMH的关系。

### 1.2 为什么提出H2

Hcy关联小血管损伤，WMH体现既有白质负担。由此提出明确的“易损背景”假设：**相同幅度的Hcy升高，在WMH较重者中对应更强的复发风险关联。** 这是由既有结构和预后证据推导的新假设，而不是已有试验结论。它的关键检验是Hcy×WMH交互，不是分别在高、低WMH组内寻找显著性。

如果H2得到支持，研究将给出相同Hcy差异在不同WMH水平下对应多少绝对复发风险差，形成联合风险分层依据。如果只有WMH主效应，临床价值落在影像的基础风险信息；Hcy与WMH交互假设则按预设判定记录。

### 1.3 为什么提出H3和H4

Shi等的卒中研究分析了急性期及恢复期Hcy，提示恢复期测量与后续复发有关。[M3_HCY] H3把这一问题改写为可检验的增量问题：在基线Hcy非线性、WMH及其交互已进入模型后，三个月Hcy是否仍有正向条件关联。

Hcy/叶酸与卒中不良功能结局、WMH与功能结局分别已有研究；脑体积也与卒中功能结局有关。[FOLATE_OUTCOME; WMH_FUNCTION; BRAIN_VOLUME] H4检验两种背景因素的联合关联是否也出现在功能结局上，并明确要求调整本次急性损伤。四个假设共同覆盖结构、复发、复测和功能四个问题；它们彼此支持研究叙事，统计检验互不作为启动门槛。

## 2. 四条可检验假设

|顺序|明确研究假设|核心检验|预期方向|
|---|---|---|---|
|H1 结构|调整后，Hcy 15对10 μmol/L对应更高WMH负担|两个固定Hcy值的非线性对比Δ1|Δ1>0|
|H2 复发，主要|WMH越重，基线Hcy与一年缺血性复发的正向关联越强|Cox中的H×W，一自由度|βHW>0|
|H3 更新|基线Hcy和WMH相近时，三个月Hcy 15对10对应更高后续复发风险|三个月Hcy非线性对比Δ3|Δ3>0|
|H4 功能|急性损伤等相近时，WMH越重，Hcy与12个月较差mRS的正向关联越强|有序logistic中的H×W，一自由度|δHW>0|

研究估计调整后的关联。统计检验用于判断数据是否支持上述命题；因果机制由这些结果提供间接证据。

## 3. 协变量如何确定：证据指定角色，统计准则选择集合

### 3.1 可执行选择规则

采用VanderWeele提出的修正析取原因准则。[VANDERWEELE] 对已提供的候选字段，按研究发生顺序判断：它是否为Hcy或结局的背景原因，或是否测量了二者的共同原因；保留满足规则者，排除已知工具变量及暴露后的结构后果。ICV和采血时点作为测量设计项另行加入。代码`adjustment.py`直接执行这一规则，并输出选择记录。

**角色指定的工作假设固定如下：**年龄、性别和长期生活方式位于当前Hcy测量之前；B12/叶酸状态和肾功能影响当前Hcy；既往血管病及卒中史影响当前脑结构和再次事件。本次NIHSS及急性梗死负担承接先前血管/脑损伤的影响；本次急性状态也可能影响随后测得的Hcy。因果图中的Hcy代谢背景指长期状态，基线Hcy是其近似测量，不能把采血结果画成本次卒中的原因。TOAST描述已发生的本次事件机制。因此主结构/复发模型调整背景集合，功能模型另对急性状态进行条件化。[HORD; HORD_METAB; NOMAS; RECURRENCE; OVERADJUST]

这里的统计依据是有理论条件的选择准则，临床依据是每项角色的证据。单凭样本内相关性无法分辨混杂因素和中间变量。文献报告关联也不等于已证明因果箭头；本研究将角色假设逐项冻结，避免分析结果反过来决定调整集合。主模型估计C0标准化的背景关联，并不宣称C0是充分因果调整集；H2的急性扩展专门量化加入本次病情后的条件关联。

### 3.2 固定调整集合

|集合|组成|进入模型|
|---|---|---|
|C0：10项背景因素|年龄、性别、B12、叶酸、CYSC、吸烟、饮酒、高血压、糖尿病、既往卒中|H1、H2、H4|
|D0：2项测量因素|真实ICV；基线采血距发病天数|H1、H2、H4|
|C3：恢复期背景|C0中B12、叶酸、CYSC替换为三个月同期值；其余保留|H3|
|H3基线信息|基线Hcy的两个RCS项、基线Hcy×WMH、WMH两个RCS项、ICV、两次采血时点|H3的M0和M1均有|
|A：急性状态|入院NIHSS、TOAST、急性梗死体积|H2敏感性；H4|
|P：既往功能|卒中前mRS|H4|
|T1补充|GM119灰质总体积|H4合格T1子样本|

**既往卒中在v2进入C0。** 它是实际提供的既往事件字段，既往损伤可以影响当前结构，且既往事件是复发背景的一部分。[RECURRENCE] 将它只留在临床扩展模型缺乏充分理由，因此本版改正。

CYSC作为肾功能代理；肌酐作为同一角色的替代敏感性分析，二者不在主模型中重复进入。B12与叶酸分别反映不同代谢环节，二者同时保留。急性状态模型回答“在本次卒中情况相近时”的条件问题，其估计量与主复发模型不同，因此单独报告。[B12WMH; HORD_METAB; OVERADJUST]

未纳入的其他代谢物和SNP没有被指定为本研究的必要背景原因，故不进入主调整集合。协变量全表见`covariate_rationale.csv`；表中列出原字段、角色、路径、转换、文献键及程序选择标记。

## 4. 人群、影像和时间

来源：成年缺血性卒中患者，临床–影像ID明确映射。主WMH来自急性梗死处理后的全掩膜体积，单位mL；原始WMH用于敏感性分析。ICV使用DLMUSE702；GM119单独形成T1子样本。影像主队列使用WMH及ICV的实际质量条件。

H1使用基线血液、WMH和ICV合格人群，不以随访结局完整为入组条件。H2以发病为时间尺度，从实际基线采血日进入风险集；采血前或同日已发生复发者不进入该风险集。退出时间为缺血性事件、死亡、最后有效随访与365天的最早者。时间约定为(entry,exit]。

按2026-09-16的研究者决定，发病到MRI的天数不再提取、不作为调整因素或入组条件。研究将本次住院合格MRI视为基线结构表型；由于不使用扫描日期，无法由本分析核定扫描与早期复发的先后。该时间信息边界在报告中说明。

**H2保留入组后0–90天复发。** H3是单独的三个月更新队列：实际采血日仍存活且此前无缺血性事件者，从采血日起随访至发病365天。H3的入组限制不会删除H2中的早期事件。采血日之前或同日事件不计为复测后的事件。

终点只用`y1_is`及`y1_is_dd`，死亡和有效随访字段遵循已冻结字典；不重新拼接另一套复发终点。H4使用12个月mRS；已确认死亡记6，失访保持缺失。急性损伤质量不合格仅影响H4或急性扩展。

## 5. H1：结构关联

模型：`E[log(1+WMH)] = α + fH(log2 Hcy) + fA(age) + C0其余项 + D0`。

Hcy和年龄采用10%、50%、90%三个节点的限制性立方样条，每个含线性和一个非线性项。固定10和15须均在观测Hcy的第5–95百分位范围内；超出时报告该对比不可估计并保留连续曲线。B12、叶酸和肾功能按log2；分类变量使用哑变量。HC3标准误处理异方差。

估计Δ1=fH(log2 15)−fH(log2 10)，使用完整系数协方差c′Vc计算方差，再按Rubin规则汇总。Δ1>0且Holm校正P<0.05，支持H1。报告Δ1及95%CI、exp(Δ1)；后者是`1+WMH`的比值，必须保留“加1后体积”含义。

原始WMH重复同一模型，检查结果对病灶校正是否敏感。主要结论以校正后WMH为准，不按P值选择影像版本。

## 6. H2：基线Hcy与WMH的复发交互

定义H=log2(Hcy)−均值，W=[log(1+WMH)−均值]/SD。原因别Cox：

`log λIS(t) = log λ0IS(t) + fH(H) + fW(W) + βHW·H·W + fA(age) + C0其余项 + D0`。

H、W和年龄主效应均为三个节点RCS，交互固定1自由度。采用PHReg、延迟入组及Efron并列处理。唯一主要检验H0:βHW=0；双侧α=0.05。βHW>0且P<0.05支持H2，βHW<0且P<0.05则支持相反方向。

exp(βHW)为Hcy每翻倍的HR在W相差1SD时的比值。在Hcy从10增至15时，相应HR之比为exp[βHW·log2(1.5)·ΔW]。主效应中的单个线性系数不能替代整个非线性Hcy效应。

分别拟合死亡原因别Cox，联合两种原因风险计算累计发生概率。给出W第25/50/75百分位、Hcy观测支持范围内的风险曲线；重点报告Hcy15对10的风险差RD和95%CI。交互检验发生在相对风险尺度；不同基础风险也可以产生不同RD，应分别呈现。

绝对风险在实际分析样本的背景协变量与进入时间分布上标准化，解释为基线采血时仍处于风险集者至发病365天的平均风险。每份插补内按患者bootstrap估计不确定性，合并插补内和插补间方差。生产配置50份插补、每份200次bootstrap；参考点之外曲线为点估计。

固定敏感性分析：完整病例；原始WMH；肌酐替代CYSC；加急性状态A；入组后至90天与90–365天。敏感性分析不形成额外确认性检验家族，不据其显著性改换主结论。

## 7. H3：三个月Hcy的条件风险更新

M0固定包含：基线Hcy两个RCS项、W两个RCS项、基线Hcy×W、C3、ICV、基线采血距发病天数及实际三个月采血距发病天数。

M1=M0+三个月log2(Hcy)的两个RCS项。M0和M1使用同一批患者、同一套插补后的协变量。保留基线Hcy的非线性和交互，可以避免把漏调的基线关系误归于复测值。

主要对比Δ3=f3(log2 15)−f3(log2 10)，报告条件HR=exp(Δ3)及95%CI。Δ3>0且Holm校正P<0.05支持H3。该检验明确针对“15对10风险更高”命题；两个样条项的总体显著性并不自动证明这一方向。

分别建立两种竞争事件模型，输出M0/M1的拟合风险分布与个体风险变化文件，以及不同WMH水平的M3Hcy连续风险曲线。本次H3不再增加M3Hcy×WMH检验，保持问题集中在复测信息。模型内风险更新是解释性结果；临床预测性能评价需要独立校准和验证。

## 8. H4：一年功能结局

使用有序logistic模型：`logit P(mRS>k) = αk + fH(H) + fW(W) + δHW·H·W + fA(age) + C0其余项 + D0 + A + P`，k=0,…,5。

较大系数对应更差mRS。δHW>0且Holm校正P<0.05支持H4，报告共同OR之比exp(δHW)。T1合格子样本在同一模型加入GM119每100mL，报告其条件关联及Hcy×WMH估计。

诊断包括模型收敛、各阈值二项模型的交互估计及比例优势一致性。如共同OR假设明显不合适，H4共同OR结论标为“模型假设不满足”，展示阈值估计作为探索性补充。全模型秩不足或不收敛时保留失败原因；不自动删除预设背景变量以得到显著结果。

## 9. 多重比较、缺失和判定

H2为唯一主要假设。H1、H3、H4作为固定三项次要家族进行Holm校正；即使某项不可估计，也保持三项家族，计算时该项P置1并单独标记不可估计。所有P值双侧。方向正且P达到阈值→支持预设假设；方向负且达到阈值→支持相反方向；否则→未获得支持。逐项95%CI反映估计精度，不能替代次要家族校正后的判定。正向交互严格表示Hcy对应的HR或OR随WMH增加而升高；各WMH水平下HR或OR是否大于1，另由条件对比及风险曲线报告。

协变量散在缺失采用miceforest链式多重插补，生产配置50份。连续正偏态标志物在对数尺度插补，分类变量保持字典编码。插补预测因子包含观测Hcy/WMH及乘积、事件类型和考虑延迟入组的累计风险信息；M3附加基线Hcy×WMH。交互由每份完成数据的暴露确定，回归矩阵按冻结定义重建。

患者ID、观测Hcy、核心WMH、关键事件时间与整项未参加的子研究保持实际观测状态。B12等整列不可用时模型明确报告该前提缺失。采用MAR工作假设，并以完整病例和有限模拟检查数值表现。[WHITE; MICOX; MICE]

## 10. 按假设排列的技术路线与项目

`白名单SAS→CSV与元数据 → ID精确连接/WMH与ICV质控 → 各假设病例清单 → 固定协变量集合 → H1结构 → H2复发 → H3复测更新 → H4功能 → 假设判定总表`。

代码：`src/wmh_hcy/hypotheses/h01_structure.py`至`h04_function.py`；阅读入口：`analyses/01_structure`至`04_function`；输出也使用同一顺序。共享模块实现提取、队列、插补、回归和绝对风险，各假设入口调用相应定义。

每次运行保存配置、协变量选择记录、患者清单、模型矩阵转换、每份插补信息、系数、完整对比结果及诊断。研究级`hypothesis_summary.csv`固定包含H1–H4四行；单独运行某假设也不会缩小多重比较家族。

## 11. 交付和结果审阅

第一张主表：四条假设的估计量、95%CI、原始/校正P和判定。H1附结构对比；H2附绝对风险曲线、参考点RD与敏感性森林图；H3附复测条件HR及M0/M1风险分布；H4附共同OR之比及阈值诊断。全部病例数量由真实接入后生成。

研究价值的核心判断为H2交互大小及精度；H1回答结构对应，H3回答风险信息更新，H4回答功能后果。真实数据分析前冻结假设、角色表、模型和参考值，运行后按同一顺序审阅。

## 参考文献

各节方括号文献键与`references.md`及PPT备注逐项对应。



**[NOMAS]** Wright CB, Paik MC, Brown TR, et al. Total homocysteine is associated with white matter hyperintensity volume: the Northern Manhattan Study. Stroke. 2005. [https://doi.org/10.1161/01.STR.0000165923.02318.22](https://doi.org/10.1161/01.STR.0000165923.02318.22)

**[SMART]** Kloppenborg RP, Geerlings MI, Visseren FL, et al. Homocysteine and progression of generalized small-vessel disease. Neurology. 2014. [https://doi.org/10.1212/WNL.0000000000000168](https://doi.org/10.1212/WNL.0000000000000168)

**[B12WMH]** de Lau LML, Smith AD, Refsum H, et al. Plasma vitamin B12 status and cerebral white-matter lesions. Journal of Neurology, Neurosurgery &amp; Psychiatry. 2009. [https://doi.org/10.1136/jnnp.2008.149286](https://doi.org/10.1136/jnnp.2008.149286)

**[NIH]** National Institutes of Health Office of Dietary Supplements. Vitamin B12 Fact Sheet for Health Professionals. NIH Office of Dietary Supplements. 2026. [https://ods.od.nih.gov/factsheets/VitaminB12-HealthProfessional/](https://ods.od.nih.gov/factsheets/VitaminB12-HealthProfessional/)

**[HORD]** Nygård O, Refsum H, Ueland PM, Vollset SE. Major lifestyle determinants of plasma total homocysteine distribution: the Hordaland Homocysteine Study. American Journal of Clinical Nutrition. 1998. [https://pubmed.ncbi.nlm.nih.gov/9459374/](https://pubmed.ncbi.nlm.nih.gov/9459374/)

**[HORD_METAB]** Refsum H, Nurk E, Smith AD, et al. The Hordaland Homocysteine Study: A Community-Based Study of Homocysteine, Its Determinants, and Associations with Disease. Journal of Nutrition. 2006. [https://doi.org/10.1093/jn/136.6.1731S](https://doi.org/10.1093/jn/136.6.1731S)

**[RECURRENCE]** Weimar C, Diener HC, Alberts MJ, et al. The Essen stroke risk score predicts recurrent cardiovascular events: a validation within the REduction of Atherothrombosis for Continued Health (REACH) registry. Stroke. 2009. [https://doi.org/10.1161/STROKEAHA.108.521419](https://doi.org/10.1161/STROKEAHA.108.521419)

**[M3_HCY]** Shi Z, Liu S, Guan Y, et al. Changes in total homocysteine levels after acute stroke and recurrence of stroke. Scientific Reports. 2018. [https://doi.org/10.1038/s41598-018-25398-5](https://doi.org/10.1038/s41598-018-25398-5)

**[FOLATE_OUTCOME]** Shi M, Zheng J, Liu Y, et al. Folate, Homocysteine, and Adverse Outcomes After Ischemic Stroke. Journal of the American Heart Association. 2024. [https://doi.org/10.1161/JAHA.124.036527](https://doi.org/10.1161/JAHA.124.036527)

**[WMH_FUNCTION]** Griessenauer CJ, et al. Effects of White Matter Hyperintensities on 90-Day Functional Outcome after Large Vessel and Non-Large Vessel Stroke. Cerebrovascular Diseases. 2020. [https://pubmed.ncbi.nlm.nih.gov/32694259/](https://pubmed.ncbi.nlm.nih.gov/32694259/)

**[BRAIN_VOLUME]** Schirmer MD, et al. Brain Volume: An Important Determinant of Functional Outcome After Acute Ischemic Stroke. Mayo Clinic Proceedings. 2020. [https://pubmed.ncbi.nlm.nih.gov/32370856/](https://pubmed.ncbi.nlm.nih.gov/32370856/)

**[CNSR]** Wang Y, Jing J, Meng X, et al. The Third China National Stroke Registry (CNSR-III) for patients with acute ischaemic stroke or transient ischaemic attack: design, rationale and baseline patient characteristics. Stroke and Vascular Neurology. 2019. [https://doi.org/10.1136/svn-2019-000242](https://doi.org/10.1136/svn-2019-000242)

**[VANDERWEELE]** VanderWeele TJ. Principles of confounder selection. European Journal of Epidemiology. 2019. [https://doi.org/10.1007/s10654-019-00494-6](https://doi.org/10.1007/s10654-019-00494-6)

**[OVERADJUST]** Schisterman EF, Cole SR, Platt RW. Overadjustment bias and unnecessary adjustment in epidemiologic studies. Epidemiology. 2009. [https://doi.org/10.1097/EDE.0b013e3181a819a1](https://doi.org/10.1097/EDE.0b013e3181a819a1)

**[WHITE]** White IR, Royston P. Imputing missing covariate values for the Cox model. Statistics in Medicine. 2009. [https://doi.org/10.1002/sim.3618](https://doi.org/10.1002/sim.3618)

**[MICOX]** Bonneville EF, Resche-Rigon M, Schetelig J, et al. Multiple imputation for cause-specific Cox models: Assessing methods for estimation and prediction. Statistical Methods in Medical Research. 2022. [https://doi.org/10.1177/09622802221102623](https://doi.org/10.1177/09622802221102623)

**[CUTOFF]** Altman DG, Royston P. The cost of dichotomising continuous variables. BMJ. 2006. [https://doi.org/10.1136/bmj.332.7549.1080](https://doi.org/10.1136/bmj.332.7549.1080)

**[PHREG]**  statsmodels PHReg: proportional hazards regression, entry and Efron ties. Official documentation. 2026. [https://www.statsmodels.org/stable/generated/statsmodels.duration.hazard_regression.PHReg.html](https://www.statsmodels.org/stable/generated/statsmodels.duration.hazard_regression.PHReg.html)

**[ORDERED]**  statsmodels OrderedModel: ordinal regression. Official documentation. 2026. [https://www.statsmodels.org/stable/generated/statsmodels.miscmodels.ordinal_model.OrderedModel.html](https://www.statsmodels.org/stable/generated/statsmodels.miscmodels.ordinal_model.OrderedModel.html)

**[MICE]**  miceforest: multiple imputation with LightGBM. Official documentation. 2026. [https://github.com/AnotherSamWilson/miceforest](https://github.com/AnotherSamWilson/miceforest)

**[PYREAD]**  pyreadstat: SAS files, metadata and special missing values. Official documentation. 2026. [https://github.com/Roche/pyreadstat](https://github.com/Roche/pyreadstat)
