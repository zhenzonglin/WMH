"""Prespecified fields and models, not data-driven variable selection."""
from __future__ import annotations

from dataclasses import dataclass, field, replace

from ..fields import FIELDS as LEGACY_FIELDS

# (canonical name, type, permitted codes). No Hcy/B12 requirement is inherited.
SOURCES = {k: v for k, v in LEGACY_FIELDS.items() if k in {
    "code_n", "AGE", "GENDER", "D_DIAG", "H_SMK", "H_DRINK", "H_HYPT", "H_DIAB",
    "H_STROKE", "H_MRS", "A_NIHSS", "IMG_C_TOAST", "BSL_CYSC", "ONSET_D", "I_BLDSAMP_DT",
    "F3_DATE", "F12_DATE", "F3_DEATH", "F6_DEATH", "F12_DEATH", "D_DEATH",
    "D_DEATH_D", "F3_DEATH_D", "F6_DEATH_D", "F12_DEATH_D", "F12_MRS",
    "Y5_IS", "Y5_IS_DD", "y1_is", "y1_is_dd", "y2_IS", "y2_IS_dd",
    "y3_IS", "y3_IS_dd", "y4_IS", "y4_IS_dd", "m24_mrs", "m36_mrs", "m48_mrs", "m60_mrs",
}}
SOURCES.update({
    "EDUC": ("education", "category", [1, 2, 3, 4, 5]),
    "BMI": ("bmi", "continuous", None),
    "H_CHD": ("chd", "category", [0, 1]),
    "H_HD": ("heart_disease_gate", "continuous", None),
    "H_CHD_TP": ("chd_type_present", "presence", None),
    "MH_HYPT": ("prior_bp_med", "category", [1, 2]),
    "DM_HYPT": ("discharge_bp_med", "category", [1, 2]),
    "DM_HYPT_ACEI": ("discharge_acei", "category", [0, 1]),
    "DM_HYPT_ARB": ("discharge_arb", "category", [0, 1]),
    "MH_LL": ("prior_lipid_med", "category", [1, 2]),
    "MH_LL_STT": ("prior_statin", "category", [0, 1]),
    "A_SBP": ("sbp0", "continuous", None),
    "F3_MRS": ("mrs3", "score", list(range(6))),
    "D_MRS": ("discharge_mrs", "score", list(range(6))),
    "IMG_ICAS": ("icas", "category", [1, 2]),
    "CEC": ("cec", "continuous", None),
    "BSL_HDL": ("hdl", "continuous", None),
    "BSL_LDL": ("ldl", "continuous", None),
    "BSL_TG": ("tg", "continuous", None),
    "Apo_AI": ("apo_ai", "continuous", None),
    "BSL_UACR": ("uacr0", "continuous", None),
    "M03_UACR": ("uacr3", "continuous", None),
    "M03_CYSC": ("cysc3", "continuous", None),
})
for month in (3, 12):
    for arm in ("L", "R"):
        for component in ("SBP", "DBP"):
            SOURCES[f"F{month}_{arm}{component}"] = (
                f"{arm.lower()}{component.lower()}{month}", "continuous", None)
for year in (2, 3, 4, 5):
    SOURCES[f"F{year}Y_DEATH"] = (f"death{12*year}", "category", [1, 2])

CATEGORIES = {alias: codes for alias, kind, codes in SOURCES.values()
              if kind == "category" or alias in {"pre_mrs", "mrs3"}}
CATEGORIES.update(pre_mrs=list(range(6)), mrs3=list(range(6)), mrs12=list(range(7)),
                  raas=[0, 1], albuminuria=[0, 1, 2, 3])
UNKNOWN = {"H_HYPT": [98], "H_DIAB": [98], "H_STROKE": [98], "EDUC": [98], "IMG_ICAS": [3]}
UNITS = {"cec": "%", "uacr0": "mg/mmol", "uacr3": "mg/mmol", "cysc": "mg/L", "cysc3": "mg/L",
         "ldl": "mmol/L", "hdl": "mmol/L", "tg": "mmol/L", "apo_ai": "g/L",
         "sbp0": "mmHg", "lsbp3": "mmHg", "rsbp3": "mmHg", "bmi": "kg/m2"}
COMMON = ("age", "sex", "smoking", "drinking", "hypertension", "diabetes", "prior_stroke")
METABOLIC = COMMON + ("bmi", "education", "cysc", "ldl", "tg", "prior_statin",
                      "sample_day", "lesion_ml", "icv_ml")
COVARIATES = {
    "recovery": COMMON + ("education", "pre_mrs", "nihss", "toast", "lesion_ml", "mrs3", "icv_ml"),
    "bp": COMMON + ("bmi", "cysc", "sbp0", "chd", "toast", "icas", "prior_bp_med",
                    "discharge_bp_med", "mrs3", "icv_ml"),
    "cec": METABOLIC + ("hdl",),
    "kidney": COMMON + ("bmi", "education", "cysc3", "sbp3", "raas", "mrs3", "icv_ml"),
}
TITLES = {"recovery": "01 早期功能独立后的远期失能", "bp": "02 恢复期血压与WMH",
          "cec": "05 HDL功能与灰质结构",
          "kidney": "06 脑肾微血管损伤"}
FOLDERS = {k: v.split()[0] + "_" + k for k, v in TITLES.items()}
STUDIES = tuple(TITLES)
IMAGE_COLUMNS = {"wmh_ml", "wmh_raw_ml", "gm119_ml", "icv_ml", "lesion_ml"}
LOG_COLUMNS = {"tg", "cysc", "cysc3"}


@dataclass
class ModelSpec:
    study: str
    name: str = "primary"
    family: str = "multinomial"
    outcome: str = "state60"
    exposures: tuple[str, ...] = ()
    covariates: tuple[str, ...] = ()
    primary: tuple[str, ...] = ()
    splines: tuple[str, ...] = ("age",)
    interactions: tuple[tuple[str, str], ...] = ()
    tier: str = "sensitivity"
    transforms: dict = field(default_factory=dict)

    @property
    def predictors(self):
        return tuple(dict.fromkeys(self.exposures + self.covariates))

    def variant(self, name, **kwargs):
        return replace(self, name=name, **{"tier": "sensitivity", **kwargs})


def primary_spec(study: str) -> ModelSpec:
    spec = ModelSpec(study=study, covariates=COVARIATES[study], tier="primary")
    if study == "recovery":
        return replace(spec, exposures=("wmh_ml", "gm119_ml"),
                       primary=("dependent:wmh_ml", "dependent:gm119_ml"))
    if study == "bp":
        return replace(spec, family="cox", outcome="event_type", exposures=("sbp3", "wmh_ml"),
                       splines=("age", "sbp3", "wmh_ml"),
                       interactions=(("sbp3", "wmh_ml"), ("sbp3_rcs", "wmh_ml")),
                       primary=("sbp3_x_wmh_ml", "sbp3_rcs_x_wmh_ml"))
    if study == "cec":
        return replace(spec, family="ols", outcome="gm119_ml", exposures=("cec",), primary=("cec",))
    return replace(spec, exposures=("wmh_ml", "albuminuria"),
                   interactions=tuple((f"albuminuria_{i}", "wmh_ml") for i in (1, 2, 3)),
                   primary=("dependent:albuminuria_3_x_wmh_ml",))


def roles(study: str) -> list[dict]:
    rows = []
    reasons = {
        "age": "年龄可同时影响脑结构、血管代谢状态与预后", "sex": "性别与头颅体积、代谢和卒中背景有关",
        "smoking": "既往吸烟可影响血管损伤、脂质功能和预后", "drinking": "饮酒史可共同影响代谢和血压背景",
        "hypertension": "既往血压损伤同时影响脑与肾表型及预后", "diabetes": "糖尿病是系统性微血管及代谢背景",
        "prior_stroke": "既往卒中反映基线脑损伤与后续事件易感性", "education": "教育与脑储备及健康行为背景有关",
        "bmi": "体型与脂质、血压和肾脏指标相关", "cysc": "基线滤过功能影响血液标志物及血管结局",
        "cysc3": "恢复期滤过功能与UACR及功能预后有关", "ldl": "控制常规胆固醇负担的条件关联",
        "tg": "控制常规代谢性脂质背景", "hdl": "固定HDL胆固醇浓度以解释CEC的额外关联",
        "prior_statin": "卒中前用药可改变基线脂质状态", "sample_day": "控制急性期采血时点差异",
        "lesion_ml": "条件化本次卒中损伤量，不估计背景表型的总因果效应", "icv_ml": "独立控制头颅大小",
        "pre_mrs": "控制卒中前功能状态", "nihss": "给定入院卒中严重程度", "toast": "给定卒中病因构成",
        "mrs3": "给定恢复期功能状态，是条件预后问题而非总效应混杂因素",
        "sbp0": "控制发病入院时血压背景", "sbp3": "恢复期实测血压；研究02为暴露，06为条件调整",
        "chd": "既往冠心病可关联血压处置及血管结局", "icas": "颅内狭窄关联灌注状态和复发风险",
        "prior_bp_med": "卒中前降压使用属于访视血压前的管理背景", "discharge_bp_med": "出院降压使用先于恢复期血压测量",
        "raas": "出院ACEI或ARB先于恢复期UACR并可能影响其测量", "wmh_ml": "背景白质损伤表型",
        "gm119_ml": "背景灰质结构指标，单次测量不代表萎缩速率",
        "cec": "基线HDL胆固醇外排能力", "albuminuria": "两次实际UACR的四种组合；主要检验持续升高相对两次均低的关联是否随WMH变化",
    }
    for name in primary_spec(study).predictors:
        role = ("exposure_or_imaging" if name in primary_spec(study).exposures
                else "measurement" if name in {"icv_ml", "sample_day"}
                else "conditional_prognosis" if name in {"mrs3", "pre_mrs", "nihss", "toast", "lesion_ml"}
                else "prespecified_background")
        rows.append({"study": study, "variable": name, "role": role,
                     "source": "H_CHD;H_HD;H_CHD_TP" if name == "chd" else (
                         ";".join(s for s, v in SOURCES.items() if v[0] == name) or "derived_or_imaging"),
                     "reason": reasons.get(name, "方案预设"),
                     "timing": "month3" if name in {"mrs3", "sbp3", "cysc3", "albuminuria"} else "baseline_or_pre_stroke",
                     "method_basis": "VanderWeele 2019; clinical temporal assumptions; see SAP study-specific references",
                     "selection": "protocol_fixed; no p-value/stepwise/ML selection"})
    return rows
