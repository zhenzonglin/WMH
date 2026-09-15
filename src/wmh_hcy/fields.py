"""Frozen clinical whitelist. No unavailable center or supplement variables."""
from __future__ import annotations

# source: (analysis name, unit/type, allowed codes); unknown codes are not guessed.
FIELDS = {
    "code_n": ("patient_id", "id", None),
    "AGE": ("age", "years", None),
    "GENDER": ("sex", "category", [1, 2]),
    "D_DIAG": ("diagnosis", "category", [1, 2]),
    "BSL_HCY": ("hcy", "umol/L", None),
    "BSL_B12": ("b12", "pmol/L", None),
    "BSL_B9": ("folate", "nmol/L", None),
    "BSL_CYSC": ("cysc", "mg/L", None),
    "BSL_Cr": ("creatinine", "umol/L", None),
    "H_SMK": ("smoking", "category", [1, 2, 3, 4]),
    "H_DRINK": ("drinking", "category", [1, 2, 3, 4]),
    "H_HYPT": ("hypertension", "category", [1, 2]),
    "H_DIAB": ("diabetes", "category", [1, 2]),
    "H_STROKE": ("prior_stroke", "category", [1, 2]),
    "A_NIHSS": ("nihss", "score", None),
    "IMG_C_TOAST": ("toast", "category", [1, 2, 3, 4, 5]),
    "H_MRS": ("pre_mrs", "score", None),
    "F12_MRS": ("mrs12", "score", None),
    "ONSET_D": ("onset_date", "date", None),
    "I_BLDSAMP_DT": ("sample_date", "datetime", None),
    "F3_BLDSAMP_D": ("sample3_date", "date", None),
    "M03_HCY": ("hcy3", "umol/L", None),
    "M03_B12": ("b123", "pmol/L", None),
    "M03_B9": ("folate3", "nmol/L", None),
    "M03_CYSC": ("cysc3", "mg/L", None),
    "y1_is": ("is_event", "category", [0, 1]),
    "y1_is_dd": ("is_day", "days", None),
    "D_DEATH": ("hospital_death", "category", [1, 2]),
    "D_DEATH_D": ("hospital_death_date", "date", None),
}
for _v in [3, 6, 12]:
    FIELDS.update({
        f"F{_v}_DATE": (f"visit{_v}_date", "date", None),
        f"F{_v}_DEATH": (f"death{_v}", "category", [1, 2]),
        f"F{_v}_DEATH_D": (f"death{_v}_date", "date", None),
    })

CORE_CATEGORIES = ["sex", "smoking", "drinking", "hypertension", "diabetes", "prior_stroke"]
CORE_CONTINUOUS = ["age", "b12", "folate", "cysc", "icv_ml", "sample_day"]
CORE = CORE_CATEGORIES + CORE_CONTINUOUS
EXTENDED = ["nihss", "toast", "lesion_ml"]
CATEGORIES = {v[0]: v[2] for v in FIELDS.values() if v[1] == "category"}

# 98 is unknown only for the original fields whose dictionary explicitly says so.
UNKNOWN_98 = {"H_HYPT", "H_DIAB", "H_STROKE"}
