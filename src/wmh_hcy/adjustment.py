"""Executable, outcome-independent modified disjunctive cause selection registry."""
from __future__ import annotations

# Evidence supports candidate role assignment; the statistical rule is applied
# to those assignments. It does not discover causal direction from patient P values.
RULE_REFERENCE = "VANDERWEELE"
BACKGROUND = [
    ("age", "AGE", True, True, False, "HORD;NOMAS;RECURRENCE"),
    ("sex", "GENDER", True, True, False, "HORD;NOMAS"),
    ("b12", "BSL_B12", True, True, False, "B12WMH;NIH"),
    ("folate", "BSL_B9", True, True, False, "HORD;FOLATE_OUTCOME"),
    ("cysc", "BSL_CYSC", False, False, True, "HORD;SMART"),
    ("smoking", "H_SMK", True, True, False, "HORD;RECURRENCE"),
    ("drinking", "H_DRINK", True, True, False, "HORD;FOLATE_OUTCOME"),
    ("hypertension", "H_HYPT", False, True, False, "NOMAS;RECURRENCE"),
    ("diabetes", "H_DIAB", False, True, False, "NOMAS;RECURRENCE"),
    ("prior_stroke", "H_STROKE", False, True, False, "RECURRENCE"),
]


def selection_registry() -> list[dict]:
    rows = []
    for name, source, exposure, outcome, proxy, references in BACKGROUND:
        rows.append({"variable": name, "source": source, "exposure_parent": exposure,
                     "outcome_parent": outcome, "shared_cause_proxy": proxy,
                     "known_instrument": False, "structural_descendant": False,
                     "selection": "modified_disjunctive", "references": references})
    for name, source, role in [
        ("icv_ml", "DLMUSE702", "measurement_size"),
        ("sample_day", "I_BLDSAMP_DT - ONSET_D", "measurement_time"),
        ("nihss", "A_NIHSS", "acute_state"),
        ("toast", "IMG_C_TOAST", "index_mechanism"),
        ("lesion_ml", "acute_lesion_mask", "acute_state"),
        ("pre_mrs", "H_MRS", "prior_function"),
        ("gm119_ml", "GM119", "structural_extension"),
    ]:
        rows.append({"variable": name, "source": source, "exposure_parent": False,
                     "outcome_parent": False, "shared_cause_proxy": False,
                     "known_instrument": False, "structural_descendant": role in {"acute_state", "index_mechanism"},
                     "selection": role, "references": "OVERADJUST;WMH_FUNCTION;BRAIN_VOLUME"})
    return rows


def select_background(rows: list[dict] | None = None) -> list[str]:
    rows = selection_registry() if rows is None else rows
    return [r["variable"] for r in rows if r["selection"] == "modified_disjunctive"
            and not r["known_instrument"] and not r["structural_descendant"]
            and (r["exposure_parent"] or r["outcome_parent"] or r["shared_cause_proxy"])]


def adjustment_columns(kind: str = "main") -> list[str]:
    columns = select_background() + ["icv_ml", "sample_day"]
    if kind == "month3":
        replace = {"b12": "b123", "folate": "folate3", "cysc": "cysc3"}
        columns = [replace.get(c, c) for c in columns]
    if kind == "creatinine":
        columns = ["creatinine" if c == "cysc" else c for c in columns]
    if kind in {"extended", "functional", "functional_t1"}:
        columns += ["nihss", "toast", "lesion_ml"]
    if kind.startswith("functional"):
        columns += ["pre_mrs"]
    return columns
