"""Independent immutable runs; original Hcy and recurrence pointers are untouched."""
from __future__ import annotations

import copy
import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from ..common import DataError, dump_json, outdir, read_csv, read_json, record_run, resolve, sha256
from ..imaging import read_imaging
from ..sas_extract import extract, inventory, resolve_owners
from . import CONTRACT
from .analyses import run_study
from .data import build_cohort, chd_rule_summary, cohort_audit, derive, read_clinical
from .design import StudyDesign
from .registry import FOLDERS, SOURCES, STUDIES, primary_spec, roles


def study_sources(study):
    selected = set(primary_spec(study).predictors) | {"patient_id", "age", "diagnosis", "hospital_death"}
    if study != "bp":
        selected |= {f"mrs{m}" for m in (3, 12, 24, 36, 48, 60)}
        selected |= {f"death{m}" for m in (3, 6, 12, 24, 36, 48, 60)}
    if study == "recovery":
        selected |= {"discharge_mrs"}
    elif study == "bp":
        selected |= {"onset_date", "visit3_date", "visit12_date", "mrs12", "hospital_death_date", "death3", "death6", "death12",
                     "death3_date", "death6_date", "death12_date", "is_event", "is_day", "heart_disease_gate", "chd_type_present"}
        selected |= {f"y{y}_is_{suffix}" for y in (2, 3, 4, 5) for suffix in ("event", "day")}
    elif study == "cec":
        selected |= {"apo_ai", "prior_lipid_med", "onset_date", "sample_date", "pre_mrs", "nihss", "toast"}
    elif study == "kidney":
        selected |= {"uacr0", "uacr3", "cysc", "sbp0", "discharge_bp_med", "discharge_acei", "discharge_arb", "nihss", "toast"}
    if study in {"bp", "kidney"}:
        selected |= {f"{arm}{part}{m}" for arm in ("l", "r") for part in ("sbp", "dbp")
                     for m in ((3, 12) if study == "bp" else (3,))}
    return {s: value for s, value in SOURCES.items() if value[0] in selected}


def config_hash(cfg):
    clean = {k: v for k, v in cfg.items() if not k.startswith("_")}
    return hashlib.sha256(json.dumps(clean, sort_keys=True).encode()).hexdigest()


def study_root(cfg, study):
    return outdir(cfg) / "studies" / FOLDERS[study]


def code_hashes():
    root = Path(__file__).parent.parent
    return {str(p.relative_to(root)): sha256(p) for p in root.rglob("*.py")}


def explicit_unit_conflicts(inv, owner, fields):
    """Detect an explicitly incompatible assay label; never infer or convert units."""
    issues = []
    for source, path in owner.items():
        if source not in fields:
            continue
        alias = fields[source][0]
        record = next((r for r in inv if r["path"] == path), None)
        if not record:
            continue
        label = str(record.get("labels", {}).get(record["field_columns"][source], "") or "")
        compact = label.lower().replace(" ", "").replace("μ", "u").replace("µ", "u")
        mismatch = False
        if alias in {"uacr0", "uacr3"}:
            mismatch = bool(re.search(r"mg/g(?![a-z])|mg/mg|ug/mg", compact))
        if mismatch:
            issues.append({"source": source, "label": label, "reason": "Explicit label differs from frozen unit contract"})
    return issues


def settings(cfg):
    # Explicit study overrides are optional. Legacy bootstrap and Hcy settings are never read.
    chosen = {k: cfg.get("analysis", {}).get(k, value) for k, value in
              {"imputations": 50, "mice_iterations": 10, "mice_threads": 2, "seed": 20260917}.items()}
    chosen.update({k: v for k, v in cfg.get("studies", {}).items() if k in chosen})
    return chosen


def prepare(cfg, study):
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    folder = study_root(cfg, study) / "runs" / stamp
    folder.mkdir(parents=True, exist_ok=False)
    local = copy.deepcopy(cfg)
    local["_out"] = str(folder)
    local["study_contract"] = CONTRACT
    dump_json(folder / "config_snapshot.json", local)
    state = {"study": study, "mode": cfg["mode"], "run": stamp, "status": "PREPARING", "contract": CONTRACT,
                 "config_hash": config_hash(cfg), "path": str(folder), "statistics": settings(cfg)}
    state["source_code_sha256"] = code_hashes()
    dump_json(folder / "status.json", state)
    try:
        fields = study_sources(study)
        inv = []
        if cfg["inputs"].get("clinical_csv"):
            source = resolve(cfg, cfg["inputs"]["clinical_csv"])
            names = read_csv(source).columns
            owner = {s: str(source) for s in fields if s.casefold() in {n.casefold() for n in names}}
            dump_json(folder / "source_files.json", [{"path": str(source), "sha256": sha256(source)}])
        else:
            inv = inventory(local, fields=fields)
            owner, issues = resolve_owners(inv, local, fields=fields)
            dump_json(folder / "source_inventory.json", inv)
            dump_json(folder / "source_issues.json", issues)
            if issues:
                raise DataError("Conflicting variable owners; see source_issues.json and set variable_sources")
            extract(local, fields=fields)
        clinical, fields_audit = read_clinical(local)
        fields_audit = fields_audit.loc[fields_audit.source.isin(fields)].copy()
        unit_issues = explicit_unit_conflicts(inv, owner, fields)
        fields_audit["file"] = fields_audit.source.map(owner).fillna("")
        fields_audit.loc[fields_audit.source.eq("code_n"), "file"] = "exact ID in each selected source"
        fields_audit.to_csv(folder / "field_audit.csv", index=False)
        images = read_imaging(local, include_gm=study != "bp")
        master = derive(clinical.merge(images, on="patient_id", how="left", validate="one_to_one"))
        master.to_csv(folder / "master.csv", index=False)
        data, exclusions, flow = build_cohort(master, study)
        data.to_csv(folder / "eligible.csv", index=False)
        exclusions.to_csv(folder / "exclusions.csv", index=False)
        flow.to_csv(folder / "cohort_flow.csv", index=False)
        pd.DataFrame(roles(study)).to_csv(folder / "covariate_roles.csv", index=False)
        audit = cohort_audit(data, study)
        audit.update(clinical_n=len(clinical), mapped_image_n=len(images),
                     exact_id_intersection=int(clinical.patient_id.isin(images.patient_id).sum()),
                     fields_absent=fields_audit.loc[~fields_audit.present, "source"].tolist(),
                     invalid_fields=fields_audit.loc[fields_audit.invalid.gt(0), ["source", "invalid"]].to_dict("records"),
                     medication_conflicts=int(master.medication_conflict.sum()),
                     bp3_conflicts=int(master.bp3_conflict.sum()),
                     state60_conflicts=int(master.state60_conflict.sum()),
                     unit_contract="uacr mg/mmol; CEC percent; blood pressure mmHg; image mL",
                     counts_are_this_study_only=True)
        audit["unit_conflicts"] = unit_issues
        audit["source_groups"] = {
            Path(path).name: ", ".join(group.source) for path, group in fields_audit.loc[
                fields_audit.present & fields_audit.source.ne("code_n")].groupby("file")}
        audit["flow_excluded"] = {r.step: int(r.excluded_here) for r in flow.itertuples() if r.excluded_here}
        if study == "bp" and len(data):
            audit["entry_days_quantiles"] = data.entry.quantile([0, .25, .5, .75, 1]).to_dict()
        if study == "bp":
            audit["chd_rules"] = {"all_clinical": chd_rule_summary(master), "eligible": chd_rule_summary(data)}
            dump_json(folder / "chd_rule_audit.json", audit["chd_rules"])
            columns = ["patient_id", "chd_recorded", "heart_disease_gate", "chd_type_present", "chd", "chd_origin"]
            master.loc[master.chd_rule_conflict, columns].to_csv(folder / "chd_conflicts.csv", index=False)
        primary_aliases = set(primary_spec(study).predictors)
        invalid_covariates = fields_audit.loc[fields_audit.canonical.isin(primary_aliases) & fields_audit.invalid.gt(0), "source"].tolist()
        audit["invalid_primary_covariates_require_review"] = invalid_covariates
        try:
            design = StudyDesign.freeze(data, primary_spec(study))
            filled = data[list(design.spec.predictors)].copy()
            for c in filled:
                filled[c] = filled[c].fillna(filled[c].mode().iloc[0] if design.coding[c]["kind"] == "category" else filled[c].median())
            parameters = design.transform(filled).shape[1]
            audit["model_parameters"] = parameters*(2 if design.spec.family == "multinomial" else 1)
            dump_json(folder / "frozen_primary_design.json", design.coding)
            # This fill is used only to count columns, never written as completed patient data.
        except (DataError, IndexError) as exc:
            audit["primary_design_blocker"] = str(exc)
        audit["status"] = "REVIEW_REQUIRED" if (
            not len(data) or audit["covariates_entirely_missing"] or "primary_design_blocker" in audit
            or invalid_covariates or unit_issues
            or (study == "bp" and audit["chd_rules"]["eligible"]["conflicts"] > 0)) else "PREPARED"
        state.update(status=audit["status"], audit=audit,
                     cohort_sha256=sha256(folder / "eligible.csv"), master_sha256=sha256(folder / "master.csv"))
        dump_json(folder / "audit.json", audit)
        record_run(local, "wmh-study-prepare", {"study": study, "contract": CONTRACT})
    except (DataError, OSError, ValueError, KeyError) as exc:
        state.update(status="INPUTS_REQUIRED", error=str(exc))
    dump_json(folder / "status.json", state)
    dump_json(study_root(cfg, study) / "latest_prepared.json", {"path": str(folder), "run": stamp})
    return state


def load_prepared(cfg, study):
    pointer = read_json(study_root(cfg, study) / "latest_prepared.json")
    if not pointer:
        return None
    folder = Path(pointer["path"])
    state = read_json(folder / "status.json")
    if state.get("status") != "PREPARED" or state.get("config_hash") != config_hash(cfg):
        return None
    if state.get("source_code_sha256") != code_hashes():
        return None
    for file, key in (("eligible.csv", "cohort_sha256"), ("master.csv", "master_sha256")):
        if sha256(folder / file) != state[key]:
            raise DataError("Prepared data were changed; rerun prepare to create an auditable snapshot")
    return state


def run(cfg, study, through="prepare"):
    state = load_prepared(cfg, study) if through == "report" else None
    if state is None:
        state = prepare(cfg, study)
    if through == "prepare" or state["status"] != "PREPARED":
        return state
    folder = Path(state["path"])
    data = pd.read_csv(folder / "eligible.csv", dtype={"patient_id": str})
    master = pd.read_csv(folder / "master.csv", dtype={"patient_id": str})
    state["status"] = "ANALYSING"
    dump_json(folder / "status.json", state)
    results = run_study(data, master, study, folder, settings(cfg))
    from .reporting import report
    report(folder, state, results)
    state.update(status="COMPLETED" if results.iloc[0].status == "ESTIMATED" else "PRIMARY_NOT_ESTIMABLE",
                 finished_utc=datetime.now(UTC).isoformat())
    dump_json(folder / "status.json", state)
    dump_json(study_root(cfg, study) / "latest_results.json", {"path": str(folder), "run": state["run"]})
    return state


def read_results(cfg):
    rows = []
    for study in STUDIES:
        pointer = read_json(study_root(cfg, study) / "latest_results.json")
        row = {"study": study, "status": "NOT_RUN", "p": float("nan")}
        if pointer:
            folder = Path(pointer["path"])
            result = read_json(folder / "primary/result.json")
            version = read_json(folder / "status.json").get("contract")
            if version == CONTRACT:
                row.update(result, path=str(folder), run=pointer["run"], contract=version)
            else:
                row.update(status="PREVIOUS_VERSION", previous_contract=version, previous_path=str(folder))
        rows.append(row)
    table = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    table["p_holm_four"] = multipletests(table.p.fillna(1), method="holm")[1]
    table.loc[table.status.ne("ESTIMATED"), "p_holm_four"] = float("nan")
    return table
