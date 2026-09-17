"""Prespecified analyses, one primary test per study and explicitly named extensions."""
from __future__ import annotations

import traceback
from dataclasses import asdict

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.special import expit, logit
from statsmodels.stats.multitest import multipletests

from ..common import DataError, dump_json
from ..imputation import pool_scalar
from .data import build_cohort, volume_ok
from .design import StudyDesign, split_time
from .imputation import impute
from .models import fit, state_probabilities
from .pooling import coefficients, contrast, terms_test
from .registry import COMMON, ModelSpec, primary_spec


def clinical_probabilities(datasets, design, fits):
    """Standardization over the observed analysis cohort; pointwise conditional CIs."""
    spec, grid = design.spec, []
    if spec.study == "recovery":
        for w in np.quantile(datasets[0].wmh_ml, [.25, .5, .75]):
            for g in np.quantile(datasets[0].gm119_ml, [.25, .5, .75]):
                grid.append({"wmh_ml": float(w), "gm119_ml": float(g)})
    elif spec.study == "kidney":
        grid = [{"albuminuria": float(v)} for v in design.coding["albuminuria"]["levels"]]
    else:
        name = spec.exposures[0]
        grid = [{name: float(v)} for v in np.quantile(datasets[0][name], [.25, .5, .75])]
    rows = []
    for changes in grid:
        values = []
        for d, f in zip(datasets, fits, strict=True):
            new = d.assign(**changes)
            values.append(state_probabilities(f, design.transform(new)))
        for s, label in enumerate(("independent", "dependent", "dead")):
            q = np.clip([v[s][0] for v in values], 1e-8, 1-1e-8)
            pooled = pool_scalar(q, [v[s][1] for v in values])
            probability = pooled["estimate"]
            margin = stats.t.ppf(.975, pooled["df"])*pooled["se"]/(probability*(1-probability))
            rows.append({**changes, "state": label, "probability": probability,
                         "lower": float(expit(logit(probability)-margin)), "upper": float(expit(logit(probability)+margin)),
                         "ci": "mean probability; Rubin variance; logit-delta pointwise CI; covariate distribution fixed"})
    return pd.DataFrame(rows)


def bp_contrasts(datasets, design, fits):
    observed = datasets[0]
    rows = []
    # Local support uses the WMH stratum around each display percentile.
    quantiles = np.quantile(observed.wmh_ml, [.125, .375, .625, .875])
    for index, percentile in enumerate((25, 50, 75)):
        w = float(observed.wmh_ml.quantile(percentile/100))
        local = observed.loc[observed.wmh_ml.between(quantiles[index], quantiles[index+1]), "sbp3"]
        for sbp in (120, 130, 150):
            supported = len(local) >= 20 and local.min() <= min(sbp, 140) and local.max() >= max(sbp, 140)
            row = {"wmh_percentile": percentile, "wmh_ml": w, "sbp": sbp, "reference_sbp": 140,
                       "local_n": len(local), "support_min": float(local.min()), "support_max": float(local.max())}
            if not supported:
                rows.append({**row, "status": "OUTSIDE_OBSERVED_SUPPORT"})
                continue
            vectors = [design.contrast(d, {"sbp3": sbp, "wmh_ml": w}, {"sbp3": 140, "wmh_ml": w}) for d in datasets]
            r = contrast(fits, vectors)
            rows.append({**row, **r, "HR": np.exp(r["estimate"]), "HR_lower": np.exp(r["lower"]),
                         "HR_upper": np.exp(r["upper"]), "status": "ESTIMATED"})
    return pd.DataFrame(rows)


def observation_weights(data, design, outcome):
    """Prespecified logistic observation model; no future predictors or p selection."""
    indicator = data[outcome].notna().astype(int).to_numpy()
    if indicator.min() == 1:
        return np.ones(len(data)), {"method": "all_outcomes_observed", "ess": len(data)}, None
    x = design.transform(data).to_numpy()
    model = sm.Logit(indicator, x).fit(method="newton", maxiter=150, disp=False)
    if not model.mle_retvals.get("converged", False) or not np.isfinite(model.params).all():
        raise DataError("Outcome-observation model failed")
    probability = np.asarray(model.predict(x))
    if probability.min() <= 0 or probability.max() >= 1:
        raise DataError("Observation model separation/positivity failure")
    raw = 1/probability[indicator == 1]
    # No truncation by default; expose the price of positivity violations.
    diag = {"method": "untruncated_inverse_observation_probability", "minimum_probability": float(probability.min()),
                "weight_max": float(raw.max()), "weight_p99": float(np.quantile(raw, .99)),
                "ess": float(raw.sum()**2/np.sum(raw**2)),
                "variance": "stacked sandwich includes estimated observation weights"}
    return raw, diag, (x, indicator, probability)


def run_one(data, spec, directory, settings, inherited=None, complete_case=False, ipw=False, time_split=False):
    path = directory / spec.name
    path.mkdir(parents=True, exist_ok=True)
    dump_json(path / "model_definition.json", asdict(spec))
    result = {"analysis": spec.name, "study": spec.study, "tier": spec.tier, "family": spec.family, "status": "NOT_ESTIMABLE"}
    try:
        working = data.copy()
        # Eligibility never comes from imputation of exposures or imaging.
        needed = list(spec.exposures)
        needed += [c for c in spec.covariates if c.endswith("_ml")]
        if not ipw:
            needed += [spec.outcome]
        if complete_case:
            needed += list(spec.covariates)
        working = working.dropna(subset=list(dict.fromkeys(needed))).reset_index(drop=True)
        if len(working) < 10:
            raise DataError("Fewer than 10 eligible observed outcomes; model not fitted")
        design = StudyDesign.freeze(data, spec, inherited)
        dump_json(path / "frozen_design.json", design.coding)
        result.update(n=len(working), outcome_unknown=int(working[spec.outcome].isna().sum()))
        datasets, mi = ([working], {"method": "complete_case", "m": 1}) if complete_case else impute(working, design, settings)
        dump_json(path / "imputation.json", mi)
        fits, model_data, weighting = [], [], []
        for index, d in enumerate(datasets):
            weight, observation = None, None
            if ipw:
                weight, diag, observation = observation_weights(d, design, spec.outcome)
                weighting.append(diag)
                d = d.loc[d[spec.outcome].notna()].copy()
            if time_split:
                d = split_time(d)
            fits.append(fit(d, design, weights=weight, observation=observation))
            model_data.append(d)
            if (index+1) % 10 == 0:
                print(f"  {spec.name}: fitted {index+1}/{len(datasets)}", flush=True)
        coefficients(fits, spec.family).to_csv(path / "coefficients.csv", index=False)
        np.savez_compressed(path / "pooled_inputs.npz", terms=np.asarray(fits[0].terms),
                            params=np.asarray([f.params for f in fits]),
                            covariance=np.asarray([f.covariance for f in fits]))
        dump_json(path / "fit_diagnostics.json", [f.diagnostics for f in fits])
        if weighting:
            dump_json(path / "observation_weights.json", weighting)
        test = terms_test(fits, spec.primary)
        result.update(status="ESTIMATED", n=len(model_data[0].patient_id.unique()),
                      parameters=len(fits[0].params), imputations=len(fits), **test)
        if len(spec.primary) == 1:
            row = coefficients(fits, spec.family).set_index("term").loc[spec.primary[0]].to_dict()
            result.update({k: row[k] for k in ("estimate", "lower", "upper", "ratio", "ratio_lower", "ratio_upper") if k in row})
        model_data[0][["patient_id", spec.outcome] + (["entry", "exit"] if spec.family == "cox" else [])].to_csv(path / "model_membership.csv", index=False)
        if time_split:
            base_terms = [t.removesuffix("_late") for t in spec.primary]
            early = terms_test(fits, base_terms)
            vectors = []
            for t in base_terms:
                vector = np.zeros(len(fits[0].terms))
                vector[fits[0].terms.index(t)] = 1
                vector[fits[0].terms.index(t+"_late")] = 1
                vectors.append(vector)
            a = np.array(vectors)
            from .pooling import joint_test
            late = joint_test([a @ f.params for f in fits], [a @ f.covariance @ a.T for f in fits])
            dump_json(path / "period_tests.json", {"early_joint": early, "late_joint": late, "difference_joint": test,
                                                  "difference_test_is_primary_for_this_sensitivity": True})
            period_contrasts = []
            for phase in (0, 1):
                completed = [d.assign(late_period=phase) for d in datasets]
                curve = bp_contrasts(completed, design, fits)
                curve["period"] = "visit_to_day365" if phase == 0 else "day365_to1825"
                period_contrasts.append(curve)
            pd.concat(period_contrasts).to_csv(path / "period_clinical_contrasts.csv", index=False)
        if spec.name == "primary":
            if spec.family == "multinomial":
                clinical_probabilities(model_data, design, fits).to_csv(path / "standardized_states.csv", index=False)
            if spec.family == "cox":
                bp_contrasts(model_data, design, fits).to_csv(path / "clinical_contrasts.csv", index=False)
        return result, design
    except (DataError, ValueError, np.linalg.LinAlgError, FloatingPointError, KeyError) as exc:
        result.update(reason=str(exc))
        (path / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return result, inherited
    finally:
        dump_json(path / "result.json", result)


def functional_spec(base, name, add_image=(), covariates=None):
    cov = tuple(covariates or base.covariates)
    extra = tuple(c for c in ("pre_mrs", "nihss", "toast") if c not in cov)
    exposure = base.exposures
    return base.variant(name, family="multinomial", outcome="state60", tier="secondary",
                        covariates=cov+extra+tuple(add_image),
                        primary=("dependent:"+exposure[0],))


def run_study(data, master, study, directory, settings):
    primary = primary_spec(study)
    rows = []

    def run(d, spec, **kwargs):
        print(f"{study} / {spec.name}", flush=True)
        row, design = run_one(d, spec, directory, settings, **kwargs)
        rows.append(row)
        pd.DataFrame(rows).to_csv(directory / "results.csv", index=False)
        return design

    frozen = run(data, primary)
    # Failure of one model never causes a data-driven change to its specification.
    run(data, primary.variant("complete_case"), inherited=frozen, complete_case=True)
    if study == "recovery":
        for month in (12, 24, 36, 48):
            run(data.loc[~data[f"state{month}_conflict"]], primary.variant(f"month{month}", outcome=f"state{month}", tier="secondary"), inherited=frozen)
        run(data.loc[data.discharge_mrs.gt(data.mrs3)], primary.variant("actual_improvement"), inherited=frozen)
        raw = data.loc[volume_ok(data, "wmh_raw_ml")].copy()
        raw["wmh_ml"] = raw.wmh_raw_ml
        run(raw, primary.variant("raw_wmh"), inherited=frozen)
        binary = data.copy()
        binary["adverse60"] = binary.state60.gt(0).astype(float).where(binary.state60.notna())
        run(binary, primary.variant("dependent_or_dead", family="binary", outcome="adverse60", primary=("wmh_ml", "gm119_ml")), inherited=frozen)
        run(data, primary.variant("observation_weighted"), inherited=frozen, ipw=True)
    elif study == "bp":
        mean = data.copy()
        mean["sbp3"] = mean.sbp3_mean
        run(mean, primary.variant("mean_arms"), inherited=frozen)
        run(data.loc[data.severe_wmh.eq(1)], primary.variant("severe_wmh"), inherited=frozen)
        landmark, _, _ = build_cohort(master, "bp", month=12)
        run(landmark, primary.variant("month12_landmark"), inherited=frozen)
        run(data, primary.variant("time_varying", primary=tuple(t+"_late" for t in primary.primary)),
            inherited=frozen, time_split=True)
    elif study == "ceramide":
        run(data, primary.variant("alternative_ratio", exposures=("cer_alt",), primary=("cer_alt",), tier="secondary"), inherited=frozen)
        gm = data.loc[volume_ok(data, "gm119_ml")]
        run(gm, primary.variant("gray_matter", outcome="gm119_ml", tier="secondary"), inherited=frozen)
        run(data, primary.variant("acute_lesion", outcome="log_lesion", tier="secondary",
                                  covariates=tuple(c for c in primary.covariates if c != "lesion_ml")), inherited=frozen)
        run(data, primary.variant("separate_components", exposures=("log_c16", "log_c24"),
                                  primary=("log_c16", "log_c24"), tier="secondary"), inherited=frozen)
        valid = data.loc[~data.state60_conflict]
        run(valid, functional_spec(primary, "function60"), inherited=frozen)
        run(valid, functional_spec(primary, "function60_with_wmh", ("wmh_ml",)), inherited=frozen)
        run(valid, functional_spec(primary, "function60_observation_weighted").variant("function60_observation_weighted"), inherited=frozen, ipw=True)
    elif study == "cec":
        wmh = data.loc[volume_ok(data, "wmh_ml")]
        run(wmh, primary.variant("white_matter", outcome="log_wmh", tier="secondary"), inherited=frozen)
        run(data, primary.variant("apo_ai_adjusted", covariates=primary.covariates+("apo_ai",), tier="secondary"), inherited=frozen)
        valid = data.loc[~data.state60_conflict]
        run(valid, functional_spec(primary, "function60"), inherited=frozen)
        run(valid, functional_spec(primary, "function60_with_gm", ("gm119_ml",)), inherited=frozen)
        run(valid, functional_spec(primary, "function60_observation_weighted").variant("function60_observation_weighted"), inherited=frozen, ipw=True)
        run(data, primary.variant("cec_spline", splines=("age", "cec"), primary=("cec", "cec_rcs")), inherited=frozen)
    else:
        # Baseline-only structural question: no future UACR, BP or function in C.
        baseline = master.loc[master.age.ge(18) & master.diagnosis.eq(1) & volume_ok(master, "wmh_ml") & volume_ok(master, "icv_ml")].copy()
        baseline_spec = ModelSpec("kidney", name="baseline_structure", family="ols", outcome="log_wmh",
                                 exposures=("uacr0_log",), primary=("uacr0_log",), tier="secondary",
                                 covariates=COMMON+("bmi", "education", "cysc", "sbp0", "icv_ml"))
        run(baseline, baseline_spec, inherited=frozen)
        run(data, primary.variant("continuous_uacr3", exposures=("wmh_ml", "uacr3_log"),
                                  covariates=primary.covariates+("uacr0_log",), primary=("dependent:uacr3_log",), tier="secondary"), inherited=frozen)
        # Nested improvement is a pooled Wald test, not a claim of predictive performance.
        run(data, primary.variant("incremental_albuminuria", primary=tuple(f"dependent:albuminuria_{i}" for i in (1, 2, 3)), tier="secondary"), inherited=frozen)
        run(data, primary.variant("clinical_wmh_only", exposures=("wmh_ml",), primary=("dependent:wmh_ml",), tier="secondary"), inherited=frozen)
        subset = data.loc[volume_ok(data, "lesion_ml")].copy()
        run(subset, primary.variant("clinical_same_subset_core"), inherited=frozen)
        run(subset, primary.variant("clinical_extended", covariates=primary.covariates+("nihss", "toast", "lesion_ml")), inherited=frozen)
        run(data, primary.variant("observation_weighted"), inherited=frozen, ipw=True)
    result = pd.DataFrame(rows)
    secondary = result.tier.eq("secondary")
    if secondary.any():
        # Fixed family includes non-estimable planned tests, represented by p=1 for adjustment only.
        result.loc[secondary, "p_bh_secondary"] = multipletests(result.loc[secondary, "p"].fillna(1), method="fdr_bh")[1]
        result.loc[secondary & result.status.ne("ESTIMATED"), "p_bh_secondary"] = np.nan
    result.to_csv(directory / "results.csv", index=False)
    return result
