"""Each study imputes only its own scattered, non-image covariate missingness."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from ..common import DataError
from ..design import rcs_nonlinear
from ..imputation import nelson_aalen_increment
from .design import raw_transform
from .registry import CATEGORIES, IMAGE_COLUMNS, LOG_COLUMNS


def impute(data, design, settings, progress=print):
    spec = design.spec
    d = data.reset_index(drop=True).copy()
    for col in spec.predictors:
        if d[col].notna().sum() == 0:
            raise DataError(f"Entirely unmeasured {col}; cannot impute")
        if (col in spec.exposures or col in IMAGE_COLUMNS) and d[col].isna().any():
            raise DataError(f"Exposure/image missing: {col}; cannot impute")
    columns, bases = {}, {}
    for c in spec.predictors:
        v = d[c]
        if c in CATEGORIES:
            columns[c] = pd.Categorical(v, categories=design.coding[c]["levels"])
            for level in design.coding[c]["levels"][1:]:
                bases[f"{c}_{level:g}"] = v.eq(level).astype(float).where(v.notna())
        else:
            columns[c] = raw_transform(v, c)
            code = design.coding[c]
            bases[c] = (columns[c]-code["center"])/code["scale"]
    # Observed spline/interaction bases inform FCS. Covariate bases are recomputed
    # from each completed dataset for fitting; this remains approximate FCS.
    for c in spec.splines:
        if d[c].notna().all():
            code = design.coding[c]
            z = (raw_transform(d[c], c)-code["center"])/code["scale"]
            columns["aux_"+c+"_rcs"] = rcs_nonlinear(z.to_numpy(), code["knots"])
            bases[c+"_rcs"] = columns["aux_"+c+"_rcs"]
    for left, right in spec.interactions:
        if left not in bases or right not in bases:
            raise DataError(f"Observed interaction basis unavailable for MI: {left} x {right}")
        interaction = np.asarray(bases[left])*np.asarray(bases[right])
        if not np.isfinite(interaction).all():
            raise DataError("Interaction auxiliaries must be based on observed exposures")
        columns[f"aux_{left}_x_{right}"] = interaction
    if spec.family == "cox":
        columns.update(aux_event=d.event_type, aux_na=nelson_aalen_increment(d, 1),
                       aux_entry=d.entry, aux_exit=d.exit)
    elif spec.family in {"multinomial", "binary"}:
        for state in (0, 1, 2):
            columns[f"aux_state_{state}"] = d[spec.outcome].eq(state).astype(int)
        columns["aux_outcome_observed"] = d[spec.outcome].notna().astype(int)
    else:
        if d[spec.outcome].isna().any():
            raise DataError("Missing structural outcome cannot be imputed")
        columns["aux_outcome"] = d[spec.outcome]
    working = pd.DataFrame(columns).reset_index(drop=True)
    missing = {c: int(d[c].isna().sum()) for c in spec.covariates}
    targets = {c: [v for v in working if v != c] for c, count in missing.items() if count}
    if not targets:
        return [d], {"method": "no_missing_covariates", "m": 1, "missing": missing}
    import miceforest as mf

    m, iterations = int(settings.get("imputations", 50)), int(settings.get("mice_iterations", 10))
    if m < 2 or iterations < 1:
        raise DataError("MI requires >=2 datasets and >=1 iteration")
    kernel = mf.ImputationKernel(working, num_datasets=m, variable_schema=targets,
                                 mean_match_candidates=min(5, min(d[c].notna().sum() for c in targets)),
                                 random_state=int(settings.get("seed", 20260917)))
    traces = []
    with warnings.catch_warnings():
        # A pandas memory-layout warning is not a statistical diagnostic.
        warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning, module="miceforest.*")
        for iteration in range(iterations):
            kernel.mice(1, num_threads=int(settings.get("mice_threads", 2)), num_iterations=60, verbosity=-1)
            progress(f"  MI {iteration+1}/{iterations}; datasets={m}")
            for j in range(m):
                completed = kernel.complete_data(dataset=j)
                traces.append(dict(iteration=iteration+1, imputation=j,
                                   **{c: float(pd.to_numeric(completed.loc[working[c].isna(), c]).mean()) for c in targets}))
    datasets = []
    for j in range(m):
        filled = kernel.complete_data(dataset=j)
        result = d.copy()
        for c in targets:
            v = pd.to_numeric(filled[c]).astype(float)
            if c in LOG_COLUMNS:
                v = 2**v
            elif c == "lesion_ml":
                v = np.expm1(v)
            result[c] = result[c].where(result[c].notna(), v)
        datasets.append(result)
    return datasets, {"method": "study_specific_miceforest_approximate_FCS", "m": m,
                      "iterations": iterations, "seed": int(settings.get("seed", 20260917)),
                      "missing": missing, "predictors": targets, "transformed_mean_traces": traces,
                      "assumption": "MAR; approximate model compatibility; outcomes/exposures/images never imputed"}
