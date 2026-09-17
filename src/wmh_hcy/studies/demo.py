"""Explicitly synthetic fixtures; never used as a fallback for real analysis."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from ..common import load_config
from .registry import SOURCES


def synthetic_tables(n=1200, seed=20260917, missing=True):
    rng = np.random.default_rng(seed)
    aliases = {}
    for (name, kind, codes) in SOURCES.values():
        if kind == "id":
            aliases[name] = [f"{i:07d}" for i in range(n)]
        elif kind in {"date", "datetime"}:
            aliases[name] = [""]*n
        elif codes:
            aliases[name] = rng.choice(codes, n).astype(float)
        else:
            aliases[name] = rng.normal(10, 1, n)
    age = np.clip(rng.normal(63, 10, n), 30, 92)
    cec = rng.normal(12, 2, n)
    wmh = np.exp(2.4+.02*(age-63)+rng.normal(0, .6, n))
    gm = 550-2*(age-63)+4*(cec-12)+rng.normal(0, 35, n)
    lesion = np.exp(rng.normal(1.3, .8, n))
    mrs3 = rng.choice([0, 1, 2, 3, 4, 5], n, p=[.25, .3, .27, .09, .06, .03])
    w = (np.log1p(wmh)-np.log1p(wmh).mean())/np.log1p(wmh).std()
    pdeath = 1/(1+np.exp(-(-2+.025*(age-63)+.2*w)))
    died = rng.random(n) < pdeath
    death_time = np.where(died, rng.uniform(200, 1800, n), np.inf)
    visit3, visit12 = rng.integers(80, 110, n), rng.integers(340, 400, n)
    sbp3 = rng.normal(135, 15, n)
    recurrent = rng.exponential(4200*np.exp(-.012*(sbp3-135)-.15*w), n)
    observed = np.minimum(np.where(rng.random(n) < .15, rng.uniform(500, 1825, n), 1825), death_time)
    stop = np.minimum(recurrent, observed)
    event = (recurrent <= observed).astype(int)
    u0 = np.exp(rng.normal(.7, 1.1, n))
    u3 = np.exp(.7+.55*(np.log(u0)-.7)+rng.normal(0, .9, n))
    baseline = pd.Timestamp("2016-01-01")

    def date(days):
        return (baseline+pd.to_timedelta(days, unit="D")).strftime("%Y-%m-%d").to_numpy()

    aliases.update(age=age, diagnosis=np.ones(n), onset_date=date(np.zeros(n)), sample_date=date(rng.integers(0, 7, n)),
                   visit3_date=date(visit3), visit12_date=date(visit12), pre_mrs=rng.choice(range(6), n),
                   nihss=rng.poisson(4, n), mrs3=mrs3, discharge_mrs=np.minimum(5, mrs3+rng.integers(0, 3, n)),
                   education=rng.choice(range(1, 6), n), bmi=rng.normal(24, 3, n),
                   cysc=np.exp(rng.normal(0, .2, n)), cysc3=np.exp(rng.normal(0, .2, n)),
                   ldl=np.exp(rng.normal(1, .2, n)), hdl=np.exp(rng.normal(.2, .2, n)),
                   tg=np.exp(rng.normal(.4, .4, n)), apo_ai=rng.normal(1.3, .15, n),
                   cec=cec, uacr0=u0, uacr3=u3, hospital_death=np.full(n, 2), sbp0=sbp3+rng.normal(10, 10, n),
                   cer24=np.exp(rng.normal(5, .4, n)), cer20=np.exp(rng.normal(3, .4, n)),
                   prior_statin=rng.binomial(1, .3, n), prior_lipid_med=np.full(n, 2),
                   discharge_bp_med=np.full(n, 2))
    # Ensure medication parent has both levels, with coherent structural negatives.
    aliases["discharge_bp_med"] = rng.choice([1, 2], n, p=[.2, .8])
    aliases["discharge_acei"] = np.where(aliases["discharge_bp_med"] == 1, 0, rng.binomial(1, .35, n))
    aliases["discharge_arb"] = np.where(aliases["discharge_bp_med"] == 1, 0, rng.binomial(1, .35, n))
    aliases["cer16"] = aliases["cer24"]*np.exp(-1+.15*w+rng.normal(0, .5, n))
    aliases["cer241"] = aliases["cer24"]*np.exp(-.5+.1*w+rng.normal(0, .4, n))
    for month, visit in [(3, visit3), (12, visit12)]:
        pressure = sbp3 if month == 3 else sbp3+rng.normal(0, 8, n)
        aliases.update({f"lsbp{month}": pressure, f"rsbp{month}": pressure+rng.normal(0, 4, n),
                        f"ldbp{month}": rng.normal(78, 5, n), f"rdbp{month}": rng.normal(79, 5, n)})
    for month in (3, 6, 12, 24, 36, 48, 60):
        horizon = 365*month/12
        is_dead = death_time <= horizon
        aliases[f"death{month}"] = np.where(is_dead, 1, 2)
        if month <= 12:
            aliases[f"death{month}_date"] = np.where(is_dead, date(np.minimum(death_time, 2000)), "")
        if month >= 12:
            pdep = 1/(1+np.exp(-(-1.5+.3*w-.012*(gm-550)+.2*mrs3+.2*(u3 >= 3))))
            dependent = rng.random(n) < pdep
            score = np.where(dependent, rng.integers(3, 6, n), rng.integers(0, 3, n)).astype(float)
            score[is_dead] = 6
            if missing:
                score[(rng.random(n) < .08) & ~is_dead] = np.nan
            aliases[f"mrs{month}"] = score
    for year in range(1, 6):
        prefix = "is" if year == 1 else f"y{year}_is"
        aliases[prefix+"_event"] = (event & (stop <= 365*year)).astype(float)
        aliases[prefix+"_day"] = np.minimum(stop, 365*year)
    frame = pd.DataFrame({s: aliases[a] for s, (a, _, _) in SOURCES.items()})
    if missing:
        for source in ("BMI", "EDUC", "H_HYPT", "BSL_CYSC", "BSL_TG", "M03_CYSC"):
            frame.loc[rng.random(n) < .025, source] = np.nan
        for source in ("CEC", "BSL_Cer_16_0", "M03_UACR"):
            frame.loc[rng.random(n) < .05, source] = np.nan
    images = pd.DataFrame({"participant_id": aliases["patient_id"], "wmh_ml": wmh, "wmh_raw_ml": wmh+.05*lesion,
                               "gm119_ml": gm, "lesion_ml": lesion, "icv_ml": rng.normal(1350, 90, n)})
    return frame, images


def make_demo(output, n=1200):
    root = Path(output).resolve()
    if (root / "config/demo.yml").exists():
        raise ValueError("Synthetic destination already exists; select a new --output directory")
    (root / "config").mkdir(parents=True, exist_ok=True)
    clinical, images = synthetic_tables(n)
    clinical.to_csv(root / "clinical.csv", index=False)
    images.to_csv(root / "imaging.csv", index=False)
    config = {"mode": "synthetic", "output_dir": str(root / "outputs"),
                  "inputs": {"clinical_csv": str(root / "clinical.csv"), "imaging_csv": str(root / "imaging.csv"), "sas_globs": []},
                  "sas": {"encoding": None, "chunksize": 100000, "date_formats": {}}, "variable_sources": {},
                  "analysis": {"imputations": 2, "mice_iterations": 2, "mice_threads": 2, "seed": 20260917}}
    path = root / "config/demo.yml"
    path.write_text(yaml.safe_dump(config, allow_unicode=True), encoding="utf-8")
    return load_config(path)
