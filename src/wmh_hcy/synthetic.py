"""Synthetic, non-patient fixtures. Their event rates are not CNSR-III estimates."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .common import load_config


def synthetic_frames(n: int = 900, seed: int = 20260912, interaction: float = 0.65,
                     missing: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    ids = [f"SYN{i:06d}" for i in range(n)]
    age = np.clip(rng.normal(64, 10, n), 25, 90)
    wmh = np.exp(rng.normal(1.9, 0.85, n))
    hcy = np.exp(rng.normal(np.log(13), 0.4, n) + 0.12*(np.log1p(wmh)-2))
    b12 = np.exp(rng.normal(np.log(300), 0.35, n) - 0.15*(np.log(hcy)-np.log(13)))
    folate = np.exp(rng.normal(np.log(17), 0.4, n))
    cysc = np.exp(rng.normal(np.log(0.95), 0.22, n))
    h = np.log2(hcy)-np.log2(hcy).mean()
    w = (np.log1p(wmh)-np.log1p(wmh).mean())/np.log1p(wmh).std(ddof=1)
    stroke_time = np.ceil(rng.exponential(1/(0.0017*np.exp(0.2*h+0.3*w+interaction*h*w))))
    death_time = np.ceil(rng.exponential(1/(0.00055*np.exp(0.25*w))))
    is_event = (stroke_time <= 365) & (stroke_time <= death_time)
    dead = death_time <= 365
    onset = pd.Timestamp("2020-01-01") + pd.to_timedelta(rng.integers(0, 120, n), unit="D")
    stamp = lambda days: (onset + pd.to_timedelta(days, unit="D")).strftime("%Y-%m-%d").to_numpy()
    sample_days = rng.integers(1, 5, n)
    m3_days = rng.integers(80, 106, n)
    df = pd.DataFrame({
        "code_n": ids, "AGE": age, "GENDER": rng.integers(1, 3, n), "D_DIAG": np.ones(n, int),
        "BSL_HCY": hcy, "BSL_B12": b12, "BSL_B9": folate, "BSL_CYSC": cysc,
        "BSL_Cr": np.exp(rng.normal(np.log(75), 0.23, n)),
        "H_SMK": rng.integers(1, 5, n), "H_DRINK": rng.integers(1, 5, n),
        "H_HYPT": rng.integers(1, 3, n), "H_DIAB": rng.integers(1, 3, n),
        "H_STROKE": rng.integers(1, 3, n), "A_NIHSS": rng.poisson(5, n),
        "IMG_C_TOAST": rng.integers(1, 6, n), "H_MRS": rng.choice([0, 1, 2, 3], n),
        "ONSET_D": onset.strftime("%Y-%m-%d"), "I_BLDSAMP_DT": stamp(sample_days),
        "F3_BLDSAMP_D": stamp(m3_days),
        "M03_HCY": hcy*np.exp(rng.normal(-0.05, 0.2, n)),
        "M03_B12": b12*np.exp(rng.normal(0.05, 0.2, n)),
        "M03_B9": folate*np.exp(rng.normal(0.05, 0.2, n)),
        "M03_CYSC": cysc*np.exp(rng.normal(0, 0.1, n)),
        "y1_is": is_event.astype(int), "y1_is_dd": np.where(is_event, stroke_time, np.nan),
        "D_DEATH": np.where(dead & (death_time < 10), 1, 2),
        "D_DEATH_D": np.where(dead & (death_time < 10), stamp(death_time), ""),
    })
    latent = 0.2*h + 0.3*w + rng.logistic(size=n)
    df["F12_MRS"] = np.digitize(latent, [-1.6, -0.7, 0.1, 0.9, 1.6]).astype(float)
    df.loc[dead, "F12_MRS"] = np.nan
    for v, day in [(3, 90), (6, 180), (12, 365)]:
        df[f"F{v}_DATE"] = stamp(np.full(n, day))
        df[f"F{v}_DEATH"] = np.where(death_time <= day, 1, 2)
        df[f"F{v}_DEATH_D"] = np.where(death_time <= day, stamp(death_time), "")
    img = pd.DataFrame({"participant_id": ids, "wmh_ml": wmh,
                        "wmh_raw_ml": wmh*np.exp(rng.normal(0.04, 0.1, n)),
                        "icv_ml": np.clip(rng.normal(1450, 110, n), 1000, 2000),
                        "gm119_ml": rng.normal(610, 45, n)-1.5*(age-64)-3*w,
                        "lesion_ml": rng.exponential(6, n),
                        "wmh_source": "SYNTHETIC_GENERATOR", "icv_source": "SYNTHETIC_GENERATOR"})
    if missing:
        # MAR depends on observed age; missingness is not a real cohort characteristic.
        for col in ["BSL_B12", "BSL_B9", "BSL_CYSC", "H_SMK", "M03_B12"]:
            absent = rng.random(n) < (0.025 + 0.025*(age > 70))
            df.loc[absent, col] = np.nan
    return df, img


def create_demo(root: Path, n: int = 900) -> dict:
    folder = root / "examples/synthetic"
    folder.mkdir(parents=True, exist_ok=True)
    clinical, images = synthetic_frames(n=n)
    clinical.to_csv(folder / "clinical.csv", index=False)
    images.to_csv(folder / "imaging.csv", index=False)
    cfg = yaml.safe_load((root / "config/analysis.yml").read_text())
    cfg["mode"] = "synthetic"
    cfg["output_dir"] = "outputs/demo"
    cfg["inputs"]["clinical_csv"] = "examples/synthetic/clinical.csv"
    cfg["inputs"]["imaging_csv"] = "examples/synthetic/imaging.csv"
    cfg["analysis"].update({"imputations": 2, "mice_iterations": 3,
                            "bootstrap_per_imputation": 3, "risk_hcy_grid_points": 5})
    target = root / "config/demo.yml"
    target.write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))
    (folder / "README.txt").write_text(
        "SYNTHETIC DATA ONLY. Not patients, not CNSR-III event rates.\n"
        "Generated CSV is used to exercise the complete statistical pipeline.\n"
        "SAS decoding is independently tested on public SAS-format fixtures.\n")
    return load_config(target)
