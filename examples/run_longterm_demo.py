"""Synthetic-only 2-5 year integration demonstration; never reads patient files."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from wmh_hcy.common import load_config
from wmh_hcy.longterm import run_longterm
from wmh_hcy.synthetic import synthetic_frames


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=600)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    clinical, images = synthetic_frames(args.n, missing=True)
    rng = np.random.default_rng(20260917)
    n = len(clinical)
    late_stroke = 365 + rng.exponential(1900, n)
    stroke = clinical.y1_is_dd.where(clinical.y1_is.eq(1), late_stroke).to_numpy()
    early_death = (
        pd.to_datetime(clinical.F12_DEATH_D, errors="coerce") - pd.to_datetime(clinical.ONSET_D)
    ).dt.days.astype(float)
    death = early_death.fillna(pd.Series(365 + rng.exponential(7500, n))).to_numpy()
    lost = np.where(rng.random(n) < 0.12, rng.uniform(366, 1825, n), 3000)
    for year in [2, 3, 4, 5]:
        horizon = year * 365
        stop = np.floor(np.minimum(np.minimum(death, lost), horizon))
        event = stroke <= stop
        for endpoint in ["IS", "STROKE"]:
            source = f"Y{year}_{endpoint}" if year == 5 else f"y{year}_{endpoint}"
            clinical[source] = event.astype(int)
            clinical[source + ("_DD" if year == 5 else "_dd")] = np.where(event, stroke, stop)
        source = f"Y{year}_HS" if year == 5 else f"y{year}_HS"
        clinical[source] = 0
        clinical[source + ("_DD" if year == 5 else "_dd")] = stop
        scores = np.digitize(
            0.3 * np.log2(clinical.BSL_HCY / 13) + 0.2 * np.log1p(images.wmh_ml) + rng.logistic(size=n),
            [-1.4, -0.5, 0.2, 1.0, 1.8],
        ).astype(float)
        scores[death <= horizon] = 6
        scores[(lost < horizon) & (death > lost)] = np.nan
        clinical[f"m{year * 12}_mrs"] = scores
    data = root / "examples/synthetic/longterm"
    data.mkdir(parents=True, exist_ok=True)
    clinical.to_csv(data / "clinical.csv", index=False)
    images.to_csv(data / "imaging.csv", index=False)
    cfg = yaml.safe_load((root / "config/analysis.yml").read_text(encoding="utf-8"))
    cfg.update(mode="synthetic", output_dir="outputs/demo_longterm")
    cfg["inputs"].update(clinical_csv=str(data / "clinical.csv"), imaging_csv=str(data / "imaging.csv"))
    cfg["analysis"].update(imputations=2, mice_iterations=2, bootstrap_per_imputation=2)
    path = root / "config/demo_longterm.local.yml"
    path.write_text(yaml.safe_dump(cfg, allow_unicode=True), encoding="utf-8")
    print("SYNTHETIC SOFTWARE VERIFICATION ONLY; not CNSR-III findings", flush=True)
    result = run_longterm(load_config(path), through="report")
    print(result["status"], result["report"])
    return 0 if result["status"] == "COMPLETED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
