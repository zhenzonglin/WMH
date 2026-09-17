"""Generate explicitly synthetic five-year inputs and exercise the production CLI.

Writes only to outputs/validation_recurrence_v3/<timestamp>; never opens real inputs.
Uses reduced MI (2 datasets, 2 iterations) for a bounded software smoke test.
"""
import argparse
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from wmh_hcy.synthetic import synthetic_frames

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--n", type=int, default=1200)
args = parser.parse_args()
repo = Path(__file__).resolve().parents[1]
root = repo / "outputs/validation_recurrence_v3" / datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
(root / "config").mkdir(parents=True)
clinical, images = synthetic_frames(args.n, seed=8311, missing=True)
source_names = ["code_n", "AGE", "GENDER", "D_DIAG", "BSL_HCY", "BSL_B12", "BSL_B9",
                "BSL_CYSC", "H_SMK", "H_DRINK", "H_HYPT", "H_DIAB", "H_STROKE",
                "ONSET_D", "I_BLDSAMP_DT", "A_NIHSS", "IMG_C_TOAST"]
clinical = clinical[source_names].copy()
rng = np.random.default_rng(9228)
h = np.log2(clinical.BSL_HCY)
h -= h.mean()
w = np.log1p(images.wmh_ml)
w = (w-w.mean())/w.std(ddof=1)
entry = (pd.to_datetime(clinical.I_BLDSAMP_DT)-pd.to_datetime(clinical.ONSET_D)).dt.days
time = entry + np.ceil(rng.exponential(1/(.0012*np.exp(.2*h+.15*w+.3*h*w))))
censor = np.where(rng.uniform(size=args.n) < .25, rng.integers(100, 1800, args.n), 1825)
stop = np.minimum(time, censor)
event = time <= censor
for year in [1, 2, 3, 4, 5]:
    status = f"Y{year}_IS" if year == 5 else f"y{year}_IS" if year > 1 else "y1_is"
    day = status + ("_DD" if year == 5 else "_dd")
    clinical[status] = (event & (time <= 365*year)).astype(int)
    clinical[day] = np.minimum(stop, 365*year)
clinical.to_csv(root / "clinical_SYNTHETIC.csv", index=False)
images.to_csv(root / "imaging_SYNTHETIC.csv", index=False)
cfg = yaml.safe_load((repo / "config/analysis.yml").read_text(encoding="utf-8"))
cfg.update(mode="synthetic", output_dir="outputs")
cfg["inputs"].update(clinical_csv=str(root / "clinical_SYNTHETIC.csv"),
                      imaging_csv=str(root / "imaging_SYNTHETIC.csv"))
cfg["analysis"].update(imputations=2, mice_iterations=2)
path = root / "config/validation.local.yml"
path.write_text(yaml.safe_dump(cfg, allow_unicode=True), encoding="utf-8")
for stage in ["prepare", "report"]:
    subprocess.run([sys.executable, "-m", "wmh_hcy.cli", "recurrence", "--config", str(path),
                    "--through", stage], check=True)
pointer = json.loads((root / "outputs/recurrence_v3/latest_results.json").read_text())
print("SYNTHETIC VALIDATION REPORT:", Path(pointer["path"]) / "report.html")
