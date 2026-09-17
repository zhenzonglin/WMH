"""One-screen, read-only aggregate summary; never reads patient-level rows."""
import argparse
from pathlib import Path

import pandas as pd

from wmh_hcy.common import read_json
from wmh_hcy.workstation import load_workstation

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--config")
parser.add_argument("--results", type=Path)
args = parser.parse_args()
cfg = load_workstation(args.config, "recurrence")
pointer = read_json(Path(cfg["_out"]) / "recurrence_v3/latest_run.json")
root = args.results or (Path(pointer["path"]) if pointer else None)
if root is None:
    parser.error("No recurrence v3 run found. Run wmh-hcy recurrence --through prepare first.")
state = read_json(root / "status.json")
print(f"RECURRENCE V3 | mode={state.get('mode')} | status={state.get('status')} | run={root.name}")
print("READ ONLY: aggregate files only; no patient IDs, no fitting.")
print("Primary=60 months; sensitivity=3,6,12,24,36,48; all from Y5_IS/Y5_IS_DD.")
if state.get("reason"):
    print("FAILURE:", state["reason"])
for name, cols in [
    ("cohort_counts", ["month", "n", "events", "entry_at_or_after_cutoff", "censored_before_cutoff"]),
    ("horizon_results", ["month", "status", "HR", "HR_lower", "HR_upper", "p", "p_holm_6"]),
    ("robustness_results", ["analysis", "n", "status", "HR", "HR_lower", "HR_upper", "p", "reason"]),
]:
    path = root / f"{name}.csv"
    if path.is_file():
        print(name.upper())
        frame = pd.read_csv(path)
        print(frame[[c for c in cols if c in frame]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
print("Report:", root / "report.html")
