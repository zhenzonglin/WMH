"""One entry point. Real input failures never trigger a synthetic fallback."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .common import DataError, dump_json, outdir
from .workstation import default_config, load_workstation


def main() -> int:
    parser = argparse.ArgumentParser(description="CNSR-III Hcy–WMH, Python-only observational analysis")
    parser.add_argument("command", choices=["configure", "audit", "image-audit", "run", "doctor", "extract", "prepare", "analyse", "report", "demo", "longterm"])
    parser.add_argument("--config", default=None, help="Defaults to workstation.local.yml when present, otherwise analysis.yml")
    parser.add_argument("--sas-dir", help="SAS directory, recursively scanned; configure command")
    parser.add_argument("--sustain-dir", help="Original SuStaIn project or derivatives directory; configure command")
    parser.add_argument("--id-map", help="Exact participant_id,code_n CSV mapping; configure command")
    parser.add_argument("--qc-csv", help="Legacy argument, ignored; image review is not an eligibility condition")
    parser.add_argument("--through", choices=["audit", "prepare", "analyse", "report"], default="prepare")
    parser.add_argument("--n", type=int, default=900, help="Synthetic sample size for demo only")
    parser.add_argument("--hypothesis", choices=["H1", "H2", "H3", "H4"], help="Run one hypothesis; default runs H1-H4")
    parser.add_argument("--years", type=int, nargs="+", choices=[2, 3, 4, 5],
                        default=[2, 3, 4, 5], help="Cumulative horizons for longterm; default 2 3 4 5")
    args = parser.parse_args()
    cfg = None
    try:
        if args.command == "configure":
            from .workstation import configure
            result = configure(args.config or default_config(args.command), args.sas_dir, args.sustain_dir, args.id_map, args.qc_csv)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0
        cfg = load_workstation(args.config, args.command)
        if args.command == "longterm":
            from .longterm import run_longterm
            result = run_longterm(cfg, args.years, args.through, args.hypothesis)
            print(result["status"], result.get("report", result["result_dir"]))
            return 2 if result["status"] == "COMPLETED_WITH_MODEL_FAILURES" else 0
        if args.command == "audit":
            from .audit import clinical_audit
            result = clinical_audit(cfg)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0 if result["status"] == "READY_FOR_EXTRACTION" else 2
        if args.command == "image-audit":
            from .workstation import image_audit
            result = image_audit(cfg)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0 if result["status"] == "READY_FOR_PREPARE" else 2
        if args.command == "run":
            from .workstation import run_pipeline
            result = run_pipeline(cfg, args.through, args.hypothesis)
            print(result["status"], outdir(cfg) / "pipeline_status.json")
            return 0 if result["status"] == "COMPLETED" else 2
        if args.command == "doctor":
            from .sas_extract import doctor
            result = doctor(cfg)
            print(result["status"])
            print("Missing required inputs:", result["missing_required_sources"])
            return 0 if result["status"] == "INPUTS_PRESENT" else 2
        if args.command == "extract":
            from .sas_extract import extract
            print(extract(cfg))
        elif args.command == "prepare":
            from .cohorts import prepare
            print(prepare(cfg))
        elif args.command == "analyse":
            from .analysis import analyse
            result = analyse(cfg, only=args.hypothesis)
            print(result["result_dir"])
        elif args.command == "report":
            from .reporting import report
            print(report(cfg))
        elif args.command == "demo":
            from .analysis import analyse
            from .cohorts import prepare
            from .reporting import report
            from .synthetic import create_demo
            cfg = create_demo(Path(cfg["_root"]), args.n)
            print("SYNTHETIC DATA ONLY. Reduced repetitions for software validation.", flush=True)
            print(prepare(cfg), flush=True)
            analyse(cfg)
            print(report(cfg))
        return 0
    except (DataError, FileNotFoundError) as exc:
        if cfg is not None:
            dump_json(outdir(cfg) / "logs" / f"{args.command}_error.json", {"error": str(exc), "mode": cfg["mode"]})
        print(f"INPUT/MODEL ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
