"""wmh-study: one-command preparation, analysis and reports; explicit synthetic mode."""
from __future__ import annotations

import argparse

from ..common import DataError
from ..workstation import load_workstation
from .registry import STUDIES


def main():
    parser = argparse.ArgumentParser(description="CNSR-III imaging studies; recovery is the default")
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("start", "audit", "run", "summary", "diagnose", "recovery-path"):
        p = sub.add_parser(command)
        p.add_argument("--config")
        if command not in {"summary", "diagnose", "recovery-path"}:
            p.add_argument("--study", choices=["all", *STUDIES], default="recovery")
        if command == "run":
            p.add_argument("--through", choices=["prepare", "report"], default="prepare")
        if command == "recovery-path":
            p.add_argument("--through", choices=["prepare", "analyse"], default="analyse")
        if command == "summary":
            p.add_argument("--page", type=int, choices=[1, 2], default=1)
        if command == "diagnose":
            p.add_argument("--page", type=int, choices=[1, 2, 3], default=1)
    p = sub.add_parser("demo")
    p.add_argument("--output", default="outputs/synthetic/studies_demo")
    p.add_argument("--n", type=int, default=1200)
    p.add_argument("--through", choices=["prepare", "report"], default="report")
    p.add_argument("--study", choices=["all", *STUDIES], default="recovery")
    args = parser.parse_args()
    try:
        from .reporting import print_audit, summary
        from .runner import run
        if args.command == "demo":
            from .demo import make_demo
            cfg = make_demo(args.output, n=args.n)
            for study in (STUDIES if args.study == "all" else (args.study,)):
                print_audit(run(cfg, study, args.through))
            if args.study == "all":
                summary(cfg)
        else:
            cfg = load_workstation(args.config, "wmh-study")
            if args.command == "recovery-path":
                from .recovery_path import run as run_recovery_path
                state = run_recovery_path(cfg, args.through)
                print_audit(state)
                if state["status"] not in {"PREPARED", "COMPLETED"}:
                    parser.exit(2, "Recovery path review required; see run audit and status.json\n")
            elif args.command == "diagnose":
                from .diagnostics import diagnose
                print("\n".join(diagnose(cfg, args.page)))
            elif args.command == "summary":
                summary(cfg, args.page)
            else:
                studies = STUDIES if args.study == "all" else [args.study]
                through = "report" if args.command == "start" else (
                    "prepare" if args.command == "audit" else args.through)
                failures = []
                attempts = {}
                for index, study in enumerate(studies, 1):
                    print(f"[{index}/{len(studies)}] {study} | {cfg['mode']} | through={through}", flush=True)
                    try:
                        state = run(cfg, study, through, fresh=args.command == "start")
                    except (DataError, OSError, ValueError) as exc:
                        print(f"{study}: FAILED: {exc}", flush=True)
                        attempts[study] = {"status": "FAILED", "error": str(exc)}
                        failures.append(study)
                        continue
                    attempts[study] = state
                    print_audit(state)
                    if state["status"] in {"INPUTS_REQUIRED", "REVIEW_REQUIRED", "PRIMARY_NOT_ESTIMABLE"}:
                        failures.append(study)
                if through == "report":
                    summary(cfg, attempts=attempts)
                if failures:
                    parser.exit(2, "Review required: " + ", ".join(failures) + "\n")
    except (DataError, OSError, ValueError) as exc:
        parser.exit(2, f"Stopped: {exc}\n")


if __name__ == "__main__":
    main()
