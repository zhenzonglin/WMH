"""Standalone publication entry point; no analysis command is invoked."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wmh_hcy.recovery_publication.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
