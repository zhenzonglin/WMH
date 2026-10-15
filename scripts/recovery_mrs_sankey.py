"""Independent descriptive mRS Sankey; approved publication figures stay untouched."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wmh_hcy.recovery_publication.sankey_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
