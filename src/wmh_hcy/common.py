"""Configuration, typed failures and reproducible file outputs."""
from __future__ import annotations

import hashlib
import json
import platform
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


class DataError(ValueError):
    """An input/model condition requiring an explicit correction, never silent fallback."""


def dump_json(path: Path, obj: object) -> None:
    def convert(x: object) -> object:
        if isinstance(x, (np.integer, np.floating, np.bool_)):
            return x.item()
        if isinstance(x, (Path, datetime)):
            return str(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
        raise TypeError(type(x).__name__)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=convert), encoding="utf-8")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}


def load_config(path: str | Path) -> dict:
    p = Path(path).resolve()
    cfg = yaml.safe_load(p.read_text(encoding="utf-8"))
    root = p.parent.parent
    cfg["_root"] = str(root)
    cfg["_config"] = str(p)
    cfg["_out"] = str(resolve(cfg, cfg.get("output_dir", "outputs/real")))
    if cfg.get("mode") not in {"real", "synthetic"}:
        raise DataError("mode must be real or synthetic")
    if "site_code" in cfg.get("variable_sources", {}):
        raise DataError("Center fields are not part of this protocol")
    return cfg


def resolve(cfg: dict, value: str | Path) -> Path:
    p = Path(value).expanduser()
    return p if p.is_absolute() else Path(cfg["_root"]) / p


def outdir(cfg: dict) -> Path:
    p = Path(cfg["_out"])
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_csv(path: Path) -> pd.DataFrame:
    # Character IDs and literal NA strings are preserved at ingestion.
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")


def unique_ids(frame: pd.DataFrame, column: str, source: str) -> None:
    if column not in frame:
        raise DataError(f"{source}: required ID column {column} absent")
    if frame[column].isna().any() or frame[column].astype(str).str.strip().eq("").any():
        raise DataError(f"{source}: empty IDs")
    if frame[column].duplicated().any():
        raise DataError(f"{source}: duplicate IDs; explicit session/source selection required")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def record_run(cfg: dict, command: str, extra: dict | None = None) -> None:
    info = {
        "utc": datetime.now(UTC).isoformat(), "command": command,
        "mode": cfg["mode"], "config": cfg,
        "python": platform.python_version(),
        "packages": {n: version(n) for n in ["wmh-hcy", "numpy", "pandas", "statsmodels",
                                             "pyreadstat", "miceforest", "nibabel"]},
        **(extra or {}),
    }
    dump_json(outdir(cfg) / "logs" / f"{command}.json", info)


def finite_positive(frame: pd.DataFrame, cols: list[str]) -> pd.Series:
    ok = pd.Series(True, index=frame.index)
    for c in cols:
        ok &= np.isfinite(frame[c]) & (frame[c] > 0)
    return ok
