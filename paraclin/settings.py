"""Configuration loading and path resolution.

Resolves every path in config.yaml against the repository root so the app can be
launched from anywhere. Loaded once at import time via ``get_settings()``.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import yaml

# Repository root = one level up from this package (paraclin/paraclin/settings.py).
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "config.yaml"


@dataclass(frozen=True)
class Settings:
    results_root: Path
    default_build: str
    database: Path
    audit_log: Path
    igv_genome: str
    reference_fasta: Path | None

    @property
    def db_url(self) -> str:
        return f"sqlite:///{self.database}"


def _resolve(base: Path, value: str) -> Path:
    p = Path(value).expanduser()
    return p if p.is_absolute() else (base / p)


@lru_cache(maxsize=1)
def get_settings(config_path: str | os.PathLike | None = None) -> Settings:
    path = Path(config_path) if config_path else DEFAULT_CONFIG
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}

    ref = raw.get("reference_fasta")
    return Settings(
        results_root=_resolve(REPO_ROOT, raw.get("results_root", "sample_data")),
        default_build=str(raw.get("default_build", "38")),
        database=_resolve(REPO_ROOT, raw.get("database", "paraclin.db")),
        audit_log=_resolve(REPO_ROOT, raw.get("audit_log", "audit.log")),
        igv_genome=raw.get("igv_genome", "hg38"),
        reference_fasta=_resolve(REPO_ROOT, ref) if ref else None,
    )
