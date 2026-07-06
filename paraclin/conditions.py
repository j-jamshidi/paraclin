"""Condition registry.

A "condition" is a Paraphase gene locus presented in clinical terms. The genomic
facts (member genes, display locus/`realign_region`, reference build) come
straight from the vendored Paraphase config so there is a single source of truth;
this module only adds the clinical layer (disease name + which interpreter runs).

Paraphase encodes the build in each sample's JSON ``phase_region`` prefix
(``38`` / ``19`` / ``chm13``), so the same gene can be resolved per-build.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import yaml

from .settings import get_settings

# Clinical metadata for the conditions we interpret. Genes absent here are still
# viewable (raw calls + IGV) but carry no derived clinical statement.
CONDITION_META: dict[str, dict] = {
    "smn1": {
        "disease": "Spinal Muscular Atrophy (SMA)",
        "interpreter": "smn1",
        "summary": "SMN1/SMN2 copy number, carrier / affected status, and 2+0 "
        "silent-carrier assessment.",
    },
    "f8": {
        "disease": "Hemophilia A (F8)",
        "interpreter": "f8",
        "summary": "F8 intron-22 inversion and exon1-22 deletion detection.",
    },
}


@dataclass(frozen=True)
class Condition:
    gene: str                 # Paraphase gene key, e.g. "smn1"
    build: str                # "38" / "19" / "chm13"
    genes: list[str]          # member gene symbols
    realign_region: str       # display locus, e.g. "chr5:70890000-71100000"
    disease: str | None       # clinical display name, if interpreted
    interpreter: str | None   # interpreter id, if any
    summary: str = ""
    has_interpreter: bool = field(default=False)

    @property
    def chrom(self) -> str:
        return self.realign_region.split(":", 1)[0]


@lru_cache(maxsize=8)
def _load_build_config(build: str) -> dict:
    settings = get_settings()
    cfg_path: Path = settings.paraphase_repo / "paraphase" / "data" / build / "config.yaml"
    if not cfg_path.exists():
        raise FileNotFoundError(f"Paraphase config not found for build {build}: {cfg_path}")
    with open(cfg_path) as fh:
        return yaml.safe_load(fh) or {}


def get_condition(gene: str, build: str) -> Condition | None:
    """Return the Condition for a gene in a build, or None if the gene is unknown
    to Paraphase's config for that build."""
    cfg = _load_build_config(build)
    entry = cfg.get(gene)
    if entry is None:
        return None
    meta = CONDITION_META.get(gene, {})
    genes = [g.strip() for g in str(entry.get("genes", gene)).split(",") if g.strip()]
    return Condition(
        gene=gene,
        build=build,
        genes=genes,
        realign_region=entry.get("realign_region", ""),
        disease=meta.get("disease"),
        interpreter=meta.get("interpreter"),
        summary=meta.get("summary", ""),
        has_interpreter=bool(meta.get("interpreter")),
    )


def list_conditions(build: str) -> list[Condition]:
    """All interpreted conditions available for a build (currently smn1, f8)."""
    out = []
    for gene in CONDITION_META:
        cond = get_condition(gene, build)
        if cond:
            out.append(cond)
    return out
