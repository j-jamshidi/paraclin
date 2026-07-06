"""Condition registry.

A "condition" is a Paraphase gene locus presented in clinical terms. paraclin only
adds the clinical layer (disease name + which interpreter runs); the genomic facts
it needs — the display locus, the member gene symbols and the reference build —
are read from each sample's own ``*.paraphase.json`` (fields ``phase_region`` and
``genes_in_region``). This keeps paraclin self-contained: it needs Paraphase
*output*, not a Paraphase installation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

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
    region: str               # display locus, e.g. "chr5:70917100-70961220"
    disease: str | None       # clinical display name, if interpreted
    interpreter: str | None   # interpreter id, if any
    summary: str = ""
    has_interpreter: bool = field(default=False)

    @property
    def chrom(self) -> str:
        return self.region.split(":", 1)[0] if self.region else ""


def parse_phase_region(phase_region: str | None) -> tuple[str | None, str]:
    """Split Paraphase's ``phase_region`` ("<build>:chrN:start-end") into
    (build, "chrN:start-end"). Returns (None, "") if it can't be parsed."""
    if isinstance(phase_region, str) and phase_region.count(":") == 2:
        build, region = phase_region.split(":", 1)
        return build, region
    return None, ""


def interpreted_genes() -> list[str]:
    """Genes paraclin has an interpreter for (currently smn1, f8)."""
    return list(CONDITION_META)


def condition_from_record(gene: str, rec: dict, default_build: str = "38") -> Condition:
    """Build a Condition from a sample's gene record in the Paraphase JSON."""
    meta = CONDITION_META.get(gene, {})
    build, region = parse_phase_region(rec.get("phase_region"))
    genes = [g.strip() for g in str(rec.get("genes_in_region", gene)).split(",") if g.strip()]
    return Condition(
        gene=gene,
        build=build or default_build,
        genes=genes or [gene],
        region=region,
        disease=meta.get("disease"),
        interpreter=meta.get("interpreter"),
        summary=meta.get("summary", ""),
        has_interpreter=bool(meta.get("interpreter")),
    )
