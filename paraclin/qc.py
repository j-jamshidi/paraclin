"""QC gating.

Extracts depth context from a Paraphase gene record so no clinical call is ever
presented without it. This is a READ-DEPTH gate only — it does not assess mapping
quality, phasing confidence, or call correctness.

Thresholds are conservative defaults for ONT and are configurable in config.yaml
(`min_region_depth`, `min_genome_depth`); tune per validated assay. A sample fails
(-> REVIEW) if a depth value is below its threshold OR is not reported.
"""

from __future__ import annotations

from .settings import get_settings


def build_qc(rec: dict) -> dict:
    settings = get_settings()
    min_region = settings.min_region_depth
    min_genome = settings.min_genome_depth

    region = rec.get("region_depth") or {}
    region_median = region.get("median") if isinstance(region, dict) else None
    genome_depth = rec.get("genome_depth")

    flags: list[str] = []
    if region_median is None:
        flags.append("Region depth not reported.")
    elif region_median < min_region:
        flags.append(f"Low region depth ({region_median}x < {min_region:g}x).")
    if genome_depth is None:
        flags.append("Genome depth not reported.")
    elif genome_depth < min_genome:
        flags.append(f"Low genome depth ({genome_depth}x < {min_genome:g}x).")

    return {
        "region_depth_median": region_median,
        "region_depth_p80": region.get("percentile80") if isinstance(region, dict) else None,
        "genome_depth": genome_depth,
        "min_region_depth": min_region,
        "min_genome_depth": min_genome,
        "flags": flags,
        "pass": len(flags) == 0,
    }
