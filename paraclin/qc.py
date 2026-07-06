"""QC gating.

Extracts depth/quality context from a Paraphase gene record so no clinical call
is ever presented without it. Thresholds are conservative defaults for HiFi;
tune per validated assay.
"""

from __future__ import annotations

# Minimum region median depth below which calls are flagged for review.
MIN_REGION_DEPTH = 20.0
MIN_GENOME_DEPTH = 15.0


def build_qc(rec: dict) -> dict:
    region = rec.get("region_depth") or {}
    region_median = region.get("median") if isinstance(region, dict) else None
    genome_depth = rec.get("genome_depth")

    flags: list[str] = []
    if region_median is not None and region_median < MIN_REGION_DEPTH:
        flags.append(f"Low region depth ({region_median}x < {MIN_REGION_DEPTH}x).")
    if genome_depth is not None and genome_depth < MIN_GENOME_DEPTH:
        flags.append(f"Low genome depth ({genome_depth}x < {MIN_GENOME_DEPTH}x).")

    return {
        "region_depth_median": region_median,
        "region_depth_p80": region.get("percentile80") if isinstance(region, dict) else None,
        "genome_depth": genome_depth,
        "min_region_depth": MIN_REGION_DEPTH,
        "flags": flags,
        "pass": len(flags) == 0,
    }
