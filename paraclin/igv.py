"""IGV desktop session (.xml) generation and portable bundles.

Reproduces Paraphase's recommended view: load the ``.paraphase.bam``, group reads
by the ``HP`` tag and color alignments by the ``YC`` tag, at the condition's
display locus. Also builds a zip bundle (BAM+BAI+VCF+session) whose session uses
relative paths so it opens anywhere.
"""

from __future__ import annotations

import io
import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.dom import minidom

from .conditions import Condition


def build_session_xml(
    sample_id: str,
    bam_path: str,
    vcf_paths: list[str],
    condition: Condition,
    igv_genome: str = "hg38",
    relative: bool = False,
) -> str:
    """Return an IGV session XML string.

    ``relative=True`` writes bare filenames (for a portable zip); otherwise
    absolute paths (for opening in place on this workstation).
    """
    def ref(p: str) -> str:
        return os.path.basename(p) if relative else str(Path(p).resolve())

    session = ET.Element(
        "Session",
        {
            "genome": igv_genome,
            "hasGeneTrack": "true",
            "hasSequenceTrack": "true",
            "locus": condition.realign_region,
            "version": "8",
        },
    )

    resources = ET.SubElement(session, "Resources")
    ET.SubElement(resources, "Resource", {"path": ref(bam_path), "type": "bam"})
    for v in vcf_paths:
        ET.SubElement(resources, "Resource", {"path": ref(v), "type": "vcf"})

    panel = ET.SubElement(session, "Panel", {"name": "DataPanel"})
    bam_id = ref(bam_path)
    ET.SubElement(
        panel, "Track",
        {"clazz": "org.broad.igv.sam.CoverageTrack",
         "id": f"{bam_id}_coverage", "name": f"{sample_id} coverage", "autoScale": "true"},
    )
    align = ET.SubElement(
        panel, "Track",
        {"clazz": "org.broad.igv.sam.AlignmentTrack",
         "id": bam_id, "name": f"{sample_id} — {condition.gene} (HP/YC phased)"},
    )
    # The two options that reproduce the Paraphase recommendation.
    ET.SubElement(
        align, "RenderOptions",
        {"colorOption": "TAG", "colorByTag": "YC",
         "groupByOption": "TAG", "groupByTag": "HP"},
    )
    for v in vcf_paths:
        vpanel = ET.SubElement(session, "Panel", {"name": f"VariantPanel_{os.path.basename(v)}"})
        ET.SubElement(
            vpanel, "Track",
            {"clazz": "org.broad.igv.variant.VariantTrack", "id": ref(v),
             "name": os.path.basename(v)},
        )

    rough = ET.tostring(session, encoding="unicode")
    return minidom.parseString(rough).toprettyxml(indent="  ")


def build_bundle_zip(
    sample_id: str,
    bam_path: str,
    bai_path: str | None,
    vcf_paths: list[str],
    condition: Condition,
    igv_genome: str = "hg38",
) -> bytes:
    """Zip BAM(+BAI)+VCF(s)+relative-path session.xml for portable sharing."""
    session_xml = build_session_xml(
        sample_id, bam_path, vcf_paths, condition, igv_genome, relative=True
    )
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(bam_path, os.path.basename(bam_path))
        if bai_path and os.path.exists(bai_path):
            zf.write(bai_path, os.path.basename(bai_path))
        for v in vcf_paths:
            zf.write(v, os.path.basename(v))
        zf.writestr(f"{sample_id}_{condition.gene}_igv_session.xml", session_xml)
    return buf.getvalue()
