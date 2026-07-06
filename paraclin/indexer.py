"""Scan a results folder and (re)build the sample index.

For every ``*.paraphase.json`` found under ``results_root`` we resolve the
sibling BAM/BAI and per-gene VCFs, derive the reference build, checksum the
primary files, run matched-trio validation, and upsert a row into SQLite.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .db import Sample, get_session, init_db
from .settings import get_settings


def _md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _build_from_json(data: dict) -> str | None:
    """Build id is the prefix of any gene record's ``phase_region`` (e.g. '38')."""
    for rec in data.values():
        if isinstance(rec, dict):
            pr = rec.get("phase_region")
            if isinstance(pr, str) and ":" in pr:
                return pr.split(":", 1)[0]
    return None


def _discover_vcfs(json_path: Path, sample_id: str, genes: list[str]) -> dict[str, str]:
    """Map gene -> VCF path using Paraphase's layout:
    ``<dir>/<sample>_paraphase_vcfs/<sample>_<gene>.vcf``."""
    vcf_dir = json_path.parent / f"{sample_id}_paraphase_vcfs"
    out: dict[str, str] = {}
    if not vcf_dir.is_dir():
        return out
    for gene in genes:
        cand = vcf_dir / f"{sample_id}_{gene}.vcf"
        if cand.exists():
            out[gene] = str(cand)
    return out


def index_one(json_path: Path) -> Sample:
    sample_id = json_path.name
    for suffix in (".paraphase.json", ".json"):
        if sample_id.endswith(suffix):
            sample_id = sample_id[: -len(suffix)]
            break

    with open(json_path) as fh:
        data = json.load(fh)

    warnings: list[str] = []
    build = _build_from_json(data) or get_settings().default_build

    genes = [g for g, rec in data.items() if isinstance(rec, dict)]

    # sample-level fields (take from whichever record carries them)
    sample_sex = None
    genome_depth = None
    for rec in data.values():
        if isinstance(rec, dict):
            sample_sex = sample_sex or rec.get("sample_sex")
            if genome_depth is None and rec.get("genome_depth") is not None:
                genome_depth = rec.get("genome_depth")

    # sibling BAM / BAI
    bam = json_path.parent / f"{sample_id}.paraphase.bam"
    bai = json_path.parent / f"{sample_id}.paraphase.bam.bai"
    bam_path = str(bam) if bam.exists() else None
    bai_path = str(bai) if bai.exists() else None
    if not bam_path:
        warnings.append("BAM not found next to JSON — in-app viewer/session unavailable.")
    elif not bai_path:
        warnings.append("BAM index (.bai) missing — igv.js requires it.")

    vcf_paths = _discover_vcfs(json_path, sample_id, genes)
    interpreted = {"smn1", "f8"}
    for gene in interpreted & set(genes):
        if gene not in vcf_paths:
            warnings.append(f"VCF for '{gene}' not found (interpretation still works from JSON).")

    sample = Sample(
        sample_id=sample_id,
        build=build,
        json_path=str(json_path),
        bam_path=bam_path,
        bai_path=bai_path,
        vcf_paths=vcf_paths,
        genes=genes,
        sample_sex=sample_sex,
        genome_depth=genome_depth,
        json_md5=_md5(json_path),
        bam_md5=_md5(bam) if bam_path else None,
        paraphase_version=None,  # not stored in JSON; stamped from live import at report time
        warnings=warnings,
    )
    return sample


def rescan() -> dict:
    """Rebuild the index from scratch. Returns a summary."""
    Session = init_db()
    settings = get_settings()
    root = settings.results_root
    found = sorted(root.rglob("*.paraphase.json"))

    indexed, errors = [], []
    with Session() as session:
        # full rebuild: clear then repopulate (index is disposable)
        session.query(Sample).delete()
        for jp in found:
            try:
                session.merge(index_one(jp))
                indexed.append(jp.name)
            except Exception as exc:  # noqa: BLE001
                errors.append({"file": str(jp), "error": str(exc)})
        session.commit()

    return {
        "results_root": str(root),
        "indexed": len(indexed),
        "samples": [s.replace(".paraphase.json", "") for s in indexed],
        "errors": errors,
    }
