"""FastAPI application wiring the Paraphase GUI backend together."""

from __future__ import annotations

import json
import os
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse

from . import __version__ as APP_VERSION
from . import audit, igv, indexer
from .conditions import get_condition, list_conditions
from .db import Sample, get_session
from .interpret import get_interpreter
from .qc import build_qc
from .settings import get_settings

app = FastAPI(title="Paraphase Clinical Review", version=APP_VERSION)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # local single-user; tighten to lab origin in server phase
    allow_methods=["*"],
    allow_headers=["*"],
)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _get_sample(sample_id: str) -> Sample:
    with get_session() as s:
        row = s.get(Sample, sample_id)
    if row is None:
        raise HTTPException(404, f"Sample '{sample_id}' not indexed. Run /api/rescan.")
    return row


def _load_record(row: Sample, gene: str) -> dict:
    with open(row.json_path) as fh:
        data = json.load(fh)
    if gene not in data:
        raise HTTPException(404, f"Gene '{gene}' not present in {row.sample_id}.")
    return data[gene]


def _range_response(path: str, request: Request, media_type: str) -> Response:
    """Byte-range-aware file response (required by igv.js for BAM/BAI/VCF)."""
    file_size = os.path.getsize(path)
    range_header = request.headers.get("range")
    if range_header is None:
        return FileResponse(path, media_type=media_type)

    try:
        units, rng = range_header.split("=")
        start_s, end_s = rng.split("-")
        start = int(start_s) if start_s else 0
        end = int(end_s) if end_s else file_size - 1
    except ValueError:
        raise HTTPException(416, "Invalid Range header")
    end = min(end, file_size - 1)
    length = end - start + 1

    def iterator():
        with open(path, "rb") as fh:
            fh.seek(start)
            remaining = length
            while remaining > 0:
                chunk = fh.read(min(1 << 20, remaining))
                if not chunk:
                    break
                remaining -= len(chunk)
                yield chunk

    headers = {
        "Content-Range": f"bytes {start}-{end}/{file_size}",
        "Accept-Ranges": "bytes",
        "Content-Length": str(length),
    }
    return StreamingResponse(iterator(), status_code=206, headers=headers, media_type=media_type)


# --------------------------------------------------------------------------- #
# meta / index
# --------------------------------------------------------------------------- #
@app.get("/api/health")
def health():
    s = get_settings()
    return {"app_version": APP_VERSION, "results_root": str(s.results_root),
            "igv_genome": s.igv_genome}


@app.post("/api/rescan")
def rescan():
    summary = indexer.rescan()
    audit.record("rescan", indexed=summary["indexed"])
    return summary


@app.get("/api/samples")
def samples():
    with get_session() as s:
        rows = s.query(Sample).order_by(Sample.sample_id).all()
    return {"samples": [r.to_dict() for r in rows]}


@app.get("/api/samples/{sample_id}")
def sample_detail(sample_id: str):
    row = _get_sample(sample_id)
    with open(row.json_path) as fh:
        data = json.load(fh)
    conds = []
    for c in list_conditions(row.build):
        if c.gene in row.genes:
            # A tighter, sample-specific default view than the whole realign_region:
            # Paraphase's per-gene phase_region ("<build>:chrN:start-end").
            pr = (data.get(c.gene) or {}).get("phase_region")
            viewer_locus = c.realign_region
            if isinstance(pr, str) and pr.count(":") == 2:
                viewer_locus = pr.split(":", 1)[1]
            conds.append({"gene": c.gene, "disease": c.disease,
                          "summary": c.summary, "region": c.realign_region,
                          "viewer_locus": viewer_locus})
    d = row.to_dict()
    d["conditions"] = conds
    return d


@app.get("/api/samples/{sample_id}/raw")
def sample_raw(sample_id: str):
    row = _get_sample(sample_id)
    with open(row.json_path) as fh:
        return JSONResponse(json.load(fh))


# --------------------------------------------------------------------------- #
# interpretation
# --------------------------------------------------------------------------- #
@app.get("/api/samples/{sample_id}/{gene}/result")
def sample_result(sample_id: str, gene: str):
    row = _get_sample(sample_id)
    cond = get_condition(gene, row.build)
    if cond is None or not cond.has_interpreter:
        raise HTTPException(404, f"No interpreter for gene '{gene}'.")
    rec = _load_record(row, gene)
    q = build_qc(rec)
    interp = get_interpreter(cond.interpreter)
    result = interp.interpret(rec, row.vcf_paths.get(gene), q)

    audit.record("view_result", sample=sample_id, gene=gene, status=result.status)
    out = result.to_dict()
    out["provenance"] = audit.provenance_stamp(row)
    out["condition"] = {"gene": cond.gene, "disease": cond.disease,
                        "region": cond.realign_region, "build": cond.build}
    return out


# --------------------------------------------------------------------------- #
# IGV assets
# --------------------------------------------------------------------------- #
@app.get("/api/samples/{sample_id}/{gene}/session.xml")
def session_xml(sample_id: str, gene: str):
    row = _get_sample(sample_id)
    cond = get_condition(gene, row.build)
    if cond is None:
        raise HTTPException(404, f"Unknown gene '{gene}'.")
    if not row.bam_path:
        raise HTTPException(409, "No BAM available for this sample.")
    vcfs = [row.vcf_paths[gene]] if gene in row.vcf_paths else []
    xml = igv.build_session_xml(sample_id, row.bam_path, vcfs, cond,
                                get_settings().igv_genome, relative=False)
    audit.record("download_session", sample=sample_id, gene=gene)
    fname = f"{sample_id}_{gene}_igv_session.xml"
    return Response(xml, media_type="application/xml",
                    headers={"Content-Disposition": f'attachment; filename="{fname}"'})


@app.get("/api/samples/{sample_id}/{gene}/bundle.zip")
def bundle_zip(sample_id: str, gene: str):
    row = _get_sample(sample_id)
    cond = get_condition(gene, row.build)
    if cond is None:
        raise HTTPException(404, f"Unknown gene '{gene}'.")
    if not row.bam_path:
        raise HTTPException(409, "No BAM available for this sample.")
    vcfs = [row.vcf_paths[gene]] if gene in row.vcf_paths else []
    data = igv.build_bundle_zip(sample_id, row.bam_path, row.bai_path, vcfs, cond,
                                get_settings().igv_genome)
    audit.record("download_bundle", sample=sample_id, gene=gene)
    fname = f"{sample_id}_{gene}_igv_bundle.zip"
    return Response(data, media_type="application/zip",
                    headers={"Content-Disposition": f'attachment; filename="{fname}"'})


# --------------------------------------------------------------------------- #
# raw file serving for igv.js (range-aware)
# --------------------------------------------------------------------------- #
@app.get("/api/files/{sample_id}/bam")
def file_bam(sample_id: str, request: Request):
    row = _get_sample(sample_id)
    if not row.bam_path:
        raise HTTPException(404, "No BAM.")
    return _range_response(row.bam_path, request, "application/octet-stream")


@app.get("/api/files/{sample_id}/bai")
def file_bai(sample_id: str, request: Request):
    row = _get_sample(sample_id)
    if not row.bai_path:
        raise HTTPException(404, "No BAI.")
    return _range_response(row.bai_path, request, "application/octet-stream")


@app.get("/api/files/{sample_id}/vcf/{gene}")
def file_vcf(sample_id: str, gene: str, request: Request):
    row = _get_sample(sample_id)
    p = row.vcf_paths.get(gene)
    if not p or not Path(p).exists():
        raise HTTPException(404, "No VCF.")
    return _range_response(p, request, "text/plain")
