# paraclin

**A clinical review & visualization GUI for [PacBio Paraphase](https://github.com/PacificBiosciences/paraphase) long-read outputs.**

Paraphase resolves segmentally-duplicated medical genes (SMN1, F8, …) from PacBio
HiFi data, but its output — copy numbers, haplotype strings, tagged BAMs — is not
convenient to review. paraclin turns it into a clean workflow: pick a sample, pick
a condition, get a clinically-framed result with QC and provenance, view the
phased/colored alignments in-app (igv.js), and download an IGV desktop session
that reproduces Paraphase's recommended view.

Conditions supported today: **SMN1 / Spinal Muscular Atrophy** and
**F8 / Hemophilia A**. Other Paraphase genes are viewable (raw calls + IGV) and
new interpreters plug into a registry.

> ⚠️ **Disclaimers**
> - **Independent project.** Not affiliated with or endorsed by Pacific
>   Biosciences. Paraphase is proprietary PacBio software and is **not** bundled —
>   install it separately under its own license.
> - **Research use.** paraclin is a research/prototype tool. Nothing here is a
>   validated diagnostic. The SMA carrier / 2+0 silent-carrier assessment is
>   explicitly marked **experimental** in the UI.
> - **Patient data stays local.** paraclin is designed to run locally; it makes no
>   external calls with sample data. Do not commit Paraphase outputs / BAMs to
>   version control (see `.gitignore`).

## Features

- **Sample browser** — scans a results folder, indexes samples into SQLite with
  MD5 checksums and matched-trio (JSON+BAM+VCF) validation.
- **Condition selector** — filter to a single condition (SMA / Hemophilia A); the
  sample list and the selected sample scope to it.
- **Clinical result view** — a color-coded headline (affected / carrier / SV
  detected …), a QC panel (depth gating), the raw Paraphase calls kept separate
  from the derived interpretation, references, and a full provenance stamp
  (app / Paraphase / interpreter versions, build, checksums).
- **Experimental carrier tab (SMA)** — affected vs not-affected is the primary
  call; carrier / silent-carrier (2+0) status is a separate tab clearly labelled
  experimental.
- **Embedded igv.js** — reproduces Paraphase's recommended view (reads grouped by
  the `HP` tag, colored by the `YC` tag), squished, at the per-sample locus.
- **Downloads** — an IGV desktop session `.xml`, a portable `bundle.zip`
  (BAM+BAI+VCF+relative-path session), and the raw JSON.
- **Audit + auth seam** — every view/download is written to an append-only audit
  log; the auth/user hook is behind an interface so a local single-user build can
  be promoted to a multi-user server (RBAC/OIDC) without touching call sites.

## Architecture

```
paraclin/            FastAPI backend (Python)
  settings.py        config loader
  conditions.py      gene -> region/build read from each sample's JSON
  db.py indexer.py   SQLite index built by scanning results_root (checksums)
  qc.py              depth gating
  interpret/         base + smn1 + f8 (raw calls kept separate from calls)
  igv.py             IGV session.xml + portable bundle.zip
  audit.py           auth/audit seam
  app.py             API + range-aware BAM/VCF serving for igv.js
frontend/            React + TypeScript + igv.js (Vite)
scripts/             standalone CLIs (e.g. SMN1 silent-carrier caller)
```

## Requirements

- Python 3.10+
- Node.js 18+ / npm
- [Paraphase](https://github.com/PacificBiosciences/paraphase) **output** to review
  (the `*.paraphase.json` + BAM + VCFs produced by running Paraphase on your HiFi
  data). paraclin reads these output files only — it does **not** need a Paraphase
  installation at runtime.

## Install & run

paraclin runs as two local processes: a **backend** (FastAPI, port 8077) and a
**frontend** dev server (Vite, port 5199) that proxies `/api` to the backend. Use
two terminals.

### 0. Clone

```bash
git clone https://github.com/j-jamshidi/paraclin.git
cd paraclin
```

### 1. Backend (terminal 1)

```bash
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
uvicorn paraclin.app:app --port 8077
```

Leave it running. During development add `--reload` to auto-restart on code
changes. A quick check: `curl http://localhost:8077/api/health` should return JSON.

### 2. Frontend (terminal 2)

```bash
cd frontend
npm install
npm run dev
```

Open **http://localhost:5199**.

### 3. Point paraclin at your data

Put your Paraphase outputs under `results_root` (default `sample_data/`) — see
[Configuration](#configuration). Everything paraclin needs (locus, member genes,
build) is read from the outputs themselves; no Paraphase installation is required.

Then click **Rescan folder** in the UI (or `curl -X POST
http://localhost:8077/api/rescan`) to index them. Pick a sample and a condition.

### 4. Access from other computers on your network (LAN)

To let colleagues on the **same local network** open paraclin from their own
browsers, start the frontend in LAN mode:

```bash
cd frontend
npm run dev:lan          # binds 0.0.0.0 and prints a Network URL
```

Vite prints something like `Network: http://192.168.1.42:5199/`. Anyone whose
machine can reach yours on that subnet can open that URL. The backend can stay on
`localhost` — the dev server proxies `/api` (and the BAM/VCF byte-range requests)
to it on your machine, so **no client-side setup is needed**.

Find your machine's IP if Vite doesn't show it:

```bash
# macOS
ipconfig getifaddr en0
# Linux
hostname -I
# Windows
ipconfig        # look for the IPv4 Address
```

Optional env overrides (see `frontend/vite.config.ts`): `PARACLIN_PORT`,
`PARACLIN_HOST`, and `PARACLIN_BACKEND` (if the backend runs on a different host).

> ⚠️ **Security / PHI:** LAN mode exposes paraclin — and any Paraphase/patient
> data it can load — to everyone who can reach your machine on that port. Use it
> only on a trusted network, behind your institution's firewall, and **never**
> forward the port to the public internet. It has no authentication yet (planned
> for the server phase).

### Stopping

Press `Ctrl-C` in each terminal (or `pkill -f uvicorn` for the backend).

### Troubleshooting

| Symptom | Fix |
|---|---|
| Sample list is empty | put outputs under `results_root`, then click **Rescan folder** |
| igv viewer is blank | ensure the `.bam` **and** `.bam.bai` sit next to the `.paraphase.json` |
| LAN URL not reachable | check your firewall allows inbound TCP 5199; confirm both machines are on the same subnet |

## Configuration

Edit [`config.yaml`](config.yaml):

| Key | Meaning |
|---|---|
| `results_root` | folder scanned for `*.paraphase.json` (+ sibling BAM/VCFs) |
| `default_build` | `38` / `19` / `chm13` fallback when a sample omits it |
| `database`, `audit_log` | SQLite index and audit log locations |
| `igv_genome` | igv.js genome id (match your data's build) |

Put your Paraphase outputs under `results_root` (default `sample_data/`), keeping
Paraphase's layout: `<sample>.paraphase.json`, `<sample>.paraphase.bam(.bai)`, and
`<sample>_paraphase_vcfs/<sample>_<gene>.vcf`.

## API (selected)

| Endpoint | Purpose |
|---|---|
| `POST /api/rescan` | rebuild the sample index |
| `GET /api/samples` | indexed samples |
| `GET /api/samples/{id}` | sample + available conditions (with viewer locus) |
| `GET /api/samples/{id}/{gene}/result` | interpretation + QC + provenance |
| `GET /api/samples/{id}/{gene}/session.xml` | IGV desktop session (HP group / YC color) |
| `GET /api/samples/{id}/{gene}/bundle.zip` | portable BAM+BAI+VCF+session |
| `GET /api/files/{id}/bam` `.../bai` `.../vcf/{gene}` | range-served for igv.js |

## Standalone CLI

`scripts/smn1_silent_carrier.py` runs the SMN1 affected/carrier + 2+0
silent-carrier assessment directly on a Paraphase JSON (or a directory of them):

```bash
python scripts/smn1_silent_carrier.py /path/to/<sample>.paraphase.json
```

## Design notes (accreditation-oriented)

Genomic facts (locus, member genes, build) are read from each sample's own
Paraphase JSON, so paraclin stays a self-contained reader of Paraphase output; every
result carries QC and a provenance stamp; raw caller output is always shown next
to the derived clinical statement; interpretation is deterministic and offline;
every view/download is audited. Interpretation thresholds (e.g. QC depth) are
conservative defaults and must be validated per assay before any clinical use.

## Author & license

Created by **Javad Jamshidi**. Released under the [MIT License](LICENSE). paraclin
reads Paraphase output but contains none of Paraphase's code; Paraphase is
separately licensed by Pacific Biosciences.
