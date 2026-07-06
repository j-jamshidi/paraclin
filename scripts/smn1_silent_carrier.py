#!/usr/bin/env python3
"""
smn1_silent_carrier.py

For samples where Paraphase reports SMN1 copy number = 2, flag whether the two
SMN1 copies are more likely in a 1+1 configuration (one copy on each chromosome
5 -> NOT a carrier) or a 2+0 "silent carrier" configuration (both SMN1 copies in
cis on one chromosome, zero on the other -> a hidden SMA carrier).

IMPORTANT — what this can and cannot do
---------------------------------------
Copy number alone CANNOT distinguish 1+1 from 2+0: both give a total of 2 copies.
Short/long reads of a single individual cannot phase the two chromosome-5
homologues across the whole locus, so 2+0 cannot be *proven* from one sample.

What we CAN do is look for evidence of an SMN1 *duplication* allele (two SMN1
copies in cis). This tool checks TWO independent lines of evidence:

(A) Duplication-linked marker variants (Luo et al. 2014, Genet Med; also used by
    the DRAGEN SMN caller):
      Marker 1 : g.27134T>G         rs143838139  chr5:70952074  T>G           (hg38)
      Marker 2 : g.27706_27707delAT rs200800214  chr5:70952646_70952647 delAT (hg38)

(B) Paraphase haplogroup assignment (Chen et al. 2023, AJHG — the Paraphase
    paper). Each SMN1 haplotype is assigned a haplogroup; some haplogroups are
    known two-copy (duplication) alleles. This signal is BROADER than marker (A):
    the paper reports g.27134T>G is a poor predictor in African populations,
    whereas the haplogroup captures the African-specific two-copy alleles S1-8
    and S1-9d directly (see DUP_HAPLOGROUPS_* below).

If a CN=2 individual shows either signal, a 2+0 silent-carrier configuration is
more likely. Per Chen et al. 2023, however, "without pedigree information, it is
currently not possible to identify silent carriers (2+0)" from a single sample —
so these are PROBABILISTIC signals, not definitive calls.

Caveats:
  * Signals are POPULATION-SPECIFIC and tag only a subset of duplication alleles.
    A negative result does NOT rule out 2+0; it removes the known positive
    evidence and leaves a residual (reduced, non-zero) silent-carrier risk.

Output calls:
  * marker SNP or two-copy haplogroup present -> "2+0 silent carrier – LIKELY"
  * marker-carrying haplogroup only           -> "possible 2+0 – WEAK evidence"
  * nothing found                             -> "consistent with 1+1
                                                  (2+0 not excluded)"

Usage:
  python3 smn1_silent_carrier.py <paraphase.json | directory> [more ...]

If given a directory it scans for *.paraphase.json files (recursively).
"""

import sys
import os
import json
import glob

# --- Marker definitions (hg38 / GRCh38) --------------------------------------
# Each marker: label, chrom, position (1-based), ref, alt, and a position window
# used to match Paraphase's per-haplotype variant strings (which are left-aligned
# and may differ slightly in representation for the indel).
MARKERS = [
    {
        "name": "g.27134T>G",
        "rsid": "rs143838139",
        "chrom": "chr5",
        "pos": 70952074,
        "ref": "T",
        "alt": "G",
        "kind": "snv",
        "window": (70952074, 70952074),
    },
    {
        "name": "g.27706_27707delAT",
        "rsid": "rs200800214",
        "chrom": "chr5",
        "pos": 70952646,
        "ref": "TAT",
        "alt": "T",
        "kind": "del",
        # allow a couple bp of wiggle for left-alignment / representation
        "window": (70952644, 70952648),
    },
]

EXON7_C840 = 70951946  # SMN1 (C) vs SMN2 (T) determinant, for reference only

# --- Duplication-associated SMN1 haplogroups (Chen et al. 2023, AJHG) ---------
# Paraphase assigns each SMN1 haplotype a haplogroup (haplotype_details[hap]
# ['haplogroup']). The Paraphase paper shows some haplogroups are two-copy /
# duplication alleles (i.e. carried in cis on a chromosome that has 2x SMN1),
# which is the mechanism behind 2+0 silent carriers. This haplogroup signal is
# broader than the g.27134T>G SNP, which the paper reports is a poor predictor
# in African populations.
#
#   S1-8, S1-9d : African-specific TWO-COPY (duplication) SMN1 alleles. The paper
#                 states they rarely occur as singletons, so seeing them in a CN=2
#                 individual makes 2+0 substantially more likely (up to ~88.5% in
#                 the African-ancestry example given).
#   S1-9, S1-9c : additional haplogroups that carry the g.27134T>G marker.
DUP_HAPLOGROUPS_STRONG = {"S1-8", "S1-9d"}          # confirmed two-copy alleles
DUP_HAPLOGROUPS_MARKER = {"S1-8", "S1-9", "S1-9c", "S1-9d"}  # carry g.27134T>G


def _variant_pos(vstr):
    """Paraphase variant strings look like '70952074_T_G' or '70929880_TA_T'.
    Return the integer position (first field)."""
    try:
        return int(vstr.split("_", 1)[0])
    except (ValueError, IndexError):
        return None


def marker_in_haplotype(marker, variant_list):
    """Return the matching variant string if this marker is present in a
    haplotype's Paraphase variant list, else None."""
    lo, hi = marker["window"]
    for v in variant_list:
        p = _variant_pos(v)
        if p is None:
            continue
        if marker["kind"] == "snv":
            # require exact position and matching alt base
            if p == marker["pos"]:
                parts = v.split("_")
                if len(parts) >= 3 and parts[2].upper() == marker["alt"].upper():
                    return v
        else:  # deletion: match by position window + net length loss
            if lo <= p <= hi:
                parts = v.split("_")
                if len(parts) >= 3 and len(parts[1]) > len(parts[2]):
                    return v
    return None


def find_vcf(json_path, sample, gene):
    """Locate the per-gene Paraphase VCF that sits next to the JSON.
    Convention: <dir>/<sample>_paraphase_vcfs/<sample>_<gene>.vcf"""
    base = os.path.dirname(os.path.abspath(json_path))
    candidates = [
        os.path.join(base, f"{sample}_paraphase_vcfs", f"{sample}_{gene}.vcf"),
        os.path.join(base, f"{sample}_paraphase_vcfs", f"{sample}_{gene}.vcf.gz"),
    ]
    # fall back to any *_<gene>.vcf under a *_paraphase_vcfs dir
    candidates += glob.glob(os.path.join(base, "*_paraphase_vcfs",
                                         f"*_{gene}.vcf"))
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None


def scan_vcf_for_markers(vcf_path):
    """Independent cross-check straight from the per-haplotype VCF genotypes.

    Returns {marker_name: {"present": bool, "hits": [(col, "REF>ALT@POS"), ...]}}.
    A marker counts as present only if an SMN1-haplotype column (name contains
    'smn1hap') has GT==1 (carries ALT) and the row's REF/ALT match the marker.
    """
    import gzip
    opener = gzip.open if vcf_path.endswith(".gz") else open
    out = {m["name"]: {"present": False, "hits": []} for m in MARKERS}
    smn1_cols = []
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                # sample columns start at index 9
                smn1_cols = [(i, name) for i, name in enumerate(fields)
                             if i >= 9 and "smn1hap" in name]
                continue
            if not smn1_cols:
                continue
            chrom, pos, _id, ref, alt = fields[0], int(fields[1]), fields[2], \
                fields[3], fields[4]
            fmt = fields[8].split(":")
            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                continue
            for marker in MARKERS:
                if chrom != marker["chrom"]:
                    continue
                lo, hi = marker["window"]
                if not (lo <= pos <= hi):
                    continue
                if marker["kind"] == "snv":
                    if not (pos == marker["pos"]
                            and ref.upper() == marker["ref"].upper()
                            and alt.upper() == marker["alt"].upper()):
                        continue
                else:  # deletion: net length loss, ALT shorter than REF
                    if len(ref) <= len(alt):
                        continue
                for i, name in smn1_cols:
                    gt = fields[i].split(":")[gt_idx]
                    if gt == "1":
                        out[marker["name"]]["present"] = True
                        out[marker["name"]]["hits"].append(
                            (name, f"{ref}>{alt}@{pos}"))
    return out


def analyze_gene(gene, rec):
    """Analyze one gene record (e.g. rec = data['smn1'])."""
    smn1_cn = rec.get("smn1_cn")
    smn2_cn = rec.get("smn2_cn")
    hap_details = rec.get("haplotype_details", {}) or {}
    smn1_haps = rec.get("smn1_haplotypes", {}) or {}
    two_copy = rec.get("two_copy_haplotypes", []) or []

    # SMN1 haplotype names (paraphase names contain 'smn1hap')
    smn1_hap_names = [name for name in smn1_haps.values()]
    if not smn1_hap_names:
        smn1_hap_names = [h for h in hap_details if "smn1hap" in h]

    result = {
        "gene": gene,
        "smn1_cn": smn1_cn,
        "smn2_cn": smn2_cn,
        "n_smn1_haplotypes": len(smn1_hap_names),
        "smn1_haplotypes_have_two_copy": bool(two_copy),
        "markers": [],
        "any_marker_positive": False,
        "applicable": (smn1_cn == 2),
    }

    vcf_hits = rec.get("_vcf_marker_hits")  # optional, injected by analyze_file

    for marker in MARKERS:
        hits = []
        for hap in smn1_hap_names:
            vlist = hap_details.get(hap, {}).get("variants", []) or []
            m = marker_in_haplotype(marker, vlist)
            if m:
                hits.append((hap, m))
        present = len(hits) > 0
        entry = {
            "name": marker["name"],
            "rsid": marker["rsid"],
            "locus": f'{marker["chrom"]}:{marker["pos"]}',
            "present": present,          # from Paraphase JSON haplotype variants
            "hits": hits,
            "vcf_present": None,          # from independent VCF cross-check
            "vcf_hits": [],
        }
        if vcf_hits is not None:
            v = vcf_hits.get(marker["name"], {})
            entry["vcf_present"] = v.get("present", False)
            entry["vcf_hits"] = v.get("hits", [])
        # positive if EITHER source sees it
        if present or entry["vcf_present"]:
            result["any_marker_positive"] = True
        result["markers"].append(entry)

    result["vcf_checked"] = vcf_hits is not None

    # --- Haplogroup signal (Chen et al. 2023) --------------------------------
    haplogroups = []
    for hap in smn1_hap_names:
        hg = hap_details.get(hap, {}).get("haplogroup")
        haplogroups.append((hap, hg))
    result["smn1_haplogroups"] = haplogroups
    strong = [(h, hg) for h, hg in haplogroups if hg in DUP_HAPLOGROUPS_STRONG]
    marker_hg = [(h, hg) for h, hg in haplogroups
                 if hg in DUP_HAPLOGROUPS_MARKER and hg not in DUP_HAPLOGROUPS_STRONG]
    result["haplogroup_strong"] = strong
    result["haplogroup_marker"] = marker_hg

    # Interpretation ----------------------------------------------------------
    # Two independent lines of evidence for a duplication (2+0) allele:
    #   (a) the g.27134T>G / g.27706_27707delAT marker variants, and
    #   (b) Paraphase haplogroup assignment to a known two-copy allele.
    snp_pos = result["any_marker_positive"]
    hg_strong = bool(result["haplogroup_strong"])
    hg_marker = bool(result["haplogroup_marker"])

    if smn1_cn != 2:
        result["call"] = "N/A"
        result["interpretation"] = (
            f"SMN1 copy number is {smn1_cn}, not 2 — the 1+1 vs 2+0 question "
            "applies only to CN=2 individuals."
        )
    elif snp_pos or hg_strong:
        result["call"] = "2+0 SILENT CARRIER — LIKELY"
        reasons = []
        if snp_pos:
            reasons.append("a duplication-linked marker SNP is present on an "
                           "SMN1 haplotype")
        if hg_strong:
            hgs = ", ".join(sorted({hg for _, hg in result["haplogroup_strong"]}))
            reasons.append(f"an SMN1 haplotype belongs to a known two-copy "
                           f"(duplication) haplogroup ({hgs})")
        result["interpretation"] = (
            "; ".join(r[0].upper() + r[1:] if i == 0 else r
                      for i, r in enumerate(reasons)) + ". "
            "This suggests two SMN1 copies in cis (2+0), i.e. a silent SMA "
            "carrier despite a copy number of 2. Per Chen et al. 2023 this "
            "cannot be confirmed from a single sample without pedigree data — "
            "confirm with family studies / orthogonal testing and counselling."
        )
    elif hg_marker:
        result["call"] = "possible 2+0 — WEAK haplogroup evidence"
        hgs = ", ".join(sorted({hg for _, hg in result["haplogroup_marker"]}))
        result["interpretation"] = (
            f"An SMN1 haplotype is in haplogroup {hgs}, which carries the "
            "g.27134T>G marker but is not itself a confirmed two-copy allele. "
            "This raises, but does not establish, the chance of a 2+0 "
            "configuration; interpret with ancestry and pedigree in mind."
        )
    else:
        result["call"] = "no silent-carrier evidence — consistent with 1+1"
        result["interpretation"] = (
            "No duplication-linked marker SNP (rs143838139, rs200800214) and no "
            "duplication-associated SMN1 haplogroup (S1-8/S1-9/S1-9c/S1-9d) were "
            "found. This is consistent with a 1+1 configuration (not a carrier). "
            "NOTE: per Chen et al. 2023, 2+0 cannot be excluded from a single "
            "sample — these signals are population-specific and tag only a subset "
            "of duplication alleles, so a residual silent-carrier risk remains, "
            "especially outside African / Ashkenazi Jewish ancestry."
        )
    return result


def analyze_file(path):
    with open(path) as fh:
        data = json.load(fh)
    sample = os.path.basename(path)
    for suffix in (".paraphase.json", ".json"):
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
            break
    out = []
    for gene, rec in data.items():
        if not isinstance(rec, dict):
            continue
        if "smn1_cn" not in rec:  # only SMN1-type records
            continue
        vcf_path = find_vcf(path, sample, gene)
        if vcf_path:
            try:
                rec["_vcf_marker_hits"] = scan_vcf_for_markers(vcf_path)
            except Exception as e:  # noqa: BLE001
                print(f"WARN: VCF cross-check failed for {vcf_path}: {e}",
                      file=sys.stderr)
        res = analyze_gene(gene, rec)
        res["sample"] = sample
        res["file"] = path
        res["vcf_path"] = vcf_path
        out.append(res)
    return out


def collect_files(args):
    files = []
    for a in args:
        if os.path.isdir(a):
            files.extend(sorted(glob.glob(os.path.join(a, "**", "*.paraphase.json"),
                                          recursive=True)))
        elif os.path.isfile(a):
            files.append(a)
        else:
            print(f"WARN: not found: {a}", file=sys.stderr)
    return files


def print_report(results):
    for r in results:
        print("=" * 72)
        print(f"Sample : {r['sample']}   (gene record: {r['gene']})")
        print(f"SMN1 CN: {r['smn1_cn']}    SMN2 CN: {r['smn2_cn']}    "
              f"distinct SMN1 haplotypes: {r['n_smn1_haplotypes']}")
        if not r["applicable"]:
            print(f"  -> {r['interpretation']}")
            continue
        vcf_note = "with VCF cross-check" if r.get("vcf_checked") else \
            "JSON only (no VCF found)"
        print(f"  Silent-carrier markers on SMN1 haplotypes ({vcf_note}):")
        for m in r["markers"]:
            json_s = "PRESENT" if m["present"] else "absent"
            if m["vcf_present"] is None:
                vcf_s = "n/a"
            else:
                vcf_s = "PRESENT" if m["vcf_present"] else "absent"
            flag = ""
            if m["vcf_present"] is not None and m["present"] != m["vcf_present"]:
                flag = "   <-- JSON/VCF DISAGREE, review manually"
            detail = ""
            hits = m["hits"] or m["vcf_hits"]
            if hits:
                detail = "  on " + ", ".join(f"{h}({v})" for h, v in hits)
            print(f"    - {m['name']:<22} {m['rsid']:<12} {m['locus']:<16} "
                  f"JSON:{json_s:<8} VCF:{vcf_s:<8}{detail}{flag}")
        hg_str = ", ".join(f"{h}={hg}" for h, hg in r.get("smn1_haplogroups", []))
        print(f"  SMN1 haplogroups: {hg_str or 'n/a'}")
        if r.get("haplogroup_strong"):
            print("    ^ two-copy (duplication) haplogroup present — "
                  "strong 2+0 signal")
        elif r.get("haplogroup_marker"):
            print("    ^ marker-carrying haplogroup present — weak 2+0 signal")
        print(f"  CALL: {r['call']}")
        print(f"  {r['interpretation']}")
    print("=" * 72)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    files = collect_files(argv[1:])
    if not files:
        print("No paraphase JSON files found.", file=sys.stderr)
        return 1
    all_results = []
    for f in files:
        try:
            all_results.extend(analyze_file(f))
        except Exception as e:  # noqa: BLE001
            print(f"ERROR processing {f}: {e}", file=sys.stderr)
    print_report(all_results)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
