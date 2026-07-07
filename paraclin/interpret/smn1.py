"""SMN1 / Spinal Muscular Atrophy interpreter.

Turns Paraphase SMN1 copy-number + haplotype output into a clinical statement:
  * SMN1 CN 0  -> biallelic loss, SMA-affected genotype
  * SMN1 CN 1  -> carrier
  * SMN1 CN 2  -> 1+1 (not a carrier) vs 2+0 silent carrier, assessed from
                  duplication-linked markers and Paraphase haplogroups
  * SMN1 CN >=3 -> not a carrier (multiple copies)

Evidence sources (verified against dbSNP / Chen et al. 2023, AJHG):
  * Silent-carrier markers  rs143838139 (chr5:70952074 T>G),
                            rs200800214 (chr5:70952646_70952647 delAT)
  * Duplication haplogroups S1-8, S1-9d (two-copy alleles); S1-9/S1-9c carry the
                            g.27134T>G marker
  * SMN2 c.859G>C positive modifier rs121909192 (chr5:70076545 G>C)
"""

from __future__ import annotations

from .base import Interpreter, Result

VERSION = "1.1.0"  # v1.1: adds SMN1 exon 7 / SMNΔ7-8 deletion evaluation

MARKERS = [
    {"name": "g.27134T>G", "rsid": "rs143838139", "chrom": "chr5", "pos": 70952074,
     "ref": "T", "alt": "G", "kind": "snv", "window": (70952074, 70952074)},
    {"name": "g.27706_27707delAT", "rsid": "rs200800214", "chrom": "chr5", "pos": 70952646,
     "ref": "TAT", "alt": "T", "kind": "del", "window": (70952644, 70952648)},
]
DUP_HAPLOGROUPS_STRONG = {"S1-8", "S1-9d"}
DUP_HAPLOGROUPS_MARKER = {"S1-8", "S1-9", "S1-9c", "S1-9d"}

MODIFIER = {"name": "SMN2 c.859G>C", "rsid": "rs121909192", "chrom": "chr5",
            "pos": 70076545, "ref": "G", "alt": "C"}

REFERENCES = [
    "Chen et al. 2023, AJHG — Paraphase SMN1/SMN2 (PMC9943720)",
    "Luo et al. 2014, Genet Med — SMN1 duplication haplotype markers",
    "dbSNP rs143838139, rs200800214, rs121909192",
]


def _variant_pos(vstr: str):
    try:
        return int(vstr.split("_", 1)[0])
    except (ValueError, IndexError):
        return None


def _marker_in_variants(marker: dict, variant_list: list[str]) -> str | None:
    lo, hi = marker["window"]
    for v in variant_list:
        p = _variant_pos(v)
        if p is None:
            continue
        parts = v.split("_")
        if marker["kind"] == "snv":
            if p == marker["pos"] and len(parts) >= 3 and parts[2].upper() == marker["alt"].upper():
                return v
        else:
            if lo <= p <= hi and len(parts) >= 3 and len(parts[1]) > len(parts[2]):
                return v
    return None


def _scan_vcf(vcf_path: str, targets: list[dict], hap_filter: str) -> dict:
    """Independent VCF cross-check. Returns {name: [(col, "REF>ALT@POS"), ...]}.
    Only counts genotype==1 on columns whose name contains ``hap_filter``."""
    import gzip
    opener = gzip.open if vcf_path.endswith(".gz") else open
    out = {t["name"]: [] for t in targets}
    cols: list[tuple[int, str]] = []
    try:
        with opener(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                f = line.rstrip("\n").split("\t")
                if line.startswith("#CHROM"):
                    cols = [(i, n) for i, n in enumerate(f) if i >= 9 and hap_filter in n]
                    continue
                if not cols:
                    continue
                chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
                fmt = f[8].split(":")
                if "GT" not in fmt:
                    continue
                gt_idx = fmt.index("GT")
                for t in targets:
                    if chrom != t["chrom"]:
                        continue
                    if t.get("kind", "snv") == "snv":
                        if not (pos == t["pos"] and ref.upper() == t["ref"].upper()
                                and alt.upper() == t["alt"].upper()):
                            continue
                    else:
                        lo, hi = t.get("window", (t["pos"], t["pos"]))
                        if not (lo <= pos <= hi and len(ref) > len(alt)):
                            continue
                    for i, name in cols:
                        if f[i].split(":")[gt_idx] == "1":
                            out[t["name"]].append((name, f"{ref}>{alt}@{pos}"))
    except OSError:
        pass
    return out


AFFECTED_REFERENCES = [
    "Chen et al. 2023, AJHG — Paraphase SMN1/SMN2 (PMC9943720)",
    "GeneReviews: Spinal Muscular Atrophy (NBK1352)",
    "Homozygous absence of SMN1 exon 7 causes ~95% of SMA (Lefebvre 1995; Prior, GeneReviews)",
    "Hybrid SMN (SMN2ex7–SMN1ex8) / gene conversion associates with milder phenotype "
    "(Niba et al. 2020, Brain Dev)",
]


def _exon7_evidence(smn1_cn, smn2_cn, del78_cn) -> tuple[list[dict], str]:
    """Build the exon-7 / SMNΔ7-8 (exon 7-8) deletion evidence lines and, when
    SMN1 exon 7 is lost, a sentence describing the likely mechanism (physical
    deletion vs SMN1->SMN2 gene conversion)."""
    if smn1_cn is None:
        exon7_status = "indeterminate"
    elif smn1_cn == 0:
        exon7_status = "absent on both alleles (homozygous loss)"
    elif smn1_cn == 1:
        exon7_status = "1 copy (one allele lost)"
    else:
        exon7_status = f"{smn1_cn} copies present"

    evidence = [
        {"label": "SMN1 exon 7 status", "value": exon7_status,
         "note": "Exon 7 (c.840C) is the functionally critical exon; its "
                 "homozygous absence causes ~95% of SMA."},
        {"label": "SMNΔ7-8 (exon 7-8) deletion copies", "value": del78_cn,
         "note": "Paraphase count of SMN alleles with exon 7-8 physically deleted."},
    ]

    mechanism = ""
    if smn1_cn is not None and smn1_cn < 2:
        if del78_cn:
            mechanism = (
                f"Paraphase detects {del78_cn} SMNΔ7-8 deletion allele(s), consistent "
                "with the exon 7-8 deletion mechanism (the cause of ~95% of SMA).")
        else:
            mechanism = (
                "No SMNΔ7-8 deletion allele was detected, so the exon 7 loss may reflect "
                "SMN1->SMN2 gene conversion (which raises SMN2 copy number) or a hybrid "
                "rearrangement rather than a simple deletion — correlate with the SMN2 "
                "copy number. Hybrid/conversion alleles are associated with a milder phenotype.")
    return evidence, mechanism


class Smn1Interpreter(Interpreter):
    gene = "smn1"
    disease = "Spinal Muscular Atrophy (SMA)"
    version = VERSION

    def interpret(self, rec: dict, vcf_path: str | None, qc: dict) -> Result:
        # SMN2 c.859G>C positive modifier (checked on SMN2 haplotypes) — computed
        # once and shared by both the affected and carrier views.
        modifier_present = False
        if vcf_path:
            mhits = _scan_vcf(vcf_path, [MODIFIER], "smn2hap").get(MODIFIER["name"], [])
            modifier_present = bool(mhits)

        affected = self._affected_result(rec, qc, modifier_present)
        carrier = self._carrier_result(rec, vcf_path, qc)
        affected.secondary = carrier.to_dict()
        affected.secondary_tab_label = "Carrier status (experimental)"
        return affected

    # ------------------------------------------------------------------ #
    # Primary view: affected vs not-affected (diagnostic)
    # ------------------------------------------------------------------ #
    def _affected_result(self, rec: dict, qc: dict, modifier_present: bool) -> Result:
        smn1_cn = rec.get("smn1_cn")
        smn2_cn = rec.get("smn2_cn")
        del78_cn = rec.get("smn_del78_cn")
        exon7_evidence, mechanism = _exon7_evidence(smn1_cn, smn2_cn, del78_cn)
        raw = {
            "smn1_cn": smn1_cn, "smn2_cn": smn2_cn, "smn_del78_cn": del78_cn,
        }
        evidence = [
            {"label": "SMN1 copy number", "value": smn1_cn},
            {"label": "SMN2 copy number", "value": smn2_cn,
             "note": "Main severity modifier when affected."},
            *exon7_evidence,
            {"label": "SMN2 c.859G>C modifier (rs121909192)",
             "value": "present" if modifier_present else "absent",
             "note": "Positive modifier — tends to milder phenotype." if modifier_present
             else "Absent (reference G on SMN2)."},
        ]

        def base(status, headline, interpretation, caveats=None):
            return Result(
                gene=self.gene, disease=self.disease, module_version=self.version,
                headline=headline, status=status, interpretation=interpretation,
                raw=raw, evidence=evidence, references=AFFECTED_REFERENCES, qc=qc,
                caveats=caveats or [],
            )

        if smn1_cn is None:
            return base("review", "SMA — copy number ambiguous, review required",
                        "Paraphase could not unambiguously resolve the SMN1 copy number.",
                        ["Ambiguous (null) SMN1 copy number; manual review required."])

        if smn1_cn == 0:
            sev = ("With 2 SMN2 copies this typically sits at the more severe end; "
                   if smn2_cn == 2 else
                   f"With {smn2_cn} SMN2 copies the phenotype trends milder; "
                   if smn2_cn and smn2_cn >= 3 else "")
            return base(
                "affected",
                f"SMA — affected genotype (homozygous SMN1 exon 7 loss, {smn2_cn} SMN2)",
                "Homozygous absence of SMN1 exon 7 (copy number 0) — the defining molecular "
                "lesion in ~95% of SMA and a genotype consistent with being affected. "
                + (mechanism + " " if mechanism else "") + sev
                + ("The c.859G>C positive modifier is present."
                   if modifier_present else "The c.859G>C positive modifier is absent.")
                + " Requires clinical correlation and confirmation.",
            )

        return base(
            "not_affected",
            f"SMA — not affected ({smn1_cn} SMN1 exon 7 copies)",
            f"{smn1_cn} functional SMN1 (exon 7) copy(ies) present, so this is not an "
            "SMA-affected genotype. "
            + (mechanism + " " if mechanism else "")
            + "Carrier status is assessed separately (experimental tab).",
        )

    # ------------------------------------------------------------------ #
    # Secondary view: carrier / silent-carrier status (EXPERIMENTAL)
    # ------------------------------------------------------------------ #
    def _carrier_result(self, rec: dict, vcf_path: str | None, qc: dict) -> Result:
        smn1_cn = rec.get("smn1_cn")
        del78_cn = rec.get("smn_del78_cn")
        hap_details = rec.get("haplotype_details", {}) or {}
        raw = {
            "smn1_cn": smn1_cn, "smn_del78_cn": del78_cn,
            "smn1_haplotypes": list((rec.get("smn1_haplotypes", {}) or {}).values()),
        }
        evidence: list[dict] = [
            {"label": "SMN1 copy number", "value": smn1_cn},
            {"label": "SMNΔ7-8 (exon 7-8) deletion copies", "value": del78_cn,
             "note": "A deletion allele on one chromosome is the usual carrier mechanism."},
        ]
        caveats = [
            "EXPERIMENTAL — carrier / silent-carrier (2+0) assessment is for research "
            "use and must not be used for clinical carrier reporting without validation.",
        ]

        def base(status, headline, interpretation):
            return Result(
                gene=self.gene, disease=self.disease, module_version=self.version,
                headline=headline, status=status, interpretation=interpretation,
                raw=raw, evidence=evidence, references=REFERENCES, qc=qc,
                caveats=caveats, experimental=True,
            )

        if smn1_cn is None:
            return base("review", "Carrier status — ambiguous copy number",
                        "SMN1 copy number could not be resolved; carrier status "
                        "cannot be assessed.")
        if smn1_cn == 0:
            return base("not_applicable", "Carrier status — not applicable (affected)",
                        "0 functional SMN1 copies (affected genotype); carrier status "
                        "does not apply.")
        if smn1_cn == 1:
            mech = (f" Paraphase detects {del78_cn} SMNΔ7-8 deletion allele(s), consistent "
                    "with an exon 7-8 deletion on the other chromosome." if del78_cn else "")
            return base("carrier", "Carrier — 1 SMN1 exon 7 copy",
                        "One functional SMN1 (exon 7) copy — an SMA carrier." + mech
                        + " Reproductive/genetic counselling indicated; partner testing "
                        "recommended.")
        if smn1_cn >= 3:
            return base("not_carrier", f"Not a carrier — {smn1_cn} SMN1 copies",
                        f"{smn1_cn} SMN1 copies. Not a carrier by copy number.")

        # smn1_cn == 2 : the 1+1 vs 2+0 question
        marker_hits: list[str] = []
        vcf_marker = _scan_vcf(vcf_path, MARKERS, "smn1hap") if vcf_path else {}
        for marker in MARKERS:
            json_hit = any(
                _marker_in_variants(marker, hap_details.get(h, {}).get("variants", []) or [])
                for h in [n for n in hap_details if "smn1hap" in n]
            )
            present = json_hit or bool(vcf_marker.get(marker["name"]))
            if present:
                marker_hits.append(marker["name"])
            evidence.append({
                "label": f"{marker['name']} ({marker['rsid']})",
                "value": "present" if present else "absent",
                "note": f"{marker['chrom']}:{marker['pos']}",
            })

        hg_values = [hap_details.get(h, {}).get("haplogroup")
                     for h in [n for n in hap_details if "smn1hap" in n]]
        hg_values = [hg for hg in hg_values if hg]
        evidence.append({"label": "SMN1 haplogroups", "value": ", ".join(hg_values) or "none"})
        strong_hg = [hg for hg in hg_values if hg in DUP_HAPLOGROUPS_STRONG]
        marker_hg = [hg for hg in hg_values
                     if hg in DUP_HAPLOGROUPS_MARKER and hg not in DUP_HAPLOGROUPS_STRONG]

        caveats.append(
            "Per Chen et al. 2023, a 2+0 silent carrier cannot be confirmed from a "
            "single sample without pedigree data. These signals are population-specific "
            "and tag only a subset of duplication alleles, so a residual risk remains "
            "even when negative (especially outside African / Ashkenazi Jewish ancestry).")

        if marker_hits or strong_hg:
            reasons = []
            if marker_hits:
                reasons.append("a duplication-linked marker SNP is present")
            if strong_hg:
                reasons.append(f"an SMN1 haplotype is a known two-copy haplogroup "
                               f"({', '.join(strong_hg)})")
            return base(
                "silent_carrier_likely",
                "Silent carrier (2+0) likely — 2 copies in cis",
                (reasons[0][0].upper() + reasons[0][1:]
                 + ("; " + "; ".join(reasons[1:]) if len(reasons) > 1 else "")
                 + ". This suggests two SMN1 copies in cis (2+0) — a silent SMA carrier "
                 "despite a copy number of 2. Confirm with family studies / orthogonal testing."),
            )
        if marker_hg:
            return base(
                "possible", "Possible 2+0 — weak haplogroup evidence",
                f"An SMN1 haplotype is in haplogroup {', '.join(marker_hg)}, which carries "
                "the g.27134T>G marker but is not itself a confirmed two-copy allele. This "
                "raises but does not establish the chance of a 2+0 configuration.",
            )
        return base(
            "not_carrier", "Not a carrier — 2 copies, consistent with 1+1",
            "No duplication-linked marker SNP and no duplication-associated SMN1 "
            "haplogroup were found. Consistent with a 1+1 configuration (not a carrier).",
        )
