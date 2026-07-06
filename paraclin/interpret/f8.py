"""F8 / Hemophilia A interpreter.

Paraphase detects two structural variants in the int22h homology region and
reports them in the ``sv_called`` field (also emitted as SVs in the VCF):
  * ``inversion`` (int22h-1 / int22h-3) -> intron-22 inversion, ~45% of severe
    Hemophilia A.
  * ``deletion``  (int22h-1 / int22h-2) -> exon 1-22 deletion.
F8 is on chrX, so haplotype counts are sex-dependent (3 in males, 6 in females).
"""

from __future__ import annotations

from .base import Interpreter, Result

VERSION = "1.0.0"

REFERENCES = [
    "Chen et al. 2023, AJHG — Paraphase F8 intron-22 inversion (PMC9943720)",
    "Paraphase docs/F8.md",
    "GeneReviews: Hemophilia A (NBK1404)",
]


class F8Interpreter(Interpreter):
    gene = "f8"
    disease = "Hemophilia A (F8)"
    version = VERSION

    def interpret(self, rec: dict, vcf_path: str | None, qc: dict) -> Result:
        sv_called = rec.get("sv_called") or {}
        sample_sex = rec.get("sample_sex")
        total_cn = rec.get("total_cn")

        svs = {str(v).lower() for v in sv_called.values()}
        raw = {
            "sv_called": sv_called,
            "sample_sex": sample_sex,
            "total_cn": total_cn,
            "flanking_summary": rec.get("flanking_summary"),
        }
        expected_haps = 6 if sample_sex == "female" else 3 if sample_sex == "male" else None
        evidence = [
            {"label": "SV called", "value": sv_called or "none"},
            {"label": "Sample sex", "value": sample_sex or "unknown"},
            {"label": "Total haplotypes (F8 region)", "value": total_cn,
             "note": f"expected ~{expected_haps} for {sample_sex}" if expected_haps else None},
        ]
        caveats = []
        if sample_sex is None:
            caveats.append("Sample sex not reported by Paraphase; expected haplotype "
                           "count (3 male / 6 female) cannot be checked.")

        def base(status, headline, interpretation):
            return Result(
                gene=self.gene, disease=self.disease, module_version=self.version,
                headline=headline, status=status, interpretation=interpretation,
                raw=raw, evidence=evidence, references=REFERENCES, qc=qc, caveats=caveats,
            )

        if "inversion" in svs:
            return base(
                "inversion_detected", "Hemophilia A — intron-22 inversion detected",
                "Paraphase called an intron-22 inversion (recombination between int22h-1 "
                "and int22h-3), the most common cause of severe Hemophilia A. Confirm and "
                "correlate clinically.",
            )
        if "deletion" in svs:
            return base(
                "deletion_detected", "Hemophilia A — exon 1-22 deletion detected",
                "Paraphase called a deletion between int22h-1 and int22h-2, consistent "
                "with an F8 exon 1-22 deletion. Confirm and correlate clinically.",
            )
        return base(
            "no_sv", "Hemophilia A — no int22h structural variant detected",
            "No intron-22 inversion or exon 1-22 deletion was detected in the F8 int22h "
            "region. Note: this analysis targets the int22h-mediated SVs only; it does not "
            "exclude point mutations or other F8 variants causing Hemophilia A.",
        )
