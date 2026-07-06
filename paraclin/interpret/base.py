"""Interpreter plugin interface.

Each interpreter turns a Paraphase gene record (+ optional VCF) into a
``Result`` that keeps the raw caller output verbatim and separate from the
derived clinical statement — a hard requirement for the accreditation goal.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field


# Status vocabulary shared across conditions; the UI maps these to colors.
# Not every status applies to every gene.
STATUS_LEVELS = {
    "affected": "critical",
    "not_affected": "normal",
    "silent_carrier_likely": "warning",
    "carrier": "warning",
    "inversion_detected": "critical",
    "deletion_detected": "critical",
    "possible": "warning",
    "not_carrier": "normal",
    "no_sv": "normal",
    "review": "review",
    "not_applicable": "info",
}


@dataclass
class Result:
    gene: str
    disease: str
    module_version: str

    headline: str                       # one-line clinical statement
    status: str                         # key from STATUS_LEVELS
    interpretation: str                 # prose explanation

    raw: dict = field(default_factory=dict)          # verbatim paraphase calls
    evidence: list[dict] = field(default_factory=list)  # [{label, value, note}]
    references: list[str] = field(default_factory=list)
    qc: dict = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    # Optional experimental flag + a nested secondary result rendered on its own
    # tab (e.g. SMA carrier status shown separately from affected status).
    experimental: bool = False
    secondary: dict | None = None
    secondary_tab_label: str | None = None

    @property
    def level(self) -> str:
        return STATUS_LEVELS.get(self.status, "info")

    def to_dict(self) -> dict:
        d = asdict(self)
        d["level"] = self.level
        return d


class Interpreter:
    """Base class. Subclasses set ``gene``, ``disease``, ``version`` and
    implement ``interpret``."""

    gene: str = ""
    disease: str = ""
    version: str = "0.0.0"

    def interpret(self, gene_record: dict, vcf_path: str | None, qc: dict) -> Result:
        raise NotImplementedError
