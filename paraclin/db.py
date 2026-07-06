"""SQLite persistence for the sample index (SQLAlchemy 2.0 style).

Stores *references* to primary data (paths + checksums), never copies. The index
is disposable — it can always be rebuilt from the results folder via a rescan.
"""

from __future__ import annotations

from datetime import datetime, timezone

from sqlalchemy import JSON, DateTime, Float, String, create_engine
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, sessionmaker

from .settings import get_settings


class Base(DeclarativeBase):
    pass


class Sample(Base):
    __tablename__ = "samples"

    sample_id: Mapped[str] = mapped_column(String, primary_key=True)
    build: Mapped[str] = mapped_column(String, index=True)

    json_path: Mapped[str] = mapped_column(String)
    bam_path: Mapped[str | None] = mapped_column(String, nullable=True)
    bai_path: Mapped[str | None] = mapped_column(String, nullable=True)

    # gene -> vcf path, and gene -> per-gene summary (cn/sv/depth) as JSON blobs.
    vcf_paths: Mapped[dict] = mapped_column(JSON, default=dict)
    genes: Mapped[list] = mapped_column(JSON, default=list)

    sample_sex: Mapped[str | None] = mapped_column(String, nullable=True)
    genome_depth: Mapped[float | None] = mapped_column(Float, nullable=True)

    # integrity + provenance
    json_md5: Mapped[str | None] = mapped_column(String, nullable=True)
    bam_md5: Mapped[str | None] = mapped_column(String, nullable=True)
    paraphase_version: Mapped[str | None] = mapped_column(String, nullable=True)

    # QC / trio-validation warnings collected at index time.
    warnings: Mapped[list] = mapped_column(JSON, default=list)

    indexed_at: Mapped[datetime] = mapped_column(
        DateTime, default=lambda: datetime.now(timezone.utc)
    )

    def to_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "build": self.build,
            "genes": self.genes,
            "sample_sex": self.sample_sex,
            "genome_depth": self.genome_depth,
            "json_md5": (self.json_md5 or "")[:12],
            "bam_md5": (self.bam_md5 or "")[:12],
            "paraphase_version": self.paraphase_version,
            "has_bam": bool(self.bam_path),
            "warnings": self.warnings,
            "indexed_at": self.indexed_at.isoformat() if self.indexed_at else None,
        }


_engine = None
_Session = None


def init_db():
    global _engine, _Session
    if _engine is None:
        settings = get_settings()
        settings.database.parent.mkdir(parents=True, exist_ok=True)
        _engine = create_engine(settings.db_url, future=True)
        Base.metadata.create_all(_engine)
        _Session = sessionmaker(bind=_engine, future=True, expire_on_commit=False)
    return _Session


def get_session():
    if _Session is None:
        init_db()
    return _Session()
