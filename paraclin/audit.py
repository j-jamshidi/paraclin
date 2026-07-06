"""Auth + audit seam.

Local single-user phase: the current user is a fixed local identity and the audit
trail is an append-only JSONL file. The interfaces (``current_user``,
``record``) are what the rest of the app calls, so a future server phase can swap
in real OIDC auth + a DB-backed audit table without touching call sites.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

from . import __version__ as APP_VERSION
from .settings import get_settings


def current_user() -> str:
    # Server phase: derive from the authenticated session instead.
    return "local"


def record(action: str, **details) -> None:
    """Append one audit event. Never raises into the request path."""
    entry = {
        "ts": datetime.now(timezone.utc).isoformat(),
        "user": current_user(),
        "app_version": APP_VERSION,
        "action": action,
        **details,
    }
    try:
        path = get_settings().audit_log
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a") as fh:
            fh.write(json.dumps(entry) + "\n")
    except OSError:
        pass


def provenance_stamp(sample) -> dict:
    """Provenance block attached to every generated artifact / report."""
    return {
        "app_version": APP_VERSION,
        "paraphase_version": _paraphase_version(),
        "genome_build": sample.build,
        "sample_id": sample.sample_id,
        "json_md5": sample.json_md5,
        "bam_md5": sample.bam_md5,
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }


def _paraphase_version() -> str | None:
    try:
        import importlib.util
        import sys

        repo = get_settings().paraphase_repo
        if str(repo) not in sys.path:
            sys.path.insert(0, str(repo))
        spec = importlib.util.find_spec("paraphase")
        if spec:
            import paraphase  # noqa: WPS433
            return getattr(paraphase, "__version__", None)
    except Exception:  # noqa: BLE001
        return None
    return None
