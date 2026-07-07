#!/usr/bin/env bash
# One-command backend launcher for paraclin.
# Creates a local .venv, installs deps into it, and starts the API using that
# venv's own interpreter — so it never accidentally uses Anaconda's base python
# (a common cause of "ModuleNotFoundError: No module named 'fastapi'").
set -euo pipefail

cd "$(dirname "$0")"

PORT="${PORT:-8077}"
PYTHON="${PYTHON:-python3}"
VENV=".venv"

if [ ! -d "$VENV" ]; then
  echo "[paraclin] creating virtualenv in $VENV ..."
  "$PYTHON" -m venv "$VENV"
fi

VENV_PY="$VENV/bin/python"
[ -x "$VENV_PY" ] || VENV_PY="$VENV/Scripts/python.exe"   # Windows fallback

echo "[paraclin] installing dependencies ..."
"$VENV_PY" -m pip install --quiet --upgrade pip
"$VENV_PY" -m pip install --quiet -r requirements.txt

echo "[paraclin] starting API on http://localhost:$PORT (Ctrl-C to stop) ..."
exec "$VENV_PY" -m uvicorn paraclin.app:app --port "$PORT" "$@"
