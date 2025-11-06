#!/usr/bin/env bash
set -euo pipefail

# repo root (this script is in scripts/)
cd "$(dirname "$0")/.."

# Prefer python3.11 if present
PY=python3
if command -v python3.11 >/dev/null 2>&1; then
  PY=python3.11
elif command -v py >/dev/null 2>&1; then
  # Windows py launcher sometimes exists on WSL/Git Bash
  PY="py -3.11"
fi

# Create venv if missing
if [ ! -d ".venv" ]; then
  $PY -m venv .venv
fi

# Activate venv
# shellcheck disable=SC1091
source .venv/bin/activate

# Upgrade tooling and install deps
python -m pip install --upgrade pip setuptools wheel
pip install -r requirements.txt

# Safer default backend for matplotlib in headless/remote runs
export MPLBACKEND=${MPLBACKEND:-Agg}

# Run Streamlit
exec python -m streamlit run src/my_project/MAINAPP.py
