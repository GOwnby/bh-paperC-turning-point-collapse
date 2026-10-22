#!/usr/bin/env bash
# One-command verification of Paper C.   bash verify.sh [--quick|--full] [--checks-only] [--list] [--only step1,step2]
#   --quick (default)  re-runs every cheap generator and re-derives every quoted number (long jobs: from the shipped data) -- see README, "Verify"
#   --full             also re-runs the long jobs (composition survey, Lagrangian N = 400..1600, ...)
# Exit status 0 iff every check passes.  Last line: RESULT <repo>: PASSED n/n checks (quick|full)
set -u
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
PY="${PYTHON:-python3}"
command -v "$PY" >/dev/null 2>&1 || { echo "python3 not found"; exit 2; }
"$PY" -c 'import numpy' 2>/dev/null || { echo "numpy is required: pip install -r requirements.txt"; exit 2; }
exec "$PY" -B verify.py "$@"
