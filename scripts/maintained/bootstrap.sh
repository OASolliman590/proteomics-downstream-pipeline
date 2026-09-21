#!/usr/bin/env bash
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd); cd "$ROOT"; PYTHON_BIN=${PYTHON_BIN:-python3}; "$PYTHON_BIN" -m venv .venv; .venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m proteomics_pipeline --help >/dev/null
.venv/bin/python -m proteomics_pipeline --version
.venv/bin/python -m proteomics_pipeline validate --config configs/examples/example-independent.json --schema-only --json >/dev/null
R_AVAILABLE=0
if command -v Rscript >/dev/null 2>&1; then
  R_AVAILABLE=1
  mkdir -p .r-lib
  export R_LIBS_USER="$ROOT/.r-lib"
  R CMD INSTALL --library="$R_LIBS_USER" r/proteomicsCore
else
  export R_LIBS_USER="$ROOT/.r-lib"
fi
set +e
.venv/bin/python -m proteomics_pipeline doctor --json
DOCTOR_RC=$?
set -e
if [ "$R_AVAILABLE" -eq 1 ]; then
  [ "$DOCTOR_RC" -eq 0 ] || exit "$DOCTOR_RC"
else
  [ "$DOCTOR_RC" -eq 3 ] || exit "$DOCTOR_RC"
fi
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python -m pytest tests/contract/foundation tests/unit/foundation -q
.venv/bin/python scripts/maintained/check_baseline.py
if [ "$R_AVAILABLE" -eq 0 ]; then
  set +e
  .venv/bin/python scripts/maintained/write_solved_versions.py
  SOLVED_RC=$?
  set -e
  [ "$SOLVED_RC" -eq 3 ] || exit "$SOLVED_RC"
  echo 'NOT_RUN: Rscript is unavailable' >&2
  exit 3
fi
R_LIBS_USER="$ROOT/.r-lib" Rscript --vanilla scripts/maintained/test_r.R --suite foundation
.venv/bin/python scripts/maintained/write_solved_versions.py
exit "$DOCTOR_RC"
