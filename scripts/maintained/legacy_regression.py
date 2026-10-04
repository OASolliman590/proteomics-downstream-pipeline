#!/usr/bin/env python3
"""Maintainer-local legacy regression (R12 V111/V112). Runs only on the Maintainer's machine.

  legacy_regression.py verify-archive --archive /path/outside/repo/original.zip
  legacy_regression.py run --config /path/outside/repo/local-legacy.json

Copy configs/local-legacy.template.json outside the repository and fill in local paths. The archive is
hash-verified against docs/audit/archive_inventory.json without extraction; comparison output goes to a private
directory outside the repository. Exit codes: 0 completed, 2 refused, 3 NOT_RUN (archive unavailable), 5 integrity.
See docs/validation/013-release/MAINTAINER_RUNBOOK.md.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from proteomics_pipeline.errors import ProteomicsError  # noqa: E402
from proteomics_pipeline import legacy_service  # noqa: E402


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    v = sub.add_parser("verify-archive"); v.add_argument("--archive", type=Path, required=True)
    r = sub.add_parser("run"); r.add_argument("--config", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "verify-archive":
            result = legacy_service.verify_archive(args.archive)
            print(json.dumps({k: result[k] for k in ("verified", "file_count", "expected_file_count", "archive_sha256_matches")} | {"n_problems": len(result["problems"])}))
            return 0 if result["verified"] else 5
        print(json.dumps(legacy_service.run(args.config)))
        return 0
    except ProteomicsError as error:
        state = "NOT_RUN" if error.code == "E_LEGACY_UNAVAILABLE" else "REFUSED"
        print(json.dumps({"state": state, "error": {"code": error.code, "message": error.message}}))
        return error.exit_code


if __name__ == "__main__":
    sys.exit(main())
