#!/usr/bin/env python3
"""Check the committed private-study status records (packet R15d: FR-188; ADR 0011, SM46).

Every record under docs/validation/016-methods-extension-validation/R15d/status/ must contain only the allowed keys
and must not contain a file name or path. A record that does is refused (E_PRIVATE_IDENTIFIER) and the check fails.
Usage: python scripts/maintained/verify_private_status.py [--status-dir DIR]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "docs" / "validation" / "016-methods-extension-validation" / "R15d" / "status"
ALLOWED = {"schema", "slot", "design", "state", "checks", "check_counts"}
PATH_LIKE = re.compile(r"[/\\]|\.(raw|mzML|txt|tsv|csv|xlsx|gz|zip)\b", re.IGNORECASE)


def check_record(path: Path) -> list[str]:
    problems = []
    record = json.loads(path.read_text(encoding="utf-8"))
    extra = set(record) - ALLOWED
    if extra:
        problems.append(f"{path.name}: keys outside the status record: {sorted(extra)}")
    values = json.dumps({k: v for k, v in record.items() if k != "schema"})   # the schema identifier itself contains a slash
    if PATH_LIKE.search(values):
        problems.append(f"{path.name}: a file name or path appears in the status record (E_PRIVATE_IDENTIFIER)")
    return problems


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--status-dir", type=Path, default=DEFAULT)
    args = parser.parse_args(argv)
    problems = []
    records = sorted(args.status_dir.glob("*.json")) if args.status_dir.is_dir() else []
    for path in records:
        problems.extend(check_record(path))
    print(json.dumps({"records": len(records), "problems": problems, "state": "PASS" if not problems else "FAIL"}))
    return 0 if not problems else 5


if __name__ == "__main__":
    sys.exit(main())
