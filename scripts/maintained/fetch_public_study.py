#!/usr/bin/env python3
"""Download a public study only against a matching Maintainer approval (packet R15d: FR-187; ADR 0011, SM46).

The approval record must name the accession, the URL and the SHA-256 of the file. A missing or mismatching record
refuses the download before any network call (E_PUBLIC_DATA_UNAPPROVED); a file whose observed SHA-256 differs from
the approved one is removed and refused (E_MANIFEST_HASH). The destination must lie outside the repository
(E_PUBLIC_DATA_UNAPPROVED). The approval, the URL and the observed hash are returned for the record. Analysis never
calls this module; it is an explicit action.

Usage: python scripts/maintained/fetch_public_study.py --approval APPROVAL.json --accession ACC --dest DIR [--url URL]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from proteomics_pipeline.errors import ProteomicsError  # noqa: E402


def approved_entry(approval: dict, accession: str) -> dict:
    for entry in approval.get("accessions", []):
        if entry.get("accession") == accession:
            return entry
    raise ProteomicsError("E_PUBLIC_DATA_UNAPPROVED", f"accession {accession!r} is not in the approval record", "/accessions", exit_code=2)


def fetch_public(approval: dict, accession: str, dest_dir: str | Path, *, url: str | None = None, opener=urllib.request.urlopen,
                 repo_root: Path = ROOT) -> dict:
    entry = approved_entry(approval, accession)
    if approval.get("approved_by") in (None, "") or approval.get("approved_on") in (None, ""):
        raise ProteomicsError("E_PUBLIC_DATA_UNAPPROVED", "the approval record has no approver or date", "/approved_by", exit_code=2)
    if url is not None and url != entry["url"]:
        raise ProteomicsError("E_PUBLIC_DATA_UNAPPROVED", "the requested URL differs from the approved URL", "/url", exit_code=2)
    dest = Path(dest_dir).resolve()
    if repo_root.resolve() in [dest, *dest.parents]:
        raise ProteomicsError("E_PUBLIC_DATA_UNAPPROVED", "public data are stored outside the repository; the destination is inside it", "/dest", exit_code=2)
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / entry["file_name"]
    with opener(entry["url"]) as response:   # the only network call; never reached without an approval
        payload = response.read()
    observed = hashlib.sha256(payload).hexdigest()
    if observed != entry["sha256"]:
        raise ProteomicsError("E_MANIFEST_HASH", "the downloaded file differs from the approved SHA-256; the file is not kept", "/sha256", exit_code=5)
    target.write_bytes(payload)
    return {"accession": accession, "approval_ref": approval.get("approval_ref"), "url": entry["url"], "observed_sha256": observed,
            "approved_by": approval.get("approved_by"), "approved_on": approval.get("approved_on"), "bytes": len(payload)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--approval", type=Path, required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--url")
    args = parser.parse_args(argv)
    try:
        record = fetch_public(json.loads(args.approval.read_text(encoding="utf-8")), args.accession, args.dest, url=args.url)
    except ProteomicsError as error:
        print(json.dumps({"state": "REFUSED", "code": error.code, "message": error.message}))
        return 2
    print(json.dumps({"state": "COMPLETED", **{k: v for k, v in record.items() if k != "bytes"}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
