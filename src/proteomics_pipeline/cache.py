"""Stage fingerprints and resume cache (packet R11, V103).

A stage output is reused only when (a) the directory holds a COMPLETED stage
result whose every output hash verifies, and (b) its fingerprint — the
normalized stage request (inputs with hashes, parameters, plan hash, RNG),
the code manifest and the environment inventory — equals the fingerprint of
the request the resumed run would send now.  Directory existence alone is
never a cache hit; interrupted temporaries are never reused.
"""
from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .provenance import canonical_json_bytes, sha256_bytes, sha256_file

VOLATILE_KEYS = ("output_temp_dir", "run_id")


def fingerprint(request: dict, *, code_sha256: str, environment_sha256: str) -> str:
    normalized = {k: v for k, v in request.items() if k not in VOLATILE_KEYS}
    params = dict(normalized.get("parameters") or {})
    snapshot = params.get("run_snapshot")
    if isinstance(snapshot, dict):   # report stages: the run identity is volatile, the stage states are not
        params["run_snapshot"] = {k: v for k, v in snapshot.items() if k != "run_id"}
    normalized["parameters"] = params
    return sha256_bytes(canonical_json_bytes({"request": normalized, "code": code_sha256, "environment": environment_sha256}))


def outputs_verify(stage_dir: Path) -> bool:
    result_path = stage_dir / "stage-result.json"
    if not result_path.is_file():
        return False
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result.get("state") != "COMPLETED":
        return False
    for output in result.get("outputs", []):
        path = stage_dir / output["relative_path"]
        if not path.is_file() or sha256_file(path) != output["sha256"]:
            return False
    return True


def reusable(stage_dir: Path, expected: str) -> bool:
    stored = stage_dir / "fingerprint.json"
    if not stored.is_file() or not outputs_verify(stage_dir):
        return False
    return json.loads(stored.read_text(encoding="utf-8")).get("fingerprint") == expected


def record(stage_dir: Path, value: str) -> None:
    (stage_dir / "fingerprint.json").write_text(json.dumps({"fingerprint": value, "rule": "normalized request + code manifest + environment"}), encoding="utf-8")


def retire(path: Path, root: Path, label: str) -> Path:
    """Move a stale or incomplete directory into logs/ (never deleted, never reused)."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    target = root / "logs" / f"{label}-{path.name.lstrip('.')}-{stamp}"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(path), str(target))
    return target


def abandon_temporaries(root: Path) -> list[str]:
    moved = []
    for temp in sorted(root.glob(".stage-*")):
        moved.append(retire(temp, root, "abandoned-temporary").name)
    return moved
