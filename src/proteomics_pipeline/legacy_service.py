"""Maintainer-local legacy regression (packet R12, FR-111/FR-112, V111/V112; SM12).

Everything here runs only on the Maintainer's machine against the original archive and a maintained run of the
private study. Nothing private is read by default, written into the source tree or printed: the archive is
verified by hashes against the public recovered inventory (`docs/audit/archive_inventory.json`), and comparison
tables go to a Maintainer-chosen output directory outside the repository.

- `verify_archive` checks the archive file hash and every member's size and SHA-256 against the inventory without
  extracting or modifying anything (V111). A missing archive is NOT_RUN (`E_LEGACY_UNAVAILABLE`), never a guess.
- `compare_tables` aligns a historical table and a maintained table on declared keys and reports, per field,
  the absolute and relative difference against the frozen tolerances (exact keys/counts; effects abs <= 1e-10 or
  rel <= 1e-8; statistics and P values abs <= 1e-8 or rel <= 1e-6). Rows present on one side only are
  NOT_REPRODUCED; differences are preserved, never reconciled (V112).
- A maintained value must come from a completed maintained run whose stage result hashes verify. A "maintained"
  table that is byte-identical to the historical table is refused (`E_LEGACY_COPIED`): stored values are not
  evidence that R was rerun.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import zipfile
from pathlib import Path

from .errors import ProteomicsError

CAPABILITY = "legacy"   # deliberately not advertised as a pipeline stage: this is Maintainer-local tooling (scripts/maintained/legacy_regression.py)
ROOT = Path(__file__).resolve().parents[2]
INVENTORY = ROOT / "docs" / "audit" / "archive_inventory.json"
TOLERANCES = {"exact": {"abs": 0.0, "rel": 0.0}, "effect": {"abs": 1e-10, "rel": 1e-8}, "statistic": {"abs": 1e-8, "rel": 1e-6}}


class LegacyError(ProteomicsError):
    def __init__(self, code: str, message: str, pointer: str | None = None, exit_code: int = 3):
        super().__init__(code, message, pointer, exit_code=exit_code)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _outside_repository(path: Path, what: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved == ROOT or ROOT in resolved.parents:
        raise LegacyError("E_PRIVATE_OUTPUT_IN_REPO", f"the {what} must lie outside the source repository so private material never enters it", exit_code=2)
    return resolved


def verify_archive(archive: str | Path, inventory: str | Path = INVENTORY) -> dict:
    """Hash-verify the Maintainer's original archive against the recovered inventory; never extracts or writes."""
    archive = Path(archive).expanduser()
    if not archive.is_file():
        raise LegacyError("E_LEGACY_UNAVAILABLE", "the original archive is not available on this machine; V111 stays NOT_RUN (Maintainer-only)")
    expected = json.loads(Path(inventory).read_text(encoding="utf-8"))
    observed_sha = _sha256(archive)
    problems = []
    if observed_sha != expected["archive_sha256"]:
        problems.append({"member": None, "reason": "archive_sha256_mismatch"})
    members = {}
    with zipfile.ZipFile(archive) as handle:
        for info in handle.infolist():
            if info.is_dir():
                continue
            data = handle.read(info.filename)   # bytes in memory only; nothing is extracted to disk
            members[info.filename] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    for item in expected["files"]:
        found = members.get(item["path"])
        if found is None:
            problems.append({"member": item["path"], "reason": "missing"})
        elif found["sha256"] != item["sha256"] or found["bytes"] != item["bytes"]:
            problems.append({"member": item["path"], "reason": "hash_or_size_mismatch"})
    extra = sorted(set(members) - {i["path"] for i in expected["files"]})
    problems += [{"member": m, "reason": "not_in_inventory"} for m in extra]
    return {"result_type": "LegacyArchiveVerification", "verified": not problems, "file_count": len(members), "expected_file_count": expected["file_count"],
            "archive_sha256_matches": observed_sha == expected["archive_sha256"], "problems": problems, "archive_modified": False}


def _read(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _close(a: float, b: float, tol: dict) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    if math.isnan(a) or math.isnan(b):
        return False
    diff = abs(a - b)
    return diff <= tol["abs"] or (b != 0 and diff / abs(b) <= tol["rel"])


def _number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def maintained_table(run_dir: str | Path, relative: str) -> Path:
    """A maintained table, only from a completed run whose stage result lists it with a matching hash."""
    path = Path(run_dir) / relative
    result = path.parent / "stage-result.json"
    if not path.is_file() or not result.is_file():
        raise LegacyError("E_LEGACY_MAINTAINED_RUN", f"{relative} is not an output of a completed maintained run")
    stage = json.loads(result.read_text(encoding="utf-8"))
    entry = next((o for o in stage.get("outputs", []) if o["relative_path"] == path.name), None)
    if stage.get("state") != "COMPLETED" or entry is None or entry["sha256"] != _sha256(path):
        raise LegacyError("E_LEGACY_MAINTAINED_RUN", f"{relative} is not a hash-verified output of a COMPLETED stage", exit_code=5)
    return path


def compare_tables(legacy: str | Path, maintained: str | Path, *, keys: list[str], fields: dict[str, dict]) -> dict:
    """fields: maintained column -> {"legacy": historical column, "kind": exact|effect|statistic}."""
    legacy, maintained = Path(legacy), Path(maintained)
    if _sha256(legacy) == _sha256(maintained):
        raise LegacyError("E_LEGACY_COPIED", "the maintained table is byte-identical to the historical table; stored values are not evidence that R was rerun", exit_code=2)
    old = {tuple(r[k] for k in keys): r for r in _read(legacy)}
    new = {tuple(r[k] for k in keys): r for r in _read(maintained)}
    rows, counts = [], {"MATCH": 0, "DIFFERENT": 0, "NOT_REPRODUCED": 0}
    for key in sorted(set(old) | set(new)):
        if key not in old or key not in new:
            rows.append({**dict(zip(keys, key)), "field": None, "status": "NOT_REPRODUCED", "side": "maintained only" if key not in old else "legacy only"})
            counts["NOT_REPRODUCED"] += 1
            continue
        for column, spec in fields.items():
            a_raw, b_raw = new[key][column], old[key][spec["legacy"]]
            tol = TOLERANCES[spec["kind"]]
            if spec["kind"] == "exact":
                same = a_raw == b_raw; diff = None
            else:
                a, b = _number(a_raw), _number(b_raw); same = _close(a, b, tol); diff = abs(a - b) if not (math.isnan(a) or math.isnan(b)) else None
            status = "MATCH" if same else "DIFFERENT"
            counts[status] += 1
            rows.append({**dict(zip(keys, key)), "field": column, "maintained": a_raw, "legacy": b_raw, "abs_difference": diff, "tolerance": tol, "status": status})
    return {"result_type": "LegacyComparison", "keys": keys, "counts": counts, "rows": rows, "tolerances": TOLERANCES}


def run(config_path: str | Path) -> dict:
    """Maintainer-local regression from configs/local-legacy.template.json (copied outside the repository)."""
    config_path = Path(config_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    output = _outside_repository(Path(config["output_dir"]), "output directory")
    if output.exists():
        raise LegacyError("E_PATH_COLLISION", f"output directory already exists: {output.name}", exit_code=2)
    archive = verify_archive(config["archive"])
    comparisons = []
    for table in config.get("tables", []):
        legacy = Path(table["legacy_table"]).expanduser()
        if not legacy.is_file():
            comparisons.append({"id": table["id"], "state": "NOT_RUN", "reason_code": "E_LEGACY_UNAVAILABLE"}); continue
        try:
            maintained = maintained_table(config["maintained_run"], table["maintained_table"])
        except LegacyError as error:
            comparisons.append({"id": table["id"], "state": "NOT_RUN", "reason_code": error.code, "reason": error.message}); continue
        result = compare_tables(legacy, maintained, keys=table["keys"], fields=table["fields"])
        comparisons.append({"id": table["id"], "state": "COMPLETED", **result})
    output.mkdir(parents=True)
    summary = {"result_type": "LegacyRegression", "archive": archive, "comparisons": comparisons,
               "note": "Maintainer-local and private; review before any disclosure. Differences are preserved, not reconciled."}
    (output / "legacy_regression.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return {"output": str(output), "archive_verified": archive["verified"], "comparisons": [{k: c.get(k) for k in ("id", "state", "reason_code", "counts")} for c in comparisons]}

