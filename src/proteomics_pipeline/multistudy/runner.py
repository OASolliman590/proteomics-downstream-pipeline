"""Validate, run and verify one study (packet R15d: FR-185, FR-186, FR-188; ADR 0011).

Each step reports its own state: validate (manifest and input hashes), run (the declared analysis configuration, run
in a local run root that the caller chooses, never in the repository), verify (the proteomics verify checks of the
run), and the design-specific acceptance. Only a non-identifying status record can be committed (FR-188)."""
from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path

from ..errors import ProteomicsError
from . import acceptance as acceptance_checks
from .manifest import validate_manifest

STATUS_SCHEMA = "methods-extension-study-status/1.0.0"


def run_study(manifest_path: str | Path, run_root: str | Path, *, approvals: dict | None = None) -> dict:
    from .. import workflow
    from ..cli import main as cli_main
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    checks = {}
    validation = validate_manifest(manifest, manifest_path.parent, approvals=approvals)     # step 1: validate
    checks["validate"] = validation["state"]
    config = (manifest_path.parent / manifest["config"]).resolve()
    payload, code = workflow.run_command(config, Path(run_root))                             # step 2: run
    checks["run"] = payload["state"]
    if code == 0 and payload["state"] == "COMPLETED":                                        # step 3: verify
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            verify_code = cli_main(["verify", "--run", str(run_root), "--json"])
        checks["verify"] = "PASS" if verify_code == 0 else "FAIL"
    else:
        checks["verify"] = "NOT_RUN"
    acc = acceptance_checks.acceptance(run_root, manifest["design"]) if checks["run"] == "COMPLETED" else {"check": f"{manifest['design']}_endpoints", "state": "NOT_RUN"}
    checks["acceptance"] = acc["state"]
    passed = sum(1 for v in checks.values() if v in ("PASS", "COMPLETED"))
    state = "PASS" if checks["validate"] == "PASS" and checks["run"] == "COMPLETED" and checks["verify"] == "PASS" and checks["acceptance"] == "PASS" else "FAIL"
    return {"schema": STATUS_SCHEMA, "slot": manifest["slot"], "design": manifest["design"], "origin": manifest["origin"], "state": state,
            "checks": checks, "check_counts": {"passed": passed, "total": len(checks)}, "_study_id": str(manifest["study_id"]), "_files": [f["name"] for f in manifest["input"].get("files", [])]}


def status_record(result: dict) -> dict:
    """The only part of a study that may be committed: slot, design type, PASS or FAIL, and check counts (FR-188)."""
    return {"schema": result["schema"], "slot": result["slot"], "design": result["design"], "state": result["state"], "checks": result["checks"],
            "check_counts": result["check_counts"]}


def commit_status(result: dict, destination: str | Path) -> Path:
    """Write the non-identifying status record; refuse if it would carry the study identifier or any input file name (FR-188)."""
    record = status_record(result)
    text = json.dumps(record, indent=2, sort_keys=True)
    for identifier in [result["_study_id"], *result["_files"]]:
        if identifier and identifier in text:
            raise ProteomicsError("E_PRIVATE_IDENTIFIER", "the status record would carry a private identifier; the commit is refused", "/status", exit_code=5)
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text + "\n", encoding="utf-8", newline="\n")
    return path
