"""Offline GMT pathway libraries (packet R15b: FR-174 to FR-178, SM44; ADR 0011).

Only user-supplied local GMT files are read. Each library is hashed before any analysis (E_RESOURCE_HASH on a
mismatch), recorded with its source, version, species, identifier type, licence identifier and redistributable
flag, and checked against the licence rule (FR-175): KEGG and MSigDB files are never redistributed, and only the
synthetic marked fixtures may appear in the repository or in an artifact. Network locations are refused
(E_RESOURCE_DOWNLOAD) and the run executes under a socket guard, so analysis cannot open a connection.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import socket
from contextlib import contextmanager
from pathlib import Path

from .errors import ProteomicsError

CAPABILITY = "offline_gene_set_libraries"
SYNTHETIC_MARKER = "SYNTHETIC_FIXTURE"
RESTRICTED_LICENCE_PREFIXES = ("KEGG", "MSIGDB")


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["welch"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_restricted(licence_id: str) -> bool:
    return str(licence_id).upper().startswith(RESTRICTED_LICENCE_PREFIXES)


def is_synthetic(path: Path) -> bool:
    """A GMT fixture is synthetic only when every set row carries the marker in its description field."""
    rows = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return bool(rows) and all(SYNTHETIC_MARKER in (line.split("\t") + [""])[1] for line in rows)


def library_record(declaration: dict, base_dir: str | Path) -> dict:
    """Verify one declared library and return its manifest record (FR-174, FR-175). Nothing is parsed before the hash check."""
    raw = str(declaration["path"])
    pointer = f"/methods_extension/pathways/libraries/{declaration['id']}"
    if "://" in raw:
        raise ProteomicsError("E_RESOURCE_DOWNLOAD", "libraries are read from declared local files only; a network location is refused", pointer, exit_code=3)
    path = Path(raw) if Path(raw).is_absolute() else Path(base_dir) / raw
    if not path.is_file():
        raise ProteomicsError("E_RESOURCE_MISSING", f"library file for {declaration['id']!r} is missing", pointer, exit_code=2)
    digest = _sha256(path)
    if digest != declaration["sha256"]:
        raise ProteomicsError("E_RESOURCE_HASH", f"library {declaration['id']!r} differs from its declared SHA-256 (file bytes changed)", pointer, exit_code=5)
    synthetic = is_synthetic(path)
    if is_restricted(declaration["licence_id"]) and declaration.get("redistributable"):
        raise ProteomicsError("E_LIBRARY_REDISTRIBUTION", f"library {declaration['id']!r} has a restricted licence and cannot be declared redistributable", pointer, exit_code=2)
    return {"id": declaration["id"], "sha256": digest, "source": declaration["source"], "version": declaration["version"], "species": declaration["species"],
            "identifier_type": declaration["identifier_type"], "licence_id": declaration["licence_id"],
            "redistributable": bool(declaration.get("redistributable", False)), "synthetic": synthetic, "_path": str(path)}


def assert_artifact_members(members: dict[str, str | Path], libraries: list[dict]) -> list[str]:
    """Refuse any artifact that would contain a restricted or non-synthetic library (FR-175). Returns the accepted member names.

    The check runs before anything is written; a refused artifact produces no file.
    """
    by_sha = {lib["sha256"]: lib for lib in libraries}
    accepted = []
    for name, source in members.items():
        path = Path(source)
        if not path.is_file():
            raise ProteomicsError("E_RESOURCE_MISSING", f"artifact member {name!r} is missing", "/artifact", exit_code=2)
        if path.suffix.lower() != ".gmt":
            accepted.append(name)
            continue
        digest = _sha256(path)
        lib = by_sha.get(digest)
        if lib is not None and is_restricted(lib["licence_id"]):
            raise ProteomicsError("E_LIBRARY_REDISTRIBUTION", f"artifact member {name!r} is a restricted library; nothing is written", "/artifact", exit_code=2)
        if not is_synthetic(path):
            raise ProteomicsError("E_LIBRARY_REDISTRIBUTION", f"artifact member {name!r} is a non-synthetic library; only marked synthetic fixtures may be redistributed",
                                  "/artifact", exit_code=2)
        accepted.append(name)
    return accepted


def build_request(libraries: list[dict], *, base_dir: str | Path, welch_results: str | Path, background: dict, min_size: int, max_size: int, run_id: str,
                  output_temp_dir: str | Path, plan_hash: str | None = None, tie_break: str = "sm15_stable_id", ora_threshold: float = 0.05,
                  stage_id: str = "offline_gene_set_libraries", seed: int = 20261013) -> dict:
    records = [library_record(decl, base_dir) for decl in libraries]
    welch_path = Path(welch_results).resolve()
    inputs = [{"artifact_id": "welch_results", "path": str(welch_path), "sha256": _sha256(welch_path)}]
    for record in records:
        inputs.append({"artifact_id": f"library_{record['id']}", "path": record["_path"], "sha256": record["sha256"]})
    public = [{k: v for k, v in record.items() if not k.startswith("_")} for record in records]
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan_hash, "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None,
            "parameters": {"libraries": public, "background": background, "min_size": int(min_size), "max_size": int(max_size),
                           "tie_break": tie_break, "ora_threshold": float(ora_threshold), "methods": ["ora", "prerank_welch"]},
            "rng": {"seed": int(seed), "kind": "L'Ecuyer-CMRG", "threads": 1}}


@contextmanager
def no_network():
    """Analysis runs under this guard: any attempt to create a socket fails with E_RESOURCE_DOWNLOAD (FR-174)."""
    original = socket.socket

    class _Refused(socket.socket):   # pragma: no cover - exercised through the guard test
        def __init__(self, *args, **kwargs):
            raise ProteomicsError("E_RESOURCE_DOWNLOAD", "network access is not permitted during an offline analysis run", exit_code=3)

    socket.socket = _Refused
    try:
        yield
    finally:
        socket.socket = original


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    with no_network():
        return execute_stage(request, rscript=rscript_executable())
