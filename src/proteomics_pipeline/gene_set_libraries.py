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


def declared(config: dict) -> list[dict]:
    """The methods_extension.pathways libraries of a configuration, or [] when the block is absent or disabled."""
    block = (config.get("methods_extension") or {}).get("pathways") or {}
    return list(block.get("libraries", [])) if block.get("enabled") else []


def verify_declared(config: dict, base_dir: str | Path) -> list[dict]:
    """Plan-time check: every declared library is hashed and licence-checked before the plan is frozen (E_RESOURCE_HASH and others)."""
    return [library_record(decl, base_dir) for decl in declared(config)]


def gene_table(dea_zero_null: Path, features_tsv: Path, out_path: Path) -> int:
    """The measured, mapped input of the stage: one row per tested feature with its representative gene (SM15: the smallest stable
    feature ID is the last tie key); features without a gene are left out and are not part of the universe."""
    import csv
    import json as _json
    genes = {}
    with features_tsv.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                ids = _json.loads(row.get("gene_ids") or "[]")
            except ValueError:
                ids = []
            genes[row["feature_id"]] = ids[0] if ids else ""
    columns = ["feature_id", "gene_id", "contrast_id", "effect", "statistic", "q_value", "eligibility"]
    rows = []
    with dea_zero_null.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["eligibility"] != "tested":
                continue
            gene = genes.get(row["feature_id"], "")
            rows.append({"feature_id": row["feature_id"], "gene_id": gene, "contrast_id": row["contrast_id"], "effect": row["effect"],
                         "statistic": row["statistic"], "q_value": row["q_value"], "eligibility": "tested"})
    rows.sort(key=lambda r: (r["contrast_id"], r["gene_id"], r["feature_id"]))
    keep, seen = [], set()
    for row in rows:                       # one feature per gene and contrast; genes without a mapping are excluded
        key = (row["contrast_id"], row["gene_id"])
        if not row["gene_id"] or key in seen:
            continue
        seen.add(key); keep.append(row)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(keep)
    return len({r["gene_id"] for r in keep})


def run_request(config: dict, *, config_dir: str | Path, root: str | Path, run_id: str, output_temp_dir: str | Path, plan_hash: str | None) -> dict:
    """The stage request of a run: the declared libraries over the primary differential table (dea/zero_null.tsv) of the run."""
    root = Path(root)
    block = (config.get("methods_extension") or {}).get("pathways") or {}
    input_dir = root / "methods_extension" / "pathways_input"
    input_path = input_dir / "pathway_input.tsv"
    universe = gene_table(root / "dea" / "zero_null.tsv", root / "preprocessing" / "primary" / "features.tsv", input_path)
    return build_request(declared(config), base_dir=config_dir, welch_results=input_path,
                         background={"kind": "measured", "size": universe}, min_size=int(block.get("min_size", 1)),
                         max_size=int(block.get("max_size", 500)), run_id=run_id, output_temp_dir=output_temp_dir, plan_hash=plan_hash,
                         tie_break="sm15_stable_id", ora_threshold=float(block.get("ora_threshold", 0.05)))


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
