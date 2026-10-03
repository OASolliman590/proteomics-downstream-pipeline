"""Versioned local resources, mapping and gene sets (packet R07: FR-061 to FR-070, SM14/SM15).

Resources are consumed only as hashed local snapshots.  Building a snapshot
is a separate explicit preparation action (``proteomics resources
prepare``) that copies declared local source files, records their
provenance and hashes, and never runs during analysis.  Analysis verifies
every snapshot (file, SHA-256, taxonomy, kind, nested file hashes) before
the R mapping stage builds the label-independent gene matrix.
"""
from __future__ import annotations

import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .errors import ProteomicsError
from .provenance import canonical_json_bytes, sha256_file

CAPABILITY = "resources"
SNAPSHOT_KINDS = ("mapping", "gene_sets", "orthology")
SINGLE_FILE_KINDS = ("independent_score", "independent_direction")
REQUIRED_COLUMNS = {"mapping": ["source_id", "id_type", "gene_id", "gene_symbol", "taxonomy_id", "status"],
                    "gene_sets": ["set_id", "set_name", "gene_id"],
                    "orthology": ["source_gene_id", "source_taxonomy_id", "target_gene_id", "target_taxonomy_id", "evidence", "ambiguous"]}


class ResourceError(ProteomicsError):
    def __init__(self, code: str, message: str, pointer: str | None = None):
        super().__init__(code, message, pointer, exit_code=2)


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["design"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma", "statmod"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


# --------------------------------------------------------------------------- preparation (explicit action)
def prepare(manifest_path: str | Path, output_dir: str | Path) -> dict:
    """Build a local snapshot from a preparation manifest naming local source files."""
    manifest_path = Path(manifest_path).resolve()
    spec = json.loads(manifest_path.read_text(encoding="utf-8"))
    required = {"resource_id", "kind", "version", "source", "source_taxonomy_id", "target_taxonomy_id", "id_type", "terms", "files"}
    missing = sorted(required - set(spec))
    if missing:
        raise ResourceError("E_RESOURCE_MANIFEST", f"preparation manifest lacks {missing}")
    if spec["kind"] not in SNAPSHOT_KINDS:
        raise ResourceError("E_RESOURCE_MANIFEST", f"kind {spec['kind']!r} is not a snapshot kind {SNAPSHOT_KINDS}")
    if spec["source_taxonomy_id"] != spec["target_taxonomy_id"] and not spec.get("projection"):
        raise ResourceError("E_RESOURCE_TAXONOMY", "a cross-species resource needs explicit orthology projection evidence")
    output_dir = Path(output_dir)
    if output_dir.exists():
        raise ResourceError("E_PATH_COLLISION", f"snapshot output already exists: {output_dir.name}")
    files = []
    staging = output_dir.parent / f".{output_dir.name}.prepare-{os.getpid()}"
    staging.mkdir(parents=True)
    try:
        for item in spec["files"]:
            source = str(item["source"])
            if "://" in source:
                raise ResourceError("E_RESOURCE_SOURCE_UNSUPPORTED", "only local source files are supported; download outside the pipeline and record its provenance")
            path = Path(source) if Path(source).is_absolute() else manifest_path.parent / source
            if not path.is_file():
                raise ResourceError("E_RESOURCE_MISSING", f"preparation source file is missing: {path.name}")
            _check_columns(path, spec["kind"])
            target = staging / item["name"]
            shutil.copyfile(path, target)
            files.append({"name": item["name"], "sha256": sha256_file(target), "source_sha256": sha256_file(path), "role": item.get("role", spec["kind"])})
        snapshot = {"schema_version": "1.2.0", "resource_id": spec["resource_id"], "kind": spec["kind"], "version": spec["version"], "source": spec["source"],
                    "source_taxonomy_id": spec["source_taxonomy_id"], "target_taxonomy_id": spec["target_taxonomy_id"], "id_type": spec["id_type"], "terms": spec["terms"],
                    "collection_label": "ortholog_projected" if spec["source_taxonomy_id"] != spec["target_taxonomy_id"] else "native",
                    "projection": spec.get("projection"), "files": files,
                    "prepared_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"), "preparation": "explicit local preparation; no network access"}
        (staging / "manifest.json").write_bytes(canonical_json_bytes(snapshot) + b"\n")
        os.replace(staging, output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return {"snapshot": str(output_dir), "manifest_sha256": sha256_file(output_dir / "manifest.json"), "resource_id": spec["resource_id"], "files": len(files)}


def _check_columns(path: Path, kind: str) -> None:
    header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    missing = [c for c in REQUIRED_COLUMNS[kind] if c not in header]
    if missing:
        raise ResourceError("E_RESOURCE_MANIFEST", f"{path.name} lacks columns {missing} for kind {kind}")


# --------------------------------------------------------------------------- verification (analysis time)
def verify_resources(config: dict, config_dir: str | Path) -> list[dict]:
    """Verify every declared resource before analysis; nothing is fetched or refreshed."""
    config_dir = Path(config_dir)
    organism = config["organism"]["taxonomy_id"]
    verified = []
    for index, resource in enumerate(config.get("resources", [])):
        pointer = f"/resources/{index}"
        path = Path(resource["path"]) if Path(resource["path"]).is_absolute() else config_dir / resource["path"]
        if set(resource["sha256"]) == {"0"}:
            raise ResourceError("E_RESOURCE_HASH", f"resource {resource['id']!r} declares a placeholder digest", f"{pointer}/sha256")
        if not path.is_file():
            raise ResourceError("E_RESOURCE_MISSING", f"resource {resource['id']!r} has no local snapshot; run 'proteomics resources prepare' explicitly (analysis never downloads)", f"{pointer}/path")
        if sha256_file(path) != resource["sha256"]:
            raise ResourceError("E_RESOURCE_HASH", f"resource {resource['id']!r} bytes differ from the declared SHA-256", f"{pointer}/sha256")
        entry = {"resource_id": resource["id"], "kind": resource["kind"], "path": str(path.resolve()), "sha256": resource["sha256"], "files": []}
        if resource["kind"] in SNAPSHOT_KINDS:
            snapshot = json.loads(path.read_text(encoding="utf-8"))
            for key in ("kind", "version", "source", "source_taxonomy_id", "target_taxonomy_id", "terms"):
                if snapshot.get(key) != resource[key]:
                    raise ResourceError("E_RESOURCE_MANIFEST", f"resource {resource['id']!r} declares {key}={resource[key]!r} but the snapshot records {snapshot.get(key)!r}", f"{pointer}/{key}")
            if resource["kind"] in ("mapping", "gene_sets") and snapshot["target_taxonomy_id"] != organism:
                raise ResourceError("E_RESOURCE_TAXONOMY", f"resource {resource['id']!r} targets taxonomy {snapshot['target_taxonomy_id']} but the study organism is {organism}; cross-species joins are refused", f"{pointer}/target_taxonomy_id")
            if snapshot["source_taxonomy_id"] != snapshot["target_taxonomy_id"] and not snapshot.get("projection"):
                raise ResourceError("E_RESOURCE_TAXONOMY", f"resource {resource['id']!r} crosses species without orthology projection evidence", pointer)
            for item in snapshot["files"]:
                file_path = path.parent / item["name"]
                if not file_path.is_file():
                    raise ResourceError("E_RESOURCE_MISSING", f"snapshot file {item['name']} of {resource['id']!r} is missing", pointer)
                if sha256_file(file_path) != item["sha256"]:
                    raise ResourceError("E_RESOURCE_HASH", f"snapshot file {item['name']} of {resource['id']!r} differs from its manifest hash", pointer)
                entry["files"].append({"name": item["name"], "path": str(file_path.resolve()), "sha256": item["sha256"]})
            entry["collection_label"] = snapshot.get("collection_label", "native")
            entry["projection"] = snapshot.get("projection")
        verified.append(entry)
    return verified


# --------------------------------------------------------------------------- stage request
def build_request(plan: dict, *, plan_path: str | Path, config: dict, config_dir: str | Path, run_id: str, output_temp_dir: str | Path, stage_id: str = "resources") -> dict:
    from .preprocessing_service import coverage_rules  # noqa: F401  (documented dependency: R03 coverage semantics)
    plan_path = Path(plan_path).resolve(); root = plan_path.parent
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]
    for artifact_id in ("primary_matrix", "primary_observed_mask", "primary_observations", "primary_features"):
        artifact = next(a for a in plan["artifacts"] if a["artifact_id"] == artifact_id)
        inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]})
    verified = verify_resources(config, config_dir)
    resources = []
    for entry in verified:
        if entry["kind"] not in SNAPSHOT_KINDS:
            continue
        inputs.append({"artifact_id": f"resource_{entry['resource_id']}", "path": entry["path"], "sha256": entry["sha256"]})
        for item in entry["files"]:
            inputs.append({"artifact_id": f"resfile_{entry['resource_id']}__{Path(item['name']).stem}", "path": item["path"], "sha256": item["sha256"]})
        resources.append({"resource_id": entry["resource_id"], "kind": entry["kind"], "files": [Path(i["name"]).stem for i in entry["files"]],
                          "collection_label": entry.get("collection_label"), "projection": entry.get("projection")})
    pathways = config.get("pathways") or {}
    primary = next(m for m in config["models"] if m["role"] == "primary")
    if primary["engine"] != "limma":   # gene-level linear model only via an explicitly named limma sensitivity (SM16)
        sensitivity = pathways.get("linear_sensitivity_model_id")
        primary = next((m for m in config["models"] if m["id"] == sensitivity and m["engine"] == "limma"), primary)
    design_artifact = f"design_{primary['design_id']}"
    artifact = next(a for a in plan["artifacts"] if a["artifact_id"] == design_artifact)
    inputs.append({"artifact_id": design_artifact, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]})
    blocking = {b["design_id"]: b for b in plan.get("blocking", [])}
    design = next(d for d in [config["design"], *config.get("additional_designs", [])] if d["id"] == primary["design_id"])
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan["plan_hash"], "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None,
            "parameters": {"organism_taxonomy_id": config["organism"]["taxonomy_id"], "resources": resources,
                           "mapping_resource_id": pathways.get("mapping_resource_id"), "gene_set_resource_ids": list(pathways.get("gene_set_resource_ids", [])),
                           "orthology_resource_ids": [r["id"] for r in config.get("resources", []) if r["kind"] == "orthology"],
                           "min_size": pathways.get("min_size", 5), "max_size": pathways.get("max_size", 500),
                           "gene_aggregation_sensitivity": pathways.get("gene_aggregation_sensitivity", "off"),
                           "model": {"model_id": primary["id"], "design_id": primary["design_id"], "trend": primary.get("trend", True), "robust": primary.get("robust", True),
                                     "blocking_mode": design["blocking"]["mode"], "subject_column": design["blocking"].get("subject_column"),
                                     "consensus_correlation": blocking.get(primary["design_id"], {}).get("consensus_correlation")},
                           "contrasts": [c for c in plan["contrasts"] if c["design_id"] == primary["design_id"]], "ci_level": config["report"]["ci_level"]},
            "rng": {"seed": int(plan["runtime"]["seed"]), "kind": "L'Ecuyer-CMRG", "threads": int(plan["runtime"]["threads"])}}


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())
