"""Design-compatible pathway and enrichment stage (packet R08: FR-071 to FR-080, SM16-SM18).

Plan-time dispatch follows the SM16 design x method table; the R handler runs
CAMERA, ROAST/mroast, exploratory fgsea and background-aware ORA on the R07
finite gene matrix, gene model and verified gene-set memberships.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from .errors import ProteomicsError
from .provenance import sha256_file

CAPABILITY = "pathways"
METHOD_TYPES = {"camera": ["competitive_enrichment"], "roast": ["self_contained_directional", "self_contained_mixed"], "fgsea": ["preranked_gene_set"], "ora": ["ora_up", "ora_down"]}


class PathwayRefusal(ProteomicsError):
    def __init__(self, code: str, message: str, pointer: str = "/pathways"):
        super().__init__(code, message, pointer, exit_code=2)


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["resources", "limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def gene_model(config: dict) -> tuple[dict | None, str]:
    """The linear model behind CAMERA/ROAST: the primary limma model, or an explicitly named limma sensitivity."""
    primary = next(m for m in config["models"] if m["role"] == "primary")
    if primary["engine"] == "limma":
        return primary, "primary"
    sensitivity = (config.get("pathways") or {}).get("linear_sensitivity_model_id")
    model = next((m for m in config["models"] if m["id"] == sensitivity and m["engine"] == "limma"), None)
    return model, "linear_sensitivity"


def plan_checks(config: dict) -> None:
    pathways = config.get("pathways") or {}
    if not pathways.get("enabled"):
        return
    methods = pathways.get("methods", [])
    if not methods:
        raise PathwayRefusal("E_PATHWAY_METHOD", "pathways.enabled requires at least one method", "/pathways/methods")
    blocking = config["design"]["blocking"]["mode"]
    if "camera" in methods and blocking != "none":
        raise PathwayRefusal("E_CAMERA_BLOCKED_DESIGN", "CAMERA is not available for paired/subject-blocked designs in this release; request roast", "/pathways/methods")
    if "roast" in methods and config["runtime"]["execution_profile"] == "production" and pathways.get("nrot", 9999) < 9999:
        raise PathwayRefusal("E_ROTATION_RESOLUTION", "production ROAST needs nrot >= 9999 (smaller values are smoke-test only)", "/pathways/nrot")
    model, role = gene_model(config)
    if any(m in methods for m in ("camera", "roast")) and model is None:
        raise PathwayRefusal("E_PATHWAY_ENGINE_MISMATCH", "matrix pathway tests of a non-limma primary model need an explicitly named limma linear_sensitivity_model_id; they are never same-likelihood confirmation", "/pathways/linear_sensitivity_model_id")
    declared = {f["hypothesis_type"] for f in config["multiplicity_families"]}
    for method in methods:
        missing = [h for h in METHOD_TYPES[method] if h not in declared]
        if missing:
            raise PathwayRefusal("E_FAMILY_MISSING", f"method {method} needs declared pathway families for {missing}", "/multiplicity_families")
    if "ora" in methods and pathways.get("ora_threshold") is None:
        raise PathwayRefusal("E_PATHWAY_METHOD", "ORA needs a declared ora_threshold and ora_foreground", "/pathways/ora_threshold")


def _stage_output(stage_dir: Path, relative: str) -> dict:
    result = json.loads((stage_dir / "stage-result.json").read_text(encoding="utf-8"))
    output = next(o for o in result["outputs"] if o["relative_path"] == relative)
    if sha256_file(stage_dir / relative) != output["sha256"]:
        raise ProteomicsError("E_INTEGRITY", f"{relative} changed after its stage completed", exit_code=5)
    return {"path": str((stage_dir / relative).resolve()), "sha256": output["sha256"]}


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, dea_dir: str | Path, resources_dir: str | Path,
                  stage_id: str = "pathways") -> dict:
    from . import inference_service
    inference_service.guard_downstream(config, dea_dir)    # audit 2026-10-02: primary results must come from observed data
    plan_path = Path(plan_path).resolve(); root = plan_path.parent
    pathways = config["pathways"]
    model, role = gene_model(config)
    primary = next(m for m in config["models"] if m["role"] == "primary")
    design_id = (model or primary)["design_id"]
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]
    for artifact_id in ("primary_observations", f"design_{design_id}"):
        artifact = next(a for a in plan["artifacts"] if a["artifact_id"] == artifact_id)
        inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]})
    for name in ("gene_matrix_finite", "gene_set_membership", "gene_set_summary", "gene_mapping", "gene_zero_null"):
        inputs.append({"artifact_id": name, **_stage_output(Path(resources_dir), f"{name}.tsv")})
    inputs.append({"artifact_id": "dea_zero_null", **_stage_output(Path(dea_dir), "zero_null.tsv")})
    design = next(d for d in [config["design"], *config.get("additional_designs", [])] if d["id"] == design_id)
    blocking = {b["design_id"]: b for b in plan.get("blocking", [])}
    families = [{"id": f["id"], "hypothesis_type": f["hypothesis_type"], "contrast_ids": f.get("contrast_ids", []), "collection_ids": f.get("collection_ids", []), "adjustment": f["adjustment"]}
                for f in config["multiplicity_families"] if f["hypothesis_type"] in {h for hs in METHOD_TYPES.values() for h in hs}]
    contrasts = [c for c in plan["contrasts"] if c["design_id"] == design_id]
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan["plan_hash"], "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None,
            "parameters": {"methods": list(pathways["methods"]), "gene_set_resource_ids": list(pathways["gene_set_resource_ids"]), "design_id": design_id,
                           "blocking_mode": design["blocking"]["mode"], "subject_column": design["blocking"].get("subject_column"),
                           "consensus_correlation": blocking.get(design_id, {}).get("consensus_correlation"), "contrasts": contrasts, "families": families,
                           "nrot": pathways.get("nrot", 9999), "eps": pathways.get("eps", 0), "rank": pathways.get("rank", "native_zero_null"),
                           "ora_foreground": pathways.get("ora_foreground", "family_q"), "ora_threshold": pathways.get("ora_threshold", 0.05),
                           "min_size": pathways.get("min_size", 5), "max_size": pathways.get("max_size", 500), "execution_profile": config["runtime"]["execution_profile"],
                           "primary_model_id": primary["id"], "gene_model_id": f"{(model or primary)['id']}.gene_level", "gene_model_role": role},
            "rng": {"seed": int(plan["runtime"]["seed"]), "kind": "L'Ecuyer-CMRG", "threads": 1}}


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())
