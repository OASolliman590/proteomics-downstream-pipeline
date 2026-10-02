"""Core limma inference and multiplicity stage (packet R05: FR-041 to FR-050).

Builds the limma stage request strictly from a verified frozen plan: every
input artifact is referenced by the plan with its SHA-256, and the R
handler refuses any request whose plan hash or input bytes differ.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

CAPABILITY = "limma"


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["design"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma", "statmod"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def _artifact(plan: dict, artifact_id: str) -> dict | None:
    matches = [a for a in plan["artifacts"] if a["artifact_id"] == artifact_id]
    if len(matches) > 1:
        raise ValueError(f"plan artifact id {artifact_id} is ambiguous")
    return matches[0] if matches else None


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, stage_id: str = "limma") -> dict:
    from .provenance import sha256_file
    plan_path = Path(plan_path).resolve()
    root = plan_path.parent
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]

    def add(artifact_id: str, required: bool = True):
        artifact = _artifact(plan, artifact_id)
        if artifact is None:
            if required:
                raise ValueError(f"frozen plan lacks required artifact {artifact_id}")
            return False
        if artifact_id not in {i["artifact_id"] for i in inputs}:
            inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]})
        return True

    for artifact_id in ("primary_matrix", "primary_observed_mask", "primary_observations", "estimability"):
        add(artifact_id)
    sensitivity_for_model = {s["model_id"]: s["id"] for s in config["preprocessing"].get("sensitivities", [])}
    blocking = {b["design_id"]: b for b in plan.get("blocking", [])}
    eligible = {row["model_id"]: row for row in plan["engine_eligibility"]}
    models = []
    for model in config["models"]:
        if model["engine"] != "limma" or not eligible[model["id"]]["eligible"]:
            continue
        design = next(d for d in [config["design"], *config.get("additional_designs", [])] if d["id"] == model["design_id"])
        add(f"design_{model['design_id']}")
        matrix_artifact = "primary_matrix"
        if model["id"] in sensitivity_for_model:
            matrix_artifact = f"sensitivity_{sensitivity_for_model[model['id']]}_matrix"
            add(matrix_artifact)
        weights = model.get("precision_weights") or {"kind": "none"}
        weights_artifact = None
        if weights["kind"] != "none":
            weights_artifact = f"weights_{model['id']}"
            add(weights_artifact)
        block = blocking.get(model["design_id"], {})
        models.append({"model_id": model["id"], "design_id": model["design_id"], "role": model["role"], "execution_requirement": model["execution_requirement"],
                       "hypotheses": [model["hypothesis"]] + list(model.get("additional_hypotheses", [])), "effect_threshold": model.get("effect_threshold"),
                       "trend": model.get("trend", True), "robust": model.get("robust", True), "matrix_artifact": matrix_artifact,
                       "uses_observed_mask": matrix_artifact == "primary_matrix", "weights_artifact": weights_artifact, "weights_kind": weights["kind"],
                       "blocking": {"mode": design["blocking"]["mode"], "subject_column": design["blocking"].get("subject_column"),
                                    "consensus_correlation": block.get("consensus_correlation")}})
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan["plan_hash"], "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None,
            "parameters": {"models": models, "contrasts": plan["contrasts"], "families": plan["families"], "group_column": config["design"]["group_column"],
                           "omnibus_tests": [{"id": o["id"], "model_id": o["model_id"], "coefficient_names": list(o["coefficient_names"])} for o in plan.get("omnibus_tests", [])],
                           "ci_level": config["report"]["ci_level"], "influence": {"enabled": True}},
            "rng": {"seed": int(plan["runtime"]["seed"]), "kind": "L'Ecuyer-CMRG", "threads": int(plan["runtime"]["threads"])}}


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())


def model_status(stage_dir: str | Path) -> list[dict]:
    path = Path(stage_dir) / "model_status.json"
    return json.loads(path.read_text(encoding="utf-8"))["models"] if path.is_file() else []
