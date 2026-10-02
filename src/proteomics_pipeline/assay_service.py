"""Assay-qualified DEqMS and proDA stage (packet R06: FR-051 to FR-060, SM10).

Alternative engines run only as separately declared models with their own
families.  Scientific eligibility is decided from the frozen plan without
consulting installation; a missing package after an eligible request is
NOT_RUN (E_ENGINE_NOT_AVAILABLE), and nothing ever falls back to limma.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from .count_evidence import check_count_table
from .errors import ProteomicsError
from .provenance import sha256_file

CAPABILITY = "assay_engines"


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def assay_models(config: dict) -> list[dict]:
    return [m for m in config["models"] if m["engine"] in ("deqms", "proda")]


def plan_checks(config: dict, config_dir: Path) -> None:
    for model in assay_models(config):
        if model["engine"] == "deqms":
            path = Path(model["count_evidence"]); path = path if path.is_absolute() else Path(config_dir) / path
            check_count_table(path, model_id=model["id"], aggregation=model["count_aggregation"])


def _stage_output(stage_dir: Path, relative: str) -> dict:
    result = json.loads((stage_dir / "stage-result.json").read_text(encoding="utf-8"))
    output = next(o for o in result["outputs"] if o["relative_path"] == relative)
    if sha256_file(stage_dir / relative) != output["sha256"]:
        raise ProteomicsError("E_INTEGRITY", f"{relative} changed after its stage completed", exit_code=5)
    return {"path": str((stage_dir / relative).resolve()), "sha256": output["sha256"]}


def build_request(plan: dict, *, plan_path: str | Path, config: dict, config_dir: str | Path, run_id: str, output_temp_dir: str | Path, dea_dir: str | Path,
                  stage_id: str = "assay_engines") -> dict:
    plan_path = Path(plan_path).resolve(); root = plan_path.parent; config_dir = Path(config_dir)
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]
    seen = {"plan"}

    def add_plan(artifact_id):
        if artifact_id in seen:
            return
        artifact = next(a for a in plan["artifacts"] if a["artifact_id"] == artifact_id)
        inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]}); seen.add(artifact_id)

    for artifact_id in ("primary_matrix", "primary_observed_mask", "primary_observations", "estimability"):
        add_plan(artifact_id)
    inputs.append({"artifact_id": "dea_zero_null", **_stage_output(Path(dea_dir), "zero_null.tsv")})
    eligible = {row["model_id"]: row for row in plan["engine_eligibility"]}
    frozen = {m["model_id"]: m for m in plan["models"]}
    designs = {config["design"]["id"]: config["design"], **{d["id"]: d for d in config.get("additional_designs", [])}}
    models = []
    for model in assay_models(config):
        if not eligible[model["id"]]["eligible"]:
            continue
        add_plan(f"design_{model['design_id']}")
        entry = {"model_id": model["id"], "design_id": model["design_id"], "engine": model["engine"], "execution_requirement": model["execution_requirement"],
                 "hypotheses": [model["hypothesis"]] + list(model.get("additional_hypotheses", [])),
                 "blocking_mode": designs[model["design_id"]]["blocking"]["mode"], "weighted": (model.get("precision_weights") or {"kind": "none"})["kind"] != "none"}
        if model["engine"] == "deqms":
            path = Path(model["count_evidence"]); path = path if path.is_absolute() else config_dir / path
            if sha256_file(path) != frozen[model["id"]]["count_evidence_sha256"]:
                raise ProteomicsError("E_PLAN_CHANGED", f"count evidence for {model['id']!r} changed after the plan was frozen", exit_code=5)
            inputs.append({"artifact_id": f"counts_{model['id']}", "path": str(path.resolve()), "sha256": frozen[model["id"]]["count_evidence_sha256"]})
            entry.update(count_aggregation=model["count_aggregation"], count_zero_policy=model["count_zero_policy"])
        models.append(entry)
    primary = next(m for m in config["models"] if m["role"] == "primary")
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan["plan_hash"], "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None,
            "parameters": {"models": models, "contrasts": plan["contrasts"], "families": [f for f in plan["families"] if any(m["model_id"] in {x["model_id"] for x in models} for m in f["members"])],
                           "primary_model_id": primary["id"], "assay": config["assay"], "prior_imputation": config["input"]["prior_imputation"],
                           "group_column": config["design"]["group_column"]},
            "rng": {"seed": int(plan["runtime"]["seed"]), "kind": "L'Ecuyer-CMRG", "threads": int(plan["runtime"]["threads"])}}


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())


def model_status(stage_dir: str | Path) -> list[dict]:
    path = Path(stage_dir) / "model_status.json"
    return json.loads(path.read_text(encoding="utf-8"))["models"] if path.is_file() else []
