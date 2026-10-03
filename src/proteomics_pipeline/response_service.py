"""Treatment response, equivalence and independent-score stage (packet R09: FR-081 to FR-090).

Plan-time checks enforce SM19-SM24 eligibility (coherent axes, declared
families, independent direction/score resources); the R handler computes
descriptive axes with covariance, TOST, the conjunction endpoint and score
randomization from the frozen plan plus the completed R05 tables.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from .errors import ProteomicsError
from .provenance import sha256_file

CAPABILITY = "response"


class ResponseRefusal(ProteomicsError):
    def __init__(self, code: str, message: str, pointer: str | None = "/response"):
        super().__init__(code, message, pointer, exit_code=2)


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def _families(config: dict, htype: str, key: str) -> dict[str, dict]:
    out = {}
    for family in config["multiplicity_families"]:
        if family["hypothesis_type"] == htype:
            for item in family.get(key, []):
                out[item] = family
    return out


def plan_checks(config: dict, contrasts: list[dict]) -> None:
    """Reject incoherent or ineligible response requests before any fit (exit 2)."""
    response = config.get("response") or {}
    if not response.get("enabled"):
        return
    by_id = {c["contrast_id"]: c for c in contrasts}
    models = {m["id"]: m for m in config["models"]}
    mode = config.get("response_mode", "descriptive_only")
    for index, axis in enumerate(response.get("axes", [])):
        pointer = f"/response/axes/{index}"
        model = models[axis["model_id"]]
        parts = [by_id[axis[k]] for k in ("disease_contrast", "treatment_contrast", "residual_contrast")]
        if any(p["design_id"] != model["design_id"] for p in parts):
            raise ResponseRefusal("E_AXIS_INCOHERENT", f"axis {axis['id']!r} contrasts must belong to model {model['id']!r}'s design", pointer)
        d, t, r = (p["weights"] for p in parts)
        if any(abs(r[i] - (d[i] + t[i])) > 1e-12 for i in range(len(r))):
            raise ResponseRefusal("E_AXIS_INCOHERENT", f"axis {axis['id']!r}: residual weights must equal disease + treatment weights exactly", pointer)
        if mode in ("equivalence", "formal_rescue"):
            if model["engine"] != "limma":
                raise ResponseRefusal("E_EQUIVALENCE_INELIGIBLE", "equivalence/conjunction are qualified only for limma models", pointer)
            if axis["id"] not in _families(config, "equivalence", "axis_ids"):
                raise ResponseRefusal("E_FAMILY_MISSING", f"axis {axis['id']!r} needs a declared equivalence family", "/multiplicity_families")
        if mode == "formal_rescue" and axis["id"] not in _families(config, "formal_rescue", "axis_ids"):
            raise ResponseRefusal("E_FAMILY_MISSING", f"axis {axis['id']!r} needs a declared formal_rescue family", "/multiplicity_families")
    resources = {r["id"]: r for r in config.get("resources", [])}
    if mode == "formal_rescue":
        resource = resources[response["direction_resource_id"]]
        if resource["kind"] != "independent_direction":
            raise ResponseRefusal("E_RESCUE_DIRECTION_NOT_INDEPENDENT", "formal_rescue needs an independent_direction resource", "/response/direction_resource_id")
    score_families = _families(config, "independent_score_test", "score_ids")
    for index, score in enumerate(response.get("scores", [])):
        if resources[score["resource_id"]]["kind"] != "independent_score":
            raise ResponseRefusal("E_SCORE_INDEPENDENCE_UNVERIFIED", f"score {score['id']!r} needs an independent_score resource", f"/response/scores/{index}/resource_id")
        if len(by_id[score["contrast_id"]]["required_groups"]) != 2:
            raise ResponseRefusal("E_SCORE_EXCHANGEABILITY_UNSUPPORTED", f"score {score['id']!r} contrast must compare exactly two groups", f"/response/scores/{index}/contrast_id")
        if config.get("score_test", "off") != "off" and score["id"] not in score_families:
            raise ResponseRefusal("E_FAMILY_MISSING", f"score {score['id']!r} needs a declared independent_score_test family", "/multiplicity_families")
        if config.get("score_test") == "independent_monte_carlo" and not score.get("draws"):
            raise ResponseRefusal("E_SCORE_EXCHANGEABILITY_UNSUPPORTED", "Monte Carlo score tests need a declared number of draws", f"/response/scores/{index}/draws")


def _stage_output(stage_dir: Path, relative: str) -> dict:
    result = json.loads((stage_dir / "stage-result.json").read_text(encoding="utf-8"))
    output = next(o for o in result["outputs"] if o["relative_path"] == relative)
    path = stage_dir / relative
    if sha256_file(path) != output["sha256"]:
        raise ProteomicsError("E_INTEGRITY", f"{relative} changed after its stage completed", exit_code=5)
    return {"path": str(path.resolve()), "sha256": output["sha256"]}


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, dea_dir: str | Path, config_dir: str | Path,
                  stage_id: str = "response") -> dict:
    from . import inference_service
    inference_service.guard_downstream(config, dea_dir)    # audit 2026-10-02: primary results must come from observed data
    plan_path = Path(plan_path).resolve(); root = plan_path.parent; dea_dir = Path(dea_dir)
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]
    seen = {"plan"}

    def add_plan(artifact_id):
        if artifact_id in seen:
            return
        artifact = next(a for a in plan["artifacts"] if a["artifact_id"] == artifact_id)
        inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]}); seen.add(artifact_id)

    for artifact_id in ("primary_matrix", "primary_observed_mask", "primary_observations"):
        add_plan(artifact_id)
    inputs.append({"artifact_id": "dea_zero_null", **_stage_output(dea_dir, "zero_null.tsv")}); seen.add("dea_zero_null")
    response = config["response"]
    contrasts = {c["contrast_id"]: c for c in plan["contrasts"]}
    models = {m["id"]: m for m in config["models"]}
    designs = {config["design"]["id"]: config["design"], **{d["id"]: d for d in config.get("additional_designs", [])}}
    blocking = {b["design_id"]: b for b in plan.get("blocking", [])}
    eq = _families(config, "equivalence", "axis_ids"); rescue = _families(config, "formal_rescue", "axis_ids")
    axes = []
    for axis in response.get("axes", []):
        model = models[axis["model_id"]]
        add_plan(f"design_{model['design_id']}")
        moderation = f"diagnostics/{model['id']}_moderation.tsv"
        if f"moderation_{model['id']}" not in seen:
            inputs.append({"artifact_id": f"moderation_{model['id']}", **_stage_output(dea_dir, moderation)}); seen.add(f"moderation_{model['id']}")
        weights_artifact = None
        if (model.get("precision_weights") or {"kind": "none"})["kind"] == "observation":
            weights_artifact = f"weights_{model['id']}"; add_plan(weights_artifact)
        design = designs[model["design_id"]]
        axes.append({"id": axis["id"], "model_id": model["id"], "design_id": model["design_id"],
                     "disease_contrast": axis["disease_contrast"], "treatment_contrast": axis["treatment_contrast"], "residual_contrast": axis["residual_contrast"],
                     "disease_weights": contrasts[axis["disease_contrast"]]["weights"], "treatment_weights": contrasts[axis["treatment_contrast"]]["weights"],
                     "residual_weights": contrasts[axis["residual_contrast"]]["weights"],
                     "disease_groups": contrasts[axis["disease_contrast"]]["required_groups"], "treatment_groups": contrasts[axis["treatment_contrast"]]["required_groups"],
                     "weights_artifact": weights_artifact, "blocking_mode": design["blocking"]["mode"], "subject_column": design["blocking"].get("subject_column"),
                     "consensus_correlation": blocking.get(model["design_id"], {}).get("consensus_correlation"),
                     "equivalence_family": eq[axis["id"]]["id"] if axis["id"] in eq else None, "equivalence_adjustment": eq[axis["id"]]["adjustment"] if axis["id"] in eq else "BH",
                     "rescue_family": rescue[axis["id"]]["id"] if axis["id"] in rescue else None, "rescue_adjustment": rescue[axis["id"]]["adjustment"] if axis["id"] in rescue else "BH"})
    resources = {r["id"]: r for r in config.get("resources", [])}
    resource_root = Path(config_dir).resolve()

    def resource_input(resource_id, artifact_id):
        resource = resources[resource_id]
        path = Path(resource["path"]); path = path if path.is_absolute() else resource_root / path
        if not path.is_file() or sha256_file(path) != resource["sha256"]:
            raise ProteomicsError("E_RESOURCE_HASH", f"resource {resource_id!r} is missing or differs from its declared SHA-256", exit_code=2)
        if resource["kind"] == "independent_score":
            from jsonschema import Draft202012Validator
            schema = json.loads((Path(__file__).parent / "schemas" / "independent-score.schema.json").read_text(encoding="utf-8"))
            errors = list(Draft202012Validator(schema).iter_errors(json.loads(path.read_text(encoding="utf-8"))))
            if errors:
                raise ProteomicsError("E_SCORE_MANIFEST", f"score manifest {resource_id!r}: {errors[0].message}", exit_code=2)
        inputs.append({"artifact_id": artifact_id, "path": str(path.resolve()), "sha256": resource["sha256"]})

    if config.get("response_mode") == "formal_rescue":
        resource_input(response["direction_resource_id"], "direction_resource")
    score_families = _families(config, "independent_score_test", "score_ids")
    scores = []
    for score in response.get("scores", []):
        resource_input(score["resource_id"], f"resource_{score['resource_id']}")
        scores.append({"id": score["id"], "model_id": score["model_id"], "contrast_id": score["contrast_id"], "required_groups": contrasts[score["contrast_id"]]["required_groups"],
                       "resource_artifact": f"resource_{score['resource_id']}", "scheme": score["scheme"], "exchangeability_evidence": score["exchangeability_evidence"],
                       "draws": score.get("draws"), "family_id": score_families[score["id"]]["id"] if score["id"] in score_families else None})
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan["plan_hash"], "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None,
            "parameters": {"axes": axes, "scores": scores, "dmin": response.get("ri_minimum_abs_disease_effect", 0.25), "ratio_uncertainty": response.get("ratio_uncertainty", "descriptive_only"),
                           "response_mode": config.get("response_mode", "descriptive_only"), "equivalence_margin": response.get("equivalence_margin"), "alpha": response.get("alpha", 0.05),
                           "disease_margin": response.get("disease_margin", 0), "treatment_margin": response.get("treatment_margin", 0), "score_test": config.get("score_test", "off"),
                           "group_column": config["design"]["group_column"]},
            "rng": {"seed": int(plan["runtime"]["seed"]), "kind": "L'Ecuyer-CMRG", "threads": int(plan["runtime"]["threads"])}}


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())
