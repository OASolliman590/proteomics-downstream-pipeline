"""Optional Welch t-test engine stage (packet R15a: FR-168 to FR-173, SM42/SM43; ADR 0011).

The Welch engine runs only when declared under ``methods_extension.welch`` (FR-169). This module
builds the stage request strictly from the frozen plan: every input is referenced by the plan with
its SHA-256, and the R handler refuses any request whose plan hash or input bytes differ. The
contrast orientation (numerator and denominator) and the adaptation records (covariates that a
two-sample Welch test cannot honour, ADR 0010) are decided here from the declared design; the R
handler decides the paired variant from the subject structure and records every decision.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

from .design_service import encode
from .errors import ProteomicsError

CAPABILITY = "welch"


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["design"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def _artifact(plan: dict, artifact_id: str) -> dict | None:
    matches = [a for a in plan["artifacts"] if a["artifact_id"] == artifact_id]
    if len(matches) > 1:
        raise ValueError(f"plan artifact id {artifact_id} is ambiguous")
    return matches[0] if matches else None


def welch_declaration(config: dict) -> dict:
    return (config.get("methods_extension") or {}).get("welch") or {}


def orient_contrast(contrast: dict) -> tuple[str, str] | None:
    """Return (numerator, denominator) for a pairwise group contrast, or None when it is not pairwise.

    The coefficient of a group level is ``group.<encoded level>``; a pairwise contrast carries +w on
    one level and -w on the other, and zero on every other coefficient (covariates are not compared).
    """
    groups = list(contrast["required_groups"])
    if len(groups) != 2:
        return None
    names = {g: f"group.{encode(g)}" for g in groups}
    weights = {key: float(value) for key, value in contrast["weights"].items() if float(value) != 0.0}
    if set(weights) != set(names.values()):
        return None
    plus = [g for g in groups if weights[names[g]] > 0]
    minus = [g for g in groups if weights[names[g]] < 0]
    if len(plus) != 1 or len(minus) != 1 or abs(weights[names[plus[0]]] + weights[names[minus[0]]]) > 1e-12:
        return None
    return plus[0], minus[0]


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, stage_id: str = "welch") -> dict:
    from .provenance import sha256_file
    plan_path = Path(plan_path).resolve()
    root = plan_path.parent
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]

    def add(artifact_id: str) -> None:
        artifact = _artifact(plan, artifact_id)
        if artifact is None:
            raise ValueError(f"frozen plan lacks required artifact {artifact_id}")
        if artifact_id not in {i["artifact_id"] for i in inputs}:
            inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]})

    for artifact_id in ("primary_matrix", "primary_observed_mask", "primary_observations", "estimability"):
        add(artifact_id)
    declaration = welch_declaration(config)
    eligible = {row["model_id"]: row for row in plan["engine_eligibility"]}
    welch_models = [m for m in config["models"] if m["engine"] == CAPABILITY and eligible[m["id"]]["eligible"]]
    if not welch_models:
        raise ProteomicsError("E_WELCH_NOT_DECLARED", "no eligible welch model was requested", "/methods_extension/welch", exit_code=2)
    welch_ids = {m["id"] for m in welch_models}
    refused: list[dict] = []
    designs = {d["id"]: d for d in [config["design"], *config.get("additional_designs", [])]}
    # The Welch family is a declared multiplicity family of the frozen plan (E_FAMILY_MISSING otherwise); a family that mixes
    # engines would pool Welch and limma P values, so it is refused.
    declared_family = declaration.get("family")
    families = []
    for family in plan.get("families", []):
        if family.get("hypothesis_type") != "protein_zero_null":
            continue
        members = family.get("members", [])
        if not any(m["model_id"] in welch_ids for m in members):
            continue
        if declared_family is not None and family["family_id"] != declared_family:
            continue
        if any(m["model_id"] not in welch_ids for m in members):
            raise ProteomicsError("E_WELCH_FAMILY_MIXED", f"family {family['family_id']!r} mixes welch and other engines; a Welch family must hold only welch endpoints",
                                  "/multiplicity_families", exit_code=2)
        families.append(family)
    if not families:
        raise ProteomicsError("E_WELCH_NOT_DECLARED", "no declared multiplicity family holds the welch endpoints", "/multiplicity_families", exit_code=2)
    models = []
    for model in welch_models:
        design = designs[model["design_id"]]
        blocking = design.get("blocking") or {"mode": "none"}
        blocked = blocking.get("mode", "none") != "none"
        member_contrasts = [m["contrast_id"] for f in families for m in f["members"] if m["model_id"] == model["id"]]
        contrasts = []
        for contrast in config["contrasts"]:
            if contrast["design_id"] != model["design_id"] or contrast["id"] not in member_contrasts:
                continue
            orientation = orient_contrast(contrast)
            if orientation is None:
                refused.append({"model_id": model["id"], "contrast_id": contrast["id"], "item": f"contrast:{contrast['id']}",
                                "requested": "pairwise_group_contrast", "reason": "E_WELCH_CONTRAST_NOT_PAIRWISE"})
                continue
            contrasts.append({"contrast_id": contrast["id"], "numerator": orientation[0], "denominator": orientation[1]})
        if not contrasts and not any(r["model_id"] == model["id"] for r in refused):
            raise ProteomicsError("E_WELCH_CONTRAST_NOT_PAIRWISE", f"welch model {model['id']!r} has no pairwise group contrast to test", "/contrasts", exit_code=2)
        adaptations = list(eligible[model["id"]].get("adaptations", []))
        models.append({"model_id": model["id"], "design_id": model["design_id"], "role": model["role"], "execution_requirement": model["execution_requirement"],
                       "paired": str(declaration.get("paired", "auto")).lower(), "blocked": blocked,
                       "subject_column": blocking.get("subject_column") if blocked else None,
                       "contrasts": contrasts, "adaptations": adaptations,
                       "refused": [r for r in refused if r["model_id"] == model["id"]]})
    family_payload = [{"family_id": f["family_id"], "role": f["role"], "hypothesis_type": f["hypothesis_type"], "adjustment": f["adjustment"],
                       "dependence_assumption": f["dependence_assumption"], "q_cutoff": f["q_cutoff"],
                       "members": [{"model_id": m["model_id"], "contrast_id": m["contrast_id"]} for m in f["members"]]} for f in families]
    if declaration.get("contrasts"):
        planned = {(m["model_id"], m["contrast_id"]) for f in family_payload for m in f["members"]}
        for model in welch_models:
            for contrast_id in declaration["contrasts"]:
                if (model["id"], contrast_id) not in planned:
                    raise ProteomicsError("E_WELCH_NOT_DECLARED", f"declared welch contrast {contrast_id!r} is not an endpoint of the declared family", "/methods_extension/welch/contrasts", exit_code=2)
    paper = declaration.get("paper_compatible") or {}
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan["plan_hash"], "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None,
            "parameters": {"group_column": config["design"]["group_column"],
                           "ci_level": float(declaration.get("ci_level", config["report"]["ci_level"])),
                           "paper_compatible": {"enabled": bool(paper.get("enabled", False)), "raw_p_cutoff": float(paper.get("raw_p_cutoff", 0.05))},
                           "models": models, "families": family_payload},
            "rng": {"seed": int(plan["runtime"]["seed"]), "kind": "L'Ecuyer-CMRG", "threads": int(plan["runtime"]["threads"])}}


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())
