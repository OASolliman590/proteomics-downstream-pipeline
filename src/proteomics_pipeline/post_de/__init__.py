"""Post-differential analysis (Phase 4, ADR 0009, slice 015, SM31-SM41): shared declarations and request plumbing (packet R14a).

Every post-DE module is declared in ``post_de.<module>`` of the configuration, has a published eligibility
rule decided by the planner before any computation (SM41), consumes the frozen plan and completed R05
outputs read-only with verified hashes, and writes ``<run>/post_de/<module>/`` with its own
``eligibility.json``.  Nothing here is study-specific.
"""
from __future__ import annotations

import importlib
import json
from pathlib import Path

from ..errors import ProteomicsError
from ..provenance import sha256_file

MODULES = ("sets", "sensitivity", "association", "biomarker", "networks")
CAPABILITIES = {"sets": "post_de_sets", "sensitivity": "post_de_sensitivity", "association": "post_de_association",
                "biomarker": "post_de_biomarker", "networks": "post_de_networks"}
ELIGIBILITY_CAPABILITY = "post_de_eligibility"
CLAIM_LABELS = ("descriptive", "exploratory_raw_p", "in_sample", "cross_validated_nested", "fixed_panel_cv", "independently_validated", "module_level")


class PostDeRefusal(ProteomicsError):
    """A post-DE declaration that cannot be honoured (typed; exit 2)."""

    def __init__(self, code: str, message: str, pointer: str | None = "/post_de"):
        super().__init__(code, message, pointer, exit_code=2)


def block(config: dict, module: str | None = None) -> dict:
    top = config.get("post_de") or {}
    if module is None:
        return top
    return top.get(module) or {}


def enabled(config: dict, module: str) -> bool:
    return bool(block(config).get("enabled")) and bool(block(config, module).get("enabled"))


def required(config: dict, module: str) -> bool:
    sub = block(config, module)
    return (sub.get("execution_requirement") or block(config).get("execution_requirement") or "optional") == "required"


def requested_modules(config: dict) -> list[str]:
    return [m for m in MODULES if enabled(config, m)]


def module_impl(module: str):
    """The module's implementation file, or None when it is not part of this release (then NOT_RUN, never INAPPLICABLE)."""
    import importlib.util
    if importlib.util.find_spec(f"{__name__}.{module}") is None:
        return None
    return importlib.import_module(f"{__name__}.{module}")   # an import error inside a present module is a defect, not absence


# --------------------------------------------------------------------------- configuration and plan-time eligibility
def precheck(raw: dict) -> None:
    """Typed declaration refusals that must win over generic schema messages (grammar, leakage, forbidden scopes)."""
    top = raw.get("post_de") or {}
    if not isinstance(top, dict) or not top.get("enabled"):
        return
    runtime = raw.get("runtime") or {}
    if runtime.get("phase") != 4:
        raise PostDeRefusal("E_PHASE_CAPABILITY", "post-differential analysis is a Phase 4 capability; runtime.phase must be 4", "/runtime/phase")
    for module in MODULES:
        sub = top.get(module) or {}
        impl = module_impl(module)
        # D-59 (Maintainer direction 2026-10-05): a declaration refusal of an *optional* module never holds up the run; it is
        # recorded as that module's INAPPLICABLE decision by plan_checks. Only a required module rejects the configuration here.
        if isinstance(sub, dict) and sub.get("enabled") and impl is not None and hasattr(impl, "precheck") and required(raw, module):
            impl.precheck(sub, raw)


def plan_checks(config: dict, context: dict) -> dict[str, dict]:
    """Per requested module: {"state": ELIGIBLE|INAPPLICABLE, "reason_code", "reason", "subanalyses": [...], "parameters": {...}}.
    ``context`` holds observations, feature ids and the compiled design request of the frozen plan."""
    decisions = {}
    for module in requested_modules(config):
        impl = module_impl(module)
        if impl is None or not hasattr(impl, "plan_checks"):
            decisions[module] = {"state": "ELIGIBLE", "reason_code": None, "reason": None, "subanalyses": [], "implemented": False}
            continue
        try:
            if hasattr(impl, "precheck"):
                impl.precheck(block(config, module), config)
            decision = impl.plan_checks(block(config, module), config, context)
        except PostDeRefusal as error:
            if required(config, module):
                raise
            decision = {"state": "INAPPLICABLE", "reason_code": error.code, "reason": error.message, "pointer": error.pointer, "declaration_refused": True,
                        "subanalyses": [], "resolved": {}}
        decision.setdefault("subanalyses", []); decision.setdefault("adaptations", []); decision["implemented"] = True
        decisions[module] = decision
    return decisions


# --------------------------------------------------------------------------- stage requests
def _plan_artifact(plan: dict, artifact_id: str) -> dict | None:
    return next((a for a in plan["artifacts"] if a["artifact_id"] == artifact_id), None)


def stage_output(stage_dir: str | Path, relative: str) -> tuple[Path, str] | None:
    """A verified output of a completed stage: (path, sha256) if it is in the stage manifest and unchanged."""
    stage_dir = Path(stage_dir)
    result_path = stage_dir / "stage-result.json"
    if not result_path.is_file():
        return None
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result.get("state") != "COMPLETED":
        return None
    output = next((o for o in result["outputs"] if o["relative_path"] == relative), None)
    if output is None:
        return None
    path = stage_dir / relative
    if sha256_file(path) != output["sha256"]:
        raise ProteomicsError("E_INTEGRITY", f"{stage_dir.name}/{relative} changed after its stage completed", exit_code=5)
    return path, output["sha256"]


def base_inputs(plan: dict, plan_path: str | Path, *, design_ids: list[str] | None = None) -> list[dict]:
    plan_path = Path(plan_path).resolve(); root = plan_path.parent
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]
    wanted = ["primary_matrix", "primary_observed_mask", "primary_observations", "primary_features", "estimability"]
    wanted += [f"design_{d}" for d in (design_ids or [])]
    for artifact_id in wanted:
        artifact = _plan_artifact(plan, artifact_id)
        if artifact is not None:
            inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]})
    return inputs


def dea_inputs(config: dict, dea_dir: str | Path, names=("zero_null.tsv", "treat.tsv", "families.tsv", "model_status.json")) -> list[dict]:
    from .. import inference_service
    inference_service.guard_downstream(config, dea_dir)   # D-32: only observed-data primary results feed downstream analyses
    inputs = []
    for name in names:
        found = stage_output(dea_dir, name)
        if found is not None:
            inputs.append({"artifact_id": f"stage__dea_{Path(name).stem}", "path": str(found[0].resolve()), "sha256": found[1]})
    return inputs


def request(capability: str, stage_id: str, plan: dict, *, run_id: str, output_temp_dir: str | Path, inputs: list[dict], parameters: dict, seed: int | None = None) -> dict:
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": capability, "plan_hash": plan["plan_hash"], "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None, "parameters": parameters,
            "rng": {"seed": int(plan["runtime"]["seed"] if seed is None else seed), "kind": "L'Ecuyer-CMRG", "threads": int(plan["runtime"]["threads"])}}


def common_parameters(config: dict, plan: dict, decision: dict | None) -> dict:
    design = config["design"]
    return {"group_column": design["group_column"], "group_levels": list(design["group_levels"]),
            "subject_column": design["blocking"].get("subject_column") if design["blocking"]["mode"] != "none" else None,
            "blocking_mode": design["blocking"]["mode"], "figure_formats": list(config["report"].get("figure_formats", [])),
            "execution_profile": config["runtime"]["execution_profile"], "eligibility": decision or {}}


def rscript_executable() -> str:
    import os
    import shutil
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def execute_r(request_value: dict) -> dict:
    from ..runtime import execute_stage
    return execute_stage(request_value, rscript=rscript_executable())


# --------------------------------------------------------------------------- shared declaration helpers
def endpoint_index(config: dict) -> dict:
    """Resolved zero-null/TREAT endpoints of the frozen families: (family_id, contrast_id) -> list of members."""
    from .. import inference_service
    sens = {s["model_id"]: s["id"] for s in config["preprocessing"].get("sensitivities", [])}
    models = {m["id"]: m for m in config["models"]}
    index: dict[tuple[str, str], list[dict]] = {}
    for family in config.get("multiplicity_families", []):
        if family["hypothesis_type"] not in ("protein_zero_null", "protein_treat"):
            continue
        for model_id in family["model_ids"]:
            if model_id not in models or models[model_id]["engine"] != "limma":
                continue
            matrix, _ = inference_service.matrix_for_model(models[model_id], sens)
            for contrast_id in family["contrast_ids"]:
                index.setdefault((family["id"], contrast_id), []).append({"family_id": family["id"], "model_id": model_id, "contrast_id": contrast_id,
                                                                         "hypothesis_type": family["hypothesis_type"], "input_matrix": matrix})
    return index


def resolve_endpoint(config: dict, leaf: dict, pointer: str, *, index: dict | None = None, codes=("E_SETRULE_UNKNOWN_CONTRAST",)) -> dict:
    index = index if index is not None else endpoint_index(config)
    key = (leaf.get("family_id"), leaf.get("contrast_id"))
    members = index.get(key, [])
    if leaf.get("model_id"):
        members = [m for m in members if m["model_id"] == leaf["model_id"]]
    if not members:
        raise PostDeRefusal(codes[0], f"{pointer}: family {key[0]!r} has no completed endpoint for contrast {key[1]!r}", pointer)
    if len(members) > 1:
        raise PostDeRefusal(codes[0], f"{pointer}: family {key[0]!r} holds contrast {key[1]!r} for several models {[m['model_id'] for m in members]}; name model_id", pointer)
    return dict(members[0])
