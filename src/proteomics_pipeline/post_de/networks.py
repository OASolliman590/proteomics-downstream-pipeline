"""Post-DE co-abundance modules and protein interaction networks (packet R14e, SM39-SM40, FR-161-FR-165).

Co-abundance is opt-in and refused below the declared minimum of biological units (E_COABUNDANCE_SMALL_N); only
features genuinely observed in every unit enter (no imputation).  Interaction analysis reads only a hashed, versioned
local PPI snapshot prepared by the explicit ``resources prepare`` action (kind ``ppi``); connectivity nulls are drawn
from the measured, mapped universe (a genome-wide null is refused: E_NETWORK_UNIVERSE).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from . import PostDeRefusal, base_inputs, common_parameters, dea_inputs, execute_r, request, stage_output

CAPABILITY = "post_de_networks"
PREREQUISITES = ("dea",)
RULES = ("wgcna_signed", "hclust_correlation")


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma", "dynamicTreeCut"]}]


def precheck(block: dict, raw: dict) -> None:
    ppi = block.get("ppi") or {}
    co = block.get("coabundance") or {}
    if not co.get("enabled") and not ppi.get("enabled"):   # SM39: nothing runs by default; an enabled module must declare an analysis
        raise PostDeRefusal("E_NETWORKS_NOT_DECLARED", "post_de.networks is enabled but declares neither coabundance nor ppi; co-abundance and network analyses never run by default",
                            "/post_de/networks")
    if ppi.get("enabled") and ppi.get("null_universe", "measured") != "measured":
        raise PostDeRefusal("E_NETWORK_UNIVERSE", f"null_universe {ppi['null_universe']!r}: connectivity nulls are drawn only from the measured, mapped universe, never the whole snapshot or genome",
                            "/post_de/networks/ppi/null_universe")
    if co.get("enabled") and co.get("rule", "wgcna_signed") not in RULES:
        raise PostDeRefusal("E_COABUNDANCE_RULE", f"unknown co-abundance rule {co.get('rule')!r}", "/post_de/networks/coabundance/rule")


def _matrix(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows[0][1:], {r[0]: r[1:] for r in rows[1:]}


def plan_checks(block: dict, config: dict, context: dict) -> dict:
    from .. import design_service
    from ..errors import ProteomicsError
    design = config["design"]; group_column = design["group_column"]
    subject = design["blocking"].get("subject_column") if design["blocking"]["mode"] != "none" else None
    observations = context["observations"]
    sub = []; parameters = {"coabundance": None, "ppi": None}
    co = block.get("coabundance") or {}
    if co.get("enabled"):
        minimum = int(co.get("min_units", 20))
        units = {(o[subject] if subject else (o.get("biological_unit_id") or o["observation_id"])) for o in observations}
        decision = {"state": "ELIGIBLE", "rule": co.get("rule", "wgcna_signed"), "soft_threshold": co.get("soft_threshold", "auto"), "min_module_size": int(co.get("min_module_size", 10)),
                    "bootstrap": int(co.get("bootstrap", 200)), "seed": int(co.get("seed", config["runtime"]["seed"])), "deep_split": int(co.get("deep_split", 2)), "min_units": minimum,
                    "cut_height": float(co.get("cut_height", 0.9 if co.get("rule", "wgcna_signed") == "wgcna_signed" else 0.5))}
        repeated = subject is not None and len(units) < len(observations)
        columns, matrix = _matrix(Path(context["preprocessing_dir"]) / "primary" / "matrix.tsv")
        _, mask = _matrix(Path(context["preprocessing_dir"]) / "primary" / "observed_mask.tsv")
        complete = [f for f in matrix if all(m == "true" for m in mask[f]) and all(v not in ("NA", "") for v in matrix[f])]
        if repeated:   # a design-type refusal comes first: more units would not make the correlation valid
            decision.update(state="INAPPLICABLE", reason_code="E_COABUNDANCE_DESIGN_UNSUPPORTED", reason="repeated observations per subject would be treated as independent units in the correlation; co-abundance supports independent-unit designs")
        elif len(units) < minimum:
            decision.update(state="INAPPLICABLE", reason_code="E_COABUNDANCE_SMALL_N", reason=f"{len(units)} biological units; co-abundance needs at least {minimum}")
        elif len(complete) < 2 * decision["min_module_size"]:
            decision.update(state="INAPPLICABLE", reason_code="E_COABUNDANCE_SMALL_N", reason=f"{len(complete)} complete features; at least {2 * decision['min_module_size']} are needed")
        if decision["state"] != "ELIGIBLE":
            sub.append({"analysis": "coabundance", "item": "modules", "state": "INAPPLICABLE", "reason_code": decision["reason_code"], "reason": decision["reason"]})
        else:
            decision["features"] = complete
            traits = []
            for trait in co.get("traits", ["group"]):
                if trait == "group":
                    planned = design_service.replan(config, observations, design_id="module_trait_group")
                    traits.append({"trait": "group", "model_id": "module_trait_group", "family_id": "module_trait__group", "design": {k: planned["design"][k] for k in ("observation_ids", "coefficients", "matrix", "term_map")},
                                   "contrasts": planned["contrasts"]})
                else:
                    if not observations or trait not in observations[0]:
                        raise PostDeRefusal("E_DESIGN_TERM", f"module trait {trait!r} is not an observation metadata column", "/post_de/networks/coabundance/traits")
                    coef = f"continuous.{design_service.encode(trait)}"
                    extra = [{"id": f"module_trait_{trait}", "design_id": f"module_trait_{trait}", "label": f"eigengene association with {trait}", "estimand": f"eigengene change per unit of {trait}",
                              "weights": {coef: 1}, "role": "secondary", "required_groups": list(design["group_levels"])}]
                    try:
                        planned = design_service.replan(config, observations, design_id=f"module_trait_{trait}", continuous=[trait], contrast_ids=[], extra_contrasts=extra)
                    except ProteomicsError as error:
                        sub.append({"analysis": "module_trait", "item": trait, "state": "INAPPLICABLE", "reason_code": error.code, "reason": error.message}); continue
                    traits.append({"trait": trait, "model_id": f"module_trait_{trait}", "family_id": f"module_trait__{trait}", "design": {k: planned["design"][k] for k in ("observation_ids", "coefficients", "matrix", "term_map")},
                                   "contrasts": planned["contrasts"]})
            decision["trait_designs"] = traits
        parameters["coabundance"] = decision
    ppi = block.get("ppi") or {}
    if ppi.get("enabled"):
        resources = {r["id"]: r for r in config.get("resources", [])}
        snapshot, mapping = resources.get(ppi.get("snapshot_id")), resources.get(ppi.get("mapping_resource_id"))
        if snapshot is None or snapshot["kind"] != "ppi":
            raise PostDeRefusal("E_RESOURCE_MISSING", f"ppi.snapshot_id {ppi.get('snapshot_id')!r} is not a declared resource of kind ppi", "/post_de/networks/ppi/snapshot_id")
        if mapping is None or mapping["kind"] != "mapping":
            raise PostDeRefusal("E_RESOURCE_MISSING", f"ppi.mapping_resource_id {ppi.get('mapping_resource_id')!r} is not a declared mapping resource", "/post_de/networks/ppi/mapping_resource_id")
        sets = {d["id"] for d in ((config.get("post_de") or {}).get("sets") or {}).get("definitions", [])}
        unknown = [s for s in ppi.get("sets", []) if s not in sets]
        if unknown or not ppi.get("sets"):
            raise PostDeRefusal("E_SETRULE_UNKNOWN_CONTRAST", f"ppi.sets must name declared post-DE sets; unknown: {unknown}", "/post_de/networks/ppi/sets")
        parameters["ppi"] = {"state": "ELIGIBLE", "snapshot_id": snapshot["id"], "mapping_resource_id": mapping["id"], "min_score": float(ppi.get("min_score", 0)),
                             "null_draws": int(ppi.get("null_draws", 2000)), "seed": int(ppi.get("seed", config["runtime"]["seed"])), "set_ids": list(ppi["sets"]),
                             "source": snapshot["source"], "release": snapshot["version"], "species": snapshot["target_taxonomy_id"], "score_type": snapshot.get("score_type"), "snapshot_sha256": snapshot["sha256"]}
    requested = [x for x in parameters.values() if x is not None]
    if requested and all(x["state"] != "ELIGIBLE" for x in requested):
        return {"state": "INAPPLICABLE", "reason_code": sub[0]["reason_code"], "reason": sub[0]["reason"], "subanalyses": sub, "resolved": parameters}
    if not requested:   # unreachable after precheck; kept typed for direct callers
        raise PostDeRefusal("E_NETWORKS_NOT_DECLARED", "post_de.networks declares neither coabundance nor ppi", "/post_de/networks")
    return {"state": "ELIGIBLE", "reason_code": None, "reason": None, "subanalyses": sub, "resolved": parameters,
            "rule": "co-abundance only when declared, >= min units, complete observed features, independent units; PPI only from a declared hashed ppi snapshot with a mapping resource and declared sets"}


def required_refusal(decision: dict):
    sub = decision.get("subanalyses") or []
    return (sub[0]["reason_code"], sub[0]["reason"]) if sub else None


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, root: Path, decision: dict,
                  config_dir: Path | None = None, stage_id: str = CAPABILITY) -> dict:
    from .. import resources as resource_service
    from ..errors import ProteomicsError
    from .sensitivity import _read_observations
    fresh = plan_checks(config["post_de"]["networks"], config, {"observations": _read_observations(root), "preprocessing_dir": root / "preprocessing"})
    parameters = fresh["resolved"]
    inputs = base_inputs(plan, plan_path) + dea_inputs(config, root / "dea", names=("families.tsv",))
    ppi = parameters.get("ppi")
    if ppi and ppi["state"] == "ELIGIBLE":
        verified = {e["resource_id"]: e for e in resource_service.verify_resources(config, config_dir or root)}   # SM14: every snapshot re-verified; never downloaded
        for rid, key in ((ppi["snapshot_id"], "edge_file"), (ppi["mapping_resource_id"], "mapping_file")):
            entry = verified[rid]
            inputs.append({"artifact_id": f"resource_{rid}", "path": entry["path"], "sha256": entry["sha256"]})
            for item in entry["files"]:
                inputs.append({"artifact_id": f"resfile_{rid}__{Path(item['name']).stem}", "path": item["path"], "sha256": item["sha256"]})
            ppi[key] = Path(entry["files"][0]["name"]).stem
        found = stage_output(root / "post_de" / "sets", "membership.tsv")
        if found is None:
            raise ProteomicsError("E_PREREQUISITE_FAILED", "PPI connectivity needs the completed post-DE sets stage", exit_code=4)
        inputs.append({"artifact_id": "stage__sets_membership", "path": str(found[0].resolve()), "sha256": found[1]})
        with found[0].open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        ppi["sets"] = [{"id": s, "members": [r["feature_id"] for r in rows if r.get(s) == "true"]} for s in ppi["set_ids"]]
    primary = next(m for m in config["models"] if m["role"] == "primary")
    blocking = {b["design_id"]: b for b in plan.get("blocking", [])}
    parameters = {**common_parameters(config, plan, {k: v for k, v in decision.items() if k != "resolved"}), **parameters,
                  "organism_taxonomy_id": config["organism"]["taxonomy_id"], "trend": primary.get("trend", True), "robust": primary.get("robust", True),
                  "ci_level": config["report"]["ci_level"], "consensus_correlation": blocking.get(config["design"]["id"], {}).get("consensus_correlation"),
                  "coverage": {"policy": "available_case", "minimum_observed_per_group": 2, "minimum_fraction": 0.5}}
    return request(CAPABILITY, stage_id, plan, run_id=run_id, output_temp_dir=output_temp_dir, inputs=inputs, parameters=parameters)


def execute(request_value: dict) -> dict:
    return execute_r(request_value)
