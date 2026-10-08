"""Post-DE co-abundance modules and protein interaction networks (packet R14e, SM39-SM40, FR-161-FR-165).

Co-abundance is opt-in.  It adapts to the data (D-59): at the declared minimum of biological units it runs normally,
from 4 units up to it as exploratory (recorded), and below 4 units it is refused (E_COABUNDANCE_SMALL_N).  Repeated
designs use subject means or subject-effect removal.  Only features genuinely observed in every unit enter (no
imputation).  Interaction analysis reads only a hashed, versioned local PPI snapshot prepared by the explicit
``resources prepare`` action (kind ``ppi``); connectivity nulls are drawn from the measured, mapped universe (a
requested genome-wide null is refused, E_NETWORK_UNIVERSE, and the measured-universe null runs instead).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from . import PostDeRefusal, base_inputs, common_parameters, dea_inputs, execute_r, request, stage_output

CAPABILITY = "post_de_networks"
COABUNDANCE_MIN_COMPUTABLE = 4   # D-59: fewer units leave every correlation with at most one degree of freedom
COABUNDANCE_MIN_MODULE = 3
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
        columns, matrix = _matrix(Path(context["preprocessing_dir"]) / "primary" / "matrix.tsv")
        _, mask = _matrix(Path(context["preprocessing_dir"]) / "primary" / "observed_mask.tsv")
        complete = [f for f in matrix if all(m == "true" for m in mask[f]) and all(v not in ("NA", "") for v in matrix[f])]
        adaptations = decision["adaptations"] = []
        # D-59: repeated observations are handled per design instead of refused.  Group constant within subject: subject means
        # (one row per biological unit, between-subject covariation kept).  Group varying within subject (paired or mixed):
        # subject effects removed with the group term kept (limma::removeBatchEffect), units resampled as whole subjects.
        transform = "none"
        if subject is not None and len(units) < len(observations):
            groups_of = {}
            for o in observations:
                groups_of.setdefault(o[subject], set()).add(o[group_column])
            transform = "subject_means" if all(len(g) == 1 for g in groups_of.values()) else "remove_subject"
            adaptations.append({"analysis": "coabundance", "item": "repeated_observations", "requested": "observations as units", "used": transform,
                                "reason": ("group is constant within subject: modules are built on subject means (one row per subject)" if transform == "subject_means" else
                                           "group varies within subject: subject effects are removed (group term kept) before correlation; whole subjects are resampled")})
        decision.update(transform=transform, subject_column=subject)
        n_units = len(units)
        if n_units < COABUNDANCE_MIN_COMPUTABLE:
            decision.update(state="INAPPLICABLE", reason_code="E_COABUNDANCE_SMALL_N",
                            reason=f"{n_units} biological units; at least {COABUNDANCE_MIN_COMPUTABLE} are needed to compute correlations (with 3 or fewer units every correlation has at most one degree of freedom)")
        elif n_units < minimum:
            adaptations.append({"analysis": "coabundance", "item": "units", "requested": f"at least {minimum} units (recommended)", "used": f"{n_units} (exploratory)",
                                "reason": "below the recommended number of units modules are unstable; results are exploratory and the soft threshold falls back to the n-dependent default when no scale-free fit is reached"})
            decision["exploratory"] = True
        if decision["state"] == "ELIGIBLE" and len(complete) < 2 * decision["min_module_size"]:
            if len(complete) >= 2 * COABUNDANCE_MIN_MODULE:
                size = max(COABUNDANCE_MIN_MODULE, len(complete) // 2)
                adaptations.append({"analysis": "coabundance", "item": "min_module_size", "requested": decision["min_module_size"], "used": size,
                                    "reason": f"{len(complete)} complete features cannot hold two modules of the declared size"})
                decision["min_module_size"] = size
            else:
                decision.update(state="INAPPLICABLE", reason_code="E_COABUNDANCE_SMALL_N",
                                reason=f"{len(complete)} complete features; at least {2 * COABUNDANCE_MIN_MODULE} are needed for two modules of {COABUNDANCE_MIN_MODULE}")
        if decision["state"] != "ELIGIBLE":
            sub.append({"analysis": "coabundance", "item": "modules", "state": "INAPPLICABLE", "reason_code": decision["reason_code"], "reason": decision["reason"]})
        else:
            decision["features"] = complete
            import copy
            cfg = copy.deepcopy(config); cfg["design"]["blocking"] = {"mode": "none"}
            if transform == "subject_means":   # trait designs on one pseudo-observation per subject
                rows = []
                for u in sorted(units):
                    members = [o for o in observations if o[subject] == u]
                    row = {k: v for k, v in members[0].items()}; row["observation_id"] = u; row["biological_unit_id"] = u
                    for col in co.get("traits", []):
                        if col != "group":
                            vals = [float(o[col]) for o in members if o.get(col) not in (None, "", "NA")]
                            row[col] = repr(sum(vals) / len(vals)) if vals else "NA"
                    rows.append(row)
                trait_rows, categorical, blocking = rows, [], "none"
            elif transform == "remove_subject":   # paired analysis of eigengenes: subject as a fixed effect
                trait_rows, categorical, blocking = observations, [subject], "none"
            else:
                trait_rows, categorical, blocking = observations, [], "primary"
            traits = []
            for trait in co.get("traits", ["group"]):
                try:
                    if trait == "group":
                        planned = design_service.replan(cfg if blocking == "none" else config, trait_rows, design_id="module_trait_group", categorical=categorical)
                    else:
                        if not trait_rows or trait not in trait_rows[0]:
                            sub.append({"analysis": "module_trait", "item": trait, "state": "INAPPLICABLE", "reason_code": "E_DESIGN_TERM", "reason": f"{trait!r} is not an observation metadata column"}); continue
                        coef = f"continuous.{design_service.encode(trait)}"
                        extra = [{"id": f"module_trait_{trait}", "design_id": f"module_trait_{trait}", "label": f"eigengene association with {trait}", "estimand": f"eigengene change per unit of {trait}",
                                  "weights": {coef: 1}, "role": "secondary", "required_groups": list(design["group_levels"])}]
                        planned = design_service.replan(cfg if blocking == "none" else config, trait_rows, design_id=f"module_trait_{trait}", continuous=[trait], categorical=categorical,
                                                        contrast_ids=[], extra_contrasts=extra)
                except ProteomicsError as error:
                    sub.append({"analysis": "module_trait", "item": trait, "state": "INAPPLICABLE", "reason_code": error.code, "reason": error.message}); continue
                if not planned["full_rank"] or not all(c["estimable"] for c in planned["contrasts"]):
                    sub.append({"analysis": "module_trait", "item": trait, "state": "INAPPLICABLE", "reason_code": "E_MODULE_TRAIT_NONESTIMABLE",
                                "reason": f"the {trait} effect is not estimable on the {transform} rows (aliased {planned['aliased']})"}); continue
                traits.append({"trait": trait, "model_id": f"module_trait_{trait}", "family_id": f"module_trait__{trait}", "blocking": blocking,
                               "design": {k: planned["design"][k] for k in ("observation_ids", "coefficients", "matrix", "term_map")}, "contrasts": planned["contrasts"]})
            decision["trait_designs"] = traits
        parameters["coabundance"] = decision
    ppi = block.get("ppi") or {}
    if ppi.get("enabled"):
        resources = {r["id"]: r for r in config.get("resources", [])}
        snapshot, mapping = resources.get(ppi.get("snapshot_id")), resources.get(ppi.get("mapping_resource_id"))
        problem = None
        if snapshot is None or snapshot["kind"] != "ppi":
            problem = ("E_RESOURCE_MISSING", f"ppi.snapshot_id {ppi.get('snapshot_id')!r} is not a declared resource of kind ppi")
        elif mapping is None or mapping["kind"] != "mapping":
            problem = ("E_RESOURCE_MISSING", f"ppi.mapping_resource_id {ppi.get('mapping_resource_id')!r} is not a declared mapping resource")
        sets = {d["id"] for d in ((config.get("post_de") or {}).get("sets") or {}).get("definitions", [])}
        unknown = [x for x in ppi.get("sets", []) if x not in sets]
        known_sets = [x for x in ppi.get("sets", []) if x in sets]
        if problem is None and not known_sets:
            problem = ("E_SETRULE_UNKNOWN_CONTRAST", f"ppi.sets must name declared post-DE sets; unknown: {unknown}")
        if problem is None and unknown:   # D-59: unknown sets are refused one by one, the others are tested
            sub.append({"analysis": "ppi", "item": ",".join(unknown), "state": "INAPPLICABLE", "reason_code": "E_SETRULE_UNKNOWN_CONTRAST", "reason": f"undeclared post-DE set(s) {unknown}"})
        if problem is None and ppi.get("null_universe", "measured") != "measured":   # validity: a genome-wide null is refused; the measured-universe null runs instead
            sub.append({"analysis": "ppi", "item": "null_universe", "state": "INAPPLICABLE", "reason_code": "E_NETWORK_UNIVERSE",
                        "reason": f"null_universe {ppi['null_universe']!r}: nulls drawn from the whole snapshot or genome are invalid; the measured mapped universe is used"})
            (parameters.setdefault("_ppi_adaptations", [])).append({"analysis": "ppi", "item": "null_universe", "requested": ppi["null_universe"], "used": "measured",
                                                                    "reason": "connectivity nulls are drawn only from the measured, mapped universe (SM40)"})
        if problem is not None:
            sub.append({"analysis": "ppi", "item": "connectivity", "state": "INAPPLICABLE", "reason_code": problem[0], "reason": problem[1]})
            parameters["ppi"] = {"state": "INAPPLICABLE", "reason_code": problem[0], "reason": problem[1]}
    if ppi.get("enabled") and parameters["ppi"] is None:
        ppi = dict(ppi, sets=known_sets)
        parameters["ppi"] = {"state": "ELIGIBLE", "snapshot_id": snapshot["id"], "mapping_resource_id": mapping["id"], "min_score": float(ppi.get("min_score", 0)),
                             "null_draws": int(ppi.get("null_draws", 2000)), "seed": int(ppi.get("seed", config["runtime"]["seed"])), "set_ids": list(ppi["sets"]),
                             "source": snapshot["source"], "release": snapshot["version"], "species": snapshot["target_taxonomy_id"], "score_type": snapshot.get("score_type"), "snapshot_sha256": snapshot["sha256"]}
    ppi_adaptations = parameters.pop("_ppi_adaptations", [])
    requested = [x for x in parameters.values() if x is not None]
    adaptations = list((parameters["coabundance"] or {}).get("adaptations", [])) + ppi_adaptations
    if requested and all(x["state"] != "ELIGIBLE" for x in requested):
        return {"state": "INAPPLICABLE", "reason_code": sub[0]["reason_code"], "reason": sub[0]["reason"], "subanalyses": sub, "adaptations": adaptations, "resolved": parameters}
    if not requested:   # unreachable after precheck; kept typed for direct callers
        raise PostDeRefusal("E_NETWORKS_NOT_DECLARED", "post_de.networks declares neither coabundance nor ppi", "/post_de/networks")
    return {"state": "ELIGIBLE", "reason_code": None, "reason": None, "subanalyses": sub, "adaptations": adaptations, "resolved": parameters,
            "rule": "co-abundance only when declared; at least 4 units (exploratory below the declared minimum); complete observed features; repeated designs via subject means or subject-effect removal (D-59); PPI only from a declared hashed ppi snapshot with a mapping resource and declared sets"}


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
