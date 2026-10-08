"""Post-DE robustness and sensitivity (packet R14b, SM32, FR-137-FR-143).

The planner classifies every requested adjustment before any refit (exact rank and estimability through
``design_service.replan``): ``estimable_additive``, ``estimable_interaction`` or ``non_estimable`` with the
aliased coefficients or empty cells named.  Subgroups are refused (E_SUBGROUP_COVERAGE) when a required
group falls below the primary model's coverage rule; no partial fit is reported.  Matched-n resampling and
leave-one-unit-out influence act on biological units (subjects in subject-blocked designs).
"""
from __future__ import annotations

import math
from pathlib import Path

from . import PostDeRefusal, _plan_artifact, base_inputs, common_parameters, dea_inputs, execute_r, request

CAPABILITY = "post_de_sensitivity"
CRITERIA = ("family_q", "exploratory_raw_p")
DEFAULT_CRITERIA = [{"criterion": "family_q", "threshold": 0.05, "direction": "any"}]


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma", "statmod"]}]


def precheck(block: dict, raw: dict) -> None:
    for i, crit in enumerate(block.get("criteria", []) or []):
        if isinstance(crit, dict) and crit.get("criterion") not in CRITERIA:
            raise PostDeRefusal("E_SETRULE_GRAMMAR", f"post_de.sensitivity.criteria[{i}]: criterion must be one of {CRITERIA}", f"/post_de/sensitivity/criteria/{i}")
    ids = [m.get("id") for m in block.get("covariate_models", []) or []] + [s.get("id") for s in block.get("subgroups", []) or []]
    if len(set(ids)) != len(ids):
        raise PostDeRefusal("E_ID_DUPLICATE", "post_de.sensitivity model and subgroup ids must be unique", "/post_de/sensitivity")


def _number(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def covariate_specs(entries: list, observations: list[dict], pointer: str) -> list[dict]:
    specs = []
    for entry in entries:
        column, kind = (entry, None) if isinstance(entry, str) else (entry["column"], entry["type"])
        if not observations or column not in observations[0]:
            raise PostDeRefusal("E_DESIGN_TERM", f"{pointer}: covariate {column!r} is not an observation metadata column", pointer)
        values = [o.get(column) for o in observations if o.get(column) not in (None, "", "NA")]   # D-59: missing values -> complete case, recorded by the caller
        if kind is None:
            kind = "continuous" if all(_number(v) is not None for v in values) and len(set(values)) > 2 else "categorical"
        specs.append({"column": column, "type": kind, "type_source": "declared" if not isinstance(entry, str) else "inferred from values"})
    return specs


def complete_rows(rows: list[dict], specs: list[dict]) -> tuple[list[dict], list[str]]:
    """Rows with every covariate present, and the ids dropped (complete-case removal is recorded, never silent)."""
    keep = [o for o in rows if all(o.get(s["column"]) not in (None, "", "NA") for s in specs)]
    kept = {o["observation_id"] for o in keep}
    return keep, [o["observation_id"] for o in rows if o["observation_id"] not in kept]


def imbalance(specs: list[dict], observations: list[dict], group_column: str, levels: list[str], units) -> list[dict]:
    rows = []
    for spec in specs:
        column = spec["column"]
        if spec["type"] == "categorical":
            for level in sorted({o[column] for o in observations}):
                for group in levels:
                    n = len({units(o) for o in observations if o[group_column] == group and o[column] == level})
                    rows.append({"covariate": column, "level": level, "group": group, "n_units": n, "mean": None, "min": None, "max": None,
                                 "aliased_with_group": len({o[group_column] for o in observations if o[column] == level}) == 1})
        else:
            for group in levels:
                xs = [float(o[column]) for o in observations if o[group_column] == group]
                rows.append({"covariate": column, "level": "continuous", "group": group, "n_units": len({units(o) for o in observations if o[group_column] == group}),
                             "mean": sum(xs) / len(xs) if xs else None, "min": min(xs) if xs else None, "max": max(xs) if xs else None, "aliased_with_group": False})
    return rows


def _units_fn(config: dict):
    design = config["design"]
    subject = design["blocking"].get("subject_column") if design["blocking"]["mode"] != "none" else None
    return (lambda o: o[subject]) if subject else (lambda o: o.get("biological_unit_id") or o["observation_id"])


def _model_record(kind, model_id, plan, description, state="ELIGIBLE", reason_code=None, reason=None):
    record = {"kind": kind, "model_id": model_id, "state": state, "reason_code": reason_code, "reason": reason, "description": description}
    if plan is not None:
        record.update({"design": {"observation_ids": plan["design"]["observation_ids"], "coefficients": plan["design"]["coefficients"], "matrix": plan["design"]["matrix"],
                                  "term_map": plan["design"]["term_map"]},
                       "contrasts": plan["contrasts"], "rank": plan["rank"], "aliased": plan["aliased"]})
    return record


def plan_checks(block: dict, config: dict, context: dict) -> dict:
    from .. import design_service
    from ..errors import ProteomicsError
    observations = context["observations"]
    design = config["design"]
    group_column, levels = design["group_column"], list(design["group_levels"])
    units = _units_fn(config)
    primary = next(m for m in config["models"] if m["role"] == "primary")
    adaptations = []
    if primary["engine"] != "limma":   # D-59: sensitivity models are limma models; a DEqMS/proDA primary does not block them
        adaptations.append({"analysis": "engine", "item": primary["id"], "requested": primary["engine"], "used": "limma",
                            "reason": "sensitivity models are fitted with limma (trend/robust settings of the primary model)"})
    coverage = primary["coverage"]
    min_units = max(2, int(coverage.get("minimum_observed_per_group", 2)))
    primary_contrasts = [c for c in config["contrasts"] if c["design_id"] == design["id"]]
    models, classification, imbalance_rows, sub = [], [], [], []
    all_specs = []
    for index, cm in enumerate(block.get("covariate_models", [])):
        pointer = f"/post_de/sensitivity/covariate_models/{index}"
        try:
            specs = covariate_specs(cm["add_covariates"], observations, pointer)
        except PostDeRefusal as error:   # a covariate that is not a metadata column refuses this model only
            sub.append({"analysis": "covariate_model", "item": cm["id"], "state": "INAPPLICABLE", "reason_code": error.code, "reason": error.message}); continue
        all_specs += [s for s in specs if s["column"] not in {x["column"] for x in all_specs}]
        cats = [s["column"] for s in specs if s["type"] == "categorical"]; conts = [s["column"] for s in specs if s["type"] == "continuous"]
        rows, dropped = complete_rows(observations, specs)
        if dropped:
            adaptations.append({"analysis": "covariate_model", "item": cm["id"], "requested": "all observations", "used": f"complete case ({len(rows)} observations)",
                                "reason": f"covariate value missing for {len(dropped)} observation(s): {dropped[:10]}"})
        try:
            additive = design_service.replan(config, rows, design_id=cm["id"], categorical=cats, continuous=conts)
        except ProteomicsError as error:
            additive = None; detail = error.message
        if additive is not None and additive["full_rank"] and all(c["estimable"] for c in additive["contrasts"]):
            classification.append({"model_id": cm["id"], "kind": "covariate", "class": "estimable_additive", "reason_code": None, "detail": f"rank {additive['rank']} of {len(additive['design']['coefficients'])}"})
            models.append(_model_record("covariate", cm["id"], additive, f"primary design + {', '.join(s['column'] + ' (' + s['type'] + ')' for s in specs)}"))
        else:
            detail = detail if additive is None else f"aliased coefficients {additive['aliased']}" if not additive["full_rank"] else "a planned contrast is not estimable"
            classification.append({"model_id": cm["id"], "kind": "covariate", "class": "non_estimable", "reason_code": "E_SENSITIVITY_NONESTIMABLE", "detail": detail})
            models.append(_model_record("covariate", cm["id"], None, "covariate-adjusted model", "INAPPLICABLE", "E_SENSITIVITY_NONESTIMABLE", detail))
            sub.append({"analysis": "covariate_model", "item": cm["id"], "state": "INAPPLICABLE", "reason_code": "E_SENSITIVITY_NONESTIMABLE", "reason": detail})
        if cm.get("interaction_with_group"):
            iid = f"{cm['id']}__interaction"
            empty = [f"{g}/{lvl}" for s in specs if s["type"] == "categorical" for lvl in sorted({o[s["column"]] for o in observations}) for g in levels
                     if not any(o[group_column] == g and o[s["column"]] == lvl for o in observations)]
            interaction = None
            if not empty:
                try:
                    probe = design_service.replan(config, rows, design_id=iid, categorical=cats, continuous=conts,
                                                  interactions=[[group_column, s["column"]] for s in specs], contrast_ids=[])
                    names = [t["coefficient"] for t in probe["design"]["term_map"] if t["term"] == "interaction"]
                    extra = [{"id": f"{iid}:{n}", "design_id": iid, "label": f"interaction {n}", "estimand": "group x covariate interaction (difference of effects)",
                              "weights": {n: 1}, "role": "secondary", "required_groups": list(probe["group_levels"])} for n in names]
                    interaction = design_service.replan(config, rows, design_id=iid, categorical=cats, continuous=conts,
                                                        interactions=[[group_column, s["column"]] for s in specs], contrast_ids=[], extra_contrasts=extra)
                except ProteomicsError as error:
                    empty = [error.message]
            if interaction is not None and interaction["full_rank"] and all(c["estimable"] for c in interaction["contrasts"]):
                classification.append({"model_id": iid, "kind": "interaction", "class": "estimable_interaction", "reason_code": None, "detail": f"{len(interaction['contrasts'])} interaction coefficient(s)"})
                models.append(_model_record("interaction", iid, interaction, f"formal group x {', '.join(s['column'] for s in specs)} interaction (own family)"))
            else:
                detail = (f"empty group x covariate cell(s): {', '.join(empty)}" if empty else f"aliased coefficients {interaction['aliased']}")
                classification.append({"model_id": iid, "kind": "interaction", "class": "non_estimable", "reason_code": "E_SENSITIVITY_NONESTIMABLE", "detail": detail})
                models.append(_model_record("interaction", iid, None, "formal interaction", "INAPPLICABLE", "E_SENSITIVITY_NONESTIMABLE", detail))
                sub.append({"analysis": "interaction", "item": iid, "state": "INAPPLICABLE", "reason_code": "E_SENSITIVITY_NONESTIMABLE", "reason": detail})
    imbalance_rows = imbalance(all_specs, observations, group_column, levels, units)
    for index, sg in enumerate(block.get("subgroups", [])):
        column, allowed = sg["filter"]["column"], set(sg["filter"]["in"])
        if not observations or column not in observations[0]:
            raise PostDeRefusal("E_DESIGN_TERM", f"subgroup {sg['id']!r}: column {column!r} is not an observation metadata column", f"/post_de/sensitivity/subgroups/{index}")
        subset = [o for o in observations if o.get(column) in allowed]
        short = []
        for contrast in primary_contrasts:
            for g in contrast["required_groups"]:
                n = len({units(o) for o in subset if o[group_column] == g})
                if n < min_units and f"{g} ({n})" not in short:
                    short.append(f"{g} ({n})")
        if short:
            reason = f"subgroup {column} in {sorted(allowed)}: required group(s) below the coverage rule of {min_units} units: {', '.join(short)}"
            models.append(_model_record("subgroup", sg["id"], None, reason, "INAPPLICABLE", "E_SUBGROUP_COVERAGE", reason))
            sub.append({"analysis": "subgroup", "item": sg["id"], "state": "INAPPLICABLE", "reason_code": "E_SUBGROUP_COVERAGE", "reason": reason})
            continue
        planned = design_service.replan(config, observations, design_id=sg["id"], observation_ids=[o["observation_id"] for o in subset])
        if not planned["full_rank"] or not all(c["estimable"] for c in planned["contrasts"]):
            reason = f"subgroup design is not estimable (aliased {planned['aliased']})"
            models.append(_model_record("subgroup", sg["id"], None, reason, "INAPPLICABLE", "E_SENSITIVITY_NONESTIMABLE", reason))
            sub.append({"analysis": "subgroup", "item": sg["id"], "state": "INAPPLICABLE", "reason_code": "E_SENSITIVITY_NONESTIMABLE", "reason": reason}); continue
        classification.append({"model_id": sg["id"], "kind": "subgroup", "class": "estimable_additive", "reason_code": None, "detail": f"{len(subset)} observations"})
        models.append(_model_record("subgroup", sg["id"], planned, f"subgroup {column} in {sorted(allowed)} (re-planned on the subset)"))
    matched = None
    if block.get("matched_n"):
        mn = block["matched_n"]
        matched = {"targets": dict(mn["targets"]), "draws": int(mn.get("draws", 500)), "seed": int(mn.get("seed", config["runtime"]["seed"])), "state": "ELIGIBLE"}
        varying = design["blocking"]["mode"] != "none" and any(len({o[group_column] for o in observations if units(o) == u}) > 1 for u in {units(o) for o in observations})
        problems = [f"{g}: not a declared group" for g in mn["targets"] if g not in levels]
        problems += [f"{g}: target {t} is not below its {len({units(o) for o in observations if o[group_column] == g})} units" for g, t in mn["targets"].items()
                     if g in levels and t >= len({units(o) for o in observations if o[group_column] == g})]
        if varying:
            matched.update(state="INAPPLICABLE", reason_code="E_SENSITIVITY_DESIGN_UNSUPPORTED", reason="matched-n subsampling of groups is undefined when group varies within subject")
        elif problems:
            matched.update(state="INAPPLICABLE", reason_code="E_SUBGROUP_COVERAGE", reason="; ".join(problems))
        if matched["state"] != "ELIGIBLE":
            sub.append({"analysis": "matched_n", "item": "matched_n", "state": "INAPPLICABLE", "reason_code": matched["reason_code"], "reason": matched["reason"]})
    criteria = [dict({"direction": "any"}, **c) for c in block.get("criteria", DEFAULT_CRITERIA)]
    eligible = [m for m in models if m["state"] == "ELIGIBLE"]
    nothing = not eligible and not (matched and matched["state"] == "ELIGIBLE") and not block.get("influence")
    parameters = {"models": models, "matched_n": matched, "influence": bool(block.get("influence", False)), "criteria": criteria,
                  "imbalance": imbalance_rows, "classification": classification, "adaptations": adaptations}
    if nothing:
        first = next((s for s in sub), {"reason_code": "E_SENSITIVITY_NONESTIMABLE", "reason": "no sensitivity analysis was declared"})
        return {"state": "INAPPLICABLE", "reason_code": first["reason_code"], "reason": first["reason"], "subanalyses": sub, "adaptations": adaptations, "resolved": parameters}
    return {"state": "ELIGIBLE", "reason_code": None, "reason": None, "subanalyses": sub, "adaptations": adaptations, "resolved": parameters, "classification": classification,
            "rule": "covariate/subgroup models need full rank and estimable contrasts (exact); subgroups need every required group at the coverage minimum; matched-n needs between-unit groups"}


def required_refusal(decision: dict) -> tuple[str, str] | None:
    """A required sensitivity module rejects the plan when any declared analysis is ineligible (V137 negative: before any fit)."""
    sub = decision.get("subanalyses") or []
    return (sub[0]["reason_code"], sub[0]["reason"]) if sub else None


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, root: Path, decision: dict,
                  config_dir: Path | None = None, stage_id: str = CAPABILITY) -> dict:
    import json
    observations = _read_observations(root)
    design_request = json.loads((root / "designs" / "stage-request.json").read_text(encoding="utf-8"))
    fresh = plan_checks(config["post_de"]["sensitivity"], config, {"observations": observations})
    primary = next(m for m in config["models"] if m["role"] == "primary")
    compiled = next(d for d in design_request["parameters"]["designs"] if d["design_id"] == config["design"]["id"])
    blocking = {b["design_id"]: b for b in plan.get("blocking", [])}
    inputs = base_inputs(plan, plan_path, design_ids=[config["design"]["id"]])
    lineage = _plan_artifact(plan, "observation_lineage")
    if lineage is not None:
        inputs.append({"artifact_id": "observation_lineage", "path": str((Path(plan_path).resolve().parent / lineage["relative_path"]).resolve()), "sha256": lineage["sha256"]})
    inputs += dea_inputs(config, root / "dea")
    parameters = {**common_parameters(config, plan, {k: v for k, v in decision.items() if k != "resolved"}), **fresh["resolved"],
                  "primary_model_id": primary["id"], "trend": primary.get("trend", True), "robust": primary.get("robust", True), "ci_level": config["report"]["ci_level"],
                  "coverage": {"policy": primary["coverage"]["policy"], "minimum_observed_per_group": primary["coverage"].get("minimum_observed_per_group", 2),
                               "minimum_fraction": primary["coverage"].get("minimum_fraction", 0.5)},
                  "consensus_correlation": blocking.get(config["design"]["id"], {}).get("consensus_correlation"),
                  "primary_design": {k: compiled[k] for k in ("observation_ids", "coefficients", "matrix", "term_map")},
                  "primary_contrasts": [c for c in design_request["parameters"]["contrasts"] if c["design_id"] == config["design"]["id"]]}
    return request(CAPABILITY, stage_id, plan, run_id=run_id, output_temp_dir=output_temp_dir, inputs=inputs, parameters=parameters)


def _read_observations(root: Path) -> list[dict]:
    import csv
    with (root / "preprocessing" / "primary" / "observations.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def execute(request_value: dict) -> dict:
    return execute_r(request_value)
