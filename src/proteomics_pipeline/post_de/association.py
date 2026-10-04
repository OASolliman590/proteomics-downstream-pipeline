"""Post-DE protein-phenotype association (packet R14c, SM33, FR-144-FR-148).

Proteins are the response; each declared numeric or ordinal phenotype is a frozen design term (re-planned through
``design_service.replan``), adjusted for declared covariates and, by default, group, and fitted with the primary
engine in one secondary family per phenotype (and per group for within-group analyses).  Phenotypes are never
imputed: units with a missing phenotype are excluded for that phenotype and n is reported.  A pooled association
is refused (E_PHENOTYPE_ALIASED) when the phenotype is aliased with group, unless within-group analysis is declared.
"""
from __future__ import annotations

import math
from pathlib import Path

from . import PostDeRefusal, base_inputs, common_parameters, dea_inputs, execute_r, request
from .sensitivity import _read_observations, _units_fn, covariate_specs

CAPABILITY = "post_de_association"
PREREQUISITES = ("dea",)
MIN_UNITS = 4


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma", "statmod"]}]


def precheck(block: dict, raw: dict) -> None:
    for i, ph in enumerate(block.get("phenotypes", []) or []):
        if not isinstance(ph, dict):
            continue
        if ph.get("missing", "complete_case") != "complete_case":
            raise PostDeRefusal("E_PHENOTYPE_IMPUTATION", f"phenotype {ph.get('column')!r}: missing values are handled by complete-case exclusion only; "
                                f"{ph['missing']!r} would impute the phenotype", f"/post_de/association/phenotypes/{i}/missing")
        if ph.get("type") == "ordinal" and not ph.get("levels"):
            raise PostDeRefusal("E_PHENOTYPE_ORDINAL_LEVELS", f"ordinal phenotype {ph.get('column')!r} needs its ordered levels", f"/post_de/association/phenotypes/{i}/levels")


def _score(value, phenotype):
    if value in (None, "", "NA"):
        return None
    if phenotype["type"] == "ordinal":
        levels = list(phenotype["levels"])
        return float(levels.index(value) + 1) if value in levels else None
    try:
        x = float(value)
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def plan_checks(block: dict, config: dict, context: dict) -> dict:
    from .. import design_service
    from ..errors import ProteomicsError
    observations = context["observations"]
    design = config["design"]
    group_column, levels = design["group_column"], list(design["group_levels"])
    units = _units_fn(config)
    primary = next(m for m in config["models"] if m["role"] == "primary")
    if primary["engine"] != "limma":
        return {"state": "INAPPLICABLE", "reason_code": "E_SENSITIVITY_ENGINE", "reason": "association models use the limma primary engine"}
    method = block.get("method", "model")
    analyses, sub = [], []
    for index, phenotype in enumerate(block.get("phenotypes", [])):
        pointer = f"/post_de/association/phenotypes/{index}"
        column = phenotype["column"]
        if not observations or column not in observations[0]:
            raise PostDeRefusal("E_DESIGN_TERM", f"{pointer}: phenotype {column!r} is not an observation metadata column", pointer)
        adjust_group = phenotype.get("adjust_for_group", True)
        scope = phenotype.get("scope", "pooled")
        complete = [dict(o, **{column: _score(o.get(column), phenotype)}) for o in observations if _score(o.get(column), phenotype) is not None]
        n_missing = len(observations) - len(complete)
        specs = covariate_specs(phenotype.get("adjust_for", []), complete, pointer) if complete else []
        cats = [s["column"] for s in specs if s["type"] == "categorical"]; conts = [s["column"] for s in specs if s["type"] == "continuous"]
        coef = f"continuous.{design_service.encode(column)}"
        by_group = {g: [o for o in complete if o[group_column] == g] for g in levels}
        present = [g for g in levels if by_group[g]]
        separated = len(present) > 1 and all(max(o[column] for o in by_group[a]) < min(o[column] for o in by_group[b]) or max(o[column] for o in by_group[b]) < min(o[column] for o in by_group[a])
                                             for i, a in enumerate(present) for b in present[i + 1:])

        def build(analysis_id, rows, group_term, analysis_scope, group=None):
            if len({units(o) for o in rows}) < MIN_UNITS:
                return None, ("E_PHENOTYPE_SMALL_N", f"{len({units(o) for o in rows})} complete units (minimum {MIN_UNITS})")
            if len({o[column] for o in rows}) < 2:
                return None, ("E_PHENOTYPE_ALIASED", "the phenotype is constant on the analysed units")
            groups_here = [g for g in levels if any(o[group_column] == g for o in rows)]
            extra = [{"id": analysis_id, "design_id": analysis_id, "label": f"association with {column}", "estimand": f"log2 abundance change per unit of {column}",
                      "weights": {coef: 1}, "role": "secondary", "required_groups": groups_here}]
            try:
                plan = design_service.replan(config, rows, design_id=analysis_id, group_levels=groups_here, continuous=[column] + conts, categorical=cats,
                                             contrast_ids=[], extra_contrasts=extra, group_term=group_term)
            except ProteomicsError as error:
                return None, (error.code, error.message)
            if not plan["full_rank"] or not plan["contrasts"][0]["estimable"]:
                return None, ("E_PHENOTYPE_ALIASED", f"the phenotype is aliased with {'group/' if group_term else ''}covariates on the analysed units (aliased {plan['aliased']})")
            return {"analysis_id": analysis_id, "phenotype": column, "phenotype_type": phenotype["type"], "scope": analysis_scope, "group": group,
                    "adjust_for": [s["column"] for s in specs], "adjust_for_group": bool(group_term), "n_units": len({units(o) for o in rows}), "n_missing": n_missing,
                    "family_id": f"phenotype__{column}" + (f"__within_{group}" if group else ""), "coefficient": coef,
                    "values": {o["observation_id"]: o[column] for o in rows},
                    "design": {k: plan["design"][k] for k in ("observation_ids", "coefficients", "matrix", "term_map")}, "contrasts": plan["contrasts"], "rank": plan["rank"],
                    "separated_by_group": separated}, None

        if scope in ("pooled", "both"):
            reason = None
            if len(present) < 2 and len(levels) > 1:
                reason = ("E_PHENOTYPE_ALIASED", f"phenotype {column!r} is measured in one group only ({present}); declare scope within_group")
            if reason is None:
                analysis, reason = build(f"assoc_{column}", complete, adjust_group, "pooled")
                if reason and reason[0] == "E_PHENOTYPE_ALIASED":
                    reason = (reason[0], reason[1] + "; declare scope within_group")
            if reason is None:
                analyses.append(analysis)
            else:
                sub.append({"analysis": "pooled", "item": column, "state": "INAPPLICABLE", "reason_code": reason[0], "reason": reason[1]})
        if scope in ("within_group", "both"):
            for g in levels:
                analysis, reason = build(f"assoc_{column}__{g}", by_group[g], False, "within_group", g)
                if reason is None:
                    analyses.append(analysis)
                else:
                    sub.append({"analysis": "within_group", "item": f"{column}/{g}", "state": "INAPPLICABLE", "reason_code": reason[0], "reason": reason[1]})
    parameters = {"analyses": analyses, "method": method, "correlation": block.get("correlation", "pearson"),
                  "permutations": int(block.get("permutations", 9999)), "seed": int(block.get("seed", config["runtime"]["seed"]))}
    if not analyses:
        first = sub[0] if sub else {"reason_code": "E_PHENOTYPE_ALIASED", "reason": "no phenotype was declared"}
        return {"state": "INAPPLICABLE", "reason_code": first["reason_code"], "reason": first["reason"], "subanalyses": sub, "resolved": parameters}
    return {"state": "ELIGIBLE", "reason_code": None, "reason": None, "subanalyses": sub, "resolved": parameters,
            "rule": "complete-case units per phenotype (>= 4, non-constant); pooled analysis refused when the phenotype is aliased with group (one group only, or determined by group); one family per phenotype"}


def required_refusal(decision: dict):
    sub = [s for s in decision.get("subanalyses") or [] if s["analysis"] == "pooled"]
    return (sub[0]["reason_code"], sub[0]["reason"]) if sub else None


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, root: Path, decision: dict,
                  config_dir: Path | None = None, stage_id: str = CAPABILITY) -> dict:
    fresh = plan_checks(config["post_de"]["association"], config, {"observations": _read_observations(root)})
    primary = next(m for m in config["models"] if m["role"] == "primary")
    blocking = {b["design_id"]: b for b in plan.get("blocking", [])}
    inputs = base_inputs(plan, plan_path) + dea_inputs(config, root / "dea", names=("families.tsv",))
    parameters = {**common_parameters(config, plan, {k: v for k, v in decision.items() if k != "resolved"}), **fresh["resolved"],
                  "primary_model_id": primary["id"], "trend": primary.get("trend", True), "robust": primary.get("robust", True), "ci_level": config["report"]["ci_level"],
                  "coverage": {"policy": primary["coverage"]["policy"], "minimum_observed_per_group": primary["coverage"].get("minimum_observed_per_group", 2),
                               "minimum_fraction": primary["coverage"].get("minimum_fraction", 0.5)},
                  "consensus_correlation": blocking.get(config["design"]["id"], {}).get("consensus_correlation")}
    return request(CAPABILITY, stage_id, plan, run_id=run_id, output_temp_dir=output_temp_dir, inputs=inputs, parameters=parameters)


def execute(request_value: dict) -> dict:
    return execute_r(request_value)
