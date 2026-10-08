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
from .sensitivity import _read_observations, _units_fn, complete_rows, covariate_specs

CAPABILITY = "post_de_association"
PREREQUISITES = ("dea",)
MIN_UNITS = 4


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma", "statmod"]}]


def precheck(block: dict, raw: dict) -> None:
    """No declaration of this module blocks a run (D-59): an imputation request is replaced by complete-case analysis and an
    ordinal phenotype without levels is refused on its own; both are recorded by plan_checks."""
    return None


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
    method = block.get("method", "model")
    requested_corr = block.get("correlation", "pearson")
    analyses, sub, adaptations = [], [], []
    if primary["engine"] != "limma":   # D-59: association models are limma models; a DEqMS/proDA primary does not block them
        adaptations.append({"analysis": "engine", "item": primary["id"], "requested": primary["engine"], "used": "limma",
                            "reason": "association models are fitted with limma (trend/robust settings of the primary model)"})
    for index, phenotype in enumerate(block.get("phenotypes", [])):
        pointer = f"/post_de/association/phenotypes/{index}"
        column = phenotype["column"]
        if not observations or column not in observations[0]:
            sub.append({"analysis": "phenotype", "item": column, "state": "INAPPLICABLE", "reason_code": "E_DESIGN_TERM", "reason": f"{column!r} is not an observation metadata column"})
            continue
        if phenotype.get("type") == "ordinal" and not phenotype.get("levels"):
            sub.append({"analysis": "phenotype", "item": column, "state": "INAPPLICABLE", "reason_code": "E_PHENOTYPE_ORDINAL_LEVELS",
                        "reason": "an ordinal phenotype needs its ordered levels to be scored; this phenotype is not analysed"})
            continue
        if phenotype.get("missing", "complete_case") != "complete_case":   # SM33 forbids phenotype imputation: refuse that piece, run complete case
            sub.append({"analysis": "imputation", "item": column, "state": "INAPPLICABLE", "reason_code": "E_PHENOTYPE_IMPUTATION",
                        "reason": f"{phenotype['missing']!r} would impute the phenotype; phenotypes are never imputed"})
            adaptations.append({"analysis": "phenotype", "item": column, "requested": f"missing: {phenotype['missing']}", "used": "complete case",
                                "reason": "phenotype imputation is not valid (SM33); units with a missing phenotype are excluded and n is reported"})
        adjust_group = phenotype.get("adjust_for_group", True)
        scope = phenotype.get("scope", "pooled")
        complete = [dict(o, **{column: _score(o.get(column), phenotype)}) for o in observations if _score(o.get(column), phenotype) is not None]
        n_missing = len(observations) - len(complete)
        try:
            specs = covariate_specs(phenotype.get("adjust_for", []), complete, pointer) if complete else []
        except PostDeRefusal as error:
            sub.append({"analysis": "phenotype", "item": column, "state": "INAPPLICABLE", "reason_code": error.code, "reason": error.message}); continue
        complete, dropped = complete_rows(complete, specs)
        if dropped:
            adaptations.append({"analysis": "phenotype", "item": column, "requested": "all units with the phenotype", "used": f"complete case for the covariates ({len(complete)} observations)",
                                "reason": f"adjustment covariate missing for {len(dropped)} observation(s): {dropped[:10]}"})
            n_missing = len(observations) - len(complete)
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
            if not group_term and analysis_scope == "pooled" and len(groups_here) > 1:
                # Review 2026-10-05 M2 (D-55): aliasing with group is checked on [phenotype + group] even when group is not fitted;
                # a phenotype determined by group stays full rank without the group term but must not be pooled (SM33, D-47).
                try:
                    check = design_service.replan(config, rows, design_id=f"{analysis_id}__alias_check", group_levels=groups_here, continuous=[column] + conts,
                                                  categorical=cats, contrast_ids=[], extra_contrasts=[dict(extra[0], id=f"{analysis_id}__alias_check",
                                                                                                             design_id=f"{analysis_id}__alias_check")], group_term=True)
                except ProteomicsError as error:   # the same terms passed the first replan, so a failure here is the group term (rank)
                    return None, ("E_PHENOTYPE_ALIASED", f"the phenotype is determined by group on the analysed units ({error.code}: {error.message})")
                if not check["full_rank"] or not check["contrasts"][0]["estimable"]:
                    return None, ("E_PHENOTYPE_ALIASED", f"the phenotype is determined by group on the analysed units (aliased {check['aliased']}), "
                                                         "although the group term is not fitted")
            adjusted = bool(group_term) or bool(specs)
            used = requested_corr if (requested_corr == "partial" or not adjusted) else f"{requested_corr}_partial"
            if method == "correlation" and used != requested_corr:   # D-59: the partial form of the requested correlation (never unadjusted values labelled adjusted)
                adaptations.append({"analysis": analysis_id, "item": column, "requested": requested_corr, "used": used,
                                    "reason": "an adjustment is declared (adjust_for_group defaults to true for pooled analyses); the partial form of the requested correlation is computed"})
            return {"analysis_id": analysis_id, "phenotype": column, "phenotype_type": phenotype["type"], "scope": analysis_scope, "group": group,
                    "correlation_requested": requested_corr, "correlation_used": used,
                    "adjust_for": [s["column"] for s in specs], "adjust_for_group": bool(group_term), "n_units": len({units(o) for o in rows}), "n_missing": n_missing,
                    "family_id": f"phenotype__{column}" + (f"__within_{group}" if group else ""), "coefficient": coef,
                    "values": {o["observation_id"]: o[column] for o in rows},
                    "design": {k: plan["design"][k] for k in ("observation_ids", "coefficients", "matrix", "term_map")}, "contrasts": plan["contrasts"], "rank": plan["rank"],
                    "separated_by_group": separated}, None

        within = scope in ("within_group", "both")
        if scope in ("pooled", "both"):
            reason = None
            if len(present) < 2 and len(levels) > 1:
                reason = ("E_PHENOTYPE_ALIASED", f"phenotype {column!r} is measured in one group only ({present})")
            if reason is None:
                analysis, reason = build(f"assoc_{column}", complete, adjust_group, "pooled")
            if reason is None:
                analyses.append(analysis)
            else:
                sub.append({"analysis": "pooled", "item": column, "state": "INAPPLICABLE", "reason_code": reason[0], "reason": reason[1]})
                if reason[0] == "E_PHENOTYPE_ALIASED" and not within:   # D-59: the pooled estimate is invalid; the within-group analyses are valid and run instead
                    within = True
                    adaptations.append({"analysis": f"assoc_{column}", "item": column, "requested": "scope: pooled", "used": "scope: within_group",
                                        "reason": "the phenotype is aliased with group, so the pooled estimate is refused and the within-group analyses are run instead"})
        if within:
            for g in levels:
                analysis, reason = build(f"assoc_{column}__{g}", by_group[g], False, "within_group", g)
                if reason is None:
                    analyses.append(analysis)
                else:
                    sub.append({"analysis": "within_group", "item": f"{column}/{g}", "state": "INAPPLICABLE", "reason_code": reason[0], "reason": reason[1]})
    parameters = {"analyses": analyses, "method": method, "correlation": requested_corr, "adaptations": adaptations,
                  "permutations": int(block.get("permutations", 9999)), "seed": int(block.get("seed", config["runtime"]["seed"]))}
    if not analyses:
        first = sub[0] if sub else {"reason_code": "E_PHENOTYPE_ALIASED", "reason": "no phenotype was declared"}
        return {"state": "INAPPLICABLE", "reason_code": first["reason_code"], "reason": "; ".join(f"{x['reason_code']}: {x['reason']}" for x in sub) or first["reason"],
                "subanalyses": sub, "adaptations": adaptations, "resolved": parameters}
    return {"state": "ELIGIBLE", "reason_code": None, "reason": None, "subanalyses": sub, "adaptations": adaptations, "resolved": parameters,
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
