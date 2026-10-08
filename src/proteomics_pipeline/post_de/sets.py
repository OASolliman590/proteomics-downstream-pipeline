"""Post-DE result structure: declarative set algebra, exact regions, direction-aware overlap, concordance (packet R14a, SM31, FR-131-FR-136).

The rule grammar is finite and validated before any computation: a rule is a JSON tree whose inner nodes are
``{"op": union|intersect|difference|complement_within_tested, "args": [...]}`` and whose leaves name one completed
endpoint ``{"family_id", "contrast_id", ["model_id"], "criterion": family_q|treat|exploratory_raw_p, "threshold", ["direction"]}``.
Strings, unknown operators and unknown keys are refused (E_SETRULE_GRAMMAR); no expression is ever evaluated.
"""
from __future__ import annotations

from pathlib import Path

from . import PostDeRefusal, base_inputs, common_parameters, dea_inputs, endpoint_index, execute_r, request, resolve_endpoint

CAPABILITY = "post_de_sets"
OPS = {"union": (2, None), "intersect": (2, None), "difference": (2, 2), "complement_within_tested": (1, 1)}
CRITERIA = {"family_q": "protein_zero_null", "treat": "protein_treat", "exploratory_raw_p": "protein_zero_null"}
LEAF_KEYS = {"family_id", "contrast_id", "model_id", "criterion", "threshold", "direction"}
MAX_UPSET_ALL_REGIONS = 12


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma"]}]


def _leaf(raw: dict, config: dict, index: dict, pointer: str) -> dict:
    extra = set(raw) - LEAF_KEYS
    if extra:
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: unknown leaf key(s) {sorted(extra)}", pointer)
    for key in ("family_id", "contrast_id", "criterion", "threshold"):
        if key not in raw:
            raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: a leaf needs {key}", pointer)
    if raw["criterion"] not in CRITERIA:
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: criterion {raw['criterion']!r} is not one of {sorted(CRITERIA)}", pointer)
    threshold = raw["threshold"]
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not 0 < threshold <= 1:
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: threshold must be a number in (0, 1]", pointer)
    if raw.get("direction", "any") not in ("up", "down", "any"):
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: direction must be up, down or any", pointer)
    if raw["contrast_id"] not in {c["id"] for c in config["contrasts"]}:
        raise PostDeRefusal("E_SETRULE_UNKNOWN_CONTRAST", f"{pointer}: contrast {raw['contrast_id']!r} is not declared", pointer)
    if raw["family_id"] not in {f["id"] for f in config.get("multiplicity_families", [])}:
        raise PostDeRefusal("E_SETRULE_UNKNOWN_CONTRAST", f"{pointer}: family {raw['family_id']!r} is not declared", pointer)
    endpoint = resolve_endpoint(config, raw, pointer, index=index)
    if CRITERIA[raw["criterion"]] != endpoint["hypothesis_type"]:
        raise PostDeRefusal("E_SETRULE_MIXED_NULL", f"{pointer}: criterion {raw['criterion']} cannot be applied to a {endpoint['hypothesis_type']} family", pointer)
    return {**endpoint, "criterion": raw["criterion"], "threshold": float(threshold), "direction": raw.get("direction", "any")}


def parse_rule(raw, config: dict, index: dict, pointer: str) -> dict:
    """Validate one rule tree against the frozen declarations; returns the resolved tree (leaves carry model and hypothesis)."""
    if not isinstance(raw, dict):
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: a rule must be a JSON object of the finite grammar; free-form expressions such as {str(raw)[:40]!r} are never evaluated", pointer)
    if "op" not in raw:
        return _leaf(raw, config, index, pointer)
    if set(raw) - {"op", "args"}:
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: an operator node holds only op and args", pointer)
    op = raw["op"]
    if op not in OPS:
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: unknown operator {op!r}; allowed {sorted(OPS)}", pointer)
    args = raw.get("args")
    low, high = OPS[op]
    if not isinstance(args, list) or len(args) < low or (high is not None and len(args) > high):
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: {op} takes {low}{'' if high == low else ('+' if high is None else f'-{high}')} argument(s)", pointer)
    return {"op": op, "args": [parse_rule(a, config, index, f"{pointer}/args/{i}") for i, a in enumerate(args)]}


def leaves(rule: dict) -> list[dict]:
    return [rule] if "op" not in rule else [leaf for arg in rule["args"] for leaf in leaves(arg)]


def resolve(block: dict, config: dict) -> dict:
    """D-59: an invalid definition, Venn/UpSet reference or concordance pair is refused on its own (recorded in `refused`);
    every valid set is still built."""
    index = endpoint_index(config)
    definitions, refused = [], []
    seen = set()
    for i, d in enumerate(block.get("definitions", [])):
        pointer = f"/post_de/sets/definitions/{i}/rule"
        try:
            if d.get("id") in seen:
                raise PostDeRefusal("E_ID_DUPLICATE", f"set id {d.get('id')!r} is declared twice; the later declaration is refused", f"/post_de/sets/definitions/{i}/id")
            rule = parse_rule(d.get("rule"), config, index, pointer)
            found = leaves(rule)
            if len({leaf["hypothesis_type"] for leaf in found}) > 1:
                raise PostDeRefusal("E_SETRULE_MIXED_NULL", f"{pointer}: rule {d['id']!r} mixes zero-null and TREAT memberships", pointer)
            if len({leaf["input_matrix"] for leaf in found}) > 1:
                raise PostDeRefusal("E_SETRULE_MIXED_MATRIX", f"{pointer}: rule {d['id']!r} combines families fitted on different matrices {sorted({l['input_matrix'] for l in found})}", pointer)
        except PostDeRefusal as error:
            refused.append({"analysis": "set", "item": str(d.get("id")), "state": "INAPPLICABLE", "reason_code": error.code, "reason": error.message}); continue
        seen.add(d["id"])
        definitions.append({"id": d["id"], "rule": rule})
    known = {d["id"] for d in definitions}
    views = {}
    for key in ("venn", "upset"):
        unknown = [x for x in block.get(key, []) if x not in known]
        if unknown:
            refused.append({"analysis": key, "item": ",".join(map(str, unknown)), "state": "INAPPLICABLE", "reason_code": "E_SETRULE_GRAMMAR",
                            "reason": f"post_de.sets.{key} names undefined or refused sets {unknown}; they are left out of the {key}"})
        views[key] = [x for x in block.get(key, []) if x in known]
    pairs = []
    for i, pair in enumerate(block.get("concordance_pairs", [])):
        criterion = pair.get("criterion", "family_q"); threshold = pair.get("threshold", 0.05)
        try:
            sides = {side: _leaf({**pair[side], "criterion": criterion, "threshold": threshold}, config, index, f"/post_de/sets/concordance_pairs/{i}/{side}") for side in ("left", "right")}
        except PostDeRefusal as error:
            refused.append({"analysis": "concordance", "item": pair.get("id", f"pair{i + 1}"), "state": "INAPPLICABLE", "reason_code": error.code, "reason": error.message}); continue
        pairs.append({"id": pair.get("id", f"pair{i + 1}"), **sides})
    return {"definitions": definitions, "venn": views["venn"], "upset": views["upset"], "concordance_pairs": pairs,
            "overlap_test": bool(block.get("overlap_test", False)), "refused": refused}


def precheck(block: dict, raw: dict) -> None:
    """Grammar errors of a *required* sets module are refused before the schema runs.  For an optional module each invalid
    definition is refused on its own in plan_checks (D-59)."""
    if (block.get("execution_requirement") or (raw.get("post_de") or {}).get("execution_requirement") or "optional") != "required":
        return
    for i, d in enumerate(block.get("definitions", []) or []):
        rule = d.get("rule") if isinstance(d, dict) else None
        _grammar_only(rule, f"/post_de/sets/definitions/{i}/rule")


def _grammar_only(raw, pointer):
    if not isinstance(raw, dict):
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: a rule must be a JSON object of the finite grammar; free-form expressions are never evaluated", pointer)
    if "op" in raw:
        if raw["op"] not in OPS or set(raw) - {"op", "args"} or not isinstance(raw.get("args"), list):
            raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: unknown operator or malformed node {str(raw)[:60]!r}", pointer)
        for i, a in enumerate(raw["args"]):
            _grammar_only(a, f"{pointer}/args/{i}")
    elif set(raw) - LEAF_KEYS:
        raise PostDeRefusal("E_SETRULE_GRAMMAR", f"{pointer}: unknown leaf key(s) {sorted(set(raw) - LEAF_KEYS)}", pointer)


def plan_checks(block: dict, config: dict, context: dict) -> dict:
    resolved = resolve(block, config)
    sub = list(resolved["refused"])
    if len(resolved["venn"]) > 3:
        sub.append({"analysis": "venn", "item": ",".join(resolved["venn"]), "state": "INAPPLICABLE", "reason_code": "E_VENN_K",
                    "reason": f"a circle Venn needs k <= 3 sets; {len(resolved['venn'])} were requested (UpSet is still produced)"})
    if not resolved["definitions"]:
        first = sub[0] if sub else {"reason_code": "E_SETRULE_GRAMMAR", "reason": "no set definitions were declared"}
        return {"state": "INAPPLICABLE", "reason_code": first["reason_code"], "reason": first["reason"], "subanalyses": sub, "resolved": resolved}
    return {"state": "ELIGIBLE", "reason_code": None, "reason": None, "subanalyses": sub, "resolved": resolved,
            "rule": "eligible for every design with completed zero-null or TREAT families; overlap P only for unit-disjoint contrasts"}


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, root: Path, decision: dict,
                  config_dir: Path | None = None, stage_id: str = CAPABILITY) -> dict:
    resolved = resolve(config["post_de"]["sets"], config)
    inputs = base_inputs(plan, plan_path) + dea_inputs(config, root / "dea")
    parameters = {**common_parameters(config, plan, {k: v for k, v in decision.items() if k != "resolved"}), **resolved}
    return request(CAPABILITY, stage_id, plan, run_id=run_id, output_temp_dir=output_temp_dir, inputs=inputs, parameters=parameters)


def execute(request_value: dict) -> dict:
    return execute_r(request_value)


def required_refusal(decision: dict):
    """A required sets module whose declared set or pair is invalid rejects the plan (explicit strictness, D-59)."""
    sub = [s for s in decision.get("subanalyses") or [] if s["reason_code"] != "E_VENN_K"]
    return (sub[0]["reason_code"], sub[0]["reason"]) if sub else None
