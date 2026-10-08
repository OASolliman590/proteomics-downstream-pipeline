"""Post-DE eligibility and dependency report (packet R14f, SM41, FR-166).

A Python stage that runs after the declared post-DE modules.  It reads the frozen plan (the planner's eligibility
decision for every requested module, recorded before any computation), the actual stage states and each completed
module's verified ``eligibility.json``, and writes ``post_de/eligibility/eligibility.json`` (per module: requiredness,
plan eligibility, actual state, typed reason, claim label, sub-analyses and input hashes) and
``dependency_summary.tsv`` (each primary discovery linked to its set, robustness, sensitivity, association,
biomarker and network outcomes).  NOT_RUN (missing software or implementation) and INAPPLICABLE (scientifically
ineligible) are kept distinct; a refused module has no values.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from ..provenance import sha256_file

CAPABILITY = "post_de_eligibility"
MODULES = ("sets", "sensitivity", "association", "biomarker", "networks")
STATE_MEANING = {"COMPLETED": "executed; results shown with their claim labels", "INAPPLICABLE": "scientifically ineligible for this design or data (typed reason); no values",
                 "NOT_RUN": "not executed because software or an implementation is unavailable; not a scientific refusal; no values",
                 "FAILED": "executed and failed (typed reason); no values", "NOT_REQUESTED": "not declared"}


def capabilities():
    return [{"id": CAPABILITY, "dependencies": [], "required_r_packages": []}]


def _table(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _verified_dir(root: Path, stage: dict) -> Path | None:
    from ..reporting.assemble import verify_stage
    if stage.get("state") != "COMPLETED":
        return None
    ok, _ = verify_stage(root, stage)
    return (root / stage["result_path"]).parent if ok else None


def records(root: Path, stages: list[dict]) -> list[dict]:
    plan = json.loads((root / "plan.json").read_text(encoding="utf-8"))
    by_stage = {s["stage_id"]: s for s in stages}
    out = []
    for entry in plan["capability_plan"]:
        capability = entry["capability"]
        if not capability.startswith("post_de_") or capability == CAPABILITY:
            continue
        module = capability[len("post_de_"):]
        decision = entry.get("post_de") or {}
        stage = by_stage.get(capability, {"state": "NOT_RUN", "reason_code": "E_CAPABILITY_NOT_IMPLEMENTED", "message": "no stage was executed"})
        directory = _verified_dir(root, stage)
        module_record = json.loads((directory / "eligibility.json").read_text(encoding="utf-8")) if directory and (directory / "eligibility.json").is_file() else {}
        out.append({"module": module, "capability": capability, "required": bool(entry.get("required")), "implemented": bool(entry.get("implemented")),
                    "plan_eligibility": decision.get("state", "ELIGIBLE" if entry.get("scientific_eligibility") == "eligible" else "INAPPLICABLE"),
                    "plan_reason_code": decision.get("reason_code"), "plan_reason": decision.get("reason"), "plan_reasons": decision.get("reasons", []),
                    "eligibility_rule": decision.get("rule"), "subanalyses": decision.get("subanalyses", []),
                    "adaptations": decision.get("adaptations", []),   # D-59: requested vs used, with the reason
                    "state": stage.get("state"), "reason_code": stage.get("reason_code"), "reason": stage.get("message") or decision.get("reason"),
                    "state_meaning": STATE_MEANING.get(stage.get("state"), ""),
                    "claim_label": module_record.get("claim_label") if stage.get("state") == "COMPLETED" else None,
                    "inference": module_record.get("inference") if stage.get("state") == "COMPLETED" else None,
                    "refusals": module_record.get("refusals", []), "input_hashes": module_record.get("input_hashes", []),
                    "outputs": f"post_de/{module}" if directory else None})
    return out


def dependency_summary(root: Path, stages: list[dict]) -> list[dict]:
    """Each primary discovery (primary family, q <= its cutoff) with what it depends on in the post-DE results."""
    by_stage = {s["stage_id"]: s for s in stages}
    dirs = {m: _verified_dir(root, by_stage.get(f"post_de_{m}", {})) for m in MODULES}
    families = {f["family_id"]: float(f["q_cutoff"]) for f in _table(root / "dea" / "families.tsv") if f.get("role") == "primary"}
    discoveries = [r for r in _table(root / "dea" / "zero_null.tsv") if r["family_id"] in families and r["eligibility"] == "tested" and r["q_value"] not in ("NA", "")
                   and float(r["q_value"]) <= families[r["family_id"]]]
    membership = {r["feature_id"]: r for r in _table(dirs["sets"] / "membership.tsv")} if dirs["sets"] else {}
    robust = {(r["contrast_id"], r["feature_id"]): r for r in _table(dirs["sensitivity"] / "robustness_summary.tsv")} if dirs["sensitivity"] else {}
    association = {}
    if dirs["association"]:
        for path in sorted(dirs["association"].glob("association_*.tsv")):
            if path.name == "association_summary.tsv":
                continue
            for r in _table(path):
                if r["eligibility"] == "tested" and r["q_value"] not in ("NA", "") and float(r["q_value"]) < 0.05:
                    association.setdefault(r["feature_id"], []).append(f"{r['phenotype']}:{r['scope']}{'/' + r['group'] if r['group'] not in ('NA', '') else ''}")
    selection = {r["feature_id"]: r["selection_frequency"] for r in _table(dirs["biomarker"] / "selection_stability.tsv")} if dirs["biomarker"] else {}
    single = {r["feature_id"]: r["auc"] for r in _table(dirs["biomarker"] / "single_feature_auc.tsv")} if dirs["biomarker"] else {}
    modules = {r["feature_id"]: r["module"] for r in _table(dirs["networks"] / "modules.tsv")} if dirs["networks"] else {}
    state = lambda m: (by_stage.get(f"post_de_{m}") or {}).get("state", "NOT_REQUESTED")
    rows = []
    for d in discoveries:
        f, c = d["feature_id"], d["contrast_id"]
        sets = [k for k, v in membership.get(f, {}).items() if v == "true"]
        r = robust.get((c, f), {})
        rows.append({"contrast_id": c, "feature_id": f, "family_id": d["family_id"], "effect": d["effect"], "q_value": d["q_value"],
                     "sets": ";".join(sets) if dirs["sets"] else state("sets"),
                     "robustness_fraction": r.get("robustness_fraction", state("sensitivity") if not dirs["sensitivity"] else "NA"),
                     "sign_stability": r.get("sign_stability", "NA"), "influence_flag": r.get("influence_flag", "NA"),
                     "sensitivity_model_fraction": r.get("sensitivity_model_fraction", "NA"), "matched_n_selection_frequency": r.get("matched_n_selection_frequency", "NA"),
                     "phenotype_associations": ";".join(association.get(f, [])) if dirs["association"] else state("association"),
                     "biomarker_selection_frequency": selection.get(f, "NA") if dirs["biomarker"] else state("biomarker"),
                     "single_feature_auc_in_sample": single.get(f, "NA") if dirs["biomarker"] else state("biomarker"),
                     "coabundance_module": modules.get(f, "NA") if dirs["networks"] else state("networks"),
                     "claim_label": "descriptive"})
    return rows


def execute(request: dict) -> dict:
    from ..errors import ProteomicsError
    from ..runtime import stage_result, utc_now, validate_stage_request
    validate_stage_request(request)
    out = Path(request["output_temp_dir"]); out.mkdir(parents=True, exist_ok=False)
    started = utc_now()
    root = Path(request["parameters"]["run_root"]); stages = request["parameters"]["stages"]
    try:
        recs = records(root, stages)
        deps = dependency_summary(root, stages)
    except (OSError, KeyError, ValueError, ProteomicsError) as error:
        result = stage_result(request["run_id"], request["stage_id"], CAPABILITY, "FAILED", plan_hash=request["plan_hash"], exit_code=4, reason_code="E_POST_DE_REPORT",
                              message=f"{type(error).__name__}: {error}")
        result["started_at"] = started; return result
    (out / "eligibility.json").write_bytes((json.dumps({"result_type": "PostDeEligibility", "state_meanings": STATE_MEANING, "modules": recs}, indent=2, ensure_ascii=False) + "\n").encode("utf-8"))
    header = ["contrast_id", "feature_id", "family_id", "effect", "q_value", "sets", "robustness_fraction", "sign_stability", "influence_flag", "sensitivity_model_fraction",
              "matched_n_selection_frequency", "phenotype_associations", "biomarker_selection_frequency", "single_feature_auc_in_sample", "coabundance_module", "claim_label"]
    with (out / "dependency_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(deps)
    result = stage_result(request["run_id"], request["stage_id"], CAPABILITY, "COMPLETED", plan_hash=request["plan_hash"], exit_code=0,
                          message=f"post-DE eligibility for {len(recs)} module(s); {len(deps)} primary discovery row(s) in the dependency summary")
    result["started_at"] = started
    result["outputs"] = [{"artifact_id": "post_de_eligibility", "relative_path": "eligibility.json", "sha256": sha256_file(out / "eligibility.json"), "result_type": "PostDeEligibility"},
                         {"artifact_id": "post_de_dependency_summary", "relative_path": "dependency_summary.tsv", "sha256": sha256_file(out / "dependency_summary.tsv"), "result_type": "DependencySummary"}]
    return result
