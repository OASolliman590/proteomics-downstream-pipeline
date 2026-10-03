"""Typed ReportData assembler (packet R10a, SM25, V091/V094).

Only values from stage outputs whose SHA-256 manifests verify are
included.  A stage that did not run, failed or cannot be verified is shown
with its actual state and reason; its numbers are absent (null), never
replaced by zeros or healthy defaults.
"""
from __future__ import annotations

import csv
import importlib
import json
from pathlib import Path

from ..provenance import sha256_file

SCHEMA_VERSION = "1.2.0"
# Fixed optional section providers (same discipline as the R01 capability map).
# Phase 2 providers are registered here by Maintainer amendment A-2026-10-01-05 (fixed map; absent modules are NOT_RUN).
SECTION_PROVIDERS = {"permanova": ("proteomics_pipeline.reporting.permanova_data", "permanova_section"),
                     "response": ("proteomics_pipeline.reporting.response_data", "response_section"),
                     "pathways": ("proteomics_pipeline.reporting.pathway_data", "pathway_section"),
                     "assay_engines": ("proteomics_pipeline.reporting.model_data", "model_section"),
                     # Amendment A-2026-10-01-14 (ADR 0009): post-DE module sections (R14a-R14e) and the integrated eligibility section (R14f)
                     "post_de_sets": ("proteomics_pipeline.reporting.post_de_sets_data", "post_de_sets_section"),
                     "post_de_sensitivity": ("proteomics_pipeline.reporting.post_de_sensitivity_data", "post_de_sensitivity_section"),
                     "post_de_association": ("proteomics_pipeline.reporting.post_de_association_data", "post_de_association_section"),
                     "post_de_biomarker": ("proteomics_pipeline.reporting.post_de_biomarker_data", "post_de_biomarker_section"),
                     "post_de_networks": ("proteomics_pipeline.reporting.post_de_networks_data", "post_de_networks_section"),
                     "post_de_eligibility": ("proteomics_pipeline.reporting.post_de_data", "post_de_section")}
PHASE_LABEL = "Phase 1 thin offline report (R10a). This is not the full R10b report and not a v1.0-defensible release."


def _table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify_stage(root: Path, stage: dict) -> tuple[bool, str | None]:
    """Recompute every output hash of a stage result; return (verified, problem)."""
    if not stage.get("result_path"):
        return False, "no stage result"
    result_path = root / stage["result_path"]
    if not result_path.is_file():
        return False, "stage result missing"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    for output in result.get("outputs", []):
        target = result_path.parent / output["relative_path"]
        if not target.is_file():
            return False, f"missing output {output['relative_path']}"
        if sha256_file(target) != output["sha256"]:
            return False, f"hash mismatch {output['relative_path']}"
    if result.get("state") != stage.get("state"):
        return False, "stage state differs from its result file"
    return True, None


def _stage_map(snapshot: dict, root: Path) -> dict[str, dict]:
    stages = {}
    for stage in snapshot["stages"]:
        verified, problem = verify_stage(root, stage) if stage["state"] == "COMPLETED" else (False, None)
        entry = {k: stage.get(k) for k in ("stage_id", "capability", "required", "state", "reason_code", "message", "result_path")}
        entry["verified"] = verified
        entry["verification_problem"] = problem
        if stage["state"] == "COMPLETED" and not verified:
            entry["display_state"] = "UNVERIFIED"
        else:
            entry["display_state"] = stage["state"]
        stages[stage["stage_id"]] = entry
    return stages


def _usable(stages: dict, stage_id: str) -> bool:
    stage = stages.get(stage_id)
    return bool(stage and stage["state"] == "COMPLETED" and stage["verified"])


def _relative(root: Path, stages: dict, stage_id: str) -> Path:
    return (root / stages[stage_id]["result_path"]).parent


def inputs_section(root: Path, stages: dict) -> dict:
    if not _usable(stages, "intake"):
        return {"state": stages.get("intake", {}).get("display_state", "NOT_RUN"), "values": None}
    directory = _relative(root, stages, "intake")
    report = json.loads((directory / "intake_report.json").read_text(encoding="utf-8"))
    observations = _table(directory / "observations.tsv")
    lineage = _table(directory / "observation_lineage.tsv")
    prefix = directory.relative_to(root).as_posix()
    return {"state": "COMPLETED", "source": f"{prefix}/intake_report.json",
            "values": {k: report[k] for k in ("grain", "assay", "input_format", "scale", "n_features", "n_source_observations", "n_canonical_observations", "n_subjects",
                                               "technical_replicates", "groups", "missingness", "declared", "provenance_gaps", "asserted_unknown_facts", "explicit_zero_missing_cells")},
            "canonical_observation_ids": [row["observation_id"] for row in observations],
            "lineage_source": f"{prefix}/observation_lineage.tsv", "n_lineage_rows": len(lineage),
            "sources": {"sample_n": f"{prefix}/sample_n.tsv", "observations": f"{prefix}/observations.tsv", "missingness": f"{prefix}/missingness_by_observation.tsv"}}


def qc_section(root: Path, stages: dict) -> dict:
    if not _usable(stages, "preprocessing"):
        return {"state": stages.get("preprocessing", {}).get("display_state", "NOT_RUN"), "values": None}
    from .qc_data import qc_section as build
    directory = _relative(root, stages, "preprocessing")
    section = build(directory, prefix=directory.relative_to(root).as_posix())
    section["state"] = "COMPLETED"
    return section


def design_section(root: Path, stages: dict) -> dict:
    if "design" not in stages:
        return {"state": "NOT_REQUESTED", "values": None}
    if not _usable(stages, "design"):
        return {"state": stages["design"]["display_state"], "values": None}
    directory = _relative(root, stages, "design")
    diagnostics = json.loads((directory / "design_diagnostics.json").read_text(encoding="utf-8"))
    prefix = directory.relative_to(root).as_posix()
    return {"state": "COMPLETED", "source": f"{prefix}/design_diagnostics.json",
            "designs": [{k: d[k] for k in ("design_id", "n_observations", "n_coefficients", "rank", "full_rank", "coefficients", "blocking_mode")} for d in diagnostics["designs"]],
            "contrasts": [{k: c[k] for k in ("contrast_id", "design_id", "role", "required_groups", "estimable", "units_by_required_group", "subjects_by_required_group")} for c in diagnostics["contrasts"]],
            "blocking": diagnostics["blocking"]}


def dea_section(root: Path, stages: dict, plan: dict | None) -> dict:
    if "limma" not in stages:
        return {"state": "NOT_REQUESTED" if plan and plan.get("scope") == "qc_only" else "NOT_RUN", "values": None,
                "note": "No differential analysis was executed; this is not evidence of no differences."}
    if not _usable(stages, "limma"):
        return {"state": stages["limma"]["display_state"], "reason_code": stages["limma"]["reason_code"], "values": None,
                "note": "The required differential stage did not complete; no differential results are reported."}
    directory = _relative(root, stages, "limma")
    prefix = directory.relative_to(root).as_posix()
    families = _table(directory / "families.tsv")
    cutoff = {f["family_id"]: float(f["q_cutoff"]) for f in families}
    tables = {}
    contrasts = {}
    for name, result_type in (("zero_null", "ProteinZeroNullResult"), ("treat", "ProteinTreatResult"), ("omnibus", "ProteinOmnibusResult")):
        path = directory / f"{name}.tsv"
        if not path.is_file():
            continue
        rows = _table(path)
        tables[name] = {"source": f"{prefix}/{name}.tsv", "n_rows": len(rows), "result_type": result_type}
        for row in rows:
            key = (row["model_id"], row["contrast_id"], row["hypothesis_type"])
            entry = contrasts.setdefault(key, {"model_id": row["model_id"], "contrast_id": row["contrast_id"], "hypothesis_type": row["hypothesis_type"],
                                               "role": row["role"], "family_id": row["family_id"] if row["family_id"] != "NA" else None,
                                               "n_planned": 0, "n_tested": 0, "n_excluded": 0, "n_nonestimable": 0, "n_numerical_failure": 0, "n_q_at_or_below_cutoff": 0})
            entry["n_planned"] += 1
            state = row["eligibility"]
            entry[{"tested": "n_tested", "excluded": "n_excluded", "nonestimable": "n_nonestimable", "numerical_failure": "n_numerical_failure"}.get(state, "n_excluded")] += 1
            if state == "tested" and row["q_value"] not in ("NA", "") and entry["family_id"] in cutoff and float(row["q_value"]) <= cutoff[entry["family_id"]]:
                entry["n_q_at_or_below_cutoff"] += 1
    status = json.loads((directory / "model_status.json").read_text(encoding="utf-8"))
    return {"state": "COMPLETED", "families": families, "families_source": f"{prefix}/families.tsv", "tables": tables,
            "contrast_summaries": list(contrasts.values()), "model_status": status["models"], "settings": status["settings"],
            "influence_sources": sorted(f"{prefix}/influence/{p.name}" for p in (directory / "influence").glob("*.tsv")) if (directory / "influence").is_dir() else [],
            "note": "Counts of q-values at or below each family's prespecified cutoff are reported for all planned contrasts, including secondary ones; zero is a real completed result, not evidence of no biological effect."}


def optional_sections(root: Path, stages: dict) -> dict:
    out = {}
    for name, (module_name, function_name) in SECTION_PROVIDERS.items():
        if name not in stages:
            continue
        if not _usable(stages, name):
            out[name] = {"state": stages[name]["display_state"], "reason_code": stages[name]["reason_code"], "values": None}
            continue
        try:
            module = importlib.import_module(module_name)
            provider = getattr(module, function_name)
        except (ImportError, AttributeError):
            out[name] = {"state": "NOT_RUN", "reason_code": "E_REPORT_SECTION_NOT_IMPLEMENTED", "values": None}
            continue
        directory = _relative(root, stages, name)
        out[name] = provider(directory, prefix=directory.relative_to(root).as_posix())
    return out


def limitations(inputs: dict, plan: dict | None) -> list[str]:
    items = ["Phase 1 scope: QC and limma differential abundance only; pathway, response, alternative engines and full publication reporting are not part of this report.",
             "Statistical significance is reported against prespecified families; a non-significant result is not evidence of equivalence or of no biological effect.",
             "Results describe the declared protein-level input; upstream identification and quantification were not re-evaluated."]
    if inputs.get("values"):
        for gap in inputs["values"]["provenance_gaps"]:
            items.append(f"Unknown: {gap['field']} — {gap['consequence']}.")
    return items


def assemble(run_root: str | Path, snapshot: dict) -> dict:
    root = Path(run_root)
    stages = _stage_map(snapshot, root)
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8")) if plan_path.is_file() else None
    inputs = inputs_section(root, stages)
    warnings_path = root / "warnings.json"
    data = {
        "schema_version": SCHEMA_VERSION,
        "result_type": "ReportData",
        "run": {"run_id": snapshot["run_id"], "plan_hash": snapshot["plan_hash"], "state": snapshot["state"], "exit_code": snapshot["exit_code"],
                "reason_code": snapshot["reason_code"], "requested_phase": snapshot["requested_phase"], "phase_label": PHASE_LABEL,
                "implemented_capabilities": snapshot["implemented_capabilities"], "scope": plan["scope"] if plan else None},
        "stages": list(stages.values()),
        "inputs": inputs,
        "qc": qc_section(root, stages),
        "design": design_section(root, stages),
        "dea": dea_section(root, stages, plan),
        "sections": optional_sections(root, stages),
        "plan": None if plan is None else {"plan_hash": plan["plan_hash"], "source": "plan.json", "families": plan["families"], "models": [
            {k: m.get(k) for k in ("model_id", "design_id", "engine", "role", "execution_requirement")} for m in plan["models"]], "engine_eligibility": plan["engine_eligibility"],
            "environment": plan["environment"], "code_sha256": plan["code_sha256"]},
        "warnings": json.loads(warnings_path.read_text(encoding="utf-8")) if warnings_path.is_file() else [],
        "limitations": limitations(inputs, plan),
    }
    return data


def validate_report_data(data: dict) -> None:
    from jsonschema import Draft202012Validator
    schema = json.loads((Path(__file__).resolve().parents[1] / "schemas" / "report-data.schema.json").read_text(encoding="utf-8"))
    errors = sorted(Draft202012Validator(schema).iter_errors(data), key=lambda e: list(e.absolute_path))
    if errors:
        from ..errors import ProteomicsError
        raise ProteomicsError("E_REPORT_DATA_SCHEMA", f"{errors[0].message} at /{'/'.join(map(str, errors[0].absolute_path))}", exit_code=4)
