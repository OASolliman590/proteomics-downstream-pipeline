"""Multivariate PERMANOVA / PERMDISP stage (packet R13, SM26-SM30, operator-authorized amendment).

Builds the stage request from the frozen plan (primary matrix, observed
mask, observations) and, for DEP-derived feature sets, the completed R05
differential table.  Plan-time checks refuse requests that cannot be
answered (insufficient permutation resolution, covariates in subject-blocked
designs, mixed or unbalanced subject blocking, unknown panel members, and
empty interaction cells when the analysis is required; for an optional
analysis only the interaction term is refused, in R).
"""
from __future__ import annotations

import csv
import json
import os
import shutil
from pathlib import Path

from .errors import ProteomicsError
from .provenance import sha256_file

CAPABILITY = "permanova"
DEFAULTS = {"distance": {"metric": "euclidean", "scaling": "zscore"}, "permutations": 9999, "pairwise": True, "adjustment": "holm", "alpha": 0.05,
            "covariates": [], "interaction": False, "feature_sets": [{"id": "all", "kind": "all_complete"}], "random_sets": 2000}


class PermanovaRefusal(ProteomicsError):
    def __init__(self, code: str, message: str, pointer: str | None = "/multivariate/permanova"):
        super().__init__(code, message, pointer, exit_code=2)


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["design"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "vegan", "permute"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def resolved_settings(config: dict) -> dict:
    block = dict(DEFAULTS)
    block.update((config.get("multivariate") or {}).get("permanova") or {})
    block["distance"] = {**DEFAULTS["distance"], **block.get("distance", {})}
    block.setdefault("seed", config["runtime"]["seed"])
    return block


def adapted(config: dict, observations: list[dict], feature_ids: list[str]) -> tuple[dict, list[dict], list[tuple[str, str]]]:
    """D-59 (Maintainer direction 2026-10-05): PERMANOVA adapts to the declaration and the design instead of holding up the
    run.  Returns the settings actually used, the adaptations (requested vs used, reason) and the refusals that remain:
    only outputs that would be invalid, scoped as narrowly as possible."""
    settings = resolved_settings(config)
    design = config["design"]; group = design["group_column"]
    adaptations: list[dict] = []; refusals: list[tuple[str, str]] = []
    def adapt(item, requested, used, reason):
        adaptations.append({"analysis": "permanova", "item": item, "requested": requested, "used": used, "reason": reason})
    if config["runtime"]["execution_profile"] == "production" and settings["permutations"] < 999:
        adapt("permutations", settings["permutations"], 999, "production PERMANOVA uses at least 999 permutations")
        settings["permutations"] = 999
    known = set(feature_ids) if feature_ids is not None else None
    model_ids = {m["id"] for m in config["models"]}; contrast_ids = {c["id"] for c in config["contrasts"]}
    kept, seen = [], set()
    for item in settings["feature_sets"]:
        problem = None
        if item["id"] in seen:
            problem = "E_ID_DUPLICATE", f"feature set id {item['id']!r} is declared twice; the later declaration is not tested"
        elif item["kind"] == "declared_panel" and known is not None and [f for f in item["feature_ids"] if f not in known]:
            problem = "E_PERMANOVA_FEATURE_SET", f"panel {item['id']!r} names features not in the primary matrix: {[f for f in item['feature_ids'] if f not in known][:5]}"
        elif item["kind"] == "dep_derived" and (item["model_id"] not in model_ids or item["contrast_id"] not in contrast_ids):
            problem = "E_REFERENCE_UNKNOWN", f"DEP-derived set {item['id']!r} references an undeclared model or contrast"
        if problem:
            adapt(f"feature_set:{item['id']}", "tested", "not tested", f"{problem[0]}: {problem[1]}")
            continue
        seen.add(item["id"]); kept.append(item)
    settings["feature_sets"] = kept
    if not kept:
        refusals.append(("E_PERMANOVA_FEATURE_SET", "no declared feature set can be tested (see the adaptations)"))
    for cov in list(settings["covariates"]):
        if any(cov["column"] not in row for row in observations[:1]):
            adapt(f"covariate:{cov['column']}", "adjusted", "dropped", f"E_PERMANOVA_COVARIATE: column {cov['column']!r} is not in the observation metadata")
            settings["covariates"] = [c for c in settings["covariates"] if c is not cov]
    blocked = design["blocking"]["mode"] != "none"
    settings["row_adaptation"] = "none"
    if blocked:
        if settings["covariates"]:
            adapt("covariates", [c["column"] for c in settings["covariates"]], [], "E_PERMANOVA_DESIGN_UNSUPPORTED: covariate terms have no valid restricted permutation in a "
                  "subject-blocked design; the group term is tested with the blocked scheme and the covariate terms are refused")
            settings["covariates"] = []; settings["interaction"] = False
        subject = design["blocking"]["subject_column"]
        levels = set(design["group_levels"])
        by_subject: dict[str, list[str]] = {}
        for row in observations:
            if row.get(group) in levels:
                by_subject.setdefault(row.get(subject, ""), []).append(row[group])
        distinct = {k: len(set(v)) for k, v in by_subject.items()}
        if distinct and all(n == 1 for n in distinct.values()) and len({len(v) for v in by_subject.values()}) != 1:
            settings["row_adaptation"] = "subject_means"
            adapt("rows", "observations with whole-subject permutation", "subject means with unrestricted permutation of subjects",
                  "E_PERMANOVA_BLOCKING_UNBALANCED: subjects have unequal numbers of observations, so whole-subject row permutation is not defined; "
                  "each subject contributes its mean profile (one row per biological unit)")
        elif distinct and any(n == 1 for n in distinct.values()) and any(n > 1 for n in distinct.values()):
            varying = sorted(k for k, n in distinct.items() if n > 1)
            settings["row_adaptation"] = "varying_subjects_only"; settings["varying_subjects"] = varying
            adapt("rows", "all subjects", f"{len(varying)} subjects observed in more than one group",
                  "E_PERMANOVA_BLOCKING_MIXED: group is constant within some subjects and varies within others, so no single permutation scheme covers all; "
                  "the within-subject test uses the subjects observed in more than one group")
    required = (config.get("multivariate") or {}).get("execution_requirement", "optional") == "required"
    dropped = [a for a in adaptations if a["item"].startswith(("feature_set:", "covariate")) and a["used"] in ("not tested", "dropped", [])]
    if required and dropped:   # explicit strictness: a required analysis whose declared sets or terms cannot be tested rejects the plan
        code = dropped[0]["reason"].split(":", 1)[0]
        refusals.append((code, "required PERMANOVA: " + "; ".join(a["reason"] for a in dropped)))
    for cov in settings["covariates"]:
        if settings["interaction"] and cov["type"] == "categorical":
            levels_c = sorted({row[cov["column"]] for row in observations})
            empty = [f"{g}/{lvl}" for g in design["group_levels"] for lvl in levels_c if not any(r[group] == g and r[cov["column"]] == lvl for r in observations)]
            # audit 2026-10-02: for an optional analysis only the interaction term is refused (in R, recorded in refusals.tsv and
            # warnings) and the other results are kept; a required analysis that requests an inestimable term is still rejected.
            if empty and required:
                refusals.append(("E_PERMANOVA_INTERACTION_NONESTIMABLE", f"group x {cov['column']} has empty cell(s) {empty}; the interaction is not estimable"))
    return settings, adaptations, refusals


def plan_checks(config: dict, observations: list[dict], feature_ids: list[str]) -> list[tuple[str, str]]:
    """Return the refusals that remain after adaptation (D-59)."""
    return adapted(config, observations, feature_ids)[2]


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, dea_dir: str | Path | None = None,
                  stage_id: str = "permanova") -> dict:
    plan_path = Path(plan_path).resolve(); root = plan_path.parent
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]
    for artifact_id in ("primary_matrix", "primary_observed_mask", "primary_observations"):
        artifact = next(a for a in plan["artifacts"] if a["artifact_id"] == artifact_id)
        inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]})
    # the plan-time adaptation is recomputed from the run's frozen preprocessing outputs (a request built outside a run
    # directory, as in the service-layer guard tests, has none: membership checks are then left to the R stage)
    observations, feature_ids = read_observations(root / "preprocessing") if (root / "preprocessing" / "primary" / "observations.tsv").is_file() else ([], None)
    settings, adaptations, _ = adapted(config, observations, feature_ids)
    if any(s["kind"] == "dep_derived" for s in settings["feature_sets"]):
        if dea_dir is None:
            raise PermanovaRefusal("E_PERMANOVA_FEATURE_SET", "DEP-derived feature sets need the completed differential stage")
        from . import inference_service
        inference_service.guard_downstream(config, dea_dir)   # audit 2026-10-02: DEP-derived sets only from observed-data primary results
        result = json.loads((Path(dea_dir) / "stage-result.json").read_text(encoding="utf-8"))
        output = next(o for o in result["outputs"] if o["relative_path"] == "zero_null.tsv")
        path = Path(dea_dir) / "zero_null.tsv"
        if sha256_file(path) != output["sha256"]:
            raise ProteomicsError("E_INTEGRITY", "differential table changed after the limma stage completed", exit_code=5)
        inputs.append({"artifact_id": "dea_zero_null", "path": str(path.resolve()), "sha256": output["sha256"]})
    design = config["design"]
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan["plan_hash"], "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": None,
            "parameters": {"group_column": design["group_column"], "group_levels": list(design["group_levels"]),
                           "subject_column": design["blocking"].get("subject_column") if design["blocking"]["mode"] != "none" else None,
                           "metric": settings["distance"]["metric"], "scaling": settings["distance"]["scaling"], "permutations": int(settings["permutations"]),
                           "seed": int(settings["seed"]), "alpha": float(settings["alpha"]), "pairwise": bool(settings["pairwise"]), "adjustment": settings["adjustment"],
                           "covariates": list(settings["covariates"]), "interaction": bool(settings["interaction"]), "feature_sets": list(settings["feature_sets"]),
                           "random_sets": int(settings["random_sets"]), "execution_profile": config["runtime"]["execution_profile"],
                           "adaptations": adaptations, "row_adaptation": settings["row_adaptation"], "varying_subjects": settings.get("varying_subjects", []),
                           "figure_formats": list(config["report"].get("figure_formats", []))},
            "rng": {"seed": int(settings["seed"]), "kind": "L'Ecuyer-CMRG", "threads": int(config["runtime"]["threads"])}}


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())


def read_observations(preprocessing_dir: str | Path) -> tuple[list[dict], list[str]]:
    base = Path(preprocessing_dir) / "primary"
    with (base / "observations.tsv").open(encoding="utf-8", newline="") as handle:
        observations = list(csv.DictReader(handle, delimiter="\t"))
    with (base / "features.tsv").open(encoding="utf-8", newline="") as handle:
        features = [row["feature_id"] for row in csv.DictReader(handle, delimiter="\t")]
    return observations, features
