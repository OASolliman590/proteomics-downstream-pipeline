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


def plan_checks(config: dict, observations: list[dict], feature_ids: list[str]) -> list[tuple[str, str]]:
    """Return scientific refusals (code, message); configuration errors raise immediately."""
    settings = resolved_settings(config)
    design = config["design"]
    if config["runtime"]["execution_profile"] == "production" and settings["permutations"] < 999:
        raise PermanovaRefusal("E_PERMANOVA_RESOLUTION", "production PERMANOVA requires at least 999 permutations", "/multivariate/permanova/permutations")
    ids = [s["id"] for s in settings["feature_sets"]]
    if len(set(ids)) != len(ids):
        raise PermanovaRefusal("E_ID_DUPLICATE", "feature set ids must be unique", "/multivariate/permanova/feature_sets")
    known = set(feature_ids)
    model_ids = {m["id"] for m in config["models"]}; contrast_ids = {c["id"] for c in config["contrasts"]}
    for index, item in enumerate(settings["feature_sets"]):
        if item["kind"] == "declared_panel":
            unknown = [f for f in item["feature_ids"] if f not in known]
            if unknown:
                raise PermanovaRefusal("E_PERMANOVA_FEATURE_SET", f"panel {item['id']!r} names features not in the primary matrix: {unknown[:5]}", f"/multivariate/permanova/feature_sets/{index}")
        if item["kind"] == "dep_derived" and (item["model_id"] not in model_ids or item["contrast_id"] not in contrast_ids):
            raise PermanovaRefusal("E_REFERENCE_UNKNOWN", f"DEP-derived set {item['id']!r} references an undeclared model or contrast", f"/multivariate/permanova/feature_sets/{index}")
    refusals = []
    blocked = design["blocking"]["mode"] != "none"
    if settings["covariates"] and blocked:
        refusals.append(("E_PERMANOVA_DESIGN_UNSUPPORTED", "covariate PERMANOVA is not supported for subject-blocked designs"))
    group = design["group_column"]
    if blocked:   # audit 2026-10-02: the group term needs a valid exchangeability scheme (see pm_group_scheme in permanova.R)
        subject = design["blocking"]["subject_column"]
        levels = set(design["group_levels"])
        by_subject: dict[str, list[str]] = {}
        for row in observations:
            if row.get(group) in levels:
                by_subject.setdefault(row.get(subject, ""), []).append(row[group])
        distinct = {k: len(set(v)) for k, v in by_subject.items()}
        if distinct and all(n == 1 for n in distinct.values()):
            if len({len(v) for v in by_subject.values()}) != 1:
                refusals.append(("E_PERMANOVA_BLOCKING_UNBALANCED", f"group is constant within each {subject} but subjects have unequal numbers of observations; whole-subject permutation needs balanced subjects"))
        elif distinct and not all(n > 1 for n in distinct.values()):
            refusals.append(("E_PERMANOVA_BLOCKING_MIXED", f"group is constant within some {subject} values and varies within others; no single valid permutation scheme exists"))
    required = (config.get("multivariate") or {}).get("execution_requirement", "optional") == "required"
    for cov in settings["covariates"]:
        if any(cov["column"] not in row for row in observations[:1]):
            raise PermanovaRefusal("E_PERMANOVA_COVARIATE", f"covariate column {cov['column']!r} is not in the observation metadata", "/multivariate/permanova/covariates")
        if settings["interaction"] and cov["type"] == "categorical":
            levels = sorted({row[cov["column"]] for row in observations})
            empty = [f"{g}/{lvl}" for g in design["group_levels"] for lvl in levels if not any(r[group] == g and r[cov["column"]] == lvl for r in observations)]
            # audit 2026-10-02: for an optional analysis only the interaction term is refused (in R, recorded in refusals.tsv and
            # warnings) and the other results are kept; a required analysis that requests an inestimable term is still rejected.
            if empty and required:
                refusals.append(("E_PERMANOVA_INTERACTION_NONESTIMABLE", f"group x {cov['column']} has empty cell(s) {empty}; the interaction is not estimable"))
    return refusals


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, dea_dir: str | Path | None = None,
                  stage_id: str = "permanova") -> dict:
    plan_path = Path(plan_path).resolve(); root = plan_path.parent
    inputs = [{"artifact_id": "plan", "path": str(plan_path), "sha256": sha256_file(plan_path)}]
    for artifact_id in ("primary_matrix", "primary_observed_mask", "primary_observations"):
        artifact = next(a for a in plan["artifacts"] if a["artifact_id"] == artifact_id)
        inputs.append({"artifact_id": artifact_id, "path": str((root / artifact["relative_path"]).resolve()), "sha256": artifact["sha256"]})
    settings = resolved_settings(config)
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
