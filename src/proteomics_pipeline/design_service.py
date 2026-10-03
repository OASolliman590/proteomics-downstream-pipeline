"""Design validation, exact contrasts, blocking and the frozen AnalysisPlan (packet R04).

Design matrices are compiled from a finite declarative grammar
(semantic-validation "Terms and coefficients"): no formula or code is ever
evaluated.  The R handler ``design_stage`` checks rank/aliasing, contrast
estimability, featurewise df and repeated-measures correlation.  The plan
is frozen (SHA-256 of canonical JSON) before any model fit (SM07, V039).
"""
from __future__ import annotations

import copy
import csv
import json
import math
import os
import shutil
from pathlib import Path

from .errors import IntegrityError, ProteomicsError
from .provenance import CONTENT_HASH_RULE, canonical_json_bytes, canonical_json_sha256, content_sha256, sha256_file

CAPABILITY = "design"
SAFE = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-")


class DesignError(ProteomicsError):
    def __init__(self, code: str, message: str, pointer: str | None = None):
        super().__init__(code, message, pointer, exit_code=2)


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["preprocessing"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def encode(text: str) -> str:
    """Percent-encode with uppercase hex; ASCII letters/digits/_/- stay literal."""
    out = []
    for char in str(text):
        if char in SAFE:
            out.append(char)
        else:
            out.extend(f"%{byte:02X}" for byte in char.encode("utf-8"))
    return "".join(out)


def decode(text: str) -> str:
    raw = bytearray(); i = 0
    while i < len(text):
        if text[i] == "%":
            raw.append(int(text[i + 1:i + 3], 16)); i += 3
        else:
            raw.extend(text[i].encode("utf-8")); i += 1
    return raw.decode("utf-8")


def _read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _nonnull(value):
    return value not in (None, "", "NA")


# --------------------------------------------------------------------------- grammar
def compile_design(design: dict, observations: list[dict], *, assay: str | None = None, tmt_strategy: str | None = None) -> dict:
    """Compile one declarative design into explicit numeric columns (V031/V032)."""
    pointer = f"/design[{design['id']}]"
    ids = [row["observation_id"] for row in observations]
    columns_present = set(observations[0]) if observations else set()
    group_column = design["group_column"]
    declared_terms = {group_column: "group"}
    for item in design["continuous_covariates"]:
        declared_terms[item["column"]] = "continuous"
    for item in design["categorical_covariates"]:
        if item["column"] in declared_terms:
            raise DesignError("E_DESIGN_TERM", f"column {item['column']!r} declared twice", pointer)
        declared_terms[item["column"]] = "categorical"
    blocking = design["blocking"]
    for column in declared_terms:
        if column not in columns_present:
            raise DesignError("E_DESIGN_TERM", f"design term {column!r} is not a declared observation metadata column; no expression is evaluated", pointer)
    for interaction in design["interactions"]:
        if len(interaction) < 2 or len(set(interaction)) != len(interaction):
            raise DesignError("E_DESIGN_TERM", f"interaction {interaction!r} must name at least two distinct declared terms", pointer)
        for term in interaction:
            if term not in declared_terms:
                raise DesignError("E_DESIGN_TERM", f"interaction term {term!r} is not a declared design term; no expression is evaluated", pointer)
    levels = list(design["group_levels"])
    if len(set(levels)) != len(levels):
        raise DesignError("E_DESIGN_LEVEL", "group_levels contains duplicates", pointer)
    for row in observations:
        if row[group_column] not in levels:
            raise DesignError("E_DESIGN_LEVEL", f"observation {row['observation_id']!r} has undeclared group level {row[group_column]!r}; observations are never silently dropped", pointer)
    empty = [level for level in levels if not any(row[group_column] == level for row in observations)]
    if empty:
        raise DesignError("E_DESIGN_GROUP_EMPTY", f"declared group level(s) {empty} have no observations; planned contrasts cannot silently lose a group", pointer)

    columns: list[tuple[str, list[float], dict]] = []

    def add(name, values, info):
        columns.append((name, [float(v) for v in values], info))

    group_cols: dict[str, tuple[str, list[float]]] = {}
    if design["intercept"]:
        add("(Intercept)", [1.0] * len(ids), {"term": "(Intercept)"})
    for index, level in enumerate(levels):
        indicator = [1.0 if row[group_column] == level else 0.0 for row in observations]
        name = f"group.{encode(level)}"
        group_cols[level] = (name, indicator)
        if design["intercept"] and index == 0:
            continue
        add(name, indicator, {"term": "group", "column": group_column, "level": level})
    centers = {}
    continuous_cols = {}
    for item in design["continuous_covariates"]:
        column = item["column"]
        raw = []
        for row in observations:
            value = row[column]
            try:
                number = float(value) if _nonnull(value) else math.nan
            except ValueError:
                number = math.nan
            if not math.isfinite(number):
                raise DesignError("E_DESIGN_COVARIATE_NONFINITE", f"required continuous covariate {column!r} is missing or nonfinite for observation {row['observation_id']!r}; complete-case removal is not silent", f"{pointer}/continuous_covariates/{column}")
            raw.append(number)
        center = math.fsum(raw) / len(raw) if item["center"] else 0.0
        centers[column] = center if item["center"] else None
        values = [v - center for v in raw]
        name = f"continuous.{encode(column)}"
        continuous_cols[column] = (name, values)
        add(name, values, {"term": "continuous", "column": column, "center": centers[column]})
    factor_cols: dict[str, list[tuple[str, list[float]]]] = {}
    for item in design["categorical_covariates"]:
        column, factor_levels, reference = item["column"], list(item["levels"]), item["reference"]
        if reference not in factor_levels:
            raise DesignError("E_DESIGN_LEVEL", f"reference level {reference!r} is not among declared levels of {column!r}", f"{pointer}/categorical_covariates/{column}/reference")
        for row in observations:
            if row[column] not in factor_levels:
                raise DesignError("E_DESIGN_LEVEL", f"observation {row['observation_id']!r} has undeclared level {row[column]!r} for {column!r}", f"{pointer}/categorical_covariates/{column}")
        factor_cols[column] = []
        for level in factor_levels:
            if level == reference:
                continue
            values = [1.0 if row[column] == level else 0.0 for row in observations]
            name = f"factor.{encode(column)}.{encode(level)}"
            factor_cols[column].append((name, values))
            add(name, values, {"term": "categorical", "column": column, "level": level, "reference": reference})
    subject_column = blocking.get("subject_column")
    if blocking["mode"] in ("fixed_subject", "duplicate_correlation"):
        if not subject_column:
            raise DesignError("E_SUBJECT_MISSING", f"blocking mode {blocking['mode']} requires subject_column", f"{pointer}/blocking")
        if subject_column not in columns_present:
            raise DesignError("E_SUBJECT_MISSING", f"subject column {subject_column!r} is absent", f"{pointer}/blocking/subject_column")
        missing = [row["observation_id"] for row in observations if not _nonnull(row[subject_column])]
        if missing:
            raise DesignError("E_SUBJECT_MISSING", f"subject IDs are missing for {missing[:5]}", f"{pointer}/blocking/subject_column")
        if subject_column in declared_terms:
            raise DesignError("E_BLOCKING_CONFLICT", f"subject column {subject_column!r} cannot also be a declared design term (fixed and random encoding of the same subject)", f"{pointer}/blocking")
    elif subject_column:
        raise DesignError("E_BLOCKING_CONFLICT", "subject_column is only allowed with a blocking mode", f"{pointer}/blocking")
    if blocking["mode"] == "fixed_subject":
        subjects = sorted({row[subject_column] for row in observations})
        for subject in subjects[1:]:
            add(f"factor.{encode(subject_column)}.{encode(subject)}", [1.0 if row[subject_column] == subject else 0.0 for row in observations],
                {"term": "fixed_subject", "column": subject_column, "level": subject, "reference": subjects[0]})
    # interactions: treatment (reference-omitted) coding for every factor in the product
    def term_columns(term):
        kind = declared_terms[term]
        if kind == "group":
            return [group_cols[level] for level in levels[1:]]
        if kind == "continuous":
            return [continuous_cols[term]]
        return factor_cols[term]
    for interaction in design["interactions"]:
        products = [("", [1.0] * len(ids))]
        for term in interaction:
            products = [(f"{name}:{cname}" if name else cname, [a * b for a, b in zip(values, cvalues)]) for name, values in products for cname, cvalues in term_columns(term)]
        for name, values in products:
            add(name, values, {"term": "interaction", "components": list(interaction)})
    if assay == "tmt" and tmt_strategy == "no_bridge" and "plex_id" not in declared_terms:
        raise DesignError("E_TMT_PLEX_TERM_REQUIRED", "TMT no-bridge inference requires plex_id as a categorical design term", pointer)
    names = [name for name, _, _ in columns]
    if len(set(names)) != len(names):
        raise DesignError("E_DESIGN_TERM", "compiled coefficient names collide", pointer)
    return {"design_id": design["id"], "observation_ids": ids, "coefficients": names,
            "matrix": [[columns[j][1][i] for j in range(len(columns))] for i in range(len(ids))],
            "term_map": [{"coefficient": name, **info} for name, _, info in columns],
            "centers": centers, "group_column": group_column, "group_levels": levels,
            "blocking": {"mode": blocking["mode"], "subject_column": subject_column,
                         "subjects": [row[subject_column] for row in observations] if subject_column else None},
            "intercept": design["intercept"]}


def compile_contrasts(config: dict, compiled: dict[str, dict]) -> list[dict]:
    out = []
    for index, contrast in enumerate(config["contrasts"]):
        pointer = f"/contrasts/{index}"
        design = compiled[contrast["design_id"]]
        unknown = [key for key in contrast["weights"] if key not in design["coefficients"]]
        if unknown:
            raise DesignError("E_CONTRAST_COEFFICIENT", f"contrast {contrast['id']!r} names unknown coefficient(s) {unknown}; valid names are {design['coefficients']}", f"{pointer}/weights")
        weights = [float(contrast["weights"].get(name, 0.0)) for name in design["coefficients"]]
        if not any(weights) or any(not math.isfinite(w) for w in weights):
            raise DesignError("E_CONTRAST_COEFFICIENT", f"contrast {contrast['id']!r} needs finite nonzero weights", f"{pointer}/weights")
        if not contrast["required_groups"]:
            raise DesignError("E_CONTRAST_GROUP", f"contrast {contrast['id']!r} needs at least one required group", f"{pointer}/required_groups")
        for group in contrast["required_groups"]:
            if group not in design["group_levels"]:
                raise DesignError("E_CONTRAST_GROUP", f"required group {group!r} is not a declared level", f"{pointer}/required_groups")
        out.append({"contrast_id": contrast["id"], "design_id": contrast["design_id"], "label": contrast["label"], "estimand": contrast["estimand"],
                    "role": contrast["role"], "required_groups": list(contrast["required_groups"]), "weights": weights})
    return out


def read_weights(model: dict, config_dir: Path, observation_ids: list[str], feature_ids: list[str]) -> dict | None:
    spec = model.get("precision_weights") or {"kind": "none"}
    if spec["kind"] == "none":
        return None
    path = Path(spec["path"]) if Path(spec["path"]).is_absolute() else config_dir / spec["path"]
    rows = _read_tsv(path)
    pointer = f"/models[{model['id']}]/precision_weights"

    def check(value, where):
        try:
            number = float(value)
        except ValueError:
            number = math.nan
        if not math.isfinite(number) or number <= 0:
            raise DesignError("E_WEIGHTS_INVALID", f"weights must be finite and positive ({where})", pointer)
        return number
    if spec["kind"] == "observation":
        by_id = {}
        for row in rows:
            if row["observation_id"] in by_id:
                raise DesignError("E_ID_DUPLICATE", f"duplicate weight for {row['observation_id']!r}", pointer)
            by_id[row["observation_id"]] = check(row["weight"], row["observation_id"])
        missing = [o for o in observation_ids if o not in by_id]
        if missing:
            raise DesignError("E_WEIGHTS_MISSING", f"weights missing for {missing[:5]}; missing entries never become 1", pointer)
        return {"kind": "observation", "values": [by_id[o] for o in observation_ids], "sha256": sha256_file(path), "file_name": path.name}
    by_feature = {row["feature_id"]: row for row in rows}
    if set(by_feature) != set(feature_ids) or any(o not in rows[0] for o in observation_ids):
        raise DesignError("E_WEIGHTS_MISSING", "matrix weights must cover every feature/observation key", pointer)
    return {"kind": "matrix", "values": [[check(by_feature[f][o], f"{f}/{o}") for o in observation_ids] for f in feature_ids], "sha256": sha256_file(path), "file_name": path.name}


def eligibility_table(config: dict, *, mask_state: str, prior_imputation: str) -> list[dict]:
    """Method-design applicability (V038): independent of hit counts and installation."""
    designs = {config["design"]["id"]: config["design"], **{d["id"]: d for d in config.get("additional_designs", [])}}
    rows = []
    for model in config["models"]:
        design = designs[model["design_id"]]
        hypotheses = [model["hypothesis"]] + list(model.get("additional_hypotheses", []))
        weighted = (model.get("precision_weights") or {"kind": "none"})["kind"] != "none"
        blocked = design["blocking"]["mode"] != "none"
        reason = None
        if model["engine"] == "limma":
            if model["coverage"]["policy"] == "native_dropout":
                reason = "E_COVERAGE_POLICY_UNSUPPORTED"
        elif model["engine"] == "deqms":
            if weighted or blocked:
                reason = "E_DEQMS_DESIGN_UNSUPPORTED"
            elif any(h != "zero_null" for h in hypotheses):
                reason = "E_DEQMS_HYPOTHESIS_UNSUPPORTED"
            elif not model.get("count_evidence"):
                reason = "E_DEQMS_COUNT_EVIDENCE"
        elif model["engine"] == "proda":
            if config["assay"] not in ("lfq_dda", "lfq_dia") or prior_imputation != "none_documented":
                reason = "E_PRODA_LFQ_UNIMPUTED_REQUIRED"
            elif weighted or blocked:
                reason = "E_PRODA_DESIGN_UNSUPPORTED"
            elif any(h != "zero_null" for h in hypotheses):
                reason = "E_PRODA_HYPOTHESIS_UNSUPPORTED"
        for hypothesis in hypotheses:
            if hypothesis == "treat" and model["engine"] == "limma" and not model.get("effect_threshold"):
                reason = reason or "E_HYPOTHESIS_UNSUPPORTED"
        if model["engine"] in ("limma",) and reason is None and mask_state != "known" and config["runtime"]["scope"] == "analysis":
            reason = "E_ORIGINAL_MASK_REQUIRED"
        rows.append({"model_id": model["id"], "engine": model["engine"], "role": model["role"], "execution_requirement": model["execution_requirement"],
                     "hypotheses": hypotheses, "eligible": reason is None, "reason_code": reason,
                     "decision_basis": "frozen adapter/design table; independent of package installation and hit counts"})
    return rows


# --------------------------------------------------------------------------- stage request
def build_request(config: dict, *, config_path: str | Path, preprocessing_dir: str | Path, run_id: str, output_temp_dir: str | Path,
                  stage_id: str = "design", plan_hash: str | None = None) -> dict:
    pre = Path(preprocessing_dir)
    observations = _read_tsv(pre / "primary" / "observations.tsv")
    feature_ids = [row["feature_id"] for row in _read_tsv(pre / "primary" / "features.tsv")]
    config_dir = Path(config_path).resolve().parent
    designs = {config["design"]["id"]: config["design"], **{d["id"]: d for d in config.get("additional_designs", [])}}
    compiled = {did: compile_design(d, observations, assay=config["assay"], tmt_strategy=config["preprocessing"].get("tmt_strategy")) for did, d in designs.items()}
    contrasts = compile_contrasts(config, compiled)
    models = []
    for model in config["models"]:
        weights = read_weights(model, config_dir, compiled[model["design_id"]]["observation_ids"], feature_ids)
        count_sha = None
        if model.get("count_evidence"):   # integration amendment A-2026-10-01-05: freeze R06 count evidence bytes in the plan
            count_path = Path(model["count_evidence"]) if Path(model["count_evidence"]).is_absolute() else config_dir / model["count_evidence"]
            if not count_path.is_file():
                raise DesignError("E_DEQMS_COUNT_EVIDENCE", f"count evidence for model {model['id']!r} is missing", f"/models[{model['id']}]/count_evidence")
            count_sha = sha256_file(count_path)
        models.append({"model_id": model["id"], "design_id": model["design_id"], "engine": model["engine"], "role": model["role"],
                       "execution_requirement": model["execution_requirement"], "weights": weights, "count_evidence_sha256": count_sha,
                       "contrast_ids": [c["contrast_id"] for c in contrasts if c["design_id"] == model["design_id"]]})
    inputs = []
    for artifact_id, relative in (("primary_matrix", "primary/matrix.tsv"), ("primary_observations", "primary/observations.tsv"),
                                  ("primary_observed_mask", "primary/observed_mask.tsv"), ("coverage", "coverage/coverage.tsv")):
        path = pre / relative
        if path.is_file():
            inputs.append({"artifact_id": artifact_id, "path": str(path.resolve()), "sha256": sha256_file(path)})
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan_hash, "inputs": inputs,
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": str(Path(config_path).resolve()),
            "parameters": {"scope": config["runtime"]["scope"], "designs": list(compiled.values()), "contrasts": contrasts, "models": models,
                           "group_column": config["design"]["group_column"], "rank_tolerance": 1e-7, "estimability_tolerance": 1e-8},
            "rng": {"seed": int(config["runtime"]["seed"]), "kind": "L'Ecuyer-CMRG", "threads": int(config["runtime"]["threads"])}}


def execute(request: dict) -> dict:
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())


def interpret(design_dir: str | Path, config: dict) -> dict:
    """Turn design diagnostics into plan-level decisions (raises DesignError)."""
    diagnostics = json.loads((Path(design_dir) / "design_diagnostics.json").read_text(encoding="utf-8"))
    for design in diagnostics["designs"]:
        if not design["full_rank"]:
            raise DesignError("E_DESIGN_CONFOUNDED", f"design {design['design_id']!r} is rank deficient (rank {design['rank']} of {design['n_coefficients']}); aliased coefficients: {design['aliased']}")
    for contrast in diagnostics["contrasts"]:
        if not contrast["estimable"]:
            raise DesignError("E_CONTRAST_NONESTIMABLE", f"contrast {contrast['contrast_id']!r} is not estimable in its design")
        if contrast["direction_check"] == "mismatch":
            raise DesignError("E_CONTRAST_DIRECTION", f"contrast {contrast['contrast_id']!r} weights put {contrast['required_groups'][0]!r} on the negative side; the first required group is the numerator")
        if config["runtime"]["scope"] == "analysis":
            small = {g: n for g, n in contrast["units_by_required_group"].items() if n < 2}
            if small:
                raise DesignError("E_INSUFFICIENT_REPLICATION", f"contrast {contrast['contrast_id']!r} has fewer than two independent biological units in {small}; descriptive QC remains available with runtime.scope=qc_only")
    return diagnostics


# --------------------------------------------------------------------------- frozen plan
def _families(config: dict, contrasts: list[dict], eligibility: list[dict]) -> list[dict]:
    """Resolve every inferential endpoint to exactly one declared family (SM12)."""
    models = {m["id"]: m for m in config["models"]}
    contrast_by_id = {c["contrast_id"]: c for c in contrasts}
    omnibus_ids = {o["id"] for o in config.get("omnibus_tests", [])}
    endpoints = {}
    for model in config["models"]:
        for hypothesis in [model["hypothesis"]] + list(model.get("additional_hypotheses", [])):
            htype = {"zero_null": "protein_zero_null", "treat": "protein_treat"}[hypothesis]
            for contrast in contrasts:
                if contrast["design_id"] == model["design_id"]:
                    endpoints[(model["id"], contrast["contrast_id"], htype)] = None
    for omnibus in config.get("omnibus_tests", []):
        endpoints[(omnibus["model_id"], omnibus["id"], "protein_omnibus")] = None
    resolved = []
    for index, family in enumerate(config["multiplicity_families"]):
        if family["hypothesis_type"] not in ("protein_zero_null", "protein_treat", "protein_omnibus"):
            continue
        members = []
        for model_id in family["model_ids"]:
            for contrast_id in family["contrast_ids"]:
                key = (model_id, contrast_id, family["hypothesis_type"])
                if key not in endpoints:
                    raise DesignError("E_FAMILY_ENDPOINT", f"family {family['id']!r} names endpoint {key} that the plan does not request", f"/multiplicity_families/{index}")
                if endpoints[key] is not None:
                    raise DesignError("E_FAMILY_OVERLAP", f"endpoint {key} is in families {endpoints[key]!r} and {family['id']!r}", f"/multiplicity_families/{index}")
                endpoints[key] = family["id"]
                members.append({"model_id": model_id, "contrast_id": contrast_id})
                model = models[model_id]
                is_primary = (model["role"] == "primary" and contrast_id not in omnibus_ids and contrast_by_id[contrast_id]["role"] == "primary"
                              and family["hypothesis_type"] == {"zero_null": "protein_zero_null", "treat": "protein_treat"}[config["primary_hypothesis"]])
                if family["role"] == "primary" and not is_primary:
                    raise DesignError("E_FAMILY_ENDPOINT", f"primary family {family['id']!r} may only contain primary-model, primary-hypothesis, primary-contrast endpoints; {key} is secondary", f"/multiplicity_families/{index}")
        resolved.append({"family_id": family["id"], "hypothesis_type": family["hypothesis_type"], "role": family["role"], "adjustment": family["adjustment"],
                         "denominator": family["denominator"], "q_cutoff": family["q_cutoff"], "members": members,
                         "dependence_assumption": "BH: independence or positive regression dependence" if family["adjustment"] == "BH" else "BY: arbitrary dependence"})
    eligible_models = {row["model_id"] for row in eligibility if row["eligible"]}
    unassigned = [key for key, family in endpoints.items() if family is None and key[0] in eligible_models]
    if unassigned:
        raise DesignError("E_FAMILY_MISSING", f"inferential endpoints without a declared family: {unassigned[:5]}; secondary contrasts must not disappear", "/multiplicity_families")
    return resolved


def plan_config(config: dict) -> dict:
    """Location-independent copy of the resolved config (paths replaced by file names)."""
    value = copy.deepcopy(config)
    for key in ("matrix", "observations", "features", "source_provenance", "original_observed_mask", "prior_imputed_mask", "mapping"):
        if value["input"].get(key) and value["input"][key] != "from_mapping":
            value["input"][key] = Path(value["input"][key]).name
    for model in value["models"]:
        weights = model.get("precision_weights") or {}
        if weights.get("path"):
            weights["path"] = Path(weights["path"]).name
        if model.get("count_evidence"):   # D-42: no location-dependent path enters the plan
            model["count_evidence"] = Path(model["count_evidence"]).name
    for resource in value.get("resources", []):
        resource["path"] = Path(resource["path"]).name
    cohort = ((value.get("post_de") or {}).get("biomarker") or {}).get("validation_cohort")
    if cohort:   # A-2026-10-01-14: declared external files enter the plan by name; their content hashes are recorded by R14d
        for key in ("matrix", "metadata"):
            if cohort.get(key):
                cohort[key] = Path(cohort[key]).name
    return value


def build_plan(config: dict, *, plan_root: str | Path, stage_dirs: dict[str, Path], sources: list[dict], environment: dict, code_sha256: str,
               capability_plan: list[dict]) -> dict:
    """Assemble the AnalysisPlan from verified planning artifacts; returns the plan with its hash."""
    plan_root = Path(plan_root).resolve()
    artifacts = []
    for stage, directory in stage_dirs.items():
        result = json.loads((Path(directory) / "stage-result.json").read_text(encoding="utf-8"))
        if result["state"] != "COMPLETED":
            raise IntegrityError(f"planning stage {stage} is not COMPLETED")
        for output in result["outputs"]:
            path = Path(directory) / output["relative_path"]
            if sha256_file(path) != output["sha256"]:
                raise IntegrityError(f"planning artifact changed: {stage}/{output['relative_path']}")
            artifacts.append({"stage": stage, "artifact_id": output["artifact_id"], "relative_path": path.resolve().relative_to(plan_root).as_posix(),
                              "sha256": output["sha256"], "content_sha256": content_sha256(path), "result_type": output["result_type"]})
    design_dir = stage_dirs.get("design")
    diagnostics = json.loads((Path(design_dir) / "design_diagnostics.json").read_text(encoding="utf-8")) if design_dir else None
    request = json.loads((Path(design_dir) / "stage-request.json").read_text(encoding="utf-8")) if design_dir else None
    contrasts = request["parameters"]["contrasts"] if request else []
    manifest = json.loads((Path(stage_dirs["intake"]) / "manifest.json").read_text(encoding="utf-8"))
    eligibility = eligibility_table(config, mask_state=manifest["original_observed_mask"]["state"], prior_imputation=config["input"]["prior_imputation"])
    plan = {
        "schema_version": "1.2.0",
        "result_type": "AnalysisPlan",
        "scope": config["runtime"]["scope"],
        "requested_phase": config["runtime"]["phase"],
        "config": plan_config(config),
        "sources": [{"artifact_id": s["artifact_id"], "file_name": Path(s["path"]).name, "sha256": s["sha256"], "hash_rule": CONTENT_HASH_RULE} for s in sources],
        "artifacts": artifacts,
        "scale": manifest["scale"],
        "original_observed_mask": manifest["original_observed_mask"],
        "designs": diagnostics["designs"] if diagnostics else [],
        "contrasts": contrasts,
        "models": [dict(m, weights_sha256=(m["weights"] or {}).get("sha256"), weights=None) for m in request["parameters"]["models"]] if request else [],
        "engine_eligibility": eligibility,
        "blocking": diagnostics["blocking"] if diagnostics else [],
        "families": _families(config, contrasts, eligibility) if config["runtime"]["scope"] == "analysis" else [],
        "omnibus_tests": config.get("omnibus_tests", []),
        "capability_plan": capability_plan,
        "environment": environment,
        "code_sha256": code_sha256,
        "runtime": {"seed": config["runtime"]["seed"], "threads": config["runtime"]["threads"], "rng_kind": "L'Ecuyer-CMRG",
                    "execution_profile": config["runtime"]["execution_profile"]},
        "score_testing": "disabled" if config.get("score_test", "off") == "off" else config["score_test"],
    }
    plan["plan_hash"] = semantic_plan_hash(plan)
    plan["integrity_sha256"] = plan_integrity_hash(plan)
    return plan


# D-42 (2026-10-03, review follow-up to CI run 37087074374): the plan hash identifies the *analysis*, not the machine.
# It covers configuration, LF-normalised source and planning-artifact content hashes, design, contrasts, families,
# capability plan, code identity and runtime settings.  The solved environment and the R session-info artifacts are
# provenance: they stay in plan.json, outside the plan hash, protected by integrity_sha256 (checked by verify_plan).
PROVENANCE_ONLY_KEYS = ("environment", "integrity_sha256")
PROVENANCE_RESULT_TYPES = ("session_info",)


def semantic_view(plan: dict) -> dict:
    view = {k: v for k, v in plan.items() if k not in PROVENANCE_ONLY_KEYS}
    view["artifacts"] = [{"stage": a["stage"], "artifact_id": a["artifact_id"], "relative_path": a["relative_path"], "result_type": a["result_type"],
                          "content_sha256": a["content_sha256"]} for a in plan["artifacts"] if a["result_type"] not in PROVENANCE_RESULT_TYPES]
    return view


def semantic_plan_hash(plan: dict) -> str:
    return canonical_json_sha256(semantic_view(plan))


def plan_integrity_hash(plan: dict) -> str:
    """Hash of the complete plan record (including provenance), excluding only the two hash fields."""
    return canonical_json_sha256({k: v for k, v in plan.items() if k != "integrity_sha256"})


def write_plan(plan: dict, path: str | Path) -> None:
    path = Path(path)
    if path.exists():
        raise IntegrityError(f"plan already exists: {path}", code="E_PATH_COLLISION")
    from .paths import atomic_write_bytes
    atomic_write_bytes(path, json.dumps(plan, indent=2, ensure_ascii=False, sort_keys=True, allow_nan=False).encode("utf-8") + b"\n")


def verify_plan(path: str | Path, *, expected_hash: str | None = None) -> dict:
    """Recompute the plan hash and every referenced artifact hash (V039 negative)."""
    path = Path(path)
    if not path.is_file():
        raise IntegrityError("no frozen plan exists; a fit request requires a valid plan", code="E_PLAN_REQUIRED")
    plan = json.loads(path.read_text(encoding="utf-8"))
    if "integrity_sha256" not in plan or any("content_sha256" not in a for a in plan.get("artifacts", [])):
        raise IntegrityError("plan predates the platform-independent plan hash (D-42); re-plan with this release", code="E_PLAN_CHANGED")
    recomputed = semantic_plan_hash(plan)
    if plan.get("plan_hash") != recomputed or (expected_hash is not None and expected_hash != recomputed):
        raise IntegrityError("plan hash does not match its content or the requested plan", code="E_PLAN_CHANGED")
    if plan.get("integrity_sha256") != plan_integrity_hash(plan):
        raise IntegrityError("plan provenance (environment or session information) changed after freeze", code="E_PLAN_CHANGED")
    root = path.parent
    for artifact in plan["artifacts"]:
        target = root / artifact["relative_path"]
        if not target.is_file() or sha256_file(target) != artifact["sha256"]:
            raise IntegrityError(f"planned artifact changed after freeze: {artifact['relative_path']}", code="E_PLAN_CHANGED")
    return plan


def plan_bytes(plan: dict) -> bytes:
    return canonical_json_bytes(plan)


# --------------------------------------------------------------------------- re-plan API for named secondary designs
# Maintainer amendment A-2026-10-01-14 (ADR 0009): R14b/R14c derive named secondary designs (covariate-adjusted,
# subgroup, phenotype) through the same finite grammar (compile_design/compile_contrasts) and exact rank and
# estimability checks, without editing the frozen primary design.  The R stages re-verify rank with design_rank.
def exact_rank(matrix: list[list[float]]) -> tuple[int, list[int]]:
    """Exact rank by rational Gaussian elimination; returns (rank, pivot column indices)."""
    from fractions import Fraction
    rows = [[Fraction(v) for v in row] for row in matrix]
    n_cols = len(rows[0]) if rows else 0
    pivots, r = [], 0
    for c in range(n_cols):
        pivot = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        lead = rows[r][c]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                factor = rows[i][c] / lead
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[r])]
        pivots.append(c); r += 1
        if r == len(rows):
            break
    return r, pivots


def replan(config: dict, observations: list[dict], *, design_id: str, group_levels: list[str] | None = None, continuous: list[str] = (),
           categorical: list[str] = (), interactions: list[list[str]] = (), observation_ids: list[str] | None = None, contrast_ids: list[str] | None = None,
           extra_contrasts: list[dict] = ()) -> dict:
    """Compile a named secondary design from the frozen primary design plus declared terms on an optional observation subset.

    Returns the compiled design, its compiled contrasts (planned contrasts whose required groups are all present, plus
    ``extra_contrasts``), exact rank/aliasing and exact contrast estimability.  Raises DesignError for grammar errors.
    """
    primary = config["design"]
    rows = [o for o in observations if observation_ids is None or o["observation_id"] in set(observation_ids)]
    levels = list(group_levels) if group_levels is not None else [l for l in primary["group_levels"] if any(o[primary["group_column"]] == l for o in rows)]
    rows = [o for o in rows if o[primary["group_column"]] in levels]
    spec = copy.deepcopy(primary)
    spec.update({"id": design_id, "group_levels": levels, "interactions": [list(i) for i in interactions]})
    spec["continuous_covariates"] = [dict(c) for c in primary["continuous_covariates"]] + [{"column": c, "center": True} for c in continuous]
    spec["categorical_covariates"] = [dict(c) for c in primary["categorical_covariates"]]
    for column in categorical:
        present = sorted({o[column] for o in rows if _nonnull(o.get(column))})
        missing = [o["observation_id"] for o in rows if not _nonnull(o.get(column))]
        if missing:
            raise DesignError("E_DESIGN_COVARIATE_NONFINITE", f"covariate {column!r} is missing for {missing[:5]}; complete-case removal is not silent", f"/design[{design_id}]")
        spec["categorical_covariates"].append({"column": column, "levels": present, "reference": present[0]})
    compiled = compile_design(spec, rows, assay=config["assay"], tmt_strategy=config["preprocessing"].get("tmt_strategy"))
    wanted = [c for c in config["contrasts"] if c["design_id"] == primary["id"] and (contrast_ids is None or c["id"] in contrast_ids)]
    contrasts = []
    for contrast in wanted:
        if not all(g in levels for g in contrast["required_groups"]):
            continue
        contrasts.append(dict(contrast, design_id=design_id))
    contrasts += [dict(c, design_id=design_id) for c in extra_contrasts]
    compiled_contrasts = compile_contrasts({"contrasts": contrasts}, {design_id: compiled})
    rank, pivots = exact_rank(compiled["matrix"])
    names = compiled["coefficients"]
    aliased = [names[j] for j in range(len(names)) if j not in pivots]
    for contrast in compiled_contrasts:
        augmented = compiled["matrix"] + [contrast["weights"]]
        contrast["estimable"] = exact_rank(augmented)[0] == rank
    return {"design": compiled, "contrasts": compiled_contrasts, "rank": rank, "full_rank": rank == len(names), "aliased": aliased,
            "observation_ids": [o["observation_id"] for o in rows], "group_levels": levels,
            "rank_method": "exact rational Gaussian elimination (re-verified by design_rank in R)"}
