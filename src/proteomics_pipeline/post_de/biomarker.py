"""Post-DE biomarker discrimination evaluation (packet R14d, SM34-SM38, FR-149-FR-160).

Binary declared contrasts only.  The planner refuses evaluation before any computation when a class has fewer than
the declared minimum of biological units (E_BIOMARKER_SMALL_N), when cross-validation folds cannot contain both
classes (E_BIOMARKER_FOLDS), when the declared work exceeds the compute guard (E_BIOMARKER_COMPUTE_GUARD), and when a
declaration would leak held-out information (E_BIOMARKER_LEAKAGE, E_BIOMARKER_GROUP_LEAKAGE,
E_BIOMARKER_THRESHOLD_LEAKAGE, E_BIOMARKER_PERMUTATION_SCOPE, E_VALIDATION_RETUNED, E_SCORE_SELECTION_OVERLAP).
The R stage enforces leakage safety structurally (fit/apply transforms audited per fold).
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path

from . import PostDeRefusal, base_inputs, block as post_de_block, common_parameters, dea_inputs, execute_r, request, stage_output
from ..provenance import sha256_file

CAPABILITY = "post_de_biomarker"
PREREQUISITES = ("dea",)
CLASSIFIERS = ("penalized_logistic", "svm_linear", "svm_polynomial", "svm_radial", "random_forest")
DEFAULTS = {"min_units_per_class": 5, "imputation": "none", "threshold_rule": "youden_train", "classifier": "penalized_logistic", "calibration_bins": 10}
COMPUTE_GUARD_HOURS = 6.0
# Reference cost model (seconds per classifier fit on the reference machine; see D-48): base cost plus a term for
# ranking every candidate feature on the training data.
COST_BASE_SECONDS, COST_PER_CELL_SECONDS = 0.02, 2e-6


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["limma"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore", "limma", "glmnet"]}]


def precheck(block: dict, raw: dict) -> None:
    rule = str(block.get("threshold_rule", DEFAULTS["threshold_rule"]))
    if rule != "youden_train" and not re.fullmatch(r"fixed:-?[0-9]+(\.[0-9]+)?", rule):
        raise PostDeRefusal("E_BIOMARKER_THRESHOLD_LEAKAGE", f"threshold rule {rule!r} is not chosen inside training folds; allowed: youden_train or fixed:<value>", "/post_de/biomarker/threshold_rule")
    scope = (block.get("permutation") or {}).get("scope", "whole_procedure")
    if scope != "whole_procedure":
        raise PostDeRefusal("E_BIOMARKER_PERMUTATION_SCOPE", f"permutation scope {scope!r}: the permutation test must rerun the whole nested procedure, including selection and tuning", "/post_de/biomarker/permutation/scope")
    imputation = block.get("imputation", "none")
    if imputation not in ("none", "train_median"):
        raise PostDeRefusal("E_BIOMARKER_LEAKAGE", f"imputation {imputation!r} would be fitted outside the training folds; allowed: none or train_median (fitted per training fold)", "/post_de/biomarker/imputation")
    cohort = block.get("validation_cohort") or {}
    if cohort.get("retune"):
        raise PostDeRefusal("E_VALIDATION_RETUNED", "the validation cohort is evaluated once with the locked discovery model; re-tuning on it is forbidden", "/post_de/biomarker/validation_cohort/retune")
    design = raw.get("design") or {}
    blocked = (design.get("blocking") or {}).get("mode", "none") != "none"
    if blocked and (block.get("cv") or {}).get("group_by_subject") is False:
        raise PostDeRefusal("E_BIOMARKER_GROUP_LEAKAGE", "a subject-blocked design must group folds by subject; group_by_subject=false would leak subjects across folds", "/post_de/biomarker/cv/group_by_subject")
    if not block.get("contrast"):
        raise PostDeRefusal("E_BIOMARKER_CONTRAST", "post_de.biomarker.contrast must name a declared binary contrast", "/post_de/biomarker/contrast")


def _read_matrix(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows[0][1:], {r[0]: r[1:] for r in rows[1:]}


def resolved(block: dict, config: dict) -> dict:
    cv = block.get("cv") or {}
    outer = cv.get("outer") or {}
    selection = block.get("selection") or {}
    single = block.get("single_feature") or {}
    perm = block.get("permutation") or {}
    return {"contrast": block["contrast"], "min_units_per_class": int(block.get("min_units_per_class", DEFAULTS["min_units_per_class"])),
            "imputation": block.get("imputation", "none"), "threshold_rule": block.get("threshold_rule", "youden_train"),
            "classifier": block.get("classifier", "penalized_logistic"), "classifier_params": block.get("classifier_params") or {},
            "calibration_bins": int(block.get("calibration_bins", 10)), "seed": int(block.get("seed", config["runtime"]["seed"])),
            "cv": {"scheme": outer.get("scheme", "repeated_stratified_kfold"), "k": int(outer.get("k", 5)), "repeats": int(outer.get("repeats", 100)),
                   "inner_k": int((cv.get("inner") or {}).get("k", 5)), "group_by_subject": cv.get("group_by_subject", "auto")},
            "selection": {"method": selection.get("method", "top_k_auc"), "k": list(selection.get("k", [5, 10])), "alpha": selection.get("alpha", 0.5)},
            "single_feature": {"enabled": single.get("enabled", True), "features": single.get("features", "all_tested"), "direction": single.get("direction", "prespecified"),
                               "expected_direction": single.get("expected_direction", "higher_in_positive"), "ci": list(single.get("ci", ["delong", "bootstrap"])),
                               "bootstrap": int(single.get("bootstrap", 2000))},
            "permutation": {"enabled": perm.get("enabled", True), "B": int(perm.get("B", 1000)), "seed": int(perm.get("seed", block.get("seed", config["runtime"]["seed"]))),
                            "compute_guard_hours": float(perm.get("compute_guard_hours", COMPUTE_GUARD_HOURS)), "scope": "whole_procedure"},
            "fixed_panels": list(block.get("fixed_panels", [])), "validation_cohort": block.get("validation_cohort")}


def compute_estimate(r: dict, n_obs: int, n_candidates: int, threads: int, n_panels: int) -> float:
    """Hours for the declared nested CV x permutation work (reference cost model, D-48)."""
    method = r["selection"]["method"]
    grid = len(r["selection"]["k"]) if method in ("top_k_auc", "top_k_t") else 10
    repeats = 1 if r["cv"]["scheme"] == "loocv" else r["cv"]["repeats"]
    outer = n_obs if r["cv"]["scheme"] == "loocv" else r["cv"]["k"]
    jobs = 1 + (r["permutation"]["B"] if r["permutation"]["enabled"] else 0) + 2 * n_panels
    fits = jobs * repeats * outer * (r["cv"]["inner_k"] * grid + 1)
    seconds = fits * (COST_BASE_SECONDS + COST_PER_CELL_SECONDS * n_obs * n_candidates)
    return seconds / 3600 / max(1, threads)


def plan_checks(block: dict, config: dict, context: dict) -> dict:
    r = resolved(block, config)
    design = config["design"]
    group_column = design["group_column"]
    contrast = next((c for c in config["contrasts"] if c["id"] == r["contrast"]), None)
    if contrast is None or len(contrast["required_groups"]) != 2 or contrast["design_id"] != design["id"]:
        raise PostDeRefusal("E_BIOMARKER_CONTRAST", f"biomarker evaluation needs a declared binary contrast of the primary design; {r['contrast']!r} is not one", "/post_de/biomarker/contrast")
    if r["classifier"] not in CLASSIFIERS:
        raise PostDeRefusal("E_BIOMARKER_CLASSIFIER", f"unknown classifier {r['classifier']!r}", "/post_de/biomarker/classifier")
    positive, negative = contrast["required_groups"]
    subject = design["blocking"].get("subject_column") if design["blocking"]["mode"] != "none" else None
    observations = [o for o in context["observations"] if o[group_column] in (positive, negative)]
    unit = (lambda o: o[subject]) if subject else (lambda o: o.get("biological_unit_id") or o["observation_id"])
    classes: dict[str, set] = {}
    for o in observations:
        classes.setdefault(unit(o), set()).add(o[group_column])
    mixed = any(len(v) > 1 for v in classes.values())
    n_pos = sum(1 for v in classes.values() if positive in v); n_neg = sum(1 for v in classes.values() if negative in v)
    reasons = []
    minimum = r["min_units_per_class"]
    if min(n_pos, n_neg) < minimum:
        reasons.append(("E_BIOMARKER_SMALL_N", f"{positive}: {n_pos}, {negative}: {n_neg} biological units; the declared minimum per class is {minimum}"))
    minority = min(n_pos, n_neg)
    k = len(classes) if r["cv"]["scheme"] == "loocv" else r["cv"]["k"]
    if r["cv"]["scheme"] != "loocv" and minority < k:
        reasons.append(("E_BIOMARKER_FOLDS", f"{k}-fold CV with {minority} minority units: some folds cannot contain both classes"))
    train_minority = minority - (1 if r["cv"]["scheme"] == "loocv" else math.ceil(minority / max(1, k)))
    if train_minority < r["cv"]["inner_k"]:
        reasons.append(("E_BIOMARKER_FOLDS", f"inner {r['cv']['inner_k']}-fold CV with {train_minority} minority units per training fold: inner folds cannot contain both classes"))
    # candidate universe: features genuinely observed in every analysed observation (label-free planning restriction)
    pre = Path(context["preprocessing_dir"]) / "primary"
    columns, matrix = _read_matrix(pre / "matrix.tsv")
    _, mask = _read_matrix(pre / "observed_mask.tsv")
    index = [columns.index(o["observation_id"]) for o in observations]
    def usable(f):
        observed = [mask[f][i] == "true" and matrix[f][i] not in ("NA", "") for i in index]
        return all(observed) if r["imputation"] == "none" else sum(observed) >= 0.5 * len(observed)
    candidates = [f for f in matrix if usable(f)]
    if len(candidates) < 2:
        reasons.append(("E_BIOMARKER_FEATURES", f"{len(candidates)} candidate feature(s) observed in the analysed units; at least 2 are needed"))
    sel = r["selection"]
    if sel["method"] in ("top_k_auc", "top_k_t"):
        sel["k"] = [k_ for k_ in sel["k"] if k_ <= len(candidates)]
        if not sel["k"]:
            reasons.append(("E_BIOMARKER_SELECTION", "every declared k exceeds the candidate universe"))
    sub = []
    sf = dict(r["single_feature"]); sf["state"] = "ELIGIBLE"
    if sf["enabled"] and mixed:
        sf.update(state="INAPPLICABLE", reason_code="E_BIOMARKER_PAIRED_AUC", reason="class varies within subject: unit-level single-feature AUC with independent-unit CIs is undefined")
        sub.append({"analysis": "single_feature", "item": "auc", "state": "INAPPLICABLE", "reason_code": sf["reason_code"], "reason": sf["reason"]})
    if sf["features"] != "all_tested":
        sets = post_de_block(config, "sets")
        if sf["features"] not in {d["id"] for d in sets.get("definitions", [])} or not sets.get("enabled"):
            raise PostDeRefusal("E_BIOMARKER_FEATURES", f"single_feature.features {sf['features']!r} is neither all_tested nor a declared post-DE set", "/post_de/biomarker/single_feature/features")
    r["single_feature"] = sf
    panels = []
    for panel in r["fixed_panels"]:
        missing = [f for f in panel["feature_ids"] if f not in candidates]
        if missing:
            sub.append({"analysis": "fixed_panel", "item": panel["id"], "state": "INAPPLICABLE", "reason_code": "E_BIOMARKER_PANEL",
                        "reason": f"panel members not observed in every analysed unit: {missing[:5]}"})
        else:
            panels.append(panel)
    r["fixed_panels"] = panels
    threads = int(config["runtime"]["threads"])
    hours = compute_estimate(r, len(observations), len(candidates), threads, len(panels))
    if hours > r["permutation"]["compute_guard_hours"]:
        reasons.append(("E_BIOMARKER_COMPUTE_GUARD", f"estimated {hours:.2f} h for nested CV x {r['permutation']['B']} permutations exceeds the declared guard of "
                                                     f"{r['permutation']['compute_guard_hours']} h; reduce repeats or B explicitly (work is never reduced silently)"))
    validation = None
    cohort = r["validation_cohort"]
    if cohort:
        base_dir = Path(context["config_dir"])
        paths = {key: (Path(cohort[key]) if Path(cohort[key]).is_absolute() else base_dir / cohort[key]) for key in ("matrix", "metadata")}
        for key, path in paths.items():
            if not path.is_file():
                raise PostDeRefusal("E_VALIDATION_INPUT", f"validation {key} file is missing: {path.name}", f"/post_de/biomarker/validation_cohort/{key}")
        with paths["metadata"].open(encoding="utf-8", newline="") as handle:
            meta = list(csv.DictReader(handle, delimiter="\t"))
        id_key = next(iter(meta[0])) if meta else "observation_id"
        discovery_ids = {o["observation_id"] for o in context["observations"]}
        discovery_units = {o.get("biological_unit_id") for o in context["observations"]} - {None, "", "NA"}
        discovery_subjects = {o.get("subject_id") for o in context["observations"]} | ({o.get(subject) for o in context["observations"]} if subject else set())
        discovery_subjects -= {None, "", "NA"}
        shared = {m[id_key] for m in meta} & discovery_ids
        shared |= {m.get(cohort.get("unit_column", "biological_unit_id")) for m in meta} & discovery_units
        shared |= {m.get(cohort.get("subject_column", "subject_id")) for m in meta} & discovery_subjects
        shared -= {None, "", "NA"}
        validation = {"class_column": cohort["class_column"], "matrix_path": str(paths["matrix"].resolve()), "metadata_path": str(paths["metadata"].resolve()),
                      "matrix_sha256": sha256_file(paths["matrix"]), "metadata_sha256": sha256_file(paths["metadata"]), "n_samples": len(meta)}
        if shared:
            reason = f"the validation cohort shares {len(shared)} subject(s)/unit(s) with the discovery data ({sorted(shared)[:5]}); SM23 requires disjoint cohorts"
            sub.append({"analysis": "validation", "item": "cohort", "state": "INAPPLICABLE", "reason_code": "E_SCORE_SELECTION_OVERLAP", "reason": reason})
            validation = None
        elif cohort["class_column"] not in (meta[0] if meta else {}):
            raise PostDeRefusal("E_VALIDATION_INPUT", f"validation metadata lacks class column {cohort['class_column']!r}", "/post_de/biomarker/validation_cohort/class_column")
    parameters = {**{k: v for k, v in r.items() if k != "validation_cohort"}, "positive": positive, "negative": negative,
                  "class_levels": [g for g in design["group_levels"] if g in (positive, negative)], "candidate_features": candidates,
                  "validation": validation, "threads": threads, "compute_estimate_hours": round(hours, 6),
                  "units": {"n_positive": n_pos, "n_negative": n_neg, "unit_level": f"subject:{subject}" if subject else "biological_unit", "class_varies_within_unit": mixed}}
    if reasons:
        return {"state": "INAPPLICABLE", "reason_code": reasons[0][0], "reason": "; ".join(f"{c}: {m}" for c, m in reasons),
                "reasons": [{"reason_code": c, "reason": m} for c, m in reasons], "subanalyses": sub, "resolved": parameters}
    return {"state": "ELIGIBLE", "reason_code": None, "reason": None, "subanalyses": sub, "resolved": parameters, "compute_estimate_hours": round(hours, 6),
            "rule": "binary declared contrast; >= min units per class; outer and inner folds must contain both classes; compute guard; leakage-safe declarations only"}


def required_refusal(decision: dict):
    sub = [s for s in decision.get("subanalyses") or [] if s["analysis"] == "validation"]
    return (sub[0]["reason_code"], sub[0]["reason"]) if sub else None


def build_request(plan: dict, *, plan_path: str | Path, config: dict, run_id: str, output_temp_dir: str | Path, root: Path, decision: dict,
                  config_dir: Path | None = None, stage_id: str = CAPABILITY) -> dict:
    from .sensitivity import _read_observations
    fresh = plan_checks(config["post_de"]["biomarker"], config, {"observations": _read_observations(root), "preprocessing_dir": root / "preprocessing",
                                                                   "config_dir": config_dir or root})
    parameters = fresh["resolved"]
    inputs = base_inputs(plan, plan_path) + dea_inputs(config, root / "dea", names=("families.tsv",))
    sf = parameters["single_feature"]
    if sf["features"] == "all_tested":
        sf["features"] = list(parameters["candidate_features"])
    else:
        found = stage_output(root / "post_de" / "sets", "membership.tsv")
        if found is None:
            from ..errors import ProteomicsError
            raise ProteomicsError("E_PREREQUISITE_FAILED", f"single_feature.features names post-DE set {sf['features']!r}, but the sets stage did not complete", exit_code=4)
        inputs.append({"artifact_id": "stage__sets_membership", "path": str(found[0].resolve()), "sha256": found[1]})
        with found[0].open(encoding="utf-8", newline="") as handle:
            members = [row["feature_id"] for row in csv.DictReader(handle, delimiter="\t") if row.get(sf["features"]) == "true"]
        sf["features"] = [f for f in members if f in set(parameters["candidate_features"])]
    validation = parameters.get("validation")
    if validation:
        inputs.append({"artifact_id": "input__validation_matrix", "path": validation.pop("matrix_path"), "sha256": validation["matrix_sha256"]})
        inputs.append({"artifact_id": "input__validation_metadata", "path": validation.pop("metadata_path"), "sha256": validation["metadata_sha256"]})
    parameters = {**common_parameters(config, plan, {k: v for k, v in decision.items() if k != "resolved"}), **parameters, "contrast": parameters["contrast"]}
    return request(CAPABILITY, stage_id, plan, run_id=run_id, output_temp_dir=output_temp_dir, inputs=inputs, parameters=parameters, seed=parameters["seed"])


def execute(request_value: dict) -> dict:
    return execute_r(request_value)
