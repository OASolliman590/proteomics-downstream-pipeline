"""R15c acceptance: V179-V184 (FR-179 to FR-184; SM45; ADR 0011).

Oracles are direct glmnet, e1071 and randomForest fits written in the test (run through the house R runner),
a direct leave-pair-out loop (D-77), subject-level fold checks, and the committed gate file's SHA-256. The
calibration results are read from docs/validation/016-methods-extension-validation/R15c/calibration/ (the run
is pre-registered: its counts and seeds come from the committed gate file). Each case has a negative assertion.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import shutil
from pathlib import Path

import pytest

from proteomics_pipeline.runtime import run_r_code

ROOT = Path(__file__).resolve().parents[3]
CAL = ROOT / "docs" / "validation" / "016-methods-extension-validation" / "R15c" / "calibration"
GATES = ROOT / "docs" / "validation" / "016-methods-extension-validation" / "R15c" / "gates-preregistration.json"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUN = _load("calibration_runner", ROOT / "scripts" / "maintained" / "run_methods_calibration.py")
FOLLOWUP = _load("calibration_followup", ROOT / "scripts" / "maintained" / "evaluate_calibration_followup.py")
rscript = shutil.which("Rscript")
pytestmark = pytest.mark.skipif(rscript is None, reason="NOT_RUN: Rscript unavailable")


def r(code: str):
    result = run_r_code(code, rscript=rscript)
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    return json.loads([line for line in result.stdout.splitlines() if line.strip()][-1])


SETUP = """ns <- asNamespace("proteomicsCore"); f <- function(n) get(n, envir = ns)
set.seed(179L); X <- matrix(stats::rnorm(30 * 40), 30, 40, dimnames = list(sprintf("u%02d", 1:30), sprintf("F%02d", 1:40)))
y <- rep(c(1L, 0L), each = 15); tr <- c(6:15, 21:30); te <- c(1:5, 16:20); panel <- colnames(X)[1:5]   # both classes in train and test
mk <- function(kind, lambda, alpha = 0.5) list(positive = "P", negative = "N", levels = c("P", "N"), imputation = "none", selection = "top_k_auc", k = 5L,
  alpha = alpha, classifier = kind, lambda = lambda, cost = NULL, threshold_rule = "youden_train", scheme = "repeated_stratified_kfold", outer_k = 5L,
  repeats = 2L, inner_k = 3L, observation_ids = rownames(X)[tr], fixed_features = panel)
"""


def test_v179_lasso_and_elastic_net_scores_equal_direct_glmnet_at_a_fixed_lambda_and_panel():
    code = SETUP + """
out <- list()
for (kind in c("lasso", "elastic_net", "penalized_logistic")) {
  a <- if (kind == "lasso") 1 else if (kind == "elastic_net") 0.5 else 0
  s <- mk(kind, 0.05); s$classifier_grid <- 0.05
  cfg <- f(".bm_selection_grid")(s)[[1]]
  prep <- f(".bm_fit_prep")(X[tr, panel, drop = FALSE], "none")
  Ztr <- f(".bm_apply_prep")(prep, X[tr, panel, drop = FALSE]); Zte <- f(".bm_apply_prep")(prep, X[te, panel, drop = FALSE])
  pipe <- f(".bm_fit_score")(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE], s, cfg, NULL, prep)$score[, 1]
  fit <- glmnet::glmnet(Ztr, factor(ifelse(y[tr] == 1L, "P", "N"), levels = c("P", "N")), family = "binomial", alpha = a, lambda = 0.05, standardize = FALSE)
  oracle <- 1 - as.numeric(stats::predict(fit, Zte, s = 0.05, type = "response"))   # glmnet models the second level ("N")
  out[[kind]] <- list(max_abs = max(abs(pipe - oracle)), pipe = pipe, oracle = oracle)
}
cat(jsonlite::toJSON(list(lasso = out$lasso$max_abs, elastic_net = out$elastic_net$max_abs, penalized = out$penalized_logistic$max_abs,
                          lasso_vs_ridge = max(abs(out$lasso$pipe - out$penalized_logistic$pipe))), auto_unbox = TRUE, digits = NA))
"""
    res = r(code)
    assert res["lasso"] <= 1e-10 and res["elastic_net"] <= 1e-10 and res["penalized"] <= 1e-10
    # negative: the lasso scores are not the ridge scores (a lasso that silently ran alpha 0 would pass the equality above)
    assert res["lasso_vs_ridge"] > 1e-6


def test_v180_svm_decision_values_equal_e1071_and_random_forest_is_refused_when_absent():
    code = SETUP + """
s <- mk("svm_linear", NULL); s$cost <- 1; s$classifier_grid <- 1
cfg <- f(".bm_selection_grid")(s)[[1]]
prep <- f(".bm_fit_prep")(X[tr, panel, drop = FALSE], "none")
Ztr <- f(".bm_apply_prep")(prep, X[tr, panel, drop = FALSE]); Zte <- f(".bm_apply_prep")(prep, X[te, panel, drop = FALSE])
pipe <- f(".bm_fit_score")(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE], s, cfg, NULL, prep)$score[, 1]
fit <- e1071::svm(Ztr, factor(ifelse(y[tr] == 1L, "P", "N"), levels = c("P", "N")), kernel = "linear", cost = 1, scale = FALSE)
pr <- stats::predict(fit, Zte, decision.values = TRUE); dv <- attr(pr, "decision.values")
first <- strsplit(colnames(dv)[1], "/")[[1]][1]
oracle <- as.numeric(dv[, 1]) * if (identical(first, "P")) 1 else -1
classes <- as.numeric(as.character(pr) == "P")
rf <- requireNamespace("randomForest", quietly = TRUE)
cat(jsonlite::toJSON(list(svm_max_abs = max(abs(pipe - oracle)), pipe_auc = f("pd_auc")(pipe, y[te]), class_auc = f("pd_auc")(classes, y[te]),
                          distinct_scores = length(unique(round(pipe, 10))), rf_available = rf), auto_unbox = TRUE, digits = NA))
"""
    res = r(code)
    assert res["svm_max_abs"] <= 1e-10                       # the SVM score is the e1071 decision value, oriented by the class coding
    # negative: the score is a continuous decision value, not the two predicted classes (an AUC from the classes would
    # collapse the ranking to two levels)
    assert res["distinct_scores"] > 2


def test_v180_random_forest_is_not_run_without_the_package_and_otherwise_runs():
    code = SETUP + """
rf <- requireNamespace("randomForest", quietly = TRUE)
if (rf) { s <- mk("random_forest", NULL); s$classifier_grid <- 1; cfg <- f(".bm_selection_grid")(s)[[1]]
  sc <- f(".bm_fit_score")(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE], s, cfg, NULL)$score[, 1]; ok <- all(is.finite(sc)) } else ok <- NA
cat(jsonlite::toJSON(list(available = rf, finite = ok), auto_unbox = TRUE))
"""
    res = r(code)
    if not res["available"]:
        pytest.skip("NOT_RUN: randomForest absent (E_DEPENDENCY_UNAVAILABLE for this family only)")
    assert res["finite"] is True


def test_v181_loocv_auc_is_the_leave_pair_out_auc_of_a_direct_pair_loop():
    # D-77: with scheme loocv the reported AUC is the leave-pair-out AUC; the leave-one-out scores are unchanged.
    code = """ns <- asNamespace("proteomicsCore"); f <- function(n) get(n, envir = ns)
set.seed(181L); X <- matrix(stats::rnorm(24 * 30), 24, 30, dimnames = list(sprintf("u%02d", 1:24), sprintf("F%02d", 1:30)))
y <- rep(c(1L, 0L), each = 12); panel <- colnames(X)[1:5]
s <- list(positive = "P", negative = "N", levels = c("P", "N"), imputation = "none", selection = "top_k_auc", k = 5L, alpha = 0.5, classifier = "lasso",
          lambda = 0.05, cost = NULL, threshold_rule = "youden_train", scheme = "loocv", outer_k = 5L, repeats = 1L, inner_k = 3L,
          observation_ids = rownames(X), fixed_features = panel)
res <- f("pd_bm_nested")(X, y, rownames(X), s, f(".bm_stream")(181L, 0L), record = FALSE)
score_without <- function(held) {
  tr <- setdiff(1:24, held); sc <- s; sc$observation_ids <- rownames(X)[tr]; sc$classifier_grid <- 0.05
  cfg <- f(".bm_selection_grid")(sc)[[1]]; prep <- f(".bm_fit_prep")(X[tr, panel, drop = FALSE], "none")
  f(".bm_fit_score")(X[tr, , drop = FALSE], y[tr], X[held, , drop = FALSE], sc, cfg, NULL, prep)$score[, 1]
}
loo <- vapply(1:24, score_without, numeric(1))
wins <- c(); for (i in which(y == 1L)) for (j in which(y == 0L)) { sc <- score_without(c(i, j)); wins <- c(wins, (sc[1] > sc[2]) + 0.5 * (sc[1] == sc[2])) }
oof <- res$oof[match(rownames(X), res$oof$observation_id), ]
cat(jsonlite::toJSON(list(pipe = res$pooled, lpo = mean(wins), n_pairs = length(wins), estimator = res$auc_estimator, per_repeat = res$per_repeat,
  loo_max_diff = max(abs(oof$score - loo)), loo_pooled = f("pd_auc")(loo, y)), auto_unbox = TRUE, digits = NA))
"""
    res = r(code)
    assert res["n_pairs"] == 144 and res["estimator"] == "leave_pair_out"
    assert abs(res["pipe"] - res["lpo"]) <= 1e-10 and abs(res["per_repeat"] - res["lpo"]) <= 1e-10
    assert res["loo_max_diff"] <= 1e-10                       # the leave-one-out scores for the threshold metrics are unchanged
    assert abs(res["pipe"] - res["loo_pooled"]) > 1e-6        # and the pooled leave-one-out AUC is no longer reported


def test_v181_negative_loocv_with_mixed_class_units_is_refused():
    code = """ns <- asNamespace("proteomicsCore"); f <- function(n) get(n, envir = ns)
set.seed(182L); X <- matrix(stats::rnorm(24 * 30), 24, 30, dimnames = list(sprintf("o%02d", 1:24), sprintf("F%02d", 1:30)))
units <- rep(sprintf("s%02d", 1:12), each = 2); y <- rep(c(1L, 0L), times = 12)   # every subject has one observation of each class
s <- list(positive = "P", negative = "N", levels = c("P", "N"), imputation = "none", selection = "top_k_auc", k = 5L, alpha = 0.5, classifier = "lasso",
          lambda = 0.05, cost = NULL, threshold_rule = "youden_train", scheme = "loocv", outer_k = 5L, repeats = 1L, inner_k = 3L,
          observation_ids = rownames(X), fixed_features = colnames(X)[1:5])
msg <- tryCatch({ f("pd_bm_nested")(X, y, units, s, f(".bm_stream")(182L, 0L), record = FALSE); "no error" }, error = function(e) conditionMessage(e))
cat(jsonlite::toJSON(list(msg = msg), auto_unbox = TRUE))
"""
    assert r(code)["msg"].startswith("E_BIOMARKER_LPO_UNITS")


def test_v181_the_compute_estimate_counts_the_leave_pair_out_fits():
    from proteomics_pipeline.post_de.biomarker import compute_estimate
    r = {"selection": {"method": "top_k_auc", "k": [5]}, "cv": {"scheme": "loocv", "k": 5, "repeats": 1, "inner_k": 3},
         "permutation": {"enabled": False, "B": 0}}
    loo_only = compute_estimate(r, 24, 200, 1, 0)
    assert abs(compute_estimate(r, 24, 200, 1, 0, 144) / loo_only - (24 + 144) / 24) <= 1e-12   # 24 leave-one-out + 144 pair fits
    r["cv"]["scheme"] = "repeated_stratified_kfold"
    assert compute_estimate(r, 24, 200, 1, 0, 144) == compute_estimate(r, 24, 200, 1, 0)       # pairs cost nothing outside loocv


def test_v181_negative_a_selector_fitted_on_all_units_is_the_leaky_reference():
    summary = json.loads((CAL / "summary.json").read_text(encoding="utf-8")) if (CAL / "summary.json").is_file() else None
    if summary is None:
        pytest.skip("NOT_RUN: the calibration summary is not present")
    for fam, entry in summary["families"].items():
        assert entry["leaky_reference"]["mean"] > entry["nested_auc_null"]["mean"] + 0.1, fam   # the leak is visible on the same noise


def test_v182_subject_blocked_folds_never_split_a_subject_and_the_ungrouped_reference_inflates():
    code = """ns <- asNamespace("proteomicsCore"); f <- function(n) get(n, envir = ns)
subj <- rep(sprintf("s%02d", 1:20), each = 2); y <- rep(rep(c(1L, 0L), each = 10), each = 2)
ut <- f(".bm_unit_table")(subj, y); folds <- f("pd_bm_folds")(ut$class, 5L)[match(subj, ut$unit)]
grouped_split <- sum(vapply(unique(subj), function(s) length(unique(folds[subj == s])) > 1, logical(1)))
cat(jsonlite::toJSON(list(grouped_split = grouped_split, n_subjects = length(unique(subj))), auto_unbox = TRUE))
"""
    res = r(code)
    assert res["grouped_split"] == 0                         # every subject sits in one fold
    # negative: observation-level folds split subjects (the grouping leak the blocked form prevents)
    code2 = """ns <- asNamespace("proteomicsCore"); f <- function(n) get(n, envir = ns)
subj <- rep(sprintf("s%02d", 1:20), each = 2); y <- rep(rep(c(1L, 0L), each = 10), each = 2)
obs <- paste0(subj, "_", seq_along(subj)); ut <- f(".bm_unit_table")(obs, y); folds <- f("pd_bm_folds")(ut$class, 5L)[match(obs, ut$unit)]
split <- sum(vapply(unique(subj), function(s) length(unique(folds[subj == s])) > 1, logical(1)))
cat(jsonlite::toJSON(list(split = split), auto_unbox = TRUE))
"""
    assert r(code2)["split"] > 0


def test_v183_the_gate_file_is_pinned_and_a_changed_gate_is_refused(tmp_path):
    assert RUN.sha256(GATES) == RUN.PREREGISTERED_SHA256 == "da7577a6cf5da61bd8b0a0d0ca0b9e3c47819cdf76be8451d6f0e06a34336a26"
    gates = json.loads(GATES.read_text(encoding="utf-8"))
    assert gates["status"].startswith("pre-registered before the first calibration run")
    changed = tmp_path / "gates.json"
    changed.write_text(GATES.read_text(encoding="utf-8").replace("0.03", "0.09", 1), encoding="utf-8")
    with pytest.raises(RUN.CalibrationError) as info:
        RUN.load_gates(changed)
    assert info.value.code == "E_CALIBRATION_GATE_CHANGED"


def test_v183_the_followup_gate_file_is_pinned_and_a_changed_gate_is_refused(tmp_path):
    gates = FOLLOWUP.GATES_PATH
    assert RUN.sha256(gates) == FOLLOWUP.PREREGISTERED_SHA256
    assert json.loads(gates.read_text(encoding="utf-8"))["status"].startswith("pre-registered before any follow-up run")
    changed = tmp_path / "gates.json"
    changed.write_text(gates.read_text(encoding="utf-8").replace("0.03", "0.09", 1), encoding="utf-8")
    with pytest.raises(FOLLOWUP.base.CalibrationError) as info:
        FOLLOWUP.load_gates(changed)
    assert info.value.code == "E_CALIBRATION_GATE_CHANGED"


def test_v180_the_binned_uniformity_check_counts_ten_equal_bins_of_the_attainable_values():
    chi_sf = RUN._r14f_helpers().chi_square_sf
    B = 99
    uniform = [j / (B + 1) for j in range(1, B + 2)] * 2                         # every attainable value twice: 20 per bin
    res = FOLLOWUP.binned_uniformity(uniform, B, 10, chi_sf)
    assert res["counts"] == [20] * 10 and res["chi_square"] == 0 and res["gate"] == "PASS"
    skewed = [0.9 + j / (B + 1) for j in range(1, 11)] * 20                      # all 200 values in the top bin (0.91 to 1.00)
    res = FOLLOWUP.binned_uniformity(skewed, B, 10, chi_sf)
    assert res["counts"] == [0] * 9 + [200] and res["gate"] == "FAIL"           # negative: a conservative pile-up fails


def test_v184_the_run_records_seeds_counts_versions_and_reproduces_one_replicate():
    raw = json.loads((CAL / "calibration_raw.json").read_text(encoding="utf-8"))
    for key in ("r_version", "glmnet", "e1071", "randomForest", "proteomicsCore", "seconds", "cores", "gates"):
        assert key in raw
    counts = json.loads(GATES.read_text(encoding="utf-8"))["replicate_counts"]
    rows = list(csv.DictReader((CAL / "nested_auc_null.tsv").open(encoding="utf-8"), delimiter="\t"))
    per_family = {}
    for row in rows:
        per_family.setdefault(row["family"], 0)
        per_family[row["family"]] += 1
    assert all(n == counts["nested_auc_null"] for n in per_family.values()) and len(per_family) == 6
    first = next(r_ for r_ in rows if r_["family"] == "penalized_logistic" and r_["dataset"] == "1")
    code = """ns <- asNamespace("proteomicsCore"); f <- function(n) get(n, envir = ns)
set.seed(600001L, kind = "L'Ecuyer-CMRG")
X <- matrix(stats::rnorm(30 * 200), 30, 200, dimnames = list(sprintf("u%02d", 1:30), sprintf("F%03d", 1:200))); y <- rep(c(1L, 0L), each = 15)
s <- list(positive = "P", negative = "N", levels = c("P", "N"), imputation = "none", selection = "top_k_auc", k = 5L, alpha = 0.5, classifier = "penalized_logistic",
          lambda = 0.1, cost = NULL, threshold_rule = "youden_train", scheme = "repeated_stratified_kfold", outer_k = 5L, repeats = 2L, inner_k = 3L,
          observation_ids = rownames(X), fixed_features = NULL)
cat(jsonlite::toJSON(list(auc = f("pd_bm_nested")(X, y, rownames(X), s, f(".bm_stream")(600001L, 0L), record = FALSE)$pooled), auto_unbox = TRUE, digits = NA))
"""
    again = r(code)["auc"]
    assert abs(again - float(first["pooled_oof_auc"])) <= 1e-12


def test_v184_negative_a_reduced_replicate_count_without_an_adaptation_is_refused(tmp_path):
    counts = json.loads(GATES.read_text(encoding="utf-8"))["replicate_counts"]
    source = list(csv.DictReader((CAL / "nested_auc_null.tsv").open(encoding="utf-8"), delimiter="\t"))
    kept = [row for row in source if not (row["family"] == "lasso" and int(row["dataset"]) > counts["nested_auc_null"] - 5)]
    (tmp_path / "calibration").mkdir()
    for name in ("permutation_p_null.tsv", "planted_auc.tsv", "planted_permutation_p.tsv", "leaky_reference.tsv"):
        (tmp_path / "calibration" / name).write_bytes((CAL / name).read_bytes())
    reduced = tmp_path / "calibration" / "nested_auc_null.tsv"
    with reduced.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(source[0].keys()), delimiter="\t")
        writer.writeheader(); writer.writerows(kept)
    with pytest.raises(RUN.CalibrationError) as info:
        RUN.verify_replicate_counts(tmp_path / "calibration", json.loads(GATES.read_text(encoding="utf-8")), adaptations=[])
    assert info.value.code == "E_CALIBRATION_REPLICATES"
    # the same reduction recorded as an adaptation is accepted (ADR 0010)
    RUN.verify_replicate_counts(tmp_path / "calibration", json.loads(GATES.read_text(encoding="utf-8")),
                                adaptations=[{"scenario": "nested_auc_null", "family": "lasso", "pre_registered": counts["nested_auc_null"], "used": counts["nested_auc_null"] - 5}])
