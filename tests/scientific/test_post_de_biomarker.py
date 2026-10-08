"""R14d acceptance through the real Python -> R pipeline: V149-V160 (SM34-SM38).

Oracles are written here: pROC::roc/ci.auc (DeLong and stratified bootstrap), a deliberately leaky reference,
an independent nested-CV loop (same documented fold dealing, seeds and selector) with direct glmnet fits, direct
e1071 decision values, hand-applied locked models and Brier/bin arithmetic on the out-of-fold predictions.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import statistics
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


F = _load("post_de_factory", ROOT / "tests" / "fixtures" / "post_de" / "design_factory.py")
B = F.B


def _ready():
    if B.rscript() is None:
        return False
    try:
        return bool(B.r_json("cat(jsonlite::toJSON(all(sapply(c('glmnet','pROC','e1071'), requireNamespace, quietly=TRUE)), auto_unbox=TRUE))"))
    except RuntimeError:
        return False


pytestmark = pytest.mark.skipif(not _ready(), reason="NOT_RUN: Rscript/proteomicsCore/glmnet/pROC/e1071 unavailable")
CLAIMS = {"descriptive", "exploratory_raw_p", "in_sample", "cross_validated_nested", "fixed_panel_cv", "independently_validated", "module_level"}
SEED = 11


def biomarker_block(**overrides):
    block = {"enabled": True, "contrast": "P-N", "cv": {"outer": {"k": 5, "repeats": 3}, "inner": {"k": 3}}, "selection": {"method": "top_k_auc", "k": [5]},
             "classifier_params": {"lambda": [0.1]}, "permutation": {"B": 19}, "single_feature": {"bootstrap": 200}, "seed": SEED}
    for key, value in overrides.items():
        block[key] = value
    return {"enabled": True, "biomarker": block}


def run(tmp_path, values, obs, post, *, groups=("P", "N"), blocked=False, threads=2, expect=0):
    files = F.write_dataset(tmp_path / "data", values, obs)
    def mutate(c):
        c["runtime"]["threads"] = threads
    overrides = {"blocking": {"mode": "duplicate_correlation", "subject_column": "subject_id"}} if blocked else None
    path = F.config(tmp_path / "data", files, groups=list(groups), contrasts=[B.contrast("P-N", "P", "N")], post_de=post, design_overrides=overrides, mutate=mutate)
    payload, code = workflow.run_command(path, tmp_path / "run")
    if expect is not None:
        assert code == expect, payload
    return payload, code, tmp_path / "run"


def _leak_replaced(tmp_path, values, obs, block, code, item, used, blocked=False):
    """D-59: a leaky declaration is refused (recorded with its code) and its leakage-safe form runs; required rejects."""
    payload, _, out = run(tmp_path / "opt", values, obs, block, blocked=blocked)
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_biomarker")
    assert stage["state"] == "COMPLETED", stage
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    entry = next(c for c in plan["capability_plan"] if c["capability"] == "post_de_biomarker")["post_de"]
    assert any(x["reason_code"] == code and x["item"] == item for x in entry["subanalyses"])
    assert any(a["item"] == item and a["used"] == used for a in entry["adaptations"])
    block = json.loads(json.dumps(block)); block["biomarker"]["execution_requirement"] = "required"
    with pytest.raises(ProteomicsError) as error:
        run(tmp_path / "req", values, obs, block, blocked=blocked, expect=None)
    assert error.value.code == code
    return out


def out_dir(run_dir):
    return run_dir / "post_de" / "biomarker"


@pytest.fixture(scope="module")
def signal(tmp_path_factory):
    values, obs = F.biomarker_design()
    return run(tmp_path_factory.mktemp("signal"), values, obs, biomarker_block())[2]


NESTED_ORACLE = """
  a <- commandArgs(TRUE); Y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
  o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(Y), o$observation_id), ]
  o <- o[o$group %in% c('P', 'N'), ]; R <- as.integer(a[3]); K <- 5L; IK <- 3L; k <- 5L; lambda <- 0.1
  units <- if (a[4] == 'subject') o$subject_id else o$observation_id
  X <- t(Y[rowSums(is.na(Y[, o$observation_id])) == 0, o$observation_id]); y <- as.integer(o$group == 'P')
  deal <- function(cls, k) { fold <- integer(length(cls)); off <- 0L
    for (c in sort(unique(cls))) { u <- which(cls == c); u <- u[sample.int(length(u))]; fold[u] <- ((seq_along(u) - 1L + off) %% k) + 1L; off <- off + length(u) }; fold }
  uclass <- function(units, y) { u <- unique(units); list(u = u, c = vapply(u, function(x) as.character(unique(y[units == x])[1]), '')) }
  mwauc <- function(s, y) { r <- rank(c(s[y == 1], s[y == 0])); n1 <- sum(y == 1); (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * sum(y == 0)) }
  set.seed(as.integer(a[5]), kind = "L'Ecuyer-CMRG")
  preds <- list(); selected <- list(); per <- c()
  for (r in seq_len(R)) {
    ut <- uclass(units, y); fold <- deal(ut$c, K)[match(units, ut$u)]; sc <- rep(NA, length(y))
    for (f in sort(unique(fold))) {
      tr <- which(fold != f); te <- which(fold == f)
      it <- uclass(units[tr], y[tr]); invisible(deal(it$c, IK))          # inner folds consume the same random numbers
      Xtr <- X[tr, , drop = FALSE]; mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); keep <- sdv > sqrt(.Machine$double.eps)
      Z <- sweep(sweep(Xtr[, keep], 2, mu[keep]), 2, sdv[keep], '/'); Zte <- sweep(sweep(X[te, keep, drop = FALSE], 2, mu[keep]), 2, sdv[keep], '/')
      au <- apply(Z, 2, mwauc, y = y[tr]); top <- head(names(sort(-abs(au - 0.5), method = 'radix')), k)
      fit <- glmnet::glmnet(Z[, top], factor(ifelse(y[tr] == 1, 'P', 'N'), levels = c('P', 'N')), family = 'binomial', alpha = 0, lambda = lambda, standardize = FALSE)
      p <- as.numeric(predict(fit, Zte[, top, drop = FALSE], s = lambda, type = 'response')); if (fit$classnames[2] != 'P') p <- 1 - p
      sc[te] <- p; selected[[length(selected) + 1]] <- top
      preds[[length(preds) + 1]] <- data.frame(r = r, id = o$observation_id[te], score = p)
    }
    per <- c(per, mwauc(sc, y))
  }
  P <- do.call(rbind, preds)
  cat(jsonlite::toJSON(list(id = P$id, r = P$r, score = P$score, per = per, freq = as.list(table(unlist(selected)) / length(selected))), digits = NA))"""


# ----------------------------------------------------------------------------- V149
def test_v149_every_output_has_exactly_one_allowed_claim_label(signal):
    root = out_dir(signal)
    for path in root.rglob("*.tsv"):
        rows = F.read_tsv(path)
        if not rows:
            continue
        assert "claim_label" in rows[0], path.name
        assert {r["claim_label"] for r in rows} <= CLAIMS, path.name
        text = path.read_text(encoding="utf-8").lower()
        for token in ("cross_validated_nested", "independently_validated", "not externally validated"):
            text = text.replace(token, "")
        assert "diagnostic" not in text and "clinical utility" not in text and "validated" not in text, path.name
    validation = F.read_tsv(root / "validation.tsv")[0]
    assert validation["state"] == "not externally validated"
    data = json.loads((signal / "report" / "report_data.json").read_text(encoding="utf-8"))
    section = data["sections"]["post_de_biomarker"]
    assert section["state"] == "COMPLETED" and section["claim_label"] in CLAIMS


def test_v149_negative_forbidden_claims_and_missing_labels_fail(tmp_path):
    bad = tmp_path / "bad.tsv"; bad.write_text("feature\tclaim_label\nF1\tin_sample diagnostic marker\n", encoding="utf-8")
    good = tmp_path / "good.tsv"; good.write_text("feature\tclaim_label\nF1\tcross_validated_nested\n", encoding="utf-8")
    out = B.r_json("""a <- commandArgs(TRUE); ns <- asNamespace('proteomicsCore')
      e1 <- tryCatch({ get('pd_bm_claim_audit', ns)(a[1]); 'none' }, error = function(e) conditionMessage(e))
      e2 <- tryCatch({ get('pd_bm_claim_audit', ns)(a[2]); 'none' }, error = function(e) conditionMessage(e))
      e3 <- tryCatch({ get('.pd_check_claims', ns)(data.frame(x = 1)); 'none' }, error = function(e) conditionMessage(e))
      cat(jsonlite::toJSON(list(bad = e1, good = e2, missing = e3), auto_unbox = TRUE))""", bad, good)
    assert out["bad"].startswith("E_BIOMARKER_CLAIM") and out["good"] == "none" and out["missing"].startswith("E_BIOMARKER_CLAIM")


# ----------------------------------------------------------------------------- V150
def test_v150_single_feature_auc_and_cis_equal_proc(signal):
    rows = {r["feature_id"]: r for r in F.read_tsv(out_dir(signal) / "single_feature_auc.tsv")}
    oracle = B.r_json("""
      a <- commandArgs(TRUE); Y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
      o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(Y), o$observation_id), ]
      out <- list()
      for (f in c('F01', 'F02', 'F04', 'F05', 'F20')) {
        r <- pROC::roc(controls = Y[f, o$group == 'N'], cases = Y[f, o$group == 'P'], direction = '<', quiet = TRUE)
        d <- pROC::ci.auc(r, method = 'delong')
        set.seed(11, kind = "L'Ecuyer-CMRG"); b <- pROC::ci.auc(r, method = 'bootstrap', boot.n = 200, boot.stratified = TRUE)
        out[[f]] <- list(auc = as.numeric(r$auc), dl = d[1], du = d[3], bl = b[1], bu = b[3]) }
      cat(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA))""", signal / "preprocessing" / "primary" / "matrix.tsv", signal / "preprocessing" / "primary" / "observations.tsv")
    for f, o in oracle.items():
        r = rows[f]
        assert abs(float(r["auc"]) - o["auc"]) <= 1e-12
        assert abs(float(r["delong_lower"]) - o["dl"]) <= 1e-10 and abs(float(r["delong_upper"]) - o["du"]) <= 1e-10
        assert abs(float(r["bootstrap_lower"]) - o["bl"]) <= 1e-12 and abs(float(r["bootstrap_upper"]) - o["bu"]) <= 1e-12
        assert r["claim_label"] == "in_sample" and r["direction_optimistic"] == "false"
    assert float(rows["F04"]["auc"]) < 0.5                                       # prespecified direction: an opposite-direction feature


def test_v150_auto_direction_is_in_sample_and_flagged(tmp_path):
    values, obs = F.biomarker_design()
    post = biomarker_block(single_feature={"direction": "train_only", "bootstrap": 200}, permutation={"enabled": False})
    _, _, out = run(tmp_path, values, obs, post)
    rows = {r["feature_id"]: r for r in F.read_tsv(out_dir(out) / "single_feature_auc.tsv")}
    assert all(r["claim_label"] == "in_sample" and r["direction_optimistic"] == "true" for r in rows.values())   # never unlabelled
    assert float(rows["F04"]["auc"]) > 0.5 and rows["F04"]["direction"] == "lower_in_positive"


# ----------------------------------------------------------------------------- V151
def test_v151_minimum_size_and_fold_refusal_report_both_reasons(tmp_path):
    """D-59: 4 positive units (below the recommended 5) run as exploratory with outer folds adapted to 4; fewer than 3
    units per class are not computable and are refused (no AUC); a required module then rejects the plan."""
    values, obs = F.biomarker_design(n_pos=4, n_neg=30)
    payload, code, out = run(tmp_path / "small", values, obs, biomarker_block(permutation={"enabled": False}))
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_biomarker")
    assert stage["state"] == "COMPLETED"
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    adaptations = {a["item"]: a for a in next(c for c in plan["capability_plan"] if c["capability"] == "post_de_biomarker")["post_de"]["adaptations"]}
    assert adaptations["units_per_class"]["used"] == "4 (exploratory)" and (adaptations["outer_k"]["requested"], adaptations["outer_k"]["used"]) == (5, 4)
    cv = F.read_tsv(out_dir(out) / "cv_performance.tsv")[0]
    assert cv["outer_k"] == "4" and cv["claim_label"] == "cross_validated_nested"
    values, obs = F.biomarker_design(n_pos=2, n_neg=30)
    payload, code, out = run(tmp_path / "opt", values, obs, biomarker_block())
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_biomarker")
    assert stage["state"] == "INAPPLICABLE" and stage["reason_code"] == "E_BIOMARKER_SMALL_N"
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    assert "at least 4 units per class" in next(c for c in plan["capability_plan"] if c["capability"] == "post_de_biomarker")["reason"]
    assert not (out / "post_de" / "biomarker").exists()                          # no AUC emitted
    post = biomarker_block(); post["biomarker"]["execution_requirement"] = "required"
    payload, code, out = run(tmp_path / "req", values, obs, post, expect=2)
    assert payload["error"]["code"] == "E_BIOMARKER_SMALL_N" and not (out / "dea").exists()


# ----------------------------------------------------------------------------- V152
@pytest.fixture(scope="module")
def noise(tmp_path_factory):
    values, obs = F.biomarker_design(n_pos=20, n_neg=20, n_features=2000, up=0, down=0, seed=52)
    rows = [(f, v) for f, v in values.items()]
    pos = [i for i, o in enumerate(obs) if o["group"] == "P"]; neg = [i for i, o in enumerate(obs) if o["group"] == "N"]
    def auc(v):
        return sum((v[i] > v[j]) + 0.5 * (v[i] == v[j]) for i in pos for j in neg) / (len(pos) * len(neg))
    top = [f for f, _ in sorted(rows, key=lambda fv: -abs(auc(fv[1]) - 0.5))[:10]]                # chosen on the full data
    panels = [{"id": "same", "feature_ids": top, "provenance": "same_data"}, {"id": "indep", "feature_ids": [f"F{i:04d}" for i in range(1001, 1011)], "provenance": "independent"}]
    post = biomarker_block(cv={"outer": {"k": 5, "repeats": 5}, "inner": {"k": 3}}, selection={"method": "top_k_auc", "k": [10]}, fixed_panels=panels,
                           single_feature={"enabled": False})
    return run(tmp_path_factory.mktemp("noise"), values, obs, post, threads=4)[2], top


def test_v152_nested_estimate_near_half_while_leaky_reference_is_inflated(noise):
    out, _ = noise
    cv = F.read_tsv(out_dir(out) / "cv_performance.tsv")[0]
    repeats = [float(r["auc"]) for r in F.read_tsv(out_dir(out) / "cv_repeats.tsv")]
    mean, sd = statistics.mean(repeats), statistics.stdev(repeats)
    assert mean - 1.96 * sd <= 0.5 <= mean + 1.96 * sd and abs(float(cv["pooled_oof_auc"]) - 0.5) < 0.2
    leaky = B.r_json("""
      a <- commandArgs(TRUE); Y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
      o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(Y), o$observation_id), ]
      X <- scale(t(Y)); y <- as.integer(o$group == 'P')
      au <- apply(X, 2, function(s) { r <- rank(c(s[y == 1], s[y == 0])); (sum(r[1:20]) - 210) / 400 }); top <- names(sort(-abs(au - 0.5)))[1:10]
      set.seed(1); res <- c()
      for (r in 1:5) { fold <- sample(rep(1:5, length.out = 40)); sc <- numeric(40)
        for (f in 1:5) { tr <- fold != f; fit <- glmnet::glmnet(X[tr, top], y[tr], family = 'binomial', alpha = 0, lambda = 0.1)
          sc[!tr] <- predict(fit, X[!tr, top], s = 0.1, type = 'response') }
        rr <- rank(c(sc[y == 1], sc[y == 0])); res <- c(res, (sum(rr[1:20]) - 210) / 400) }
      cat(jsonlite::toJSON(list(mean = mean(res)), auto_unbox = TRUE, digits = NA))""", out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv")
    assert leaky["mean"] > 0.75                                                   # selection and scaling on all data, then CV
    audit = [json.loads(line) for line in (out_dir(out) / "transform_audit.jsonl").read_text(encoding="utf-8").splitlines()]
    oof = F.read_tsv(out_dir(out) / "oof_predictions.tsv")
    test_ids = {}
    for r in oof:
        test_ids.setdefault((int(r["repeat_id"]), int(r["fold"])), set()).add(r["observation_id"])
    transforms = {a["transform"] for a in audit}
    assert {"scaling", "selection", "tuning", "classifier", "threshold"} <= transforms
    assert len(audit) == 5 * 5 * len(transforms)
    for a in audit:                                                               # every transform fitted on the training fold only
        assert not set(a["fitted_on"]) & test_ids[(a["repeat_id"], a["fold"])] and set(a["applied_to"]) == test_ids[(a["repeat_id"], a["fold"])]


def test_v152_negative_a_transform_fitted_on_held_out_samples_fails(tmp_path):
    out = B.r_json("""ns <- asNamespace('proteomicsCore')
      audit <- list(list(repeat_id = 1L, fold = 1L, transform = 'scaling', fitted_on = c('a', 'b', 'c'), applied_to = c('c')))
      msg <- tryCatch({ get('pd_bm_audit_check', ns)(audit, list(`1 1` = 'c')); 'none' }, error = function(e) conditionMessage(e))
      cat(jsonlite::toJSON(msg, auto_unbox = TRUE))""")
    assert out.startswith("E_BIOMARKER_LEAKAGE")
    values, obs = F.biomarker_design()
    _leak_replaced(tmp_path, values, obs, biomarker_block(imputation="global_median", permutation={"enabled": False}), "E_BIOMARKER_LEAKAGE", "imputation", "train_median")


# ----------------------------------------------------------------------------- V153 / V154
def test_v153_v154_nested_cv_and_selection_equal_an_independent_loop(signal):
    oracle = B.r_json(NESTED_ORACLE, signal / "preprocessing" / "primary" / "matrix.tsv", signal / "preprocessing" / "primary" / "observations.tsv", "3", "observation", str(SEED))
    oof = {(int(r["repeat_id"]), r["observation_id"]): float(r["score"]) for r in F.read_tsv(out_dir(signal) / "oof_predictions.tsv")}
    assert len(oof) == len(oracle["id"])
    for i, r, s in zip(oracle["id"], oracle["r"], oracle["score"]):
        assert abs(oof[(r, i)] - s) <= 1e-8, (r, i)
    repeats = [float(r["auc"]) for r in F.read_tsv(out_dir(signal) / "cv_repeats.tsv")]
    assert all(abs(a - b) <= 1e-12 for a, b in zip(repeats, oracle["per"]))
    cv = F.read_tsv(out_dir(signal) / "cv_performance.tsv")[0]
    assert abs(float(cv["repeat_auc_mean"]) - statistics.mean(oracle["per"])) <= 1e-12 and cv["claim_label"] == "cross_validated_nested"
    freq = {r["feature_id"]: float(r["selection_frequency"]) for r in F.read_tsv(out_dir(signal) / "selection_stability.tsv")}
    for f, v in oracle["freq"].items():
        assert abs(freq[f] - (v[0] if isinstance(v, list) else v)) <= 1e-12
    planted = ["F01", "F02", "F03", "F04", "F05", "F06"]
    assert sum(freq[f] >= 0.9 for f in planted) >= 4 and sum(v < 0.2 for f, v in freq.items() if f not in planted) >= 0.9 * (len(freq) - 6)
    assert freq != {f: (1.0 if f in sorted(freq, key=lambda x: -freq[x])[:5] else 0.0) for f in freq}   # not one full-data selection


def test_v153_repeated_design_folds_never_split_a_subject(tmp_path):
    values, obs = F.biomarker_design(n_pos=8, n_neg=8, reps=2, seed=53)
    post = biomarker_block(permutation={"B": 19}, single_feature={"bootstrap": 200})
    _, _, out = run(tmp_path, values, obs, post, blocked=True)
    oof = F.read_tsv(out_dir(out) / "oof_predictions.tsv")
    folds = {}
    for r in oof:
        folds.setdefault((r["repeat_id"], r["unit_id"]), set()).add(r["fold"])
    assert all(len(v) == 1 for v in folds.values()) and {r["unit_id"] for r in oof} == {o["subject_id"] for o in obs}
    oracle = B.r_json(NESTED_ORACLE, out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv", "3", "subject", str(SEED))
    scores = {(int(r["repeat_id"]), r["observation_id"]): float(r["score"]) for r in oof}
    assert all(abs(scores[(r, i)] - s) <= 1e-8 for i, r, s in zip(oracle["id"], oracle["r"], oracle["score"]))
    cv = F.read_tsv(out_dir(out) / "cv_performance.tsv")[0]
    assert cv["grouped_by"] == "subject_id" and cv["permutation_scheme"] == "whole_subjects:subject_id"


def test_v153_negative_ungrouped_folds_in_a_repeated_design_fail(tmp_path):
    values, obs = F.biomarker_design(n_pos=8, n_neg=8, reps=2, seed=53)
    out_run = _leak_replaced(tmp_path, values, obs, biomarker_block(cv={"outer": {"k": 5, "repeats": 2}, "inner": {"k": 3}, "group_by_subject": False}, permutation={"enabled": False}),
                             "E_BIOMARKER_GROUP_LEAKAGE", "cv.group_by_subject", "true", blocked=True)
    assert F.read_tsv(out_dir(out_run) / "cv_performance.tsv")[0]["grouped_by"] == "subject_id"   # folds kept whole subjects
    out = B.r_json("""ns <- asNamespace('proteomicsCore'); X <- matrix(rnorm(8), 4, 2, dimnames = list(NULL, c('a', 'b')))
      msg <- tryCatch({ get('.bm_outer_fold', ns)(X, c(1L, 0L, 1L, 0L), c('s1', 's1', 's2', 's3'), 1:2, 2:4, list(), 1L, 1L, NULL); 'none' }, error = function(e) conditionMessage(e))
      cat(jsonlite::toJSON(msg, auto_unbox = TRUE))""")
    assert out.startswith("E_BIOMARKER_GROUP_LEAKAGE")


# ----------------------------------------------------------------------------- V155
def test_v155_auc_does_not_depend_on_class_level_order(tmp_path, signal):
    values, obs = F.biomarker_design()
    _, _, out = run(tmp_path, values, obs, biomarker_block(permutation={"enabled": False}), groups=("N", "P"))
    a = float(F.read_tsv(out_dir(signal) / "cv_performance.tsv")[0]["pooled_oof_auc"]); b = float(F.read_tsv(out_dir(out) / "cv_performance.tsv")[0]["pooled_oof_auc"])
    assert abs(a - b) <= 1e-10
    sa = {(r["repeat_id"], r["observation_id"]): float(r["score"]) for r in F.read_tsv(out_dir(signal) / "oof_predictions.tsv")}
    sb = {(r["repeat_id"], r["observation_id"]): float(r["score"]) for r in F.read_tsv(out_dir(out) / "oof_predictions.tsv")}
    assert all(abs(sa[k] - sb[k]) <= 1e-6 for k in sa)                          # scores oriented from the model's class coding


def test_v155_svm_uses_oriented_decision_values(tmp_path):
    values, obs = F.biomarker_design(shift=2.0)
    _, _, out = run(tmp_path, values, obs, biomarker_block(classifier="svm_linear", classifier_params={"cost": [1]}, permutation={"enabled": False}))
    cv = F.read_tsv(out_dir(out) / "cv_performance.tsv")[0]
    assert float(cv["pooled_oof_auc"]) > 0.9                                     # an orientation flip would give 1 - AUC
    scores = [float(r["score"]) for r in F.read_tsv(out_dir(out) / "oof_predictions.tsv")]
    assert len(set(scores)) > 2                                                   # decision values, not predicted classes
    refusals = F.read_tsv(out_dir(out) / "refusals.tsv")
    assert any(r["reason_code"] == "W_CALIBRATION_NOT_PROBABILISTIC" for r in refusals)


def test_v155_random_forest_is_not_run_when_absent(tmp_path, monkeypatch):
    lib = tmp_path / "shadow-lib"; (lib / "randomForest").mkdir(parents=True)
    (lib / "randomForest" / "DESCRIPTION").write_text("Package: randomForest\nVersion: 0.0.0\nTitle: shadow\n", encoding="utf-8")
    previous = os.environ.get("R_LIBS")
    monkeypatch.setenv("R_LIBS", os.pathsep.join([str(lib)] + ([previous] if previous else [])))
    values, obs = F.biomarker_design()
    payload, code, out = run(tmp_path, values, obs, biomarker_block(classifier="random_forest", permutation={"enabled": False}), expect=None)
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_biomarker")
    assert stage["state"] == "NOT_RUN" and stage["reason_code"] == "E_ENGINE_NOT_AVAILABLE" and payload["state"] == "PARTIAL"


# ----------------------------------------------------------------------------- V156
def test_v156_whole_procedure_permutation(signal, noise):
    s = F.read_tsv(out_dir(signal) / "cv_performance.tsv")[0]
    null = [float(r["pooled_oof_auc"]) for r in F.read_tsv(out_dir(signal) / "permutation_null.tsv")]
    k = sum(a >= float(s["pooled_oof_auc"]) - 1e-12 for a in null)
    assert len(null) == 19 and int(s["permutation_k"]) == k and abs(float(s["permutation_p"]) - (k + 1) / 20) <= 1e-15
    assert float(s["permutation_p"]) == 1 / 20                                    # signal: P at the floor 1/(B+1)
    assert s["permutation_scope"].startswith("whole_procedure") and s["permutation_scheme"] == "biological_units"
    n = F.read_tsv(out_dir(noise[0]) / "cv_performance.tsv")[0]
    assert float(n["permutation_p"]) > 0.1                                        # noise: a non-small P
    assert len({round(a, 12) for a in null}) > 1                                  # selection and tuning really rerun


def test_v156_negative_permuting_only_the_final_classifier_is_refused(tmp_path):
    values, obs = F.biomarker_design()
    out_run = _leak_replaced(tmp_path, values, obs, biomarker_block(permutation={"B": 19, "scope": "final_classifier_only"}),
                             "E_BIOMARKER_PERMUTATION_SCOPE", "permutation.scope", "whole_procedure")
    assert F.read_tsv(out_dir(out_run) / "cv_performance.tsv")[0]["permutation_scope"].startswith("whole_procedure")
    out = B.r_json("""ns <- asNamespace('proteomicsCore'); set.seed(3, kind = "L'Ecuyer-CMRG")
      units <- rep(paste0('s', 1:6), each = 2); y <- rep(c(1L, 0L), each = 6)
      ok <- all(replicate(50, { p <- get('pd_bm_permute_labels', ns)(y, units); all(tapply(p, units, function(v) length(unique(v))) == 1) }))
      cat(jsonlite::toJSON(ok, auto_unbox = TRUE))""")
    assert out is True                                                            # block-restricted: whole subjects keep one label


# ----------------------------------------------------------------------------- V157
def test_v157_fixed_panels_and_circularity(noise):
    out, top = noise
    panels = {r["panel_id"]: r for r in F.read_tsv(out_dir(out) / "fixed_panels.tsv")}
    same, indep = panels["same"], panels["indep"]
    assert same["claim_label"] == indep["claim_label"] == "fixed_panel_cv"
    assert same["selection_optimistic"] == "true" and float(same["pooled_oof_auc"]) > 0.7                       # optimistic
    assert same["nested_reselection_auc"] != "NA" and abs(float(same["nested_reselection_auc"]) - 0.5) < 0.2   # beside the nested estimate
    assert indep["selection_optimistic"] == "false" and indep["nested_reselection_auc"] == "NA"
    assert same["panel_sha256"] == hashlib.sha256("\n".join(sorted(top)).encode()).hexdigest()


# ----------------------------------------------------------------------------- V158
def test_v158_calibration_and_threshold_metrics_equal_the_oracle(signal):
    oof = F.read_tsv(out_dir(signal) / "oof_predictions.tsv")
    p = [float(r["score"]) for r in oof]; y = [int(r["label"]) for r in oof]
    cv = F.read_tsv(out_dir(signal) / "cv_performance.tsv")[0]
    assert abs(float(cv["brier"]) - sum((a - b) ** 2 for a, b in zip(p, y)) / len(p)) <= 1e-12
    cal = F.read_tsv(out_dir(signal) / "calibration.tsv")
    for row in cal:
        i = int(row["bin"]); members = [(a, b) for a, b in zip(p, y) if min(int(a * 10) + 1, 10) == i]
        assert int(row["n"]) == len(members)
        if members:
            assert abs(float(row["mean_predicted"]) - sum(a for a, _ in members) / len(members)) <= 1e-12
            assert abs(float(row["observed_fraction"]) - sum(b for _, b in members) / len(members)) <= 1e-12
    predicted = [int(r["predicted"]) for r in oof]
    assert all(pr == int(float(r["score"]) >= float(r["threshold"])) for pr, r in zip(predicted, oof))
    sens = sum(a == 1 for a, b in zip(predicted, y) if b == 1) / sum(y); spec = sum(a == 0 for a, b in zip(predicted, y) if b == 0) / (len(y) - sum(y))
    assert abs(float(cv["sensitivity"]) - sens) <= 1e-12 and abs(float(cv["specificity"]) - spec) <= 1e-12
    boot = B.r_json("""a <- commandArgs(TRUE); oof <- read.delim(a[1], colClasses = 'character', encoding = 'UTF-8')
      lab <- as.integer(oof$label); pred <- as.integer(oof$predicted); units <- oof$unit_id
      u <- sort(unique(units), method = 'radix'); cls <- vapply(u, function(x) as.character(lab[units == x][1]), '')
      set.seed(11, kind = "L'Ecuyer-CMRG"); s <- .Random.seed; for (i in 1:20) s <- parallel::nextRNGStream(s); assign('.Random.seed', s, envir = globalenv())
      b <- t(vapply(1:200, function(i) { pick <- unlist(lapply(split(u, cls), function(x) x[sample.int(length(x), length(x), replace = TRUE)]))
        rows <- unlist(lapply(pick, function(x) which(units == x))); c(mean(pred[rows][lab[rows] == 1] == 1), mean(pred[rows][lab[rows] == 0] == 0)) }, numeric(2)))
      cat(jsonlite::toJSON(list(sl = quantile(b[, 1], 0.025, names = FALSE), su = quantile(b[, 1], 0.975, names = FALSE), pl = quantile(b[, 2], 0.025, names = FALSE)), auto_unbox = TRUE, digits = NA))""",
                    out_dir(signal) / "oof_predictions.tsv")
    assert abs(float(cv["sensitivity_lower"]) - boot["sl"]) <= 1e-12 and abs(float(cv["sensitivity_upper"]) - boot["su"]) <= 1e-12 and abs(float(cv["specificity_lower"]) - boot["pl"]) <= 1e-12
    audit = [json.loads(line) for line in (out_dir(signal) / "transform_audit.jsonl").read_text(encoding="utf-8").splitlines() if '"threshold"' in line]
    assert audit and all("inner out-of-fold scores of the training fold" in a["detail"] for a in audit)


def test_v158_negative_test_fold_threshold_is_refused(tmp_path):
    values, obs = F.biomarker_design()
    out_run = _leak_replaced(tmp_path, values, obs, biomarker_block(threshold_rule="youden_test", permutation={"enabled": False}),
                             "E_BIOMARKER_THRESHOLD_LEAKAGE", "threshold_rule", "youden_train")
    assert F.read_tsv(out_dir(out_run) / "cv_performance.tsv")[0]["threshold_rule"] == "youden_train"


# ----------------------------------------------------------------------------- V159
def test_v159_locked_model_evaluated_once_on_a_disjoint_cohort(tmp_path):
    values, obs = F.biomarker_design()
    vvalues, vobs = F.biomarker_design(seed=59, prefix="V")
    cohort = F.write_validation_cohort(tmp_path / "cohort", vvalues, vobs)
    _, _, out = run(tmp_path, values, obs, biomarker_block(validation_cohort=cohort, permutation={"enabled": False}))
    row = F.read_tsv(out_dir(out) / "validation.tsv")[0]
    assert row["state"] == "evaluated_once" and row["claim_label"] == "independently_validated"
    locked_bytes = (out_dir(out) / "locked_model.json").read_bytes()
    assert hashlib.sha256(locked_bytes).hexdigest() == row["locked_model_sha256"]
    oracle = B.r_json("""a <- commandArgs(TRUE); m <- jsonlite::fromJSON(a[1]); V <- as.matrix(read.delim(a[2], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
      meta <- read.delim(a[3], colClasses = 'character', encoding = 'UTF-8'); f <- m$features
      Z <- sweep(sweep(t(V[f, meta$sample_id, drop = FALSE]), 2, unlist(m$center)[f]), 2, unlist(m$scale)[f], '/')
      b <- unlist(m$classifier$coefficients); eta <- b[1] + Z %*% b[f]; p <- 1 / (1 + exp(-eta)); if (m$classifier$sign < 0) p <- 1 - p
      r <- pROC::roc(controls = p[meta$group == 'N'], cases = p[meta$group == 'P'], direction = '<', quiet = TRUE)
      cat(jsonlite::toJSON(list(auc = as.numeric(r$auc), sens = mean(p[meta$group == 'P'] >= m$threshold)), auto_unbox = TRUE, digits = NA))""",
                      out_dir(out) / "locked_model.json", cohort["matrix"], cohort["metadata"])
    assert abs(float(row["auc"]) - oracle["auc"]) <= 1e-10 and abs(float(row["sensitivity"]) - oracle["sens"]) <= 1e-12
    assert row["matrix_sha256"] == hashlib.sha256(Path(cohort["matrix"]).read_bytes()).hexdigest()


def test_v159_negative_overlapping_cohort_and_retuning_fail(tmp_path):
    values, obs = F.biomarker_design()
    vvalues, vobs = F.biomarker_design(seed=59, prefix="V")
    for o, new in zip(vobs[:3], ("P1", "P2", "N1")):
        o["observation_id"] = new; o["biological_unit_id"] = new
    cohort = F.write_validation_cohort(tmp_path / "cohort", vvalues, vobs)
    post = biomarker_block(validation_cohort=cohort, permutation={"enabled": False}); post["biomarker"]["execution_requirement"] = "required"
    payload, code, out = run(tmp_path / "req", values, obs, post, expect=2)
    assert payload["error"]["code"] == "E_SCORE_SELECTION_OVERLAP"
    clean = F.write_validation_cohort(tmp_path / "clean", *F.biomarker_design(seed=59, prefix="V"))
    out_run = _leak_replaced(tmp_path / "retune", values, obs, biomarker_block(validation_cohort={**clean, "retune": True}, permutation={"enabled": False}),
                             "E_VALIDATION_RETUNED", "validation_cohort.retune", "false")
    assert F.read_tsv(out_dir(out_run) / "validation.tsv")[0]["state"] == "evaluated_once"


# ----------------------------------------------------------------------------- V160
def test_v160_figures_have_sources_and_the_report_shows_the_design(signal):
    eligibility = json.loads((out_dir(signal) / "eligibility.json").read_text(encoding="utf-8"))
    stems = {Path(f["file"]).stem for f in eligibility["figures"] if f["state"] == "COMPLETED"}
    assert {"roc_cv", "cv_auc_repeats", "permutation_null", "calibration", "selection_stability"} <= stems
    for stem in stems:
        assert (out_dir(signal) / "figure_sources" / f"{stem}.tsv").is_file(), stem
    roc = F.read_tsv(out_dir(signal) / "figure_sources" / "roc_cv.tsv")
    oof = F.read_tsv(out_dir(signal) / "oof_predictions.tsv")
    for point in roc[1:]:
        t = float(point["threshold"])
        assert abs(float(point["tpr"]) - sum(float(r["score"]) >= t for r in oof if r["label"] == "1") / sum(r["label"] == "1" for r in oof)) <= 1e-12
    data = json.loads((signal / "report" / "report_data.json").read_text(encoding="utf-8"))
    section = data["sections"]["post_de_biomarker"]
    text = json.dumps(section)
    for needle in ("cross_validated_nested", "in_sample", "repeated_stratified_kfold", "seed 11", "B = 19", "15 P", "15 N", "not externally validated"):
        assert needle in text, needle
    for figure in section["figures"]:
        assert (signal / figure["source"]).is_file() and (signal / figure["src"]).is_file()


def test_v156_compute_guard_refuses_instead_of_reducing_work(tmp_path):
    """D-59: default compute_policy adapt scales repeats and B down to the guard (documented minimums), recorded with the
    attainable P floor; compute_policy refuse keeps the strict refusal and never lowers B."""
    values, obs = F.biomarker_design()
    payload, code, out = run(tmp_path / "adapt", values, obs, biomarker_block(cv={"outer": {"k": 5, "repeats": 3}, "inner": {"k": 3}}, permutation={"B": 400, "compute_guard_hours": 1e-6}))
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_biomarker")
    assert stage["state"] == "COMPLETED"
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    record = next(a for a in next(c for c in plan["capability_plan"] if c["capability"] == "post_de_biomarker")["post_de"]["adaptations"] if a["item"] == "nested_cv_permutation")
    assert record["requested"].startswith("repeats 3, B 400") and record["used"].startswith("repeats 3, B 99") and "1/100" in record["reason"]
    cv = F.read_tsv(out_dir(out) / "cv_performance.tsv")[0]
    assert cv["permutation_B"] == "99" and len(F.read_tsv(out_dir(out) / "permutation_null.tsv")) == 99
    assert float(cv["permutation_p"]) >= 1 / 100
    payload, code, out = run(tmp_path / "refuse", values, obs, biomarker_block(permutation={"B": 5000, "compute_guard_hours": 0.001}, compute_policy="refuse"))
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_biomarker")
    assert stage["state"] == "INAPPLICABLE" and stage["reason_code"] == "E_BIOMARKER_COMPUTE_GUARD"
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    entry = next(c for c in plan["capability_plan"] if c["capability"] == "post_de_biomarker")
    assert "5000 permutations" in entry["reason"] and not (out / "post_de" / "biomarker").exists()   # strict mode: B is never lowered