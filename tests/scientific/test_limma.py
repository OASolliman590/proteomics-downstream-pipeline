"""R05 acceptance: V041-V050 (observed-data limma, exact contrasts, moderation, TREAT, CIs, families,
omnibus, typed tables, influence, golden suite).  Oracles are separate direct limma/GLS calls."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

from proteomics_pipeline import workflow

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")
EXAMPLE = B.EXAMPLES / "example-independent.json"
ABS_SE, REL_SE = 1e-8, 1e-6


def close(a, b, abs_tol=ABS_SE, rel_tol=REL_SE):
    return abs(a - b) <= abs_tol or abs(a - b) <= rel_tol * max(abs(a), abs(b))


def rows_by(path, *keys):
    return {tuple(r[k] for k in keys): r for r in B.read_tsv(path)}


def run_example(tmp_path, mutate=None, name="example-independent.json"):
    config = json.loads((B.EXAMPLES / name).read_text())
    for key in ("matrix", "observations", "features", "source_provenance"):
        config["input"][key] = str((B.EXAMPLES / config["input"][key]).resolve())
    if mutate:
        mutate(config)
    path = tmp_path / "analysis.json"; path.write_text(json.dumps(config))
    payload, code = workflow.run_command(path, tmp_path / "run")
    return payload, code, tmp_path / "run"


REFERENCE = """
args <- commandArgs(TRUE); m <- as.matrix(read.delim(args[1], row.names=1, check.names=FALSE, na.strings='NA'))
g <- factor(sub('[0-9]+$', '', colnames(m)), levels=c('C','U','T'))
out <- list()
for (spec in list(c('disease-control','C','U'), c('treated-disease','U','T'), c('treated-control','C','T'))) {
  X <- model.matrix(~ relevel(g, ref=spec[2])); colnames(X) <- c('int', levels(relevel(g, ref=spec[2]))[-1])
  fit <- limma::eBayes(limma::lmFit(m, X), trend=TRUE, robust=TRUE)
  k <- spec[3]
  out[[spec[1]]] <- list(effect=unname(fit$coefficients[,k]), se=unname(sqrt(fit$s2.post)*fit$stdev.unscaled[,k]), t=unname(fit$t[,k]),
                          p=unname(fit$p.value[,k]), df=unname(fit$df.total), dfres=unname(fit$df.residual), ids=rownames(m))
}
cat(jsonlite::toJSON(out, digits=NA, auto_unbox=FALSE))
"""


# ----------------------------------------------------------------------------- V041 / V045
def test_v041_v045_matches_direct_lmfit_ebayes_with_original_na(tmp_path):
    payload, code, run = run_example(tmp_path)
    assert payload["stages"][3] == {"stage_id": "limma", "state": "COMPLETED", "reason_code": None}
    reference = B.r_json(REFERENCE, str(run / "inputs" / "matrix.tsv"))
    table = rows_by(run / "dea" / "zero_null.tsv", "contrast_id", "feature_id")
    for contrast, ref in reference.items():
        for i, feature in enumerate(ref["ids"]):
            row = table[(contrast, feature)]
            assert row["eligibility"] == "tested" and row["exactness_path"] == "contrast_as_coefficient_refit"
            assert close(float(row["effect"]), ref["effect"][i], 1e-10, 1e-8)
            for field, key in (("effect_se", "se"), ("statistic", "t"), ("p_value", "p"), ("df_inference", "df")):
                assert close(float(row[field]), ref[key][i]), (contrast, feature, field)
            q = 0.975
            from statistics import NormalDist
            half = float(row["ci_upper"]) - float(row["effect"])
            # V045 oracle: qt(.975, df_total) * moderated SE, computed by R for this df
            assert row["ci_level"] == "0.95" and "qt(" in row["ci_method"]
            assert abs(half - (float(row["ci_upper"]) - float(row["ci_lower"])) / 2) <= 1e-12
            assert abs(half - NormalDist().inv_cdf(q) * float(row["effect_se"])) > 1e-6     # not a normal-quantile interval
    qt = B.r_json("cat(jsonlite::toJSON(qt(0.975, as.numeric(commandArgs(TRUE))), digits=NA))", table[("disease-control", "P01")]["df_inference"])
    row = table[("disease-control", "P01")]
    assert close(float(row["ci_upper"]), float(row["effect"]) + qt[0] * float(row["effect_se"]), 1e-10, 1e-8)
    # the two original NA cells stay unobserved: P07 (U2 missing) has n=11 and df_residual 8
    assert table[("disease-control", "P07")]["df_residual"] == "8" and json.loads(table[("disease-control", "P07")]["n_obs_by_required_group"]) == {"disease": 3, "control": 4}


def test_v041_negative_a_filling_adapter_would_be_detected(tmp_path):
    payload, code, run = run_example(tmp_path)
    filled = B.r_json("m <- as.matrix(read.delim(commandArgs(TRUE)[1], row.names=1, check.names=FALSE)); for (i in seq_len(nrow(m))) m[i, is.na(m[i,])] <- mean(m[i,], na.rm=TRUE); "
                      "g <- factor(sub('[0-9]+$','',colnames(m)), levels=c('C','U','T')); X <- model.matrix(~ g); fit <- limma::eBayes(limma::lmFit(m, X), trend=TRUE, robust=TRUE); "
                      "cat(jsonlite::toJSON(unname(fit$coefficients['P07','gU']), digits=NA))", str(run / "inputs" / "matrix.tsv"))
    ours = float(rows_by(run / "dea" / "zero_null.tsv", "contrast_id", "feature_id")[("disease-control", "P07")]["effect"])
    assert abs(ours - filled[0]) > 1e-6


# ----------------------------------------------------------------------------- V042
def _weighted_dataset(tmp_path, blocking):
    x = [1.0, 2.5, 3.1, 4.0, 2.2, 3.5, 5.1, 4.4, 6.0, 5.2]
    w = [1, 2, 0.5, 1.5, 1, 3, 0.7, 1, 2, 1]
    obs = [{"observation_id": f"o{i}", "group": "A" if i < 5 else "B", "x": str(x[i]), "subject_id": f"s{i % 5}"} for i in range(10)]
    noise = [[0.1, -0.2, 0.05, 0.12, -0.07, 0.2, -0.1, 0.03, -0.15, 0.08], [-0.05, 0.11, -0.2, 0.07, 0.02, -0.09, 0.14, -0.03, 0.06, -0.12],
             [0.3, -0.1, 0.2, -0.25, 0.05, 0.1, -0.2, 0.15, -0.05, -0.2], [0.02, 0.04, -0.06, 0.1, -0.1, 0.05, -0.02, 0.07, -0.08, -0.02]]
    values = {}
    for f in range(4):
        values[f"F{f + 1}"] = [10 + f + (0.8 if i >= 5 else 0) * (f % 2) + 0.3 * x[i] + noise[f][i] for i in range(10)]
    values["F2"][6] = None
    files = B.dataset(tmp_path, values, obs, extra_columns=("x",))
    B.write_tsv(tmp_path / "weights.tsv", [["observation_id", "weight"]] + [[f"o{i}", repr(w[i])] for i in range(10)])
    design = {"continuous_covariates": [{"column": "x", "center": True}], "blocking": blocking}
    model = {"precision_weights": {"kind": "observation", "path": "weights.tsv"}}
    return B.config(tmp_path, files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], design_overrides=design, model_overrides=model), x, w


ORACLE_GLS = """
a <- commandArgs(TRUE); m <- as.matrix(read.delim(a[1], row.names=1, check.names=FALSE)); x <- as.numeric(strsplit(a[2], ',')[[1]]); w <- as.numeric(strsplit(a[3], ',')[[1]])
rho <- as.numeric(a[4]); grp <- rep(0:1, each=5); X <- cbind(A=1-grp, B=grp, x=x-mean(x)); cc <- c(-1,1,0)
Tm <- rbind(cc, c(1,1,0)/sqrt(2), c(0,0,1)); Xs <- X %*% solve(Tm)
block <- if (is.finite(rho)) rep(0:4, 2) else NULL
fit <- if (is.null(block)) limma::lmFit(m, Xs, weights=w) else limma::lmFit(m, Xs, weights=w, block=block, correlation=rho)
eb <- limma::eBayes(fit, trend=TRUE, robust=TRUE)
unscaled <- sapply(seq_len(nrow(m)), function(i) { u <- !is.na(m[i,]); n <- sum(u); V <- diag(n); if (!is.null(block)) { b <- block[u]; V[outer(b,b,'==')] <- rho; diag(V) <- 1 }
  S <- diag(sqrt(w[u]), n) %*% solve(V) %*% diag(sqrt(w[u]), n); sqrt(drop(t(cc) %*% solve(t(X[u,]) %*% S %*% X[u,]) %*% cc)) })
cat(jsonlite::toJSON(list(effect=unname(eb$coefficients[,1]), se=unname(sqrt(eb$s2.post)*eb$stdev.unscaled[,1]), p=unname(eb$p.value[,1]),
  stdev=unname(eb$stdev.unscaled[,1]), gls=unscaled, ids=rownames(m)), digits=NA))
"""


@pytest.mark.parametrize("mode", ["none", "duplicate_correlation"])
def test_v042_weighted_sparse_blocked_contrasts_through_export(tmp_path, mode):
    blocking = {"mode": mode} if mode == "none" else {"mode": mode, "subject_column": "subject_id"}
    path, x, w = _weighted_dataset(tmp_path / "data", blocking)
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert any(s["stage_id"] == "limma" and s["state"] == "COMPLETED" for s in payload["stages"]), payload
    plan = json.loads((tmp_path / "run" / "plan.json").read_text())
    rho = plan["blocking"][0]["consensus_correlation"] if mode != "none" else float("nan")
    oracle = B.r_json(ORACLE_GLS, str(tmp_path / "run" / "inputs" / "matrix.tsv"), ",".join(map(repr, x)), ",".join(map(repr, w)), repr(rho))
    table = rows_by(tmp_path / "run" / "dea" / "zero_null.tsv", "feature_id")
    status = json.loads((tmp_path / "run" / "dea" / "model_status.json").read_text())["settings"]["limma-main"]
    assert status["weights"] == "observation" and status["exactness"]["B-A"]["weights"] is True
    assert (status["block_correlation"] is None) == (mode == "none")
    for i, feature in enumerate(oracle["ids"]):
        assert close(oracle["stdev"][i], oracle["gls"][i], 1e-9, 1e-7)                 # direct coefficient refit == independent GLS
        row = table[(feature,)]
        assert close(float(row["effect"]), oracle["effect"][i], 1e-9, 1e-7)
        assert close(float(row["effect_se"]), oracle["se"][i]) and close(float(row["p_value"]), oracle["p"][i])


# ----------------------------------------------------------------------------- V043 / V044 / V047
def _trend_dataset(tmp_path):
    obs = [{"observation_id": f"{g}{i}", "group": g} for g in ("C", "U", "T") for i in range(4)]
    values = {}
    effects = {"F1": 0.2, "F2": 0.5, "F3": 1.0}
    for f in range(18):
        name = f"F{f + 1}"
        eff = effects.get(name, 0.0)
        scale = 3.0 if name == "F18" else 0.15 + 0.02 * (f % 5)            # F18 is hypervariable
        values[name] = [6 + 0.6 * f + (eff if o["group"] == "U" else 0) + scale * (((i * 5 + f * 3) % 7) - 3) / 3 for i, o in enumerate(obs)]
    files = B.dataset(tmp_path, values, obs)
    contrasts = [B.contrast("d", "U", "C"), B.contrast("t", "T", "U"), B.contrast("r", "T", "C", role="secondary"),
                 {"id": "interaction", "design_id": "joint", "label": "t minus d", "estimand": "difference of treatment and disease contrasts",
                  "weights": {"(Intercept)": 0, "group.U": -2, "group.T": 1}, "role": "secondary", "required_groups": ["T", "U", "C"]}]
    def mutate(c):
        c["design"]["intercept"] = True
        for k in contrasts[:3]:
            num, den = k["required_groups"]; k["weights"] = {key: val for key, val in ((f"group.{num}", 1), (f"group.{den}", -1)) if key != "group.C"}
        c["contrasts"] = contrasts
        c["models"][0]["additional_hypotheses"] = ["treat"]; c["models"][0]["effect_threshold"] = 0.5
        c["omnibus_tests"] = [{"id": "omni", "model_id": "limma-main", "coefficient_names": ["group.U", "group.T"]}]
        c["multiplicity_families"] = [
            {"id": "zn-primary", "hypothesis_type": "protein_zero_null", "model_ids": ["limma-main"], "contrast_ids": ["d", "t"], "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "primary"},
            {"id": "zn-secondary", "hypothesis_type": "protein_zero_null", "model_ids": ["limma-main"], "contrast_ids": ["r", "interaction"], "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"},
            {"id": "treat-secondary", "hypothesis_type": "protein_treat", "model_ids": ["limma-main"], "contrast_ids": ["d", "t", "r", "interaction"], "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"},
            {"id": "omnibus", "hypothesis_type": "protein_omnibus", "model_ids": ["limma-main"], "contrast_ids": ["omni"], "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"}]
    return B.config(tmp_path, files, groups=["C", "U", "T"], contrasts=contrasts, mutate=mutate)


TREND_ORACLE = """
m <- as.matrix(read.delim(commandArgs(TRUE)[1], row.names=1, check.names=FALSE)); g <- factor(sub('[0-9]+$','',colnames(m)), levels=c('C','U','T'))
X <- model.matrix(~ g); fit <- limma::lmFit(m, X); eb <- limma::eBayes(fit, trend=TRUE, robust=TRUE)
tr <- limma::treat(fit, lfc=0.5, trend=TRUE, robust=TRUE); f2 <- limma::eBayes(fit[, 2:3], trend=TRUE, robust=TRUE)
cat(jsonlite::toJSON(list(ids=rownames(m), s2prior=unname(eb$s2.prior), dfprior=unname(eb$df.prior), s2post=unname(eb$s2.post), treat_p=unname(tr$p.value[,2]), treat_t=unname(tr$t[,2]),
  zn_p=unname(eb$p.value[,2]), F=unname(f2$F), F_p=unname(f2$F.p.value), int_effect=unname(fit$coefficients[,3] - 2*fit$coefficients[,2])), digits=NA))
"""


def test_v043_v044_v047_moderation_treat_and_omnibus(tmp_path):
    path = _trend_dataset(tmp_path / "data")
    payload, code = workflow.run_command(path, tmp_path / "run")
    run = tmp_path / "run"
    oracle = B.r_json(TREND_ORACLE, str(run / "inputs" / "matrix.tsv"))
    diag = rows_by(run / "dea" / "diagnostics" / "limma-main_moderation.tsv", "feature_id")
    for i, feature in enumerate(oracle["ids"]):
        d = diag[(feature,)]
        assert close(float(d["s2_prior"]), oracle["s2prior"][i]) and close(float(d["df_prior"]), oracle["dfprior"][i]) and close(float(d["s2_posterior"]), oracle["s2post"][i])
    settings = json.loads((run / "dea" / "model_status.json").read_text())["settings"]["limma-main"]
    assert settings["trend_applied"] is True and settings["robust_applied"] is True and "not robust sample regression" in settings["robust_meaning"]
    assert settings["prior_identical_across_contrasts"] is True
    zn = rows_by(run / "dea" / "zero_null.tsv", "contrast_id", "feature_id")
    tr = rows_by(run / "dea" / "treat.tsv", "contrast_id", "feature_id")
    for i, feature in enumerate(oracle["ids"]):
        assert close(float(tr[("d", feature)]["p_value"]), oracle["treat_p"][i]) and close(float(zn[("d", feature)]["p_value"]), oracle["zn_p"][i])
        assert tr[("d", feature)]["result_type"] == "ProteinTreatResult" and tr[("d", feature)]["family_id"] == "treat-secondary"
        assert zn[("d", feature)]["result_type"] == "ProteinZeroNullResult" and zn[("d", feature)]["family_id"] == "zn-primary"
        assert close(float(zn[("interaction", feature)]["effect"]), oracle["int_effect"][i], 1e-10, 1e-8)
    assert len(tr) == len(zn) == 4 * 18                                         # TREAT is not a filtered zero-null table
    for feature in ("F1", "F2", "F3"):
        assert tr[("d", feature)]["p_value"] != zn[("d", feature)]["p_value"]
        flag = tr[("d", feature)]["display_effect_filter"]
        assert flag == ("true" if abs(float(tr[("d", feature)]["effect"])) >= 0.5 else "false")
    om = rows_by(run / "dea" / "omnibus.tsv", "feature_id")
    for i, feature in enumerate(oracle["ids"]):
        row = om[(feature,)]
        assert row["hypothesis_type"] == "protein_omnibus" and row["family_id"] == "omnibus"
        assert close(float(row["statistic"]), oracle["F"][i]) and close(float(row["p_value"]), oracle["F_p"][i])
    header = B.read_tsv(run / "dea" / "zero_null.tsv")[0].keys()
    assert not any("interaction_significant" in h for h in header)                # no significance-difference claims


# ----------------------------------------------------------------------------- V046 / V048
def test_v048_complete_typed_table_with_exclusion_and_numerical_failure(tmp_path):
    obs = [{"observation_id": f"{g}{i}", "group": g} for g in ("C", "U", "T") for i in range(4)]
    values = {f"F{f + 1}": [10 + f + 0.1 * (((i * 3 + f) % 5) - 2) + (0.7 if o["group"] == "U" else 0) for i, o in enumerate(obs)] for f in range(8)}
    values["F7"] = [v if not (o["group"] == "T" and o["observation_id"] != "T0") else None for v, o in zip(values["F7"], obs)]   # 1/4 in T
    values["F8"] = [1e300 if i == 0 else (-1e300 if i == 1 else 5.0 + 0.1 * i) for i in range(12)]                            # eligible, numerically fails
    files = B.dataset(tmp_path / "data", values, obs)
    contrasts = [B.contrast("d", "U", "C"), B.contrast("t", "T", "U"), B.contrast("r", "T", "C", role="secondary")]
    path = B.config(tmp_path / "data", files, groups=["C", "U", "T"], contrasts=contrasts)
    payload, code = workflow.run_command(path, tmp_path / "run")
    rows = B.read_tsv(tmp_path / "run" / "dea" / "zero_null.tsv")
    keys = {(r["contrast_id"], r["feature_id"]) for r in rows}
    assert keys == {(c, f"F{f + 1}") for c in ("d", "t", "r") for f in range(8)}   # cartesian planned keys
    by = {(r["contrast_id"], r["feature_id"]): r for r in rows}
    excl = by[("t", "F7")]
    assert excl["eligibility"] == "excluded" and excl["reason_code"] == "coverage_below_minimum:T" and excl["p_value"] == "NA" and excl["q_value"] == "NA"
    assert by[("d", "F7")]["eligibility"] == "tested"
    fail = by[("d", "F8")]
    assert fail["eligibility"] == "numerical_failure" and fail["p_value"] == "NA" and fail["reason_code"].startswith("E_NUMERICAL")
    for r in rows:
        for field in ("engine", "hypothesis_type", "estimable", "n_obs_by_required_group", "family_id"):
            assert r[field] not in ("",)
    families = rows_by(tmp_path / "run" / "dea" / "families.tsv", "family_id")
    primary = families[("protein-primary",)]
    assert primary["completeness"] == "incomplete" and primary["n_numerical_failure"] == "2" and primary["n_planned"] == "16"
    assert int(primary["n_finite"]) == 16 - 2 - 1
    tested_p = sorted(float(by[(c, f)]["p_value"]) for c in ("d", "t") for f in [f"F{i}" for i in range(1, 9)] if by[(c, f)]["eligibility"] == "tested")
    n = len(tested_p); bh = [min(1.0, min(tested_p[j] * n / (j + 1) for j in range(i, n))) for i in range(n)]
    got = sorted(float(by[(c, f)]["q_value"]) for c in ("d", "t") for f in [f"F{i}" for i in range(1, 9)] if by[(c, f)]["eligibility"] == "tested")
    assert all(close(a, b, 1e-12, 1e-10) for a, b in zip(sorted(bh), got))        # pooled once across both primary contrasts


# ----------------------------------------------------------------------------- V049
def test_v049_whole_subject_influence_after_technical_aggregation(tmp_path):
    obs = []
    values = {f"F{f}": [] for f in range(1, 7)}
    for s in range(1, 4):
        for visit in ("pre", "post"):
            for inj in (1, 2):
                obs.append({"observation_id": f"S{s}_{visit}_i{inj}", "biological_unit_id": f"S{s}_{visit}", "subject_id": f"S{s}", "technical_replicate_id": f"i{inj}", "group": visit})
                for f in range(1, 7):
                    values[f"F{f}"].append(8 + f + 0.3 * s + (0.6 * (f % 2) if visit == "post" else 0) + 0.01 * inj + 0.05 * ((s * f + inj) % 3))
    files = B.dataset(tmp_path / "data", values, obs)
    def mutate(c):
        c["preprocessing"]["technical_replicates"] = {"mode": "aggregate", "method": "mean_log2", "biological_unit_column": "biological_unit_id"}
    path = B.config(tmp_path / "data", files, groups=["pre", "post"], contrasts=[B.contrast("post-pre", "post", "pre")],
                    design_overrides={"blocking": {"mode": "fixed_subject", "subject_column": "subject_id"}}, mutate=mutate)
    payload, code = workflow.run_command(path, tmp_path / "run")
    influence = rows_by(tmp_path / "run" / "dea" / "influence" / "limma-main_subject_influence.tsv", "subject_id")
    assert set(influence) == {("S1",), ("S2",), ("S3",)}
    for row in influence.values():
        assert row["n_omitted_observations"] == "2" and row["n_omitted_in_required_groups"] == "2" and row["state"] == "refit"
        assert "not held-out validation" in row["interpretation"]
    oracle = B.r_json("m <- as.matrix(read.delim(commandArgs(TRUE)[1], row.names=1, check.names=FALSE)); keep <- !grepl('^S1_', colnames(m)); d <- m[, keep]; "
                      "full <- rowMeans(m[, grepl('post', colnames(m))] - m[, grepl('pre', colnames(m))]); part <- rowMeans(d[, grepl('post', colnames(d))] - d[, grepl('pre', colnames(d))]); "
                      "cat(jsonlite::toJSON(list(med=median(abs(part-full)), mx=max(abs(part-full))), auto_unbox=TRUE, digits=NA))", str(tmp_path / "run" / "inputs" / "matrix.tsv"))
    assert close(float(influence[("S1",)]["max_abs_effect_change"]), oracle["mx"], 1e-10, 1e-8)
    assert close(float(influence[("S1",)]["median_abs_effect_change"]), oracle["med"], 1e-10, 1e-8)


# ----------------------------------------------------------------------------- V050
def test_v050_zero_discoveries_completed_and_shared_control_algebra(tmp_path):
    obs = [{"observation_id": f"{g}{i}", "group": g} for g in ("C", "U", "T") for i in range(4)]
    offsets = [-0.3, -0.1, 0.1, 0.3]
    values = {f"F{f + 1}": [10 + f + offsets[(i + f) % 4] for i, o in enumerate(obs)] for f in range(6)}   # equal group means, symmetric
    files = B.dataset(tmp_path / "data", values, obs)
    contrasts = [B.contrast("d", "U", "C"), B.contrast("t", "T", "U"), B.contrast("r", "T", "C", role="secondary")]
    path = B.config(tmp_path / "data", files, groups=["C", "U", "T"], contrasts=contrasts)
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert payload["state"] == "COMPLETED" and code == 0, payload
    families = rows_by(tmp_path / "run" / "dea" / "families.tsv", "family_id")
    assert families[("protein-primary",)]["rejection_count"] == "0" and families[("protein-primary",)]["completeness"] == "complete"
    cov = rows_by(tmp_path / "run" / "designs" / "designs" / "joint" / "contrast_covariance.tsv", "contrast_a", "contrast_b")
    assert float(cov[("d", "t")]["unscaled_covariance"]) == -0.25 and float(cov[("r", "r")]["unscaled_covariance"]) == 0.5
    assert not [p for p in (tmp_path / "run").rglob("*") if "score" in p.name.lower() and p.name != "pca_scores.tsv"]   # PCA scores are QC, not score tests
    for table in (tmp_path / "run" / "dea").glob("*.tsv"):
        assert "independent_score_test" not in table.read_text() and "descriptive_score" not in table.read_text()
    status = json.loads((tmp_path / "run" / "dea" / "model_status.json").read_text())
    assert status["score_p_values_emitted"] is False


# ----------------------------------------------------------------------------- V029 integration rerun after R05
def test_v029_rerun_detection_family_stays_separate_from_abundance_families(tmp_path):
    obs = [{"observation_id": f"{g}{i}", "group": g} for g in ("A", "B") for i in range(4)]
    values = {f"F{f + 1}": [10 + f + 0.1 * ((i + f) % 3) + (0.5 if o["group"] == "B" else 0) for i, o in enumerate(obs)] for f in range(4)}
    values["F1"] = [10.0, 10.2, 10.1, None, 11.0, None, None, None]
    files = B.dataset(tmp_path / "data", values, obs)
    def mutate(c):
        c["preprocessing"]["detection_enabled"] = True
        c["multiplicity_families"].append({"id": "detection-family", "hypothesis_type": "detection", "model_ids": ["limma-main"], "contrast_ids": ["B-A"],
                                           "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"})
    path = B.config(tmp_path / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], mutate=mutate)
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert code == 0, payload
    detection = B.read_tsv(tmp_path / "run" / "preprocessing" / "detection" / "detection.tsv")
    assert {r["family_id"] for r in detection} == {"detection-family"} and {r["hypothesis_type"] for r in detection} == {"detection"}
    f1 = next(r for r in detection if r["feature_id"] == "F1")
    assert close(float(f1["p_value"]), 34 / 70, 1e-12, 1e-10)                      # 3/1 vs 1/3 enumeration
    families = {r["family_id"] for r in B.read_tsv(tmp_path / "run" / "dea" / "families.tsv")}
    assert "detection-family" not in families
    for row in B.read_tsv(tmp_path / "run" / "dea" / "zero_null.tsv"):
        assert row["family_id"] != "detection-family" and row["hypothesis_type"] == "protein_zero_null"
