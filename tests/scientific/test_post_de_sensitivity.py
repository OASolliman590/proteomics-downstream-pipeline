"""R14b acceptance through the real Python -> R pipeline: V137-V143 (SM32).

Oracles are written here: cross-tabulations and design rank by qr() in R, direct limma fits
(lmFit/eBayes with the frozen trend/robust settings and BH within the family), a direct seeded
resampling loop, leave-one-unit-out refits and Python set arithmetic.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_pipeline import workflow

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


F = _load("post_de_factory", ROOT / "tests" / "fixtures" / "post_de" / "design_factory.py")
B = F.B
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")

THREE = [B.contrast("B-A", "B", "A"), B.contrast("C-A", "C", "A", role="secondary")]
THREE_FAMILIES = [F.family("protein-primary", ["B-A"], role="primary"), F.family("protein-secondary", ["C-A"])]


def sens_block(**overrides):
    block = {"enabled": True, "covariate_models": [{"id": "adj_sex", "add_covariates": [{"column": "sex", "type": "categorical"}], "interaction_with_group": True}],
             "subgroups": [{"id": "males", "filter": {"column": "sex", "in": ["M"]}}], "influence": True}
    block.update(overrides)
    return {"enabled": True, "sensitivity": block}


def three_group_run(tmp_path, variant="partial", required="optional", **overrides):
    values, obs = F.sex_three_groups(variant=variant)
    files = F.write_dataset(tmp_path / "data", values, obs, extra_columns=("sex",))
    post = sens_block(**overrides); post["sensitivity"]["execution_requirement"] = required
    path = F.config(tmp_path / "data", files, groups=["A", "B", "C"], contrasts=THREE, families=THREE_FAMILIES, post_de=post)
    payload, code = workflow.run_command(path, tmp_path / "run")
    return payload, code, tmp_path / "run"


@pytest.fixture(scope="module")
def partial(tmp_path_factory):
    payload, code, out = three_group_run(tmp_path_factory.mktemp("partial"))
    assert code == 0, payload
    return out, payload


LIMMA_ORACLE = """
  a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
  o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]
  keep <- if (nzchar(a[3])) o[[a[3]]] %in% strsplit(a[4], ',')[[1]] else rep(TRUE, nrow(o))
  y <- y[, keep, drop = FALSE]; o <- o[keep, ]
  g <- factor(o$group, levels = c('A', 'B', 'C')); g <- droplevels(g)
  X <- if (a[5] == 'sex') model.matrix(~ g + factor(o$sex)) else model.matrix(~ g)
  fit <- limma::lmFit(y, X)
  res <- list(rank = qr(X)$rank)
  for (lev in setdiff(levels(g), 'A')) {
    eb <- limma::eBayes(limma::contrasts.fit(fit, as.numeric(colnames(X) == paste0('g', lev))), trend = TRUE, robust = TRUE)
    res[[lev]] <- list(effect = unname(eb$coefficients[, 1]), t = unname(eb$t[, 1]), p = unname(eb$p.value[, 1]), feature = rownames(y))
  }
  cat(jsonlite::toJSON(res, auto_unbox = TRUE, digits = NA))"""


def _oracle(out, subset_column="", subset_levels="", covariate=""):
    return B.r_json(LIMMA_ORACLE, out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv",
                    subset_column, subset_levels, covariate)


def _rows(path, contrast=None):
    rows = F.read_tsv(path)
    return [r for r in rows if contrast is None or r["contrast_id"] == contrast]


def _bh(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i], reverse=True); n = len(ps); q = [0.0] * n; running = 1.0
    for rank, i in enumerate(order):
        running = min(running, ps[i] * n / (n - rank)); q[i] = running
    return q


# ----------------------------------------------------------------------------- V137
def test_v137_imbalance_diagnostics_and_estimability_classes(partial, tmp_path):
    out, _ = partial
    imbalance = _rows(out / "post_de" / "sensitivity" / "imbalance.tsv")
    values, obs = F.sex_three_groups(variant="partial")
    for row in imbalance:
        assert int(row["n_units"]) == sum(1 for o in obs if o["group"] == row["group"] and o["sex"] == row["level"])   # cross-tabulation oracle
    classes = {r["model_id"]: r for r in _rows(out / "post_de" / "sensitivity" / "estimability_classes.tsv")}
    assert classes["adj_sex"]["class"] == "estimable_additive" and classes["adj_sex__interaction"]["class"] == "estimable_interaction"
    # confounded: C entirely male -> additive estimable when rank allows, interaction non-estimable with the empty cell named
    payload, code, run = three_group_run(tmp_path / "conf", variant="confounded")
    assert code == 0, payload
    classes = {r["model_id"]: r for r in _rows(run / "post_de" / "sensitivity" / "estimability_classes.tsv")}
    rank = B.r_json("""a <- commandArgs(TRUE); o <- read.delim(a[1], colClasses = 'character', encoding = 'UTF-8')
      X <- model.matrix(~ factor(o$group) + factor(o$sex)); Xi <- model.matrix(~ factor(o$group) * factor(o$sex))
      cat(jsonlite::toJSON(list(additive = qr(X)$rank == ncol(X), interaction = qr(Xi)$rank == ncol(Xi)), auto_unbox = TRUE))""",
                    run / "preprocessing" / "primary" / "observations.tsv")
    assert rank == {"additive": True, "interaction": False}
    assert classes["adj_sex"]["class"] == "estimable_additive"
    assert classes["adj_sex__interaction"]["class"] == "non_estimable" and "C/F" in classes["adj_sex__interaction"]["detail"]
    assert not (run / "post_de" / "sensitivity" / "models" / "adj_sex__interaction").exists()


def test_v137_negative_required_interaction_in_confounded_design_fails_before_any_fit(tmp_path):
    payload, code, run = three_group_run(tmp_path, variant="confounded", required="required")
    assert code == 2 and payload["error"]["code"] == "E_SENSITIVITY_NONESTIMABLE" and not (run / "dea").exists()


# ----------------------------------------------------------------------------- V138
def test_v138_covariate_adjusted_model_equals_direct_limma_and_keeps_its_own_family(partial):
    out, _ = partial
    oracle = _oracle(out, covariate="sex")
    rows = {(r["contrast_id"], r["feature_id"]): r for r in _rows(out / "post_de" / "sensitivity" / "models" / "adj_sex" / "zero_null.tsv")}
    for lev, cid in (("B", "B-A"), ("C", "C-A")):
        o = oracle[lev]
        for f, effect, t, p in zip(o["feature"], o["effect"], o["t"], o["p"]):
            r = rows[(cid, f)]
            assert abs(float(r["effect"]) - effect) <= 1e-10 and abs(float(r["statistic"]) - t) <= 1e-8 and abs(float(r["p_value"]) - p) <= 1e-10
        q = _bh(o["p"])
        for f, qq in zip(o["feature"], q):
            assert abs(float(rows[(cid, f)]["q_value"]) - qq) <= 1e-10
    families = {r["family_id"] for r in rows.values()}
    assert families == {"protein-primary__adj_sex", "protein-secondary__adj_sex"} and {r["role"] for r in rows.values()} == {"sensitivity"}
    # the primary family is byte-identical to what the limma stage recorded before any post-DE stage ran
    result = json.loads((out / "dea" / "stage-result.json").read_text(encoding="utf-8"))
    recorded = next(o["sha256"] for o in result["outputs"] if o["relative_path"] == "zero_null.tsv")
    assert hashlib.sha256((out / "dea" / "zero_null.tsv").read_bytes()).hexdigest() == recorded
    assert {r["model_id"] for r in _rows(out / "dea" / "zero_null.tsv")} == {"limma-main"}         # no sensitivity row in the primary tables


# ----------------------------------------------------------------------------- V139
def test_v139_subgroup_matches_direct_fit_on_the_subset(partial):
    out, _ = partial
    oracle = _oracle(out, "sex", "M")
    rows = {(r["contrast_id"], r["feature_id"]): r for r in _rows(out / "post_de" / "sensitivity" / "models" / "males" / "zero_null.tsv")}
    for f, effect, p in zip(oracle["B"]["feature"], oracle["B"]["effect"], oracle["B"]["p"]):
        assert abs(float(rows[("B-A", f)]["effect"]) - effect) <= 1e-10 and abs(float(rows[("B-A", f)]["p_value"]) - p) <= 1e-10
    model = next(r for r in _rows(out / "post_de" / "sensitivity" / "sensitivity_models.tsv") if r["model_id"] == "males")
    assert int(model["n_observations"]) == 12                                # 2 + 6 + 4 males


def test_v139_negative_one_male_group_is_refused_without_partial_fit(tmp_path):
    payload, code, run = three_group_run(tmp_path, variant="one_male")
    assert code == 0, payload
    refusals = _rows(run / "post_de" / "sensitivity" / "refusals.tsv")
    refused = next(r for r in refusals if r["item"] == "males")
    assert refused["reason_code"] == "E_SUBGROUP_COVERAGE" and "A (1)" in refused["reason"]
    assert not (run / "post_de" / "sensitivity" / "models" / "males").exists()
    assert "males" not in {r["model_id"] for r in _rows(run / "post_de" / "sensitivity" / "comparison.tsv")}


# ----------------------------------------------------------------------------- V140
def test_v140_comparison_equals_set_and_regression_oracles(partial):
    out, _ = partial
    primary = {r["feature_id"]: r for r in _rows(out / "dea" / "zero_null.tsv", "B-A") if r["eligibility"] == "tested"}
    sens = {r["feature_id"]: r for r in _rows(out / "post_de" / "sensitivity" / "models" / "adj_sex" / "zero_null.tsv", "B-A") if r["eligibility"] == "tested"}
    a = {f for f, r in primary.items() if float(r["q_value"]) < 0.05}; b = {f for f, r in sens.items() if float(r["q_value"]) < 0.05}
    row = next(r for r in _rows(out / "post_de" / "sensitivity" / "comparison.tsv") if r["model_id"] == "adj_sex" and r["contrast_id"] == "B-A")
    split = lambda key: set(filter(None, row[key].split(";")))
    assert split("retained") == a & b and split("lost") == a - b and split("gained") == b - a
    assert abs(float(row["jaccard"]) - len(a & b) / len(a | b)) <= 1e-12
    common = sorted(set(primary) & set(sens))
    x = [float(primary[f]["effect"]) for f in common]; y = [float(sens[f]["effect"]) for f in common]
    stats = B.r_json("a <- as.numeric(strsplit(commandArgs(TRUE)[1], ',')[[1]]); b <- as.numeric(strsplit(commandArgs(TRUE)[2], ',')[[1]]); "
                     "cat(jsonlite::toJSON(list(r = cor(a, b), s = unname(coef(lm(b ~ 0 + a))[1])), auto_unbox = TRUE, digits = NA))",
                     ",".join(repr(v) for v in x), ",".join(repr(v) for v in y))
    assert abs(float(row["effect_correlation"]) - stats["r"]) <= 1e-10 and abs(float(row["attenuation_slope"]) - stats["s"]) <= 1e-10
    assert "not evidence of an interaction" in row["warning"] and row["interaction_family"] == "interaction__adj_sex__interaction"
    inter = _rows(out / "post_de" / "sensitivity" / "models" / "adj_sex__interaction" / "zero_null.tsv")
    assert {r["family_id"] for r in inter} == {"interaction__adj_sex__interaction"} and {r["result_type"] for r in inter} == {"InteractionContrastResult"}
    # negative: "significant in subgroup only" is never labelled an interaction effect
    for r in _rows(out / "post_de" / "sensitivity" / "comparison.tsv"):
        assert r["claim_label"] == "descriptive" and "interaction effect" not in r["gained"]


# ----------------------------------------------------------------------------- V141 / V142
@pytest.fixture(scope="module")
def unbalanced(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("unbalanced")
    values, obs = F.unbalanced_with_injections()
    files = F.write_dataset(tmp / "data", values, obs)
    post = {"enabled": True, "sensitivity": {"enabled": True, "matched_n": {"targets": {"A": 11}, "draws": 20, "seed": 4242}, "influence": True}}
    def mutate(c):
        c["preprocessing"]["technical_replicates"] = {"mode": "aggregate", "method": "mean_log2", "biological_unit_column": "biological_unit_id"}
    path = F.config(tmp / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], post_de=post, mutate=mutate)
    payload, code = workflow.run_command(path, tmp / "run")
    assert code == 0, payload
    return tmp / "run"


def test_v141_matched_n_resampling_equals_a_direct_seeded_loop(unbalanced):
    out = unbalanced
    draws = _rows(out / "post_de" / "sensitivity" / "matched_n.tsv")
    assert len(draws) == 20 and {r["claim_label"] for r in draws} == {"descriptive"}
    oracle = B.r_json("""
      a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
      o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]
      units <- o$biological_unit_id; groups <- o$group
      set.seed(4242, kind = "L'Ecuyer-CMRG"); counts <- integer(); hits <- list(); sel <- list()
      for (d in 1:20) {
        chosen <- c(sample(unique(units[groups == 'A']), 11), unique(units[groups != 'A']))
        keep <- units %in% chosen; g <- factor(groups[keep], levels = c('A', 'B'))
        eb <- limma::eBayes(limma::lmFit(y[, keep], model.matrix(~ g)), trend = TRUE, robust = TRUE)
        q <- p.adjust(eb$p.value[, 2], 'BH'); counts <- c(counts, sum(q < 0.05)); hits[[d]] <- rownames(y)[q < 0.05]; sel[[d]] <- sort(chosen)
      }
      cat(jsonlite::toJSON(list(counts = counts, freq = as.list(table(factor(unlist(hits), levels = rownames(y))) / 20), units = sel), auto_unbox = TRUE, digits = NA))""",
                      out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv")
    assert [int(r["n_discoveries"]) for r in sorted(draws, key=lambda r: int(r["draw"]))] == oracle["counts"]
    for r in sorted(draws, key=lambda r: int(r["draw"])):
        units = r["units"].split(";")
        assert units == oracle["units"][int(r["draw"]) - 1] and len(units) == 22 and "B1" in units   # biological units, never injections
    freq = {r["feature_id"]: float(r["selection_frequency"]) for r in _rows(out / "post_de" / "sensitivity" / "matched_n_selection.tsv")}
    assert all(abs(freq[f] - v) <= 1e-12 for f, v in oracle["freq"].items())
    settings = _rows(out / "post_de" / "sensitivity" / "matched_n_settings.tsv")[0]
    assert settings["seed"] == "4242" and settings["draws"] == "20" and settings["unit_level"] == "biological_unit"


def test_v142_unit_influence_equals_leave_one_unit_out_refits(unbalanced):
    out = unbalanced
    rows = _rows(out / "post_de" / "sensitivity" / "influence.tsv", "B-A")
    oracle = B.r_json("""
      a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
      o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]
      fit1 <- function(keep) { g <- factor(o$group[keep], levels = c('A', 'B')); eb <- limma::eBayes(limma::lmFit(y[, keep], model.matrix(~ g)), trend = TRUE, robust = TRUE)
        list(e = eb$coefficients[, 2], q = p.adjust(eb$p.value[, 2], 'BH')) }
      full <- fit1(rep(TRUE, ncol(y))); res <- list()
      for (u in unique(o$biological_unit_id)) { f <- fit1(o$biological_unit_id != u)
        res[[u]] <- list(max = max(abs(f$e - full$e)), lost = sort(setdiff(names(full$q)[full$q < 0.05], names(f$q)[f$q < 0.05])), gained = sort(setdiff(names(f$q)[f$q < 0.05], names(full$q)[full$q < 0.05]))) }
      cat(jsonlite::toJSON(res, auto_unbox = TRUE, digits = NA))""",
                      out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv")
    by_unit = {r["unit_id"]: r for r in rows}
    assert set(by_unit) == set(oracle)                                         # one row per biological unit, not per injection
    for unit, o in oracle.items():
        r = by_unit[unit]
        assert abs(float(r["max_abs_effect_change"]) - o["max"]) <= 1e-10
        assert sorted(filter(None, r["lost"].split(";"))) == sorted(o["lost"] if isinstance(o["lost"], list) else [o["lost"]])
    top = min(rows, key=lambda r: int(r["rank"]))
    assert top["unit_id"] == "B1"                                               # the planted high-leverage unit ranks first
    assert set(top["omitted_source_observations"].split(";")) == {"B1_i1", "B1_i2"} and top["n_omitted_observations"] == "1"


# ----------------------------------------------------------------------------- V143
def test_v143_robustness_summary_equals_test_fractions(partial):
    out, _ = partial
    summary = {(r["contrast_id"], r["feature_id"]): r for r in _rows(out / "post_de" / "sensitivity" / "robustness_summary.tsv")}
    assert not any("q_value" in k or "adjusted" in k for k in next(iter(summary.values())))   # no adjusted-P column
    primary = {r["feature_id"]: r for r in _rows(out / "dea" / "zero_null.tsv", "B-A") if r["eligibility"] == "tested"}
    member = lambda rows, f: f in rows and rows[f]["eligibility"] == "tested" and float(rows[f]["q_value"]) < 0.05
    analyses = []
    for model in ("adj_sex", "males"):
        analyses.append({r["feature_id"]: r for r in _rows(out / "post_de" / "sensitivity" / "models" / model / "zero_null.tsv", "B-A")})
    influence = [r for r in _rows(out / "post_de" / "sensitivity" / "influence.tsv", "B-A") if r["state"] == "refit"]
    for f, prow in primary.items():
        row = summary[("B-A", f)]
        is_member = float(prow["q_value"]) < 0.05
        meets_models = [member(a, f) for a in analyses]
        dropped = [f in set(filter(None, r["lost"].split(";"))) for r in influence]
        retained_infl = [(is_member and not d) or (not is_member and f in set(filter(None, r["gained"].split(";")))) for r, d in zip(influence, dropped)]
        meets = meets_models + retained_infl
        assert abs(float(row["robustness_fraction"]) - sum(meets) / len(meets)) <= 1e-12, f
        assert abs(float(row["sensitivity_model_fraction"]) - sum(meets_models) / 2) <= 1e-12
        assert row["influence_flag"] == ("true" if is_member and any(dropped) else "false")
        signs = [float(a[f]["effect"]) > 0 for a in analyses if f in a and a[f]["eligibility"] == "tested"]
        assert row["claim_label"] == "descriptive" and row["primary_member"] == ("true" if is_member else "false")
    # primary q values are unchanged by the robustness summary
    result = json.loads((out / "dea" / "stage-result.json").read_text(encoding="utf-8"))
    assert hashlib.sha256((out / "dea" / "zero_null.tsv").read_bytes()).hexdigest() == next(o["sha256"] for o in result["outputs"] if o["relative_path"] == "zero_null.tsv")


def test_r14b_report_section(partial):
    out, _ = partial
    data = json.loads((out / "report" / "report_data.json").read_text(encoding="utf-8"))
    section = data["sections"]["post_de_sensitivity"]
    assert section["state"] == "COMPLETED" and any("not evidence of an interaction" in s for s in section["summary"] + section["notes"])
    for figure in section["figures"]:
        assert (out / figure["source"]).is_file() and (out / figure["src"]).is_file()
