"""R14c acceptance through the real Python -> R pipeline: V144-V148 (SM33).

Oracles are written here: direct limma fits with ~ group + covariate + phenotype, cor.test() and a seeded
permutation loop over units, complete-case fits, design rank and within-group/pooled slopes, and a cell-by-cell
re-derivation of the heatmap from the association tables.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


F = _load("post_de_factory", ROOT / "tests" / "fixtures" / "post_de" / "design_factory.py")
B = F.B
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")


def run(tmp_path, *, kind="linear", missing=0, n=10, phenotype=None, method="model", correlation="pearson", required="optional", permutations=199):
    values, obs = F.phenotype_design(kind=kind, missing=missing, n=n)
    files = F.write_dataset(tmp_path / "data", values, obs, extra_columns=("age", "score"))
    ph = {"column": "score", "type": "numeric", "adjust_for": [{"column": "age", "type": "continuous"}], "adjust_for_group": True, "scope": "pooled"}
    ph.update(phenotype or {})
    post = {"enabled": True, "association": {"enabled": True, "execution_requirement": required, "phenotypes": [ph], "method": method, "correlation": correlation,
                                             "permutations": permutations, "seed": 777}}
    path = F.config(tmp_path / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], post_de=post)
    payload, code = workflow.run_command(path, tmp_path / "run")
    return payload, code, tmp_path / "run"


ORACLE = """
  a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
  o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]
  keep <- !is.na(o$score) & o$score != 'NA'; if (nzchar(a[3])) keep <- keep & o$group == a[3]
  y <- y[, keep]; o <- o[keep, ]; s <- as.numeric(o$score); age <- as.numeric(o$age)
  X <- if (a[4] == 'group') model.matrix(~ factor(o$group) + age + s) else model.matrix(~ age + s)
  eb <- limma::eBayes(limma::lmFit(y, X), trend = TRUE, robust = TRUE); j <- ncol(X)
  cat(jsonlite::toJSON(list(feature = rownames(y), effect = unname(eb$coefficients[, j]), t = unname(eb$t[, j]), p = unname(eb$p.value[, j]),
                            q = unname(p.adjust(eb$p.value[, j], 'BH')), n = ncol(y)), digits = NA))"""


def _oracle(out, group="", terms="group"):
    return B.r_json(ORACLE, out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv", group, terms)


def _assoc(out):
    return F.read_tsv(out / "post_de" / "association" / "association_score.tsv")


@pytest.fixture(scope="module")
def linear(tmp_path_factory):
    payload, code, out = run(tmp_path_factory.mktemp("linear"))
    assert code == 0, payload
    return out


# ----------------------------------------------------------------------------- V144
def test_v144_phenotype_model_equals_direct_limma_with_one_family(linear):
    oracle = _oracle(linear)
    rows = {r["feature_id"]: r for r in _assoc(linear)}
    for f, e, t, p, q in zip(oracle["feature"], oracle["effect"], oracle["t"], oracle["p"], oracle["q"]):
        r = rows[f]
        assert abs(float(r["effect"]) - e) <= 1e-10 and abs(float(r["statistic"]) - t) <= 1e-8 and abs(float(r["p_value"]) - p) <= 1e-10 and abs(float(r["q_value"]) - q) <= 1e-10
    assert {r["family_id"] for r in rows.values()} == {"phenotype__score"} and {r["adjusted_for"] for r in rows.values()} == {"group;age"}
    assert sum(float(rows[f"F{i:02d}"]["q_value"]) < 0.05 for i in range(1, 11)) == 10        # the planted features are found
    # negative: phenotype rows never enter the group-contrast families
    assert {r["family_id"] for r in F.read_tsv(linear / "dea" / "zero_null.tsv")} == {"protein-primary"}
    assert all(r["hypothesis_type"] == "protein_zero_null" and r["model_id"] == "limma-main" for r in F.read_tsv(linear / "dea" / "zero_null.tsv"))


# ----------------------------------------------------------------------------- V145
@pytest.mark.parametrize("method, adjust", [("spearman", []), ("partial", [{"column": "age", "type": "continuous"}])])
def test_v145_correlation_with_unit_permutation(tmp_path, method, adjust):
    payload, code, out = run(tmp_path, n=4, method="correlation", correlation=method, phenotype={"adjust_for": adjust, "adjust_for_group": False}, permutations=199)
    assert code == 0, payload
    rows = {r["feature_id"]: r for r in _assoc(out)}
    oracle = B.r_json("""
      a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
      o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]
      s <- as.numeric(o$score); age <- as.numeric(o$age); method <- a[3]
      set.seed(777, kind = "L'Ecuyer-CMRG"); perms <- lapply(1:199, function(b) sample.int(ncol(y)))
      res <- lapply(rownames(y), function(f) {
        v <- y[f, ]
        if (method == 'partial') { rv <- residuals(lm(v ~ age)); rs <- residuals(lm(s ~ age)); r <- cor(rv, rs); null <- sapply(perms, function(p) cor(rv, rs[p])) }
        else { r <- unname(cor.test(v, s, method = 'spearman', exact = FALSE)$estimate); null <- sapply(perms, function(p) cor(v, s[p], method = 'spearman')) }
        k <- sum(abs(null) >= abs(r) - 1e-12); list(f = f, r = r, k = k, p = (k + 1) / 200) })
      p <- sapply(res, function(x) x$p)
      cat(jsonlite::toJSON(list(f = sapply(res, function(x) x$f), r = sapply(res, function(x) x$r), k = sapply(res, function(x) x$k), p = p, q = p.adjust(p, 'BH')), digits = NA))""",
                      out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv", method)
    for f, r, k, p, q in zip(oracle["f"], oracle["r"], oracle["k"], oracle["p"], oracle["q"]):
        row = rows[f]
        assert abs(float(row["effect"]) - r) <= 1e-12 and int(row["permutation_k"]) == k
        assert abs(float(row["p_value"]) - (k + 1) / 200) <= 1e-15 and abs(float(row["q_value"]) - q) <= 1e-12
    assert {r["permutation_scheme"] for r in rows.values()} == {"biological_units"} and {r["n_units"] for r in rows.values()} == {"8"}


# ----------------------------------------------------------------------------- V146
def test_v146_missing_phenotype_is_complete_case_never_imputed(tmp_path):
    payload, code, out = run(tmp_path, missing=3)
    assert code == 0, payload
    oracle = _oracle(out)
    rows = {r["feature_id"]: r for r in _assoc(out)}
    assert oracle["n"] == [17] and {r["n_units"] for r in rows.values()} == {"17"} and {r["n_missing_phenotype"] for r in rows.values()} == {"3"}
    for f, e, p in zip(oracle["feature"], oracle["effect"], oracle["p"]):
        assert abs(float(rows[f]["effect"]) - e) <= 1e-10 and abs(float(rows[f]["p_value"]) - p) <= 1e-10
    request = json.loads((out / "post_de" / "association" / "stage-request.json").read_text(encoding="utf-8"))
    values = request["parameters"]["analyses"][0]["values"]
    assert len(values) == 17 and not {"A1", "A2", "A3"} & set(values)                 # no imputed phenotype values exist


def test_v146_negative_phenotype_imputation_is_refused(tmp_path):
    with pytest.raises(ProteomicsError) as error:
        run(tmp_path, missing=3, phenotype={"missing": "mean_impute"})
    assert error.value.code == "E_PHENOTYPE_IMPUTATION"


# ----------------------------------------------------------------------------- V147
def test_v147_one_group_phenotype_pooled_refused_within_group_runs(tmp_path):
    payload, code, out = run(tmp_path, kind="one_group", phenotype={"scope": "both"})
    assert code == 0, payload
    refusals = F.read_tsv(out / "post_de" / "association" / "refusals.tsv")
    assert any(r["analysis"] == "pooled" and r["reason_code"] == "E_PHENOTYPE_ALIASED" for r in refusals)
    rows = _assoc(out)
    assert {r["scope"] for r in rows} == {"within_group"} and {r["group"] for r in rows} == {"B"} and {r["family_id"] for r in rows} == {"phenotype__score__within_B"}
    oracle = _oracle(out, group="B", terms="none")
    by = {r["feature_id"]: r for r in rows}
    assert all(abs(float(by[f]["effect"]) - e) <= 1e-10 for f, e in zip(oracle["feature"], oracle["effect"]))


def test_v147_negative_required_pooled_association_on_aliased_phenotype_fails(tmp_path):
    payload, code, out = run(tmp_path, kind="one_group", required="required")
    assert code == 2 and payload["error"]["code"] == "E_PHENOTYPE_ALIASED" and not (out / "dea").exists()


def test_v147_simpson_pattern_reports_both_with_warning(tmp_path):
    payload, code, out = run(tmp_path, kind="simpson", phenotype={"scope": "both", "adjust_for_group": False, "adjust_for": []})
    assert code == 0, payload
    rows = _assoc(out)
    pooled = {r["feature_id"]: r for r in rows if r["scope"] == "pooled"}
    within = {(r["group"], r["feature_id"]): r for r in rows if r["scope"] == "within_group"}
    signs = B.r_json("""a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
      o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]; s <- as.numeric(o$score)
      out <- lapply(1:5, function(i) list(pooled = unname(coef(lm(y[i, ] ~ s))[2]), A = unname(coef(lm(y[i, o$group == 'A'] ~ s[o$group == 'A']))[2]), B = unname(coef(lm(y[i, o$group == 'B'] ~ s[o$group == 'B']))[2])))
      cat(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA))""", out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv")
    for i, s in enumerate(signs, start=1):
        f = f"F{i:02d}"
        assert s["pooled"] > 0 and s["A"] < 0 and s["B"] < 0                              # the fixture really is a Simpson pattern
        assert float(pooled[f]["effect"]) > 0 and float(within[("A", f)]["effect"]) < 0 and float(within[("B", f)]["effect"]) < 0
        assert pooled[f]["simpson_flag"] == "true"
    eligibility = json.loads((out / "post_de" / "association" / "eligibility.json").read_text(encoding="utf-8"))
    assert "Simpson" in eligibility["simpson_warning"]
    warnings = json.loads((out / "warnings.json").read_text(encoding="utf-8"))
    assert any(w["code"] == "W_PHENOTYPE_SIMPSON" for w in warnings)


# ----------------------------------------------------------------------------- V148
def test_v148_heatmap_source_and_report_agree_with_tables(linear):
    heat = F.read_tsv(linear / "post_de" / "association" / "heatmap_source.tsv")
    stats = {(r["feature_id"], r["analysis_id"]): r["statistic"] for r in _assoc(linear)}
    assert len(heat) == len([v for v in stats.values() if v != "NA"])
    for cell in heat:
        assert abs(float(cell["value"]) - float(stats[(cell["feature_id"], cell["analysis_id"])])) <= 1e-12
    assert sorted(int(c["row_order"]) for c in heat) == list(range(1, len(heat) + 1))
    eligibility = json.loads((linear / "post_de" / "association" / "eligibility.json").read_text(encoding="utf-8"))
    assert "figures/association_heatmap.png" in {f["file"] for f in eligibility["figures"]}
    data = json.loads((linear / "report" / "report_data.json").read_text(encoding="utf-8"))
    section = data["sections"]["post_de_association"]
    table = next(t for t in section["tables"] if t["source"].endswith("association_summary.tsv"))
    assert table["rows"][0][table["columns"].index("family")] == "phenotype__score" and table["rows"][0][table["columns"].index("n units")] == "20"
    assert table["rows"][0][table["columns"].index("adjusted for")] == "group;age"
    for figure in section["figures"]:
        assert (linear / figure["source"]).is_file() and (linear / figure["src"]).is_file()
