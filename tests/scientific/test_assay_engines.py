"""R06 acceptance: V051-V060 (count evidence, DEqMS, proDA, comparison, dependency/failure propagation)."""
from __future__ import annotations

import importlib.util
import json
import os
import random
import shutil
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.count_evidence import check_count_table
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)


def _engines():
    if B.rscript() is None:
        return False
    try:
        return bool(B.r_json("cat(jsonlite::toJSON(requireNamespace('DEqMS', quietly=TRUE) && requireNamespace('proDA', quietly=TRUE), auto_unbox=TRUE))"))
    except RuntimeError:
        return False


pytestmark = pytest.mark.skipif(not _engines(), reason="NOT_RUN: Rscript/DEqMS/proDA unavailable")
OBS = [{"observation_id": f"{g}{i}", "group": g} for g in ("A", "B") for i in range(1, 5)]


def dataset(tmp_path, *, dropout=False, counts=None, n=16, seed=5):
    rng = random.Random(seed)
    values = {}
    for f in range(n):
        shift = 0.8 if f < 4 else 0.0
        sd = 0.15 + 0.25 * (f % 4) / 3
        values[f"F{f + 1:02d}"] = [20 + f * 0.2 + (shift if o["group"] == "B" else 0) + rng.gauss(0, sd) for o in OBS]
    if dropout:
        values["F05"][1] = None; values["F06"][5] = None
        values["F07"][4:8] = [None] * 4                     # absent from group B
    files = B.dataset(tmp_path / "data", values, OBS)
    counts = counts or {f: 1 + (i % 6) for i, f in enumerate(values)}
    B.write_tsv(tmp_path / "data" / "counts.tsv", [["feature_id", "count_type", "count_value", "count_source", "count_aggregation", "pseudocount_policy"]] +
                [[f, "peptide", str(c), "search_engine_export", "protein_level", "none"] for f, c in counts.items()])
    return files, values


def family(fid, model):
    return {"id": fid, "hypothesis_type": "protein_zero_null", "model_ids": [model], "contrast_ids": ["B-A"], "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"}


def make_config(tmp_path, files, *, deqms="optional", proda=None, mutate=None):
    def m(c):
        c["runtime"]["phase"] = 2; c["assay"] = "lfq_dda"
        base = dict(c["models"][0])
        if deqms:
            d = {k: v for k, v in base.items() if k not in ("trend", "robust")}
            d.update({"id": "deqms", "engine": "deqms", "role": "sensitivity", "execution_requirement": deqms, "count_evidence": "counts.tsv", "count_aggregation": "protein_level",
                      "count_zero_policy": {"mode": "reject"}})
            c["models"].append(d); c["multiplicity_families"].append(family("deqms-family", "deqms"))
        if proda:
            pmod = {k: v for k, v in base.items() if k not in ("trend", "robust", "precision_weights")}
            pmod.update({"id": "proda", "engine": "proda", "role": "sensitivity", "execution_requirement": proda,
                         "coverage": {"policy": "native_dropout", "mask_source": "original_observed"}})
            c["models"].append(pmod); c["multiplicity_families"].append(family("proda-family", "proda"))
        if mutate:
            mutate(c)
    return B.config(tmp_path / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], mutate=m)


DEQMS_ORACLE = """
a <- commandArgs(TRUE); m <- as.matrix(read.delim(a[1], row.names=1, check.names=FALSE)); cnt <- read.delim(a[2])
g <- factor(sub('[0-9]+$','',colnames(m)), levels=c('A','B')); X <- cbind(A=as.numeric(g=='A'), B=as.numeric(g=='B'))
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(m, X), cbind(c(-1,1))))
fit$count <- setNames(cnt$count_value, cnt$feature_id)[rownames(fit$coefficients)]
f4 <- DEqMS::spectraCounteBayes(fit); r <- DEqMS::outputResult(f4, coef_col=1); r <- r[rownames(m),]
cat(jsonlite::toJSON(list(ids=rownames(m), p=r$sca.P.Value, t=r$sca.t, limma_p=r$P.Value), digits=NA))
"""


def test_v051_count_evidence_schema(tmp_path):
    files, values = dataset(tmp_path)
    assert len(check_count_table(tmp_path / "data" / "counts.tsv", model_id="deqms", aggregation="protein_level")) == 16
    rows = (tmp_path / "data" / "counts.tsv").read_text(encoding="utf-8").replace("search_engine_export", "observed_sample_count")
    (tmp_path / "proxy.tsv").write_text(rows, encoding="utf-8")
    with pytest.raises(ProteomicsError) as error:
        check_count_table(tmp_path / "proxy.tsv", model_id="deqms", aggregation="protein_level")
    assert error.value.code == "E_DEQMS_COUNT_EVIDENCE"


def test_v052_v053_deqms_matches_official_sequence(tmp_path):
    files, values = dataset(tmp_path)
    payload, code = workflow.run_command(make_config(tmp_path, files), tmp_path / "run")
    assert code == 0, payload
    oracle = B.r_json(DEQMS_ORACLE, str(tmp_path / "run" / "inputs" / "matrix.tsv"), str(tmp_path / "data" / "counts.tsv"))
    rows = {r["feature_id"]: r for r in B.read_tsv(tmp_path / "run" / "assay_engines" / "assay_results.tsv") if r["model_id"] == "deqms"}
    for i, f in enumerate(oracle["ids"]):
        r = rows[f]
        assert abs(float(r["p_value"]) - oracle["p"][i]) <= 1e-8 + 1e-6 * oracle["p"][i] and abs(float(r["statistic"]) - oracle["t"][i]) <= 1e-8
        assert abs(float(r["native_limma_p_value"]) - oracle["limma_p"][i]) <= 1e-8 + 1e-6 * oracle["limma_p"][i]
        assert r["engine"] == "deqms" and r["statistic_type"] == "deqms_sca_t" and r["family_id"] == "deqms-family"
    assert any(r["p_value"] != r["native_limma_p_value"] for r in rows.values())
    # V053: shuffled count rows give identical results; duplicate IDs fail
    lines = (tmp_path / "data" / "counts.tsv").read_text(encoding="utf-8").splitlines()
    (tmp_path / "data" / "counts.tsv").write_text("\n".join([lines[0]] + list(reversed(lines[1:]))) + "\n", encoding="utf-8")
    payload2, _ = workflow.run_command(make_config(tmp_path, files), tmp_path / "run2")
    again = {r["feature_id"]: r["p_value"] for r in B.read_tsv(tmp_path / "run2" / "assay_engines" / "assay_results.tsv") if r["model_id"] == "deqms"}
    assert again == {f: r["p_value"] for f, r in rows.items()}
    (tmp_path / "data" / "counts.tsv").write_text("\n".join(lines + [lines[1]]) + "\n", encoding="utf-8")
    payload3, code3 = workflow.run_command(make_config(tmp_path, files), tmp_path / "run3")
    status = next(s for s in payload3["stages"] if s["stage_id"] == "model.deqms")
    assert status["state"] == "FAILED" and status["reason_code"] == "E_DEQMS_COUNT_EVIDENCE" and payload3["state"] == "PARTIAL"


def test_v054_sparse_exact_route_and_design_guards(tmp_path):
    files, values = dataset(tmp_path, dropout=True)
    values_ok = {k: v for k, v in values.items() if k != "F07"}
    payload, code = workflow.run_command(make_config(tmp_path, files), tmp_path / "run")
    assert code == 0, payload
    oracle = B.r_json("""a <- commandArgs(TRUE); m <- as.matrix(read.delim(a[1], row.names=1, check.names=FALSE)); cnt <- read.delim(a[2])
      m <- m[rowSums(!is.na(m[,5:8])) >= 2 & rowSums(!is.na(m[,1:4])) >= 2,]
      g <- c(0,0,0,0,1,1,1,1); X <- cbind(A=1-g, B=g); Tm <- rbind(c(-1,1), c(1,1)/sqrt(2)); Xs <- X %*% solve(Tm)
      fit <- limma::eBayes(limma::lmFit(m, Xs)); fit$count <- setNames(cnt$count_value, cnt$feature_id)[rownames(m)]
      f4 <- DEqMS::spectraCounteBayes(fit, coef_col=1); r <- DEqMS::outputResult(f4, coef_col=1)[rownames(m),]
      cat(jsonlite::toJSON(list(ids=rownames(m), p=r$sca.P.Value), digits=NA))""", str(tmp_path / "run" / "inputs" / "matrix.tsv"), str(tmp_path / "data" / "counts.tsv"))
    rows = {r["feature_id"]: r for r in B.read_tsv(tmp_path / "run" / "assay_engines" / "assay_results.tsv") if r["model_id"] == "deqms"}
    for i, f in enumerate(oracle["ids"]):
        assert abs(float(rows[f]["p_value"]) - oracle["p"][i]) <= 1e-8 + 1e-6 * oracle["p"][i]
    assert rows["F07"]["eligibility"] != "tested"
    blocked = make_config(tmp_path / "b", dataset(tmp_path / "b")[0], deqms="required",
                          mutate=lambda c: c["design"].update({"blocking": {"mode": "fixed_subject", "subject_column": "biological_unit_id"}}))
    payload, code = workflow.run_command(blocked, tmp_path / "b" / "run")
    assert code == 2 and payload["error"]["code"] == "E_DEQMS_DESIGN_UNSUPPORTED"


PRODA_ORACLE = """
a <- commandArgs(TRUE); m <- as.matrix(read.delim(a[1], row.names=1, check.names=FALSE)); m <- m[rowSums(!is.na(m)) > 0,]
g <- c(0,0,0,0,1,1,1,1); X <- cbind(A=1-g, B=g); set.seed(as.integer(a[2]), kind="L'Ecuyer-CMRG"); fit <- proDA::proDA(m, design=X); td <- proDA::test_diff(fit, contrast=c(-1,1))
cat(jsonlite::toJSON(list(ids=td$name, p=td$pval, diff=td$diff, se=td$se, df=td$df, n_obs=td$n_obs), digits=NA))
"""


def test_v055_v056_v057_proda_native_results(tmp_path):
    files, values = dataset(tmp_path, dropout=True)
    payload, code = workflow.run_command(make_config(tmp_path, files, deqms=None, proda="optional"), tmp_path / "run")
    assert code == 0, payload
    seed = json.loads((tmp_path / "run" / "plan.json").read_text(encoding="utf-8"))["runtime"]["seed"]
    oracle = B.r_json(PRODA_ORACLE, str(tmp_path / "run" / "inputs" / "matrix.tsv"), str(seed))
    rows = {r["feature_id"]: r for r in B.read_tsv(tmp_path / "run" / "assay_engines" / "assay_results.tsv") if r["model_id"] == "proda"}
    for i, f in enumerate(oracle["ids"]):
        r = rows[f]
        assert abs(float(r["p_value"]) - oracle["p"][i]) <= 1e-6 * max(oracle["p"][i], 1e-12) + 1e-9
        assert abs(float(r["effect_se"]) - oracle["se"][i]) <= 1e-8 and abs(float(r["df_inference"]) - oracle["df"][i]) <= 1e-12
        assert r["native_n_obs"] == str(oracle["n_obs"][i])
    assert rows["F07"]["eligibility"] == "tested" and rows["F07"]["dropout_prior_dependent"] == "true" and rows["F07"]["native_n_obs"] == "4"
    limma_df = {r["feature_id"]: r["df_inference"] for r in B.read_tsv(tmp_path / "run" / "dea" / "zero_null.tsv")}
    assert rows["F01"]["df_inference"] != limma_df["F01"]                      # native df, not borrowed from limma
    tmt = make_config(tmp_path / "t", dataset(tmp_path / "t")[0], deqms=None, proda="required", mutate=lambda c: c.update({"assay": "tmt"}))
    payload, code = workflow.run_command(tmt, tmp_path / "t" / "run")
    assert code == 2 and payload["error"]["code"] == "E_PRODA_LFQ_UNIMPUTED_REQUIRED"


def test_v058_comparison_never_changes_primary(tmp_path):
    files, values = dataset(tmp_path, dropout=True)
    with_alt, _ = workflow.run_command(make_config(tmp_path, files, proda="optional"), tmp_path / "run")
    without, _ = workflow.run_command(make_config(tmp_path / "n", dataset(tmp_path / "n", dropout=True)[0], deqms=None), tmp_path / "n" / "run")
    assert (tmp_path / "run" / "dea" / "zero_null.tsv").read_text(encoding="utf-8").split("\n")[0] == (tmp_path / "n" / "run" / "dea" / "zero_null.tsv").read_text(encoding="utf-8").split("\n")[0]
    strip = lambda p: [{k: v for k, v in r.items() if k not in ("run_id", "plan_hash")} for r in B.read_tsv(p)]
    assert strip(tmp_path / "run" / "dea" / "zero_null.tsv") == strip(tmp_path / "n" / "run" / "dea" / "zero_null.tsv")
    comp = B.read_tsv(tmp_path / "run" / "assay_engines" / "engine_comparison.tsv")
    assert {r["role"] for r in comp} == {"sensitivity_comparison_not_primary"}
    proda = [r for r in comp if r["model_id"] == "proda"]
    assert any(r["in_alternative_universe"] == "true" and r["in_primary_universe"] == "false" for r in proda)   # F07: proDA-only universe disclosed
    plan = json.loads((tmp_path / "run" / "plan.json").read_text(encoding="utf-8"))
    assert next(m for m in plan["models"] if m["role"] == "primary")["engine"] == "limma"


def _lib_without(tmp_path, package):
    """Shadow one package without assuming where R libraries live (CI run 36982402784: `.r-lib` absent on runners).
    A library placed first on R_LIBS holds an invalid installation of the package, so requireNamespace() is FALSE
    while every other package still resolves through R's own .libPaths(); no symlinks or copies (portable to Windows)."""
    lib = tmp_path / "shadow-lib"; (lib / package).mkdir(parents=True)
    (lib / package / "DESCRIPTION").write_text(f"Package: {package}\nVersion: 0.0.0\nTitle: shadow of a missing package\n", encoding="utf-8")
    return lib


def _shadow_env(monkeypatch, lib, package):
    from proteomics_pipeline.runtime import run_r_code
    previous = os.environ.get("R_LIBS")
    monkeypatch.setenv("R_LIBS", os.pathsep.join([str(lib)] + ([previous] if previous else [])))
    probe = run_r_code(f"cat(requireNamespace('{package}', quietly = TRUE), requireNamespace('limma', quietly = TRUE), requireNamespace('proteomicsCore', quietly = TRUE))")
    if probe.stdout.split() != ["FALSE", "TRUE", "TRUE"]:
        # Review follow-up 2026-10-03: a setup that cannot simulate the missing package is a FAILED acceptance run, never a
        # SKIP (a skip would let V059 pass CI without executing; R and the engines are present whenever this module runs).
        pytest.fail(f"V059 setup failed: could not shadow {package} portably ({probe.stdout!r} {probe.stderr[-200:]!r})")


def test_v059_shadow_setup_failure_fails_instead_of_skipping(monkeypatch, tmp_path):
    import subprocess
    from proteomics_pipeline import runtime
    monkeypatch.setattr(runtime, "run_r_code", lambda *a, **k: subprocess.CompletedProcess([], 0, stdout="TRUE TRUE TRUE", stderr=""))
    with pytest.raises(BaseException) as outcome:
        _shadow_env(monkeypatch, _lib_without(tmp_path, "DEqMS"), "DEqMS")
    assert outcome.type is pytest.fail.Exception, outcome.type          # pytest.skip.Exception before the fix


@pytest.mark.parametrize("requirement,state,code", [("optional", "PARTIAL", 3), ("required", "FAILED", 3)])
def test_v059_missing_engine_package(tmp_path, monkeypatch, requirement, state, code):
    files, _ = dataset(tmp_path)
    _shadow_env(monkeypatch, _lib_without(tmp_path, "DEqMS"), "DEqMS")
    payload, exit_code = workflow.run_command(make_config(tmp_path, files, deqms=requirement), tmp_path / "run")
    assert payload["state"] == state and exit_code == code
    plan = json.loads((tmp_path / "run" / "plan.json").read_text(encoding="utf-8"))
    assert any(m["model_id"] == "deqms" for m in plan["models"])                 # never dropped from the frozen plan
    stage = next(s for s in payload["stages"] if s["stage_id"] in ("model.deqms", "assay_engines") and s["state"] != "COMPLETED")
    assert stage["state"] == "NOT_RUN" and stage["reason_code"] == "E_ENGINE_NOT_AVAILABLE"
    assert not (tmp_path / "run" / "assay_engines" / "assay_results.tsv").exists() or "deqms" not in (tmp_path / "run" / "assay_engines" / "assay_results.tsv").read_text(encoding="utf-8")


@pytest.mark.parametrize("requirement,state", [("optional", "PARTIAL"), ("required", "FAILED")])
def test_v059_v060_real_child_failure(tmp_path, requirement, state):
    files, values = dataset(tmp_path, counts={f"F{i:02d}": 3 for i in range(1, 17)})          # constant counts: DEqMS loess genuinely fails
    payload, exit_code = workflow.run_command(make_config(tmp_path, files, deqms=requirement), tmp_path / "run")
    assert payload["state"] == state and exit_code == 4
    bad = next(s for s in payload["stages"] if s["state"] == "FAILED")
    assert bad["reason_code"] in ("E_ENGINE_FAILED", "E_STAGE_EXCEPTION")
    if requirement == "optional":
        assert (tmp_path / "run" / "dea" / "zero_null.tsv").is_file()           # primary results preserved


def test_v051_negative_missing_count_evidence_file_is_typed():
    # Found by the 2026-10-02 full regression: a missing file previously escaped as an untyped traceback (exit 1).
    from proteomics_pipeline.count_evidence import check_count_table
    with pytest.raises(ProteomicsError) as error:
        check_count_table(ROOT / "does-not-exist.tsv", model_id="deqms", aggregation="protein_level")
    assert error.value.code == "E_DEQMS_COUNT_EVIDENCE" and error.value.exit_code == 2
