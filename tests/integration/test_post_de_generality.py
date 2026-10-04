"""R14f acceptance: V166 (eligibility and dependency report) and V167 (generality qualification matrix), SM41.

V167 runs every post-DE module through the real pipeline on every synthetic design cell of FR-167 (one shared
design factory). The oracle is the expected state table below (COMPLETED or the specific typed refusal, written
before the runs) plus, for every completed module, its own V-case oracle on that cell in reduced form.
"""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import pytest

from proteomics_pipeline import resources, workflow

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


F = _load("post_de_factory", ROOT / "tests" / "fixtures" / "post_de" / "design_factory.py")
B = F.B


def _ready():
    if B.rscript() is None:
        return False
    try:
        return bool(B.r_json("cat(jsonlite::toJSON(all(sapply(c('glmnet','pROC','e1071','dynamicTreeCut'), requireNamespace, quietly=TRUE)), auto_unbox=TRUE))"))
    except RuntimeError:
        return False


pytestmark = pytest.mark.skipif(not _ready(), reason="NOT_RUN: Rscript/proteomicsCore/post-DE R packages unavailable")
CLAIMS = {"descriptive", "exploratory_raw_p", "in_sample", "cross_validated_nested", "fixed_panel_cv", "independently_validated", "module_level"}
MODULES = ("post_de_sets", "post_de_sensitivity", "post_de_association", "post_de_biomarker", "post_de_networks")
C = "COMPLETED"
EXPECTED = {   # module -> expected state or (state, reason code); written before any run
    "two_group": dict.fromkeys(MODULES, C),
    "three_group": dict.fromkeys(MODULES, C),
    "paired": {**dict.fromkeys(MODULES, C), "post_de_networks": ("INAPPLICABLE", "E_COABUNDANCE_DESIGN_UNSUPPORTED")},
    "repeated": {**dict.fromkeys(MODULES, C), "post_de_networks": ("INAPPLICABLE", "E_COABUNDANCE_DESIGN_UNSUPPORTED")},
    "continuous_exposure": dict.fromkeys(MODULES, C),
    "small_3v3": {**dict.fromkeys(MODULES, C), "post_de_biomarker": ("INAPPLICABLE", "E_BIOMARKER_SMALL_N"), "post_de_networks": ("INAPPLICABLE", "E_COABUNDANCE_SMALL_N")},
    "unbalanced_15v80": dict.fromkeys(MODULES, C),
    "large_200": dict.fromkeys(MODULES, C),
    "human": dict.fromkeys(MODULES, C), "mouse": dict.fromkeys(MODULES, C), "rat": dict.fromkeys(MODULES, C),
}
ASSOCIATION_ORACLE = """a <- commandArgs(TRUE); Y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(Y), o$observation_id), ]
s <- as.numeric(o$score); s <- s - mean(s); g <- factor(o$group, levels = sort(unique(o$group)))
X <- if (nzchar(a[3])) { d <- as.numeric(o[[a[3]]]); model.matrix(~ g + I(d - mean(d)) + s) } else model.matrix(~ g + s)
rho <- suppressWarnings(as.numeric(a[4]))
fit <- if (is.finite(rho)) limma::lmFit(Y, X, block = o$subject_id, correlation = rho) else limma::lmFit(Y, X)
eb <- limma::eBayes(fit, trend = TRUE, robust = TRUE); j <- ncol(X)
cat(jsonlite::toJSON(list(f = rownames(Y), t = unname(eb$t[, j]), p = unname(eb$p.value[, j])), digits = NA))"""
SUBANALYSES = {"paired": {("post_de_biomarker", "single_feature", "E_BIOMARKER_PAIRED_AUC")}}


def run_cell(tmp_path, name):
    cell = F.generality_cell(name)
    files = F.write_dataset(tmp_path / "data", cell["values"], cell["observations"], species=cell["species"], extra_columns=cell["extra_columns"])
    species = name in ("human", "mouse", "rat")
    resources_entries = []
    if species:
        snap = F.write_ppi_snapshot(tmp_path / "sources", list(cell["values"]), cluster=[f"F{i:02d}" for i in range(1, 6)], species=name)
        style = F.SPECIES[name]
        for kind in ("mapping", "ppi"):
            prep = resources.prepare(snap[f"{kind}_manifest"], tmp_path / "snapshots" / kind)
            entry = {"id": "map" if kind == "mapping" else "ppi", "kind": kind, "path": str(Path(prep["snapshot"]) / "manifest.json"), "sha256": prep["manifest_sha256"],
                     "version": "synthetic-ppi-1", "source": "synthetic test fixture (not a real database)", "source_taxonomy_id": style["taxonomy_id"],
                     "target_taxonomy_id": style["taxonomy_id"], "terms": "synthetic; no redistribution restrictions"}
            if kind == "ppi":
                entry["score_type"] = "combined_score"
            resources_entries.append(entry)
    def mutate(c):
        c["runtime"]["threads"] = 2
        if resources_entries:
            c["resources"] = resources_entries
    path = F.config(tmp_path / "data", files, groups=cell["groups"], contrasts=cell["contrasts"], post_de=F.generality_post_de(cell, ppi=species),
                    design_overrides=cell["design_overrides"], mutate=mutate, species=cell["species"])
    payload, code = workflow.run_command(path, tmp_path / "run")
    return payload, code, tmp_path / "run", cell


def _auc(scores, labels):
    pos = [s for s, l in zip(scores, labels) if l == 1]; neg = [s for s, l in zip(scores, labels) if l == 0]
    return sum((a > b) + 0.5 * (a == b) for a in pos for b in neg) / (len(pos) * len(neg))


def _check_completed(out, module, cell):
    root = out / "post_de" / module.removeprefix("post_de_")
    for path in root.rglob("*.tsv"):                                               # no ineligible claim anywhere
        rows = F.read_tsv(path)
        if rows and "claim_label" in rows[0]:
            assert {r["claim_label"] for r in rows} <= CLAIMS, path
    primary = [r for r in F.read_tsv(out / "dea" / "zero_null.tsv") if r["contrast_id"] == cell["contrasts"][0]["id"] and r["model_id"] == "limma-main"]
    members = {r["feature_id"] for r in primary if r["eligibility"] == "tested" and float(r["q_value"]) < 0.05}
    if module == "post_de_sets":
        got = {r["feature_id"] for r in F.read_tsv(root / "membership.tsv") if r["primary_any"] == "true"}
        assert got == members
    elif module == "post_de_sensitivity":
        summary = [r for r in F.read_tsv(root / "robustness_summary.tsv") if r["contrast_id"] == cell["contrasts"][0]["id"]]
        assert {r["feature_id"] for r in summary if r["primary_member"] == "true"} == members
        assert all(0 <= float(r["robustness_fraction"]) <= 1 for r in summary if r["robustness_fraction"] != "NA")
    elif module == "post_de_association":
        assert json.loads((root / "eligibility.json").read_text(encoding="utf-8"))["method"] == "model"
        rows = F.read_tsv(root / "association_score.tsv")
        assert {r["family_id"] for r in rows} == {"phenotype__score"} and int(rows[0]["n_units"]) > 0
        plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
        rho = next((b["consensus_correlation"] for b in plan.get("blocking", []) if b.get("consensus_correlation") is not None), None)
        oracle = B.r_json(ASSOCIATION_ORACLE, out / "preprocessing" / "primary" / "matrix.tsv", out / "preprocessing" / "primary" / "observations.tsv",
                          "dose" if cell["name"] == "continuous_exposure" else "", "NA" if rho is None else repr(rho))
        by = {r["feature_id"]: r for r in rows}
        for f, t, pv in zip(oracle["f"], oracle["t"], oracle["p"]):                 # the R14c oracle in reduced form: direct limma on group + centred phenotype
            assert by[f]["eligibility"] == "tested" and abs(float(by[f]["statistic"]) - t) <= 1e-6 and abs(float(by[f]["p_value"]) - pv) <= 1e-8, (cell["name"], f)
        if not cell["design_overrides"] or not cell["design_overrides"].get("blocking"):   # power on the independent-unit cells with >= 24 units
            if len(cell["observations"]) >= 24:
                assert sum(float(by[f]["q_value"]) < 0.05 for f in cell["planted_phenotype"]) >= 3
    elif module == "post_de_biomarker":
        oof = F.read_tsv(root / "oof_predictions.tsv")
        cv = F.read_tsv(root / "cv_performance.tsv")[0]
        assert (cv["classifier"], cv["selection"], cv["permutation_B"]) == ("penalized_logistic", "top_k_auc", "19")   # declared method, no silent switch
        assert cv["grouped_by"] == ("subject_id" if cell["design_overrides"] and cell["design_overrides"].get("blocking") else "biological_unit")
        assert abs(_auc([float(r["score"]) for r in oof], [int(r["label"]) for r in oof]) - float(cv["pooled_oof_auc"])) <= 1e-10
        null = F.read_tsv(root / "permutation_null.tsv")
        k = sum(float(r["pooled_oof_auc"]) >= float(cv["pooled_oof_auc"]) - 1e-12 for r in null)
        assert abs(float(cv["permutation_p"]) - (k + 1) / (len(null) + 1)) <= 1e-15
    elif module == "post_de_networks":
        rows = F.read_tsv(root / "modules.tsv")
        assert {r["rule"] for r in rows} == {"wgcna_signed"}                        # the declared rule, never a silent switch
        modules = {r["feature_id"]: r["module"] for r in rows}
        planted = [modules[f] for f in cell["planted_module"]]
        assert len(set(planted)) == 1 and planted[0] != "unassigned"                # the planted latent module
        trait = {r["module"]: r for r in F.read_tsv(root / "module_trait.tsv")}
        assert {r["family_id"] for r in trait.values()} == {"module_trait__group"} and planted[0] in trait
        if cell["name"] in ("human", "mouse", "rat"):
            row = F.read_tsv(root / "connectivity.tsv")[0]
            assert abs(float(row["p_value"]) - (int(row["k"]) + 1) / (int(row["null_draws"]) + 1)) <= 1e-15
            summary = json.loads((root / "eligibility.json").read_text(encoding="utf-8"))["summary"]["ppi"]
            assert summary["species"] == F.SPECIES[cell["name"]]["taxonomy_id"] and summary["universe_size"] == 30   # every measured feature maps in its own convention


@pytest.mark.parametrize("name", F.GENERALITY_CELLS)
def test_v167_every_module_completes_or_refuses_as_declared(tmp_path, name):
    payload, code, out, cell = run_cell(tmp_path, name)
    states = {s["stage_id"]: s for s in payload["stages"]}
    for module, expected in EXPECTED[name].items():
        stage = states[module]
        if expected == C:
            assert stage["state"] == C, (name, module, stage)
            _check_completed(out, module, cell)
        else:
            assert (stage["state"], stage["reason_code"]) == expected, (name, module, stage)
            assert not (out / "post_de" / module.removeprefix("post_de_") / "stage-result.json").exists()
    assert code == 0 and payload["state"] == "COMPLETED", payload                  # no crash, no untyped failure
    for stage in payload["stages"]:                                               # every reason is a typed code; no stage exception anywhere
        assert stage["reason_code"] is None or (stage["reason_code"].startswith("E_") and stage["reason_code"] != "E_STAGE_EXCEPTION"), stage
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    found = {(c["capability"], s["analysis"], s["reason_code"]) for c in plan["capability_plan"] if c["capability"].startswith("post_de_") for s in (c.get("post_de") or {}).get("subanalyses", [])}
    assert SUBANALYSES.get(name, set()) <= found
    report = json.loads((out / "post_de" / "eligibility" / "eligibility.json").read_text(encoding="utf-8"))
    assert {m["capability"]: m["state"] for m in report["modules"]} == {m: (e if e == C else e[0]) for m, e in EXPECTED[name].items()}


# ----------------------------------------------------------------------------- V166
def test_v166_report_lists_every_module_with_distinct_states_and_dependencies(tmp_path, monkeypatch):
    lib = tmp_path / "shadow-lib"; (lib / "randomForest").mkdir(parents=True)
    (lib / "randomForest" / "DESCRIPTION").write_text("Package: randomForest\nVersion: 0.0.0\nTitle: shadow\n", encoding="utf-8")
    previous = os.environ.get("R_LIBS")
    monkeypatch.setenv("R_LIBS", os.pathsep.join([str(lib)] + ([previous] if previous else [])))
    values, obs = F.phenotype_design(kind="one_group", n=10)
    for i, o in enumerate(obs):
        o["sex"] = "M" if i % 2 else "F"
    files = F.write_dataset(tmp_path / "data", values, obs, extra_columns=("age", "score", "sex"))
    post = {"enabled": True,
            "sets": {"enabled": True, "definitions": [{"id": "primary_any", "rule": F.leaf("protein-primary", "B-A")}]},
            "sensitivity": {"enabled": True, "covariate_models": [{"id": "adj_sex", "add_covariates": [{"column": "sex", "type": "categorical"}]}], "influence": True},
            "association": {"enabled": True, "phenotypes": [{"column": "score", "type": "numeric", "scope": "pooled"}]},                 # aliased: one group only
            "biomarker": {"enabled": True, "contrast": "B-A", "classifier": "random_forest", "cv": {"outer": {"k": 5, "repeats": 2}, "inner": {"k": 3}},
                          "selection": {"method": "top_k_auc", "k": [3]}, "permutation": {"enabled": False}, "single_feature": {"bootstrap": 200}},    # package absent
            "networks": {"enabled": True, "coabundance": {"enabled": True, "min_units": 50}}}                                              # small n
    path = F.config(tmp_path / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], post_de=post)
    payload, code = workflow.run_command(path, tmp_path / "run")
    out = tmp_path / "run"
    states = {s["stage_id"]: (s["state"], s["reason_code"]) for s in payload["stages"]}
    oracle = {"post_de_sets": ("COMPLETED", None), "post_de_sensitivity": ("COMPLETED", None), "post_de_association": ("INAPPLICABLE", "E_PHENOTYPE_ALIASED"),
              "post_de_biomarker": ("NOT_RUN", "E_ENGINE_NOT_AVAILABLE"), "post_de_networks": ("INAPPLICABLE", "E_COABUNDANCE_SMALL_N")}
    assert {k: states[k] for k in oracle} == oracle and payload["state"] == "PARTIAL"
    report = json.loads((out / "post_de" / "eligibility" / "eligibility.json").read_text(encoding="utf-8"))
    by = {m["capability"]: m for m in report["modules"]}
    for capability, (state, reason) in oracle.items():
        assert by[capability]["state"] == state and by[capability]["reason_code"] == reason
        if state != "COMPLETED":
            assert by[capability]["claim_label"] is None and by[capability]["inference"] is None   # no claim for refused modules
        assert bool(by[capability]["input_hashes"]) == (state == "COMPLETED")
    assert by["post_de_biomarker"]["plan_eligibility"] == "ELIGIBLE" and by["post_de_association"]["plan_eligibility"] == "INAPPLICABLE"
    deps = F.read_tsv(out / "post_de" / "eligibility" / "dependency_summary.tsv")
    primary = {r["feature_id"] for r in F.read_tsv(out / "dea" / "zero_null.tsv") if r["family_id"] == "protein-primary" and r["eligibility"] == "tested" and float(r["q_value"]) <= 0.05}
    assert {d["feature_id"] for d in deps} == primary and primary
    robust = {r["feature_id"]: r for r in F.read_tsv(out / "post_de" / "sensitivity" / "robustness_summary.tsv")}
    for d in deps:
        assert d["robustness_fraction"] == robust[d["feature_id"]]["robustness_fraction"] and d["sets"] == "primary_any"
        assert d["phenotype_associations"] == "INAPPLICABLE" and d["biomarker_selection_frequency"] == "NOT_RUN"     # module states, never zeros
    html = (out / "report-full" / "index.html").read_text(encoding="utf-8")
    assert "state-NOT_RUN" in html and "state-INAPPLICABLE" in html
    data = json.loads((out / "report" / "report_data.json").read_text(encoding="utf-8"))
    assert data["sections"]["post_de_association"]["values"] is None and data["sections"]["post_de_biomarker"]["values"] is None
    assert data["sections"]["post_de_association"]["state"] == "INAPPLICABLE" and data["sections"]["post_de_biomarker"]["state"] == "NOT_RUN"
    rows = data["sections"]["post_de_eligibility"]["tables"][0]["rows"]
    assert {(r[0], r[3]) for r in rows} == {(c.removeprefix("post_de_"), s) for c, (s, _) in oracle.items()}
    for r in rows:                                                                # a refused module is never shown as a zero result
        if r[3] != "COMPLETED":
            assert r[6] == "none (no values)" and r[7] == "none (not run)" and r[4].startswith("E_")
    by_row = {r[0]: r for r in rows}
    assert by_row["biomarker"][8] != by_row["association"][8]                     # NOT_RUN and INAPPLICABLE carry different meanings
