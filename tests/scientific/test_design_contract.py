"""R04 acceptance: V031-V040 (design grammar, alignment, rank, blocking, contrasts, plan)."""
from __future__ import annotations

import copy
import csv
import importlib.util
import json
import sys
import math
from pathlib import Path

import pytest

from proteomics_pipeline import design_service, workflow
from proteomics_pipeline.errors import ProteomicsError
from proteomics_pipeline.provenance import canonical_json_sha256

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)

needs_r = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")


def _obs(n_per_group, groups=("A", "B"), **extra):
    rows = []
    for g in groups:
        for i in range(1, n_per_group + 1):
            row = {"observation_id": f"{g}{i}", "group": g}
            for key, fn in extra.items():
                row[key] = fn(g, i)
            rows.append(row)
    return rows


# ----------------------------------------------------------------------------- V031
def test_v031_grammar_matches_explicit_rows():
    observations = [{"observation_id": "o1", "group": "A", "age": "30", "batch": "B1"}, {"observation_id": "o2", "group": "A", "age": "40", "batch": "B2"},
                    {"observation_id": "o3", "group": "B", "age": "50", "batch": "B1"}, {"observation_id": "o4", "group": "B", "age": "60", "batch": "B2"},
                    {"observation_id": "o5", "group": "B", "age": "45", "batch": "B2"}]
    design = {"id": "d", "group_column": "group", "group_levels": ["A", "B"], "continuous_covariates": [{"column": "age", "center": True}],
              "categorical_covariates": [{"column": "batch", "levels": ["B1", "B2"], "reference": "B1"}], "interactions": [["group", "age"]],
              "intercept": True, "blocking": {"mode": "none"}}
    compiled = design_service.compile_design(design, observations)
    assert compiled["coefficients"] == ["(Intercept)", "group.B", "continuous.age", "factor.batch.B2", "group.B:continuous.age"]
    mean_age = (30 + 40 + 50 + 60 + 45) / 5
    oracle = []
    for o in observations:
        b = 1.0 if o["group"] == "B" else 0.0
        age = float(o["age"]) - mean_age
        oracle.append([1.0, b, age, 1.0 if o["batch"] == "B2" else 0.0, b * age])
    assert compiled["matrix"] == oracle
    assert compiled["centers"] == {"age": mean_age}


def test_v031_percent_encoding_is_reversible():
    assert design_service.encode("a.b c%Ω") == "a%2Eb%20c%25%CE%A9"
    assert design_service.decode(design_service.encode("a.b c%Ω")) == "a.b c%Ω"
    assert design_service.encode("ctl_1-x") == "ctl_1-x"


@pytest.mark.parametrize("bad", [["group", "system('rm -rf /')"], ["source(evil.R)", "group"], ["group", "../../etc/passwd"]])
def test_v031_negative_code_or_undeclared_terms_fail_before_r(bad, monkeypatch):
    def forbidden(*_a, **_k):
        raise AssertionError("R must not be invoked")
    monkeypatch.setattr(design_service, "execute", forbidden)
    observations = [{"observation_id": "o1", "group": "A"}, {"observation_id": "o2", "group": "B"}]
    design = {"id": "d", "group_column": "group", "group_levels": ["A", "B"], "continuous_covariates": [], "categorical_covariates": [],
              "interactions": [bad], "intercept": False, "blocking": {"mode": "none"}}
    with pytest.raises(ProteomicsError) as error:
        design_service.compile_design(design, observations)
    assert error.value.code == "E_DESIGN_TERM"
    design["interactions"] = []; design["group_column"] = "system()"
    with pytest.raises(ProteomicsError) as error:
        design_service.compile_design(design, observations)
    assert error.value.code == "E_DESIGN_TERM"


# ----------------------------------------------------------------------------- V032
def test_v032_reordering_preserves_aligned_coefficients():
    observations = [{"observation_id": f"o{i}", "group": g, "batch": b} for i, (g, b) in enumerate([("A", "x"), ("A", "y"), ("B", "x"), ("B", "y"), ("C", "x"), ("C", "y")])]
    design = {"id": "d", "group_column": "group", "group_levels": ["A", "B", "C"], "continuous_covariates": [],
              "categorical_covariates": [{"column": "batch", "levels": ["x", "y"], "reference": "x"}], "interactions": [], "intercept": True, "blocking": {"mode": "none"}}
    first = design_service.compile_design(design, observations)
    shuffled = design_service.compile_design(design, list(reversed(observations)))
    assert first["coefficients"] == shuffled["coefficients"]
    by_id_first = dict(zip(first["observation_ids"], first["matrix"]))
    by_id_second = dict(zip(shuffled["observation_ids"], shuffled["matrix"]))
    assert by_id_first == by_id_second
    assert first["term_map"][1] == {"coefficient": "group.B", "term": "group", "column": "group", "level": "B"}


def test_v032_negative_missing_reference_and_nonfinite_covariate():
    observations = [{"observation_id": "o1", "group": "A", "batch": "x", "age": "1"}, {"observation_id": "o2", "group": "B", "batch": "y", "age": "NA"}]
    design = {"id": "d", "group_column": "group", "group_levels": ["A", "B"], "continuous_covariates": [], "intercept": False, "interactions": [],
              "categorical_covariates": [{"column": "batch", "levels": ["x", "y"], "reference": "z"}], "blocking": {"mode": "none"}}
    with pytest.raises(ProteomicsError) as error:
        design_service.compile_design(design, observations)
    assert error.value.code == "E_DESIGN_LEVEL" and error.value.pointer.endswith("/reference")
    design["categorical_covariates"] = []; design["continuous_covariates"] = [{"column": "age", "center": False}]
    with pytest.raises(ProteomicsError) as error:
        design_service.compile_design(design, observations)
    assert error.value.code == "E_DESIGN_COVARIATE_NONFINITE" and "age" in error.value.pointer


# ----------------------------------------------------------------------------- shared run helper
def _values(observations, features=6, effect=lambda f, o: 0.0, missing=()):
    out = {}
    for f in range(features):
        vals = []
        for i, o in enumerate(observations):
            noise = ((i * 7 + f * 3) % 5 - 2) * 0.11 + ((i + f) % 3) * 0.05
            vals.append(None if (f, o["observation_id"]) in missing else 10 + f + effect(f, o) + noise)
        out[f"F{f + 1}"] = vals
    return out


def _plan(tmp_path, observations, contrasts, groups, **kw):
    files = B.dataset(tmp_path / "data", _values(observations, **kw.pop("values_kw", {})), observations, extra_columns=kw.pop("extra_columns", ()))
    path = B.config(tmp_path / "data", files, groups=groups, contrasts=contrasts, **kw)
    return path


# ----------------------------------------------------------------------------- V033
@needs_r
def test_v033_confounded_batch_rejected_crossed_accepted(tmp_path):
    obs = _obs(4, batch=lambda g, i: "b1" if g == "A" else "b2")
    over = {"categorical_covariates": [{"column": "batch", "levels": ["b1", "b2"], "reference": "b1"}], "intercept": True}
    contrasts = [B.contrast("B-A", "B", "A", weights={"group.B": 1})]
    path = _plan(tmp_path / "conf", obs, contrasts, ["A", "B"], design_overrides=over, extra_columns=("batch",))
    with pytest.raises(ProteomicsError) as error:
        workflow.plan_command(path, tmp_path / "conf" / "out" / "plan.json")
    assert error.value.code == "E_DESIGN_CONFOUNDED" and "factor.batch.b2" in error.value.message
    assert not (tmp_path / "conf" / "out" / "plan.json").exists()
    crossed = _obs(4, batch=lambda g, i: "b1" if i % 2 else "b2")
    path = _plan(tmp_path / "cross", crossed, contrasts, ["A", "B"], design_overrides=over, extra_columns=("batch",))
    payload, code = workflow.plan_command(path, tmp_path / "cross" / "out" / "plan.json")
    assert code == 0 and payload["fits_performed"] == 0
    diagnostics = json.loads((tmp_path / "cross" / "out" / "plan-artifacts" / "designs" / "design_diagnostics.json").read_text())
    assert diagnostics["designs"][0]["rank"] == 3 and diagnostics["designs"][0]["full_rank"] is True


@needs_r
def test_v033_rank_matches_independent_qr(tmp_path):
    code = ("X <- cbind(1, c(0,0,0,0,1,1,1,1), c(0,0,0,0,1,1,1,1)); X2 <- cbind(1, c(0,0,0,0,1,1,1,1), c(0,1,0,1,0,1,0,1)); "
            "r <- get('design_rank', asNamespace('proteomicsCore')); colnames(X) <- colnames(X2) <- c('a','b','c'); "
            "cat(jsonlite::toJSON(list(conf=r(X)$rank, conf_ref=qr(X)$rank, alias=r(X)$aliased, crossed=r(X2)$rank, ref=sum(svd(X2)$d > 1e-10)), auto_unbox=TRUE))")
    value = B.r_json(code)
    assert value["conf"] == value["conf_ref"] == 2 and value["alias"] == "c" and value["crossed"] == value["ref"] == 3


# ----------------------------------------------------------------------------- V034
@needs_r
def test_v034_paired_fixed_subject_and_repeated_correlation(tmp_path):
    obs = [{"observation_id": f"S{s}_{v}", "group": v, "subject_id": f"S{s}", "biological_unit_id": f"S{s}_{v}"} for s in range(1, 5) for v in ("before", "after")]
    contrasts = [B.contrast("after-before", "after", "before")]
    fixed = _plan(tmp_path / "fixed", obs, contrasts, ["before", "after"], design_overrides={"blocking": {"mode": "fixed_subject", "subject_column": "subject_id"}})
    workflow.plan_command(fixed, tmp_path / "fixed" / "out" / "plan.json")
    term_map = json.loads((tmp_path / "fixed" / "out" / "plan-artifacts" / "designs" / "designs" / "joint" / "term_map.json").read_text())
    assert term_map["coefficients"] == ["group.before", "group.after", "factor.subject_id.S2", "factor.subject_id.S3", "factor.subject_id.S4"]
    diag = json.loads((tmp_path / "fixed" / "out" / "plan-artifacts" / "designs" / "design_diagnostics.json").read_text())
    assert diag["contrasts"][0]["subjects_by_required_group"] == {"after": 4, "before": 4}
    repeated = _plan(tmp_path / "rep", obs, contrasts, ["before", "after"], design_overrides={"blocking": {"mode": "duplicate_correlation", "subject_column": "subject_id"}})
    workflow.plan_command(repeated, tmp_path / "rep" / "out" / "plan.json")
    plan = json.loads((tmp_path / "rep" / "out" / "plan.json").read_text())
    rho = plan["blocking"][0]["consensus_correlation"]
    oracle = B.r_json("m <- as.matrix(read.delim(commandArgs(TRUE)[1], row.names=1, check.names=FALSE)); X <- cbind(before=rep(c(1,0),4), after=rep(c(0,1),4)); "
                      "cat(jsonlite::toJSON(limma::duplicateCorrelation(m, X, block=rep(paste0('S',1:4), each=2))$consensus.correlation, digits=NA))",
                      str(tmp_path / "rep" / "data" / "abundance.tsv"))
    assert abs(rho - oracle[0]) <= 1e-10


def test_v034_negative_missing_subject_and_double_encoding():
    obs = [{"observation_id": "a", "group": "x", "subject_id": "S1"}, {"observation_id": "b", "group": "y", "subject_id": "NA"}]
    design = {"id": "d", "group_column": "group", "group_levels": ["x", "y"], "continuous_covariates": [], "categorical_covariates": [], "interactions": [],
              "intercept": False, "blocking": {"mode": "duplicate_correlation", "subject_column": "subject_id"}}
    with pytest.raises(ProteomicsError) as error:
        design_service.compile_design(design, obs)
    assert error.value.code == "E_SUBJECT_MISSING"
    obs[1]["subject_id"] = "S2"
    design["categorical_covariates"] = [{"column": "subject_id", "levels": ["S1", "S2"], "reference": "S1"}]
    with pytest.raises(ProteomicsError) as error:
        design_service.compile_design(design, obs)
    assert error.value.code == "E_BLOCKING_CONFLICT"


# ----------------------------------------------------------------------------- V035
@needs_r
def test_v035_numeric_contrasts_and_interaction_from_weights(tmp_path):
    means = {"C": 10.0, "U": 12.0, "T": 11.0}
    offsets = [-0.3, -0.1, 0.1, 0.3]
    obs = [{"observation_id": f"{g}{i}", "group": g} for g in "CUT" for i in range(4)]
    values = {"F1": [means[o["group"]] + offsets[int(o["observation_id"][1])] for o in obs], "F2": [5 + offsets[(int(o["observation_id"][1]) + 1) % 4] for o in obs]}
    files = B.dataset(tmp_path / "d", values, obs)
    contrasts = [B.contrast("d", "U", "C"), B.contrast("t", "T", "U"), B.contrast("r", "T", "C", role="secondary")]
    path = B.config(tmp_path / "d", files, groups=["C", "U", "T"], contrasts=contrasts)
    payload, code = workflow.run_command(path, tmp_path / "run")
    rows = {(r["contrast_id"], r["feature_id"]): r for r in B.read_tsv(tmp_path / "run" / "dea" / "zero_null.tsv")}
    d, t, r = (float(rows[(c, "F1")]["effect"]) for c in ("d", "t", "r"))
    assert abs(d - 2) <= 1e-10 and abs(t + 1) <= 1e-10 and abs(r - 1) <= 1e-10 and abs(r - (d + t)) <= 1e-10
    # interaction of two treatment effects computed directly from coefficient weights
    w_t1 = {"group.T": 1, "group.U": -1}; w_t2 = {"group.U": 1, "group.C": -1}
    interaction = {k: w_t1.get(k, 0) - w_t2.get(k, 0) for k in set(w_t1) | set(w_t2)}
    assert sum(v * means[k.split(".")[1]] for k, v in interaction.items()) == pytest.approx(t - d, abs=1e-12)


def test_v035_negative_reversed_weights_and_unknown_coefficient():
    compiled = {"joint": {"coefficients": ["group.C", "group.U"], "group_levels": ["C", "U"]}}
    config = {"contrasts": [{"id": "x", "design_id": "joint", "label": "", "estimand": "", "role": "primary", "required_groups": ["U", "C"], "weights": {"group.Q": 1}}]}
    with pytest.raises(ProteomicsError) as error:
        design_service.compile_contrasts(config, compiled)
    assert error.value.code == "E_CONTRAST_COEFFICIENT"


@needs_r
def test_v035_negative_reversed_direction_is_rejected(tmp_path):
    obs = _obs(3)
    path = _plan(tmp_path, obs, [B.contrast("B-A", "B", "A", weights={"group.B": -1, "group.A": 1})], ["A", "B"])
    with pytest.raises(ProteomicsError) as error:
        workflow.plan_command(path, tmp_path / "out" / "plan.json")
    assert error.value.code == "E_CONTRAST_DIRECTION"


# ----------------------------------------------------------------------------- V036
@needs_r
def test_v036_featurewise_n_df_and_named_reasons(tmp_path):
    obs = _obs(4)
    missing = {(1, "B1"), (1, "B2"), (1, "B3"), (1, "B4"), (2, "A3")}
    path = _plan(tmp_path, obs, [B.contrast("B-A", "B", "A")], ["A", "B"], values_kw={"features": 4, "missing": missing})
    workflow.plan_command(path, tmp_path / "out" / "plan.json")
    est = {r["feature_id"]: r for r in B.read_tsv(tmp_path / "out" / "plan-artifacts" / "designs" / "estimability.tsv")}
    assert est["F1"]["n_obs"] == "8" and est["F1"]["df_residual"] == "6" and est["F1"]["eligibility"] == "eligible"
    assert est["F2"]["eligibility"] == "nonestimable" and est["F2"]["reason"] == "all_missing_required_group:B"
    assert est["F3"]["n_obs"] == "7" and est["F3"]["df_residual"] == "5" and est["F3"]["n_obs_by_required_group"] == '{"B":4,"A":3}'
    assert len(est) == 4                                         # excluded rows retained


# ----------------------------------------------------------------------------- V038
def _phase2(engine, **model):
    config = json.loads((B.EXAMPLES / "example-independent.json").read_text())
    config["runtime"]["phase"] = 2
    m = dict(config["models"][0]); m.update({"id": f"{engine}-opt", "engine": engine, "role": "sensitivity", "execution_requirement": "optional"}); m.update(model)
    config["models"].append(m)
    return config


@pytest.mark.parametrize("engine,model,design,assay,code", [
    ("deqms", {"precision_weights": {"kind": "observation", "path": "w.tsv"}, "count_evidence": "c.tsv"}, None, None, "E_DEQMS_DESIGN_UNSUPPORTED"),
    ("deqms", {"count_evidence": "c.tsv"}, {"mode": "fixed_subject", "subject_column": "subject_id"}, None, "E_DEQMS_DESIGN_UNSUPPORTED"),
    ("deqms", {"count_evidence": "c.tsv", "hypothesis": "treat", "effect_threshold": 0.5}, None, None, "E_DEQMS_HYPOTHESIS_UNSUPPORTED"),
    ("deqms", {}, None, None, "E_DEQMS_COUNT_EVIDENCE"),
    ("proda", {}, None, "tmt", "E_PRODA_LFQ_UNIMPUTED_REQUIRED"),
    ("proda", {}, {"mode": "fixed_subject", "subject_column": "subject_id"}, "lfq_dda", "E_PRODA_DESIGN_UNSUPPORTED"),
    ("proda", {"hypothesis": "treat", "effect_threshold": 0.5}, None, "lfq_dda", "E_PRODA_HYPOTHESIS_UNSUPPORTED"),
])
def test_v038_adapter_eligibility_is_declarative(engine, model, design, assay, code):
    config = _phase2(engine, **model)
    if design:
        config["design"]["blocking"] = design
    if assay:
        config["assay"] = assay
    rows = {r["model_id"]: r for r in design_service.eligibility_table(config, mask_state="known", prior_imputation="none_documented")}
    assert rows["limma-main"]["eligible"] is True
    assert rows[f"{engine}-opt"]["eligible"] is False and rows[f"{engine}-opt"]["reason_code"] == code


def test_v038_negative_absent_adapter_is_not_run_never_limma(tmp_path, monkeypatch):
    # Amendment A-2026-10-01-08: R06 now implements assay_engines, so absence is simulated explicitly.
    monkeypatch.setitem(sys.modules, "proteomics_pipeline.assay_service", None)
    config = _phase2("deqms", count_evidence="c.tsv")
    requested = {r["capability"]: r for r in workflow.requested_capabilities(config)}
    assert requested["assay_engines"]["implemented"] is False and requested["assay_engines"]["required"] is False
    config["models"][1]["execution_requirement"] = "required"
    requested = workflow.requested_capabilities(config)
    with pytest.raises(ProteomicsError) as error:
        workflow._required_unavailable(requested)
    assert error.value.code == "E_CAPABILITY_NOT_IMPLEMENTED" and error.value.exit_code == 3


# ----------------------------------------------------------------------------- V039
@needs_r
def test_v039_plan_frozen_before_fit_and_meaningful_changes_change_hash(tmp_path):
    obs = _obs(4)
    contrasts = [B.contrast("B-A", "B", "A")]
    base = _plan(tmp_path / "a", obs, contrasts, ["A", "B"])
    p1, _ = workflow.plan_command(base, tmp_path / "a" / "o1" / "plan.json")
    p2, _ = workflow.plan_command(base, tmp_path / "a" / "o2" / "plan.json")
    assert p1["plan_hash"] == p2["plan_hash"] and p1["fits_performed"] == 0
    assert not list((tmp_path / "a" / "o1").rglob("zero_null.tsv"))
    plan = json.loads((tmp_path / "a" / "o1" / "plan.json").read_text())
    kinds = {a["result_type"] for a in plan["artifacts"]}
    assert {"Design", "CoverageResult", "OriginalObservedMask", "DesignDiagnostics"} <= kinds and plan["families"][0]["members"] == [{"model_id": "limma-main", "contrast_id": "B-A"}]
    weighted = _plan(tmp_path / "b", obs, [B.contrast("B-A", "B", "A", weights={"group.B": 2, "group.A": -2})], ["A", "B"])
    excluded = _plan(tmp_path / "c", obs, contrasts, ["A", "B"], mutate=lambda c: c["preprocessing"].update({"exclusions": [{"observation_id": "A1", "reason": "prespecified"}]}))
    reason = _plan(tmp_path / "d", obs, contrasts, ["A", "B"], mutate=lambda c: c["preprocessing"].update({"exclusions": [{"observation_id": "A1", "reason": "another reason"}]}))
    hashes = {workflow.plan_command(p, p.parent / "o" / "plan.json")[0]["plan_hash"] for p in (weighted, excluded, reason)}
    assert len(hashes | {p1["plan_hash"]}) == 4                  # V039 + V027 rerun on the actual AnalysisPlan
    resource = copy.deepcopy(plan["config"]); resource["resources"] = [{"id": "r", "sha256": "a" * 64}]
    changed = copy.deepcopy(resource); changed["resources"][0]["sha256"] = "b" * 64
    assert canonical_json_sha256(resource) != canonical_json_sha256(changed)


@needs_r
def test_v039_negative_fit_without_valid_plan_or_with_changed_input(tmp_path):
    from proteomics_pipeline import inference_service
    obs = _obs(4)
    path = _plan(tmp_path, obs, [B.contrast("B-A", "B", "A")], ["A", "B"])
    workflow.plan_command(path, tmp_path / "out" / "plan.json")
    plan_path = tmp_path / "out" / "plan.json"
    plan = json.loads(plan_path.read_text())
    config = workflow.load(path)
    with pytest.raises(ProteomicsError) as error:
        design_service.verify_plan(tmp_path / "missing.json")
    assert error.value.code == "E_PLAN_REQUIRED"
    request = inference_service.build_request(plan, plan_path=plan_path, config=config, run_id="neg", output_temp_dir=tmp_path / "fit1")
    request["plan_hash"] = "0" * 64
    result = inference_service.execute(request)
    assert result["state"] == "FAILED" and result["reason_code"] == "E_PLAN_CHANGED" and result["exit_code"] == 5
    matrix = tmp_path / "out" / "plan-artifacts" / "preprocessing" / "primary" / "matrix.tsv"
    request = inference_service.build_request(plan, plan_path=plan_path, config=config, run_id="neg", output_temp_dir=tmp_path / "fit2")
    from proteomics_pipeline.provenance import sha256_file
    matrix.write_text(matrix.read_text().replace("\t1", "\t2", 1))
    for item in request["inputs"]:
        if item["artifact_id"] == "primary_matrix":
            item["sha256"] = sha256_file(matrix)
    result = inference_service.execute(request)
    assert result["state"] == "FAILED" and result["reason_code"] == "E_PLAN_CHANGED"
    with pytest.raises(ProteomicsError) as error:
        design_service.verify_plan(plan_path)
    assert error.value.code == "E_PLAN_CHANGED"


@needs_r
def test_v026_negative_display_matrix_cannot_feed_the_model(tmp_path):
    from proteomics_pipeline import inference_service
    from proteomics_pipeline.provenance import sha256_file
    obs = _obs(4)
    path = _plan(tmp_path, obs, [B.contrast("B-A", "B", "A")], ["A", "B"])
    workflow.plan_command(path, tmp_path / "out" / "plan.json")
    plan_path = tmp_path / "out" / "plan.json"; plan = json.loads(plan_path.read_text())
    request = inference_service.build_request(plan, plan_path=plan_path, config=workflow.load(path), run_id="neg", output_temp_dir=tmp_path / "fit")
    display = next(a for a in plan["artifacts"] if a["artifact_id"] == "qc_pca_display_input")
    request["inputs"].append({"artifact_id": "qc_pca_display_input", "path": str(tmp_path / "out" / display["relative_path"]), "sha256": display["sha256"]})
    request["parameters"]["models"][0]["matrix_artifact"] = "qc_pca_display_input"
    result = inference_service.execute(request)
    assert result["state"] == "FAILED" and result["reason_code"] == "E_DISPLAY_MATRIX_REJECTED"


# ----------------------------------------------------------------------------- V040
@needs_r
def test_v040_adversarial_designs(tmp_path):
    unequal = [{"observation_id": f"A{i}", "group": "A"} for i in range(3)] + [{"observation_id": f"B{i}", "group": "B"} for i in range(5)]
    payload, code = workflow.plan_command(_plan(tmp_path / "u", unequal, [B.contrast("B-A", "B", "A")], ["A", "B"]), tmp_path / "u" / "o" / "plan.json")
    assert code == 0
    singleton = [{"observation_id": "A0", "group": "A"}] + [{"observation_id": f"B{i}", "group": "B"} for i in range(4)]
    path = _plan(tmp_path / "s", singleton, [B.contrast("B-A", "B", "A")], ["A", "B"])
    with pytest.raises(ProteomicsError) as error:
        workflow.plan_command(path, tmp_path / "s" / "o" / "plan.json")
    assert error.value.code == "E_INSUFFICIENT_REPLICATION"
    qc_only = _plan(tmp_path / "q", singleton, [B.contrast("B-A", "B", "A")], ["A", "B"], mutate=lambda c: c["runtime"].update({"scope": "qc_only"}))
    assert workflow.plan_command(qc_only, tmp_path / "q" / "o" / "plan.json")[1] == 0
    absent = _plan(tmp_path / "a", _obs(3), [B.contrast("B-A", "B", "A"), B.contrast("C-A", "C", "A", role="secondary")], ["A", "B", "C"])
    with pytest.raises(ProteomicsError) as error:
        workflow.plan_command(absent, tmp_path / "a" / "o" / "plan.json")
    assert error.value.code == "E_DESIGN_GROUP_EMPTY"
    imbalanced = [{"observation_id": f"S{s}_{v}", "group": v, "subject_id": f"S{s}"} for s in range(1, 5) for v in ("pre", "post") if not (s == 4 and v == "post")]
    path = _plan(tmp_path / "i", imbalanced, [B.contrast("post-pre", "post", "pre")], ["pre", "post"], design_overrides={"blocking": {"mode": "fixed_subject", "subject_column": "subject_id"}})
    assert workflow.plan_command(path, tmp_path / "i" / "o" / "plan.json")[1] == 0
    diag = json.loads((tmp_path / "i" / "o" / "plan-artifacts" / "designs" / "design_diagnostics.json").read_text())
    assert diag["contrasts"][0]["subjects_by_required_group"] == {"post": 3, "pre": 4}


@needs_r
def test_v040_reordered_ids_are_invariant(tmp_path):
    obs = _obs(4)
    first = _plan(tmp_path / "f", obs, [B.contrast("B-A", "B", "A")], ["A", "B"])
    workflow.plan_command(first, tmp_path / "f" / "o" / "plan.json")
    reordered = list(reversed(obs))
    files = B.dataset(tmp_path / "r" / "data", {k: list(reversed(v)) for k, v in _values(obs).items()}, reordered)
    path = B.config(tmp_path / "r" / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")])
    workflow.plan_command(path, tmp_path / "r" / "o" / "plan.json")
    def keyed(root):
        return {(r["feature_id"], r["contrast_id"]): (r["n_obs"], r["df_residual"], r["eligibility"]) for r in B.read_tsv(root / "plan-artifacts" / "designs" / "estimability.tsv")}
    assert keyed(tmp_path / "f" / "o") == keyed(tmp_path / "r" / "o")
