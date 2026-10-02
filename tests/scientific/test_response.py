"""R09 acceptance through the real pipeline: V081, V082, V085-V087, V089, V090 (with R-level V083/V084/V088 in test-response.R)."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_pipeline import workflow

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "configs" / "examples"
FIX = ROOT / "tests" / "fixtures" / "response"
spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")
SHA = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
AXIS = {"id": "axis1", "model_id": "limma-main", "disease_contrast": "disease-control", "treatment_contrast": "treated-disease", "residual_contrast": "treated-control"}


def resource(rid, kind, path):
    return {"id": rid, "kind": kind, "path": str(path), "sha256": SHA(path), "version": "synthetic-1", "source": "synthetic", "source_taxonomy_id": 10116,
            "target_taxonomy_id": 10116, "terms": "synthetic"}


def family(fid, htype, **sel):
    return {"id": fid, "hypothesis_type": htype, "model_ids": ["limma-main"], "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary", "contrast_ids": [], **sel}


def run(tmp_path, mutate):
    raw = json.loads((EXAMPLES / "example-independent.json").read_text())
    for key in ("matrix", "observations", "features", "source_provenance"):
        raw["input"][key] = str((EXAMPLES / raw["input"][key]).resolve())
    raw["runtime"]["phase"] = 2
    raw["response"] = {"enabled": True, "axes": [dict(AXIS)], "selection": "all_eligible", "ri_minimum_abs_disease_effect": 0.25, "ratio_uncertainty": "fieller", "alpha": 0.05}
    mutate(raw)
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "analysis.json"; path.write_text(json.dumps(raw))
    payload, code = workflow.run_command(path, tmp_path / "run")
    return payload, code, tmp_path / "run"


def test_v081_v090_axes_covariance_and_algebra(tmp_path):
    payload, code, out = run(tmp_path, lambda raw: None)
    assert code == 0, payload
    rows = B.read_tsv(out / "response" / "descriptive.tsv")
    assert len(rows) == 8 and all(r["result_type"] == "DescriptiveResponse" for r in rows)
    for r in rows:
        if r["eligibility"] != "tested":
            continue
        assert abs(float(r["r"]) - (float(r["d"]) + float(r["t"]))) <= 1e-10
        if r["feature_id"] in ("P01", "P02", "P03", "P04", "P05", "P06"):           # complete features: group means of 4 -> diag(.25)
            assert abs(float(r["unscaled_var_d"]) - 0.5) <= 1e-12 and abs(float(r["unscaled_cov_dt"]) + 0.25) <= 1e-12 and abs(float(r["unscaled_var_r"]) - 0.5) <= 1e-12
    header = set(rows[0]); assert not {"p_value", "q_value", "native_q_value", "statistic"} & header      # descriptive: physically no P/q
    warnings = json.loads((out / "warnings.json").read_text())
    assert any(w["code"] == "W_SHARED_CONTROL" for w in warnings)
    p03 = next(r for r in rows if r["feature_id"] == "P03")
    assert p03["descriptive_class"] == "near_restoration" and p03["crossed_control"] in ("true", "false")


def test_v089_near_restoration_is_not_equivalence_and_tost_family(tmp_path):
    def mutate(raw):
        raw["response_mode"] = "equivalence"; raw["response"]["equivalence_margin"] = 0.02
        raw["multiplicity_families"].append(family("equiv", "equivalence", axis_ids=["axis1"]))
    payload, code, out = run(tmp_path, mutate)
    assert code == 0, payload
    eq = {r["feature_id"]: r for r in B.read_tsv(out / "response" / "equivalence.tsv")}
    desc = {r["feature_id"]: r for r in B.read_tsv(out / "response" / "descriptive.tsv")}
    from statistics import NormalDist  # noqa: F401
    for fid, r in eq.items():
        if r["eligibility"] != "tested":
            continue
        est, se, df = float(r["residual_estimate"]), float(r["residual_se"]), float(r["df"])
        oracle = B.r_json("a <- as.numeric(commandArgs(TRUE)); cat(jsonlite::toJSON(c(pt((a[1]+a[4])/a[2], a[3], lower.tail=FALSE), pt((a[1]-a[4])/a[2], a[3])), digits=NA))", repr(est), repr(se), repr(df), "0.02")
        assert abs(float(r["p_lower"]) - oracle[0]) <= 1e-10 and abs(float(r["p_upper"]) - oracle[1]) <= 1e-10
        assert float(r["p_value"]) == max(float(r["p_lower"]), float(r["p_upper"])) and r["family_id"] == "equiv"
    near = [f for f, r in desc.items() if r["descriptive_class"] == "near_restoration"]
    assert near and any(float(eq[f]["p_value"]) >= 0.05 for f in near)          # near_restoration that fails TOST stays not-equivalent
    assert "equivalen" not in json.dumps([desc[f] for f in near]) and not {"p_value", "q_value"} & set(next(iter(desc.values())))


def test_v086_formal_rescue_with_independent_directions(tmp_path):
    def mutate(raw):
        raw["response_mode"] = "formal_rescue"
        raw["response"].update({"equivalence_margin": 0.5, "direction_resource_id": "dir", "disease_margin": 0.0, "treatment_margin": 0.0})
        raw["resources"] = [resource("dir", "independent_direction", FIX / "direction_independent.json")]
        raw["multiplicity_families"] += [family("equiv", "equivalence", axis_ids=["axis1"]), family("rescue", "formal_rescue", axis_ids=["axis1"])]
    payload, code, out = run(tmp_path, mutate)
    assert code == 0, payload
    rows = {r["feature_id"]: r for r in B.read_tsv(out / "response" / "formal_rescue.tsv")}
    assert rows["P01"]["eligibility"] == "excluded" and rows["P01"]["reason_code"] == "no_independent_direction"
    for r in rows.values():
        if r["eligibility"] == "tested":
            comps = [float(r[k]) for k in ("p_disease", "p_treatment", "p_residual_lower", "p_residual_upper")]
            assert float(r["p_value"]) == max(comps) and r["direction_resource_sha256"] == SHA(FIX / "direction_independent.json")


def test_v086_negative_direction_overlapping_tested_units(tmp_path):
    bad = json.loads((FIX / "direction_independent.json").read_text()); bad["selection_subject_ids"] = ["unit_C1"]
    (tmp_path / "dir.json").parent.mkdir(parents=True, exist_ok=True); (tmp_path / "dir.json").write_text(json.dumps(bad))
    def mutate(raw):
        raw["response_mode"] = "formal_rescue"
        raw["response"].update({"equivalence_margin": 0.5, "direction_resource_id": "dir", "disease_margin": 0.0, "treatment_margin": 0.0})
        raw["resources"] = [resource("dir", "independent_direction", tmp_path / "dir.json")]
        raw["multiplicity_families"] += [family("equiv", "equivalence", axis_ids=["axis1"]), family("rescue", "formal_rescue", axis_ids=["axis1"])]
    payload, code, out = run(tmp_path / "r", mutate)
    stage = next(s for s in payload["stages"] if s["stage_id"] == "response")
    assert stage["state"] == "FAILED" and stage["reason_code"] == "E_RESCUE_DIRECTION_NOT_INDEPENDENT" and code != 0
    assert not (out / "response").exists()


@pytest.mark.parametrize("manifest,inferential", [("score_independent.json", True), ("score_overlap.json", False)])
def test_v082_v087_v088_scores(tmp_path, manifest, inferential):
    def mutate(raw):
        raw["score_test"] = "independent_exact"
        raw["resources"] = [resource("score", "independent_score", FIX / manifest)]
        raw["response"]["scores"] = [{"id": "s1", "model_id": "limma-main", "contrast_id": "disease-control", "resource_id": "score", "scheme": "independent_labels",
                                      "exchangeability_evidence": "synthetic independent specimens"}]
        raw["multiplicity_families"].append(family("scores", "independent_score_test", score_ids=["s1"]))
    payload, code, out = run(tmp_path, mutate)
    assert code == 0, payload
    if inferential:
        test = B.read_tsv(out / "response" / "independent_score_tests.tsv")[0]
        assert test["N"] == "70" and test["observed_allocation_count"] == "1" and test["method"] == "exact"
        values = json.loads((FIX / manifest).read_text())["features"]
        m = {r["feature_id"]: r for r in B.read_tsv(out / "inputs" / "matrix.tsv")}
        groups = {"C": [f"C{i}" for i in range(1, 5)], "U": [f"U{i}" for i in range(1, 5)]}
        score = {o: sum(f["weight"] * (float(m[f["feature_id"]][o]) - f["center"]) / f["scale"] for f in values) for g in groups.values() for o in g}
        import itertools
        ids = groups["U"] + groups["C"]
        obs = sum(score[o] for o in groups["U"]) / 4 - sum(score[o] for o in groups["C"]) / 4
        k = sum(1 for comb in itertools.combinations(ids, 4) if abs(sum(score[o] for o in comb) / 4 - sum(score[o] for o in ids if o not in comb) / 4) >= abs(obs) - (1e-12 + 1e-10 * abs(obs)))
        assert int(test["k"]) == k and abs(float(test["p_value"]) - k / 70) <= 1e-12
        assert not (out / "response" / "descriptive_scores.tsv").exists()
    else:
        rows = B.read_tsv(out / "response" / "descriptive_scores.tsv")
        assert {r["score_eligibility_reason"] for r in rows} == {"E_SCORE_SELECTION_OVERLAP"} and rows[0]["hypothesis_type"] == "descriptive_score"
        assert not {"p_value", "q_value", "statistic"} & set(rows[0])
        assert not (out / "response" / "independent_score_tests.tsv").exists()


def test_v081_negative_incoherent_axis(tmp_path):
    def mutate(raw):
        raw["response"]["axes"][0]["residual_contrast"] = "disease-control"
    payload, code, out = run(tmp_path, mutate)
    assert code == 2 and payload["error"]["code"] == "E_AXIS_INCOHERENT"
