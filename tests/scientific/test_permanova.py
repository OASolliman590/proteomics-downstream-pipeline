"""R13 acceptance through the real Python -> R pipeline: V121, V122, V125, V126, V129, V130."""
from __future__ import annotations

import html.parser
import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


B = _load("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
F = _load("permanova_fixtures", ROOT / "tests" / "fixtures" / "permanova" / "make_fixtures.py")


def _has_vegan():
    if B.rscript() is None:
        return False
    try:
        return bool(B.r_json("cat(jsonlite::toJSON(requireNamespace('vegan', quietly=TRUE) && requireNamespace('permute', quietly=TRUE) && exists('permanova_stage', envir=asNamespace('proteomicsCore')), auto_unbox=TRUE))"))
    except RuntimeError:
        return False


pytestmark = pytest.mark.skipif(not _has_vegan(), reason="NOT_RUN: Rscript/proteomicsCore/vegan/permute unavailable")


def run(tmp_path, kind, *, groups=("A", "B"), permanova=None, required="optional", extra=None, contrasts=None, profile="smoke_test", seed=20261001):
    values, obs = F.matrix(kind, groups=groups, seed=seed)
    if extra:
        for o in obs:
            o.update(extra(o))
    files = B.dataset(tmp_path / "data", values, obs, extra_columns=tuple(extra(obs[0]).keys()) if extra else ())
    contrasts = contrasts or [B.contrast(f"{g}-{groups[0]}", g, groups[0], role="primary" if i == 0 else "secondary") for i, g in enumerate(groups[1:])]
    def mutate(c):
        c["runtime"]["phase"] = 2; c["runtime"]["execution_profile"] = profile
        c["multivariate"] = {"enabled": True, "execution_requirement": required, "permanova": {"permutations": 999, "random_sets": 300, **(permanova or {})}}
    path = B.config(tmp_path / "data", files, groups=list(groups), contrasts=contrasts, mutate=mutate)
    payload, code = workflow.run_command(path, tmp_path / "run")
    return payload, code, tmp_path / "run"


def _tests_table(run_dir):
    return B.read_tsv(run_dir / "permanova" / "tests.tsv")


# ----------------------------------------------------------------------------- V130 / V126
def test_v130_null_location_and_dispersion_fixtures_behave_as_constructed(tmp_path):
    _, code_null, null_run = run(tmp_path / "null", "null")
    _, code_loc, loc_run = run(tmp_path / "loc", "location", groups=("A", "B", "C"))   # 3 groups: no label-swap symmetry, so the floor is attainable
    _, code_disp, disp_run = run(tmp_path / "disp", "dispersion")
    assert code_null == code_loc == code_disp == 0
    null = next(r for r in _tests_table(null_run) if r["analysis"] == "global")
    loc = next(r for r in _tests_table(loc_run) if r["analysis"] == "global")
    disp = next(r for r in _tests_table(disp_run) if r["analysis"] == "global")
    assert float(null["p_value"]) > 0.05 and float(null["r2"]) < float(loc["r2"])
    assert float(loc["p_value"]) == float(loc["p_floor"]) == 1 / 1000 and loc["p_display"] == "< 0.001" and float(loc["r2"]) > 0.3
    assert float(disp["permdisp_p"]) < 0.05 and float(disp["p_value"]) > 0.05
    assert disp["interpretation"] == "dispersion_difference_location_not_established"
    assert loc["interpretation"] in ("location_shift", "location_and_or_dispersion")


def test_v130_identical_seed_gives_identical_tables(tmp_path):
    run(tmp_path / "a", "location", groups=("A", "B", "C"))
    run(tmp_path / "b", "location", groups=("A", "B", "C"))
    for name in ("tests.tsv", "feature_r2.tsv", "permutation_null.tsv", "ordination.tsv"):
        assert (tmp_path / "a" / "run" / "permanova" / name).read_bytes() == (tmp_path / "b" / "run" / "permanova" / name).read_bytes(), name


# ----------------------------------------------------------------------------- V121 / V128
def test_v121_v128_dep_derived_and_same_data_panel_get_circularity_context(tmp_path):
    sets = [{"id": "all", "kind": "all_complete"},
            {"id": "same", "kind": "declared_panel", "feature_ids": ["F01", "F02", "F03"], "selection_provenance": "same_data"},
            {"id": "indep", "kind": "declared_panel", "feature_ids": ["F07", "F08"], "selection_provenance": "independent"},
            {"id": "deps", "kind": "dep_derived", "model_id": "limma-main", "contrast_id": "B-A", "criterion": "family_q", "threshold": 0.05}]
    payload, code, out = run(tmp_path, "location", permanova={"feature_sets": sets})
    assert code == 0, payload
    zero = [r for r in B.read_tsv(out / "dea" / "zero_null.tsv") if r["contrast_id"] == "B-A" and r["eligibility"] == "tested" and float(r["q_value"]) <= 0.05]
    usage = [r for r in B.read_tsv(out / "permanova" / "feature_usage.tsv") if r["feature_set_id"] == "deps" and r["test"] == "global:all_groups"]
    assert sorted(r["feature_id"] for r in usage) == sorted(r["feature_id"] for r in zero) and zero
    context = {r["feature_set_id"]: r for r in B.read_tsv(out / "permanova" / "selection_context.tsv")}
    assert set(context) == {"same", "deps"} and all(r["in_sample_optimistic"] == "true" for r in context.values())
    assert int(context["same"]["null_draws"]) == 300
    random = [r for r in B.read_tsv(out / "permanova" / "random_set_null.tsv") if r["feature_set_id"] == "same"]
    assert len(random) == 300
    r2 = {r["feature_id"]: float(r["univariate_group_r2"]) for r in B.read_tsv(out / "permanova" / "feature_r2.tsv") if r["feature_set_id"] == "all"}
    assert abs(float(context["same"]["best_possible_r2"]) - sum(sorted(r2.values(), reverse=True)[:3]) / 3) <= 1e-10
    report = (out / "report" / "index.html").read_text()
    assert "in-sample R2 is optimistic" in report


def test_v121_negative_unknown_panel_member_is_refused_at_plan(tmp_path):
    sets = [{"id": "bad", "kind": "declared_panel", "feature_ids": ["F01", "NOT_A_FEATURE"], "selection_provenance": "unknown"}]
    payload, code, out = run(tmp_path, "location", permanova={"feature_sets": sets})
    assert code == 2 and payload["error"]["code"] == "E_PERMANOVA_FEATURE_SET" and not (out / "dea").exists()


# ----------------------------------------------------------------------------- V122 / V125
def test_v122_negative_production_resolution(tmp_path):
    payload, code, out = run(tmp_path, "location", permanova={"permutations": 199}, profile="production")
    assert code == 2 and payload["error"]["code"] == "E_PERMANOVA_RESOLUTION" and not (out / "dea").exists()


def test_v125_empty_interaction_cell_refused_typed(tmp_path):
    sex = lambda o: {"sex": "M" if o["group"] == "B" else ("F" if int(o["observation_id"][1:]) % 2 else "M")}
    perm = {"covariates": [{"column": "sex", "type": "categorical"}], "interaction": True}
    payload, code, out = run(tmp_path / "req", "location", permanova=perm, required="required", extra=sex)
    assert code == 2 and payload["error"]["code"] == "E_PERMANOVA_INTERACTION_NONESTIMABLE" and not (out / "dea").exists()
    payload, code, out = run(tmp_path / "opt", "location", permanova=perm, required="optional", extra=sex)
    stage = next(s for s in payload["stages"] if s["stage_id"] == "permanova")
    assert stage["state"] == "INAPPLICABLE" and stage["reason_code"] == "E_PERMANOVA_INTERACTION_NONESTIMABLE" and code == 0
    assert not (out / "permanova").exists()


# ----------------------------------------------------------------------------- V129
class _Links(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(); self.urls = []
    def handle_starttag(self, tag, attrs):
        for k, v in attrs:
            if k in ("href", "src"):
                self.urls.append(v)


def test_v129_outputs_figures_and_report_section(tmp_path):
    payload, code, out = run(tmp_path, "location", groups=("A", "B", "C"))
    assert code == 0
    result = json.loads((out / "permanova" / "stage-result.json").read_text())
    outputs = {o["relative_path"]: o for o in result["outputs"]}
    for name in ("tests.tsv", "feature_r2.tsv", "permutation_null.tsv", "ordination.tsv", "centroids.tsv", "identity_check.json", "permanova_result.json"):
        assert name in outputs
    pngs = [p for p in outputs if p.endswith(".png")]
    assert {p.split("/")[-1].split("_")[0] for p in pngs} >= {"ordination", "permdisp", "null", "r2", "feature"}
    for p in pngs:
        assert (out / "permanova" / p).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    data = json.loads((out / "report" / "report_data.json").read_text())
    section = data["sections"]["permanova"]
    assert section["state"] == "COMPLETED" and any("non-significant PERMDISP indicates a location" in s for s in section["summary"])
    for figure in section["figures"]:
        assert figure["source"] and (out / figure["source"]).is_file() and (out / figure["src"]).is_file()
    links = _Links(); page = (out / "report" / "index.html").read_text(); links.feed(page)
    for url in links.urls:
        assert not url.startswith(("http:", "https:", "//")), url
        if not url.startswith("#"):
            assert (out / "report" / url).resolve().is_file(), url
    rows = _tests_table(out)
    nulls = B.read_tsv(out / "permanova" / "permutation_null.tsv")
    for row in rows:
        observed = [n for n in nulls if n["feature_set_id"] == row["feature_set_id"] and n["test_id"].startswith(row["analysis"])]
        assert observed, row
        assert all(r["permdisp_p"] != "NA" for r in rows if r["analysis"] in ("global", "pairwise"))


def test_v129_negative_failed_permanova_is_not_shown_as_completed(tmp_path):
    sets = [{"id": "tiny", "kind": "declared_panel", "feature_ids": ["F01", "F02"], "selection_provenance": "independent"}]
    def extra(o):
        return {}
    values, obs = F.matrix("location")
    values["F02"] = [5.0] * len(obs)                          # constant -> fewer than two usable features
    files = B.dataset(tmp_path / "data", values, obs)
    def mutate(c):
        c["runtime"]["phase"] = 2
        c["multivariate"] = {"enabled": True, "execution_requirement": "optional", "permanova": {"permutations": 199, "feature_sets": sets}}
    path = B.config(tmp_path / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], mutate=mutate)
    payload, code = workflow.run_command(path, tmp_path / "run")
    stage = next(s for s in payload["stages"] if s["stage_id"] == "permanova")
    assert stage["state"] == "FAILED" and stage["reason_code"] == "E_PERMANOVA_FEATURE_SET" and payload["state"] == "PARTIAL" and code != 0
    data = json.loads((tmp_path / "run" / "report" / "report_data.json").read_text())
    assert data["sections"]["permanova"]["state"] == "FAILED" and data["sections"]["permanova"]["values"] is None
