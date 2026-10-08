"""R14a acceptance through the real Python -> R pipeline: V131-V136 (SM31).

Oracles are written here: Python set operations on the completed zero_null.tsv, brute-force region
enumeration, R's cor()/lm()/phyper() called directly, and the hand-labelled planted features of the
synthetic fixture (tests/fixtures/post_de/design_factory.py).
"""
from __future__ import annotations

import copy
import importlib.util
import itertools
import json
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.config import load_config
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


F = _load("post_de_factory", ROOT / "tests" / "fixtures" / "post_de" / "design_factory.py")
B = F.B
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")

CONTRASTS = [B.contrast("B-A", "B", "A"), B.contrast("C-A", "C", "A", role="secondary"), B.contrast("D-C", "D", "C", role="secondary")]
FAMILIES = [F.family("protein-primary", ["B-A"], role="primary"), F.family("protein-secondary", ["C-A", "D-C"]),
            F.family("treat-all", ["B-A", "C-A", "D-C"], hypothesis="protein_treat")]


def sets_block(venn=("B_any", "C_any", "DC")):
    leaf = F.leaf
    definitions = [
        {"id": "B_any", "rule": leaf("protein-primary", "B-A")},
        {"id": "C_any", "rule": leaf("protein-secondary", "C-A")},
        {"id": "DC", "rule": leaf("protein-secondary", "D-C")},
        {"id": "B_up", "rule": leaf("protein-primary", "B-A", direction="up")},
        {"id": "B_or_C", "rule": {"op": "union", "args": [leaf("protein-primary", "B-A"), leaf("protein-secondary", "C-A")]}},
        {"id": "B_and_C", "rule": {"op": "intersect", "args": [leaf("protein-primary", "B-A"), leaf("protein-secondary", "C-A")]}},
        {"id": "B_not_C", "rule": {"op": "difference", "args": [leaf("protein-primary", "B-A"), leaf("protein-secondary", "C-A")]}},
        {"id": "not_B", "rule": {"op": "complement_within_tested", "args": [leaf("protein-primary", "B-A")]}},
        {"id": "C_raw", "rule": leaf("protein-secondary", "C-A", criterion="exploratory_raw_p", threshold=0.01)},
        {"id": "B_treat", "rule": leaf("treat-all", "B-A", criterion="treat")},
    ]
    return {"enabled": True, "sets": {"enabled": True, "definitions": definitions, "venn": list(venn), "upset": ["B_any", "C_any", "DC", "B_up", "C_raw"],
                                      "concordance_pairs": [{"id": "BA_vs_CA", "left": {"family_id": "protein-primary", "contrast_id": "B-A"},
                                                             "right": {"family_id": "protein-secondary", "contrast_id": "C-A"}}],
                                      "overlap_test": True}}


def make_config(tmp_path, post_de):
    values, obs, shared = F.four_group_sets()
    files = F.write_dataset(tmp_path / "data", values, obs, shared_genes=shared)
    return F.config(tmp_path / "data", files, groups=["A", "B", "C", "D"], contrasts=CONTRASTS, families=FAMILIES, post_de=post_de,
                    model_overrides={"additional_hypotheses": ["treat"], "effect_threshold": 0.5})


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("sets")
    payload, code = workflow.run_command(make_config(tmp, sets_block()), tmp / "run")
    assert code == 0, payload
    return tmp / "run", payload


def _endpoint(out, contrast, hypothesis="protein_zero_null"):
    name = "zero_null.tsv" if hypothesis == "protein_zero_null" else "treat.tsv"
    return [r for r in F.read_tsv(out / "dea" / name) if r["contrast_id"] == contrast and r["model_id"] == "limma-main"]


def _members(rows, criterion="q_value", threshold=0.05, direction="any"):
    out = set()
    for r in rows:
        if r["eligibility"] != "tested" or r[criterion] == "NA" or not float(r[criterion]) < threshold:
            continue
        if direction == "up" and not float(r["effect"]) > 0: continue
        if direction == "down" and not float(r["effect"]) < 0: continue
        out.add(r["feature_id"])
    return out


# ----------------------------------------------------------------------------- V131
def test_v131_set_algebra_equals_python_oracle(run):
    out, payload = run
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_sets")
    assert stage["state"] == "COMPLETED"
    ba, ca, dc = _endpoint(out, "B-A"), _endpoint(out, "C-A"), _endpoint(out, "D-C")
    b, c = _members(ba), _members(ca)
    tested_b = {r["feature_id"] for r in ba if r["eligibility"] == "tested"}
    oracle = {"B_any": b, "C_any": c, "DC": _members(dc), "B_up": _members(ba, direction="up"), "B_or_C": b | c, "B_and_C": b & c, "B_not_C": b - c,
              "not_B": tested_b - b, "C_raw": _members(ca, "p_value", 0.01), "B_treat": _members(_endpoint(out, "B-A", "protein_treat"))}
    membership = F.read_tsv(out / "post_de" / "sets" / "membership.tsv")
    tested_any = {r["feature_id"] for r in ba + ca + dc if r["eligibility"] == "tested"}
    assert sorted(r["feature_id"] for r in membership) == sorted(tested_any)          # one row per tested feature
    for set_id, expected in oracle.items():
        assert {r["feature_id"] for r in membership if r[set_id] == "true"} == expected, set_id
    rows = {r["set_id"]: r for r in F.read_tsv(out / "post_de" / "sets" / "sets.tsv")}
    assert rows["C_raw"]["claim_label"] == "exploratory_raw_p" and rows["B_any"]["claim_label"] == "descriptive"
    assert set(F.SETS_PLANTED["concordant_up"] + F.SETS_PLANTED["only_b"]) <= oracle["B_any"]   # the planted signal is found


def test_v131_negative_free_form_rule_and_unknown_operator_fail_grammar(tmp_path):
    """D-59: each malformed definition is refused on its own (E_SETRULE_GRAMMAR, never evaluated); the valid sets are still
    built; a required sets module with a malformed definition rejects the plan."""
    bads = ("B-A & C-A", {"op": "xor", "args": [F.leaf("protein-primary", "B-A"), F.leaf("protein-secondary", "C-A")]},
            {"op": "union", "args": [F.leaf("protein-primary", "B-A")], "expr": "eval"})
    block = sets_block(); block["sets"]["definitions"] = [{"id": f"bad{i}", "rule": bad} for i, bad in enumerate(bads)] + [{"id": "ok", "rule": F.leaf("protein-primary", "B-A")}]
    block["sets"]["venn"] = []
    payload, code = workflow.run_command(make_config(tmp_path / "opt", block), tmp_path / "opt" / "run")
    assert code == 0, payload
    refusals = {r["item"]: r["reason_code"] for r in F.read_tsv(tmp_path / "opt" / "run" / "post_de" / "sets" / "refusals.tsv")}
    assert {refusals[f"bad{i}"] for i in range(3)} == {"E_SETRULE_GRAMMAR"}
    assert {r["set_id"] for r in F.read_tsv(tmp_path / "opt" / "run" / "post_de" / "sets" / "sets.tsv")} == {"ok"}
    block["sets"]["execution_requirement"] = "required"
    with pytest.raises(ProteomicsError) as error:
        workflow.run_command(make_config(tmp_path / "req", block), tmp_path / "req" / "run")
    assert error.value.code == "E_SETRULE_GRAMMAR"


# ----------------------------------------------------------------------------- V132
def test_v132_provenance_on_every_set(run):
    out, _ = run
    rows = {r["set_id"]: r for r in F.read_tsv(out / "post_de" / "sets" / "sets.tsv")}
    for row in rows.values():
        for key in ("rule", "criteria", "thresholds", "family_ids", "null_type", "endpoints", "input_matrix"):
            assert row[key] not in ("", "NA"), (row["set_id"], key)
    assert rows["B_treat"]["null_type"] == "treat" and rows["B_any"]["null_type"] == "zero_null"
    assert rows["B_or_C"]["family_ids"] == "protein-primary;protein-secondary" and json.loads(rows["B_or_C"]["rule"])["op"] == "union"
    assert rows["C_raw"]["criteria"] == "exploratory_raw_p" and rows["C_raw"]["thresholds"] == "0.01"


def _loaded(tmp_path, block):
    return load_config(make_config(tmp_path, block))


def test_v132_negative_mixed_null_unknown_contrast_and_mixed_matrix(tmp_path):
    """Each invalid definition is refused with its own code (D-59: one refused definition never blocks the others)."""
    from proteomics_pipeline.post_de import sets
    base = _loaded(tmp_path / "a", sets_block())
    ok = {"id": "ok", "rule": F.leaf("protein-primary", "B-A")}
    def refused(definitions, config):
        out = sets.resolve({"definitions": definitions}, config)
        assert [d["id"] for d in out["definitions"]] == ["ok"]
        return {r["item"]: r["reason_code"] for r in out["refused"]}
    mixed = {"op": "union", "args": [F.leaf("protein-primary", "B-A"), F.leaf("treat-all", "C-A", criterion="treat")]}
    assert refused([{"id": "m", "rule": mixed}, ok], base) == {"m": "E_SETRULE_MIXED_NULL"}
    assert refused([{"id": "u", "rule": F.leaf("protein-secondary", "E-A")}, ok], base) == {"u": "E_SETRULE_UNKNOWN_CONTRAST"}
    sens = copy.deepcopy(base)
    model = copy.deepcopy(sens["models"][0]); model.update({"id": "limma-sens", "role": "sensitivity", "execution_requirement": "optional"})
    model.pop("additional_hypotheses", None)
    sens["models"].append(model)
    sens["preprocessing"]["sensitivities"] = [{"id": "mind", "method": "min_deterministic", "model_id": "limma-sens"}]
    sens["multiplicity_families"].append(F.family("protein-sensitivity", ["B-A", "C-A", "D-C"], models=("limma-sens",)))
    rule = {"op": "intersect", "args": [F.leaf("protein-primary", "B-A"), F.leaf("protein-sensitivity", "B-A")]}
    assert refused([{"id": "mm", "rule": rule}, ok], sens) == {"mm": "E_SETRULE_MIXED_MATRIX"}


# ----------------------------------------------------------------------------- V133
def test_v133_exact_regions_and_gene_counts(run):
    out, _ = run
    membership = F.read_tsv(out / "post_de" / "sets" / "membership.tsv")
    regions = F.read_tsv(out / "post_de" / "sets" / "regions.tsv")
    ids = ["B_any", "C_any", "DC", "B_up", "C_raw"]
    assert len(regions) == 2 ** 5 - 1
    members = {s: {r["feature_id"] for r in membership if r[s] == "true"} for s in ids}
    union = set().union(*members.values())
    for pattern in itertools.product([True, False], repeat=5):          # brute-force oracle over the 2^k - 1 regions
        if not any(pattern):
            continue
        inside = [s for s, p in zip(ids, pattern) if p]
        expected = {f for f in union if all((f in members[s]) == p for s, p in zip(ids, pattern))}
        row = next(r for r in regions if r["defining_sets"] == ";".join(inside))
        assert int(row["n_protein_groups"]) == len(expected) and set(filter(None, row["members"].split(";"))) == expected
    assert sum(int(r["n_protein_groups"]) for r in regions) == len(union)
    genes = {r["feature_id"]: r["gene_ids"] for r in membership}
    assert genes["F01"] == genes["F02"]                                    # two protein groups share one gene
    both = next(r for r in regions if "F01" in r["members"].split(";"))
    assert "F02" in both["members"].split(";") and int(both["n_unique_genes"]) == int(both["n_protein_groups"]) - 1   # counts differ where genes are shared
    eligibility = json.loads((out / "post_de" / "sets" / "eligibility.json").read_text(encoding="utf-8"))
    files = {f["file"] for f in eligibility["figures"]}
    assert "figures/upset_all.png" in files and "figures/venn_B_any_C_any_DC.png" in files
    assert (out / "post_de" / "sets" / "figure_sources" / "venn_B_any_C_any_DC.tsv").is_file()


def test_v133_negative_venn_for_four_sets_is_refused_and_upset_still_produced(tmp_path):
    payload, code = workflow.run_command(make_config(tmp_path, sets_block(venn=("B_any", "C_any", "DC", "B_up"))), tmp_path / "run")
    assert code == 0, payload
    out = tmp_path / "run" / "post_de" / "sets"
    refusals = F.read_tsv(out / "refusals.tsv")
    assert [r["reason_code"] for r in refusals if r["analysis"] == "venn"] == ["E_VENN_K"]
    assert not list((out / "figures").glob("venn_*")) and (out / "figures" / "upset_all.png").is_file()
    plan = json.loads((tmp_path / "run" / "plan.json").read_text(encoding="utf-8"))
    entry = next(c for c in plan["capability_plan"] if c["capability"] == "post_de_sets")
    assert any(s["reason_code"] == "E_VENN_K" for s in entry["post_de"]["subanalyses"])       # decided by the planner before computation


# ----------------------------------------------------------------------------- V134
def test_v134_direction_aware_overlap_equals_hand_labels(run):
    out, _ = run
    row = next(r for r in F.read_tsv(out / "post_de" / "sets" / "direction_overlap.tsv") if (r["set_a"], r["set_b"]) == ("B_any", "C_any"))
    split = lambda key: sorted(filter(None, row[key].split(";")))
    assert split("concordant_up") == F.SETS_PLANTED["concordant_up"] and split("concordant_down") == F.SETS_PLANTED["concordant_down"]
    assert split("discordant") == F.SETS_PLANTED["discordant"]                                  # never counted as overlap
    assert split("only_a") == F.SETS_PLANTED["only_b"] and split("only_b") == F.SETS_PLANTED["only_c"]
    assert (int(row["n_concordant_up"]), int(row["n_concordant_down"]), int(row["n_discordant"])) == (2, 2, 2)


# ----------------------------------------------------------------------------- V135
def test_v135_concordance_matches_r_oracle_and_has_no_p(run):
    out, _ = run
    row = F.read_tsv(out / "post_de" / "sets" / "concordance.tsv")[0]
    assert not any("p_value" in k or k.startswith("p_") for k in row) and row["claim_label"] == "descriptive"
    oracle = B.r_json("""
      a <- commandArgs(TRUE); z <- read.delim(a[1], colClasses = 'character', encoding = 'UTF-8')
      x <- z[z$contrast_id == 'B-A' & z$eligibility == 'tested', ]; y <- z[z$contrast_id == 'C-A' & z$eligibility == 'tested', ]
      common <- intersect(x$feature_id, y$feature_id); ex <- as.numeric(x$effect[match(common, x$feature_id)]); ey <- as.numeric(y$effect[match(common, y$feature_id)])
      qx <- as.numeric(x$q_value[match(common, x$feature_id)]); qy <- as.numeric(y$q_value[match(common, y$feature_id)]); sig <- qx < 0.05 | qy < 0.05
      cat(jsonlite::toJSON(list(n = length(common), pearson = cor(ex, ey), spearman = cor(ex, ey, method = 'spearman'), slope = unname(coef(lm(ey ~ 0 + ex))[1]),
                                agree = sum(sign(ex[sig]) == sign(ey[sig])), nsig = sum(sig)), auto_unbox = TRUE, digits = NA))""", out / "dea" / "zero_null.tsv")
    assert int(row["n_common_tested"]) == oracle["n"]
    for key, column in (("pearson", "pearson"), ("spearman", "spearman"), ("slope", "origin_slope")):
        assert abs(float(row[column]) - oracle[key]) <= 1e-10, key
    assert int(row["n_sign_agree"]) == oracle["agree"] and int(row["n_significant_either"]) == oracle["nsig"]


# ----------------------------------------------------------------------------- V136
def test_v136_overlap_p_only_for_unit_disjoint_contrasts(run):
    out, _ = run
    eligibility = {(r["set_a"], r["set_b"]): r for r in F.read_tsv(out / "post_de" / "sets" / "overlap_eligibility.tsv")}
    shared = eligibility[("B_any", "C_any")]
    assert shared["reason"] == "shared_biological_units" and shared["overlap_p_eligible"] == "false" and int(shared["n_shared_units"]) == 5
    assert eligibility[("B_any", "DC")]["reason"] == "disjoint_biological_units"
    tests = {(r["set_a"], r["set_b"]): r for r in F.read_tsv(out / "post_de" / "sets" / "overlap_test.tsv")}
    assert ("B_any", "C_any") not in tests and ("B_any", "DC") in tests            # no P for the shared-control case
    row = tests[("B_any", "DC")]
    oracle = B.r_json("""
      a <- commandArgs(TRUE); z <- read.delim(a[1], colClasses = 'character', encoding = 'UTF-8')
      x <- z[z$contrast_id == 'B-A' & z$eligibility == 'tested', ]; y <- z[z$contrast_id == 'D-C' & z$eligibility == 'tested', ]
      u <- intersect(x$feature_id, y$feature_id); a1 <- intersect(x$feature_id[as.numeric(x$q_value) < 0.05], u); b1 <- intersect(y$feature_id[as.numeric(y$q_value) < 0.05], u)
      k <- length(intersect(a1, b1))
      cat(jsonlite::toJSON(list(N = length(u), k = k, p = phyper(k - 1, length(a1), length(u) - length(a1), length(b1), lower.tail = FALSE)), auto_unbox = TRUE, digits = NA))""",
                      out / "dea" / "zero_null.tsv")
    assert int(row["universe_n"]) == oracle["N"] and int(row["overlap"]) == oracle["k"] and abs(float(row["p_value"]) - oracle["p"]) <= 1e-12
    assert all(set(r) >= {"set_a", "set_b", "reason"} and "p_value" not in r for r in eligibility.values())


def test_r14a_report_section_and_eligibility_record(run):
    out, _ = run
    data = json.loads((out / "report" / "report_data.json").read_text(encoding="utf-8"))
    section = data["sections"]["post_de_sets"]
    assert section["state"] == "COMPLETED" and any("sets.tsv" in t["source"] for t in section["tables"])
    for figure in section["figures"]:
        assert (out / figure["source"]).is_file() and (out / figure["src"]).is_file()
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    entry = next(c for c in plan["capability_plan"] if c["capability"] == "post_de_sets")
    assert entry["scientific_eligibility"] == "eligible" and entry["post_de"]["state"] == "ELIGIBLE"


# ----------------------------------------------------------------------------- shared amendment A-2026-10-01-14
def test_post_de_is_phase_4_only_and_capabilities_are_registered(tmp_path):
    from proteomics_pipeline import runtime
    from proteomics_pipeline.errors import ConfigurationError
    path = make_config(tmp_path, sets_block())
    raw = json.loads(path.read_text(encoding="utf-8")); raw["runtime"]["phase"] = 3
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ProteomicsError) as error:
        workflow.run_command(path, tmp_path / "run")
    assert error.value.code == "E_PHASE_CAPABILITY" and not (tmp_path / "run").exists()
    with pytest.raises(ConfigurationError) as error:                      # the schema rule holds on its own as well
        load_config(path)
    assert error.value.code == "E_CONFIG_SCHEMA"
    ids = {c["id"] for c in runtime.capabilities()}
    assert {"post_de_sets", "post_de_sensitivity", "post_de_association", "post_de_biomarker", "post_de_networks", "post_de_eligibility"} <= ids
    assert {"post_de_sets", "post_de_sensitivity", "post_de_association", "post_de_biomarker", "post_de_networks"} <= runtime.R_BACKED_CAPABILITIES
