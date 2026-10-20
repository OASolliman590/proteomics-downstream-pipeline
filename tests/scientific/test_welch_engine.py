"""R15a acceptance through the real Python -> R pipeline: V168-V173 (SM42, SM43; ADR 0011).

Oracles are written here as direct stats::t.test(), p.adjust() and limma calls run on the stage's own
observed matrix, so no Welch statistic is compared with another Welch statistic. Every case has a
negative assertion next to its positive one. Synthetic fixtures only.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.cli import main as cli_main
from proteomics_pipeline.config import load_config
from proteomics_pipeline.errors import ConfigurationError
from proteomics_pipeline.reporting.welch_data import welch_report_data

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


F = _load("welch_factory", ROOT / "tests" / "fixtures" / "methods_extension" / "welch" / "factory.py")
B = F.B
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")

T_ORACLE = """
  a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
  o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]
  ix <- which(o$group == a[3]); iz <- which(o$group == a[4]); paired <- a[5] == 'paired'
  res <- lapply(rownames(y), function(f) {
    if (paired) {
      common <- intersect(o$subject_id[ix], o$subject_id[iz])
      x <- as.numeric(y[f, ix[match(common, o$subject_id[ix])]]); z <- as.numeric(y[f, iz[match(common, o$subject_id[iz])]])
      keep <- !is.na(x) & !is.na(z); tt <- t.test(x[keep], z[keep], paired = TRUE)
      effect <- unname(tt$estimate)
    } else {
      x <- as.numeric(y[f, ix]); z <- as.numeric(y[f, iz]); x <- x[!is.na(x)]; z <- z[!is.na(z)]
      tt <- t.test(x, z, var.equal = FALSE); effect <- unname(tt$estimate[1] - tt$estimate[2])
    }
    c(statistic = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value, effect = effect, lower = tt$conf.int[1], upper = tt$conf.int[2])
  })
  names(res) <- rownames(y)
  cat(jsonlite::toJSON(lapply(res, as.list), digits = NA, auto_unbox = TRUE))"""

BH_ORACLE = """
  a <- commandArgs(TRUE); p <- as.numeric(readLines(a[1]))
  cat(jsonlite::toJSON(p.adjust(p, method = 'BH'), digits = NA))"""

ADJUSTED_ORACLE = """
  a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
  o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]
  grp <- as.numeric(o$group == 'B'); sex <- o$sex
  out <- sapply(rownames(y), function(f) unname(coef(lm(as.numeric(y[f, ]) ~ grp + factor(sex)))['grp']))
  cat(jsonlite::toJSON(as.list(out), digits = NA, auto_unbox = TRUE))"""

LIMMA_ORACLE = """
  a <- commandArgs(TRUE); y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
  o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); o <- o[match(colnames(y), o$observation_id), ]
  design <- model.matrix(~ 0 + factor(o$group, levels = c('A', 'B'))); colnames(design) <- c('A', 'B')
  fit <- limma::lmFit(y, design); fit <- limma::contrasts.fit(fit, limma::makeContrasts(B - A, levels = design))
  eb <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  cat(jsonlite::toJSON(list(feature = rownames(y), effect = unname(eb$coefficients[, 1]), p = unname(eb$p.value[, 1])), digits = NA))"""


def _run(config_path: Path, out: Path):
    payload, code = workflow.run_command(config_path, out)
    return payload, code


def _rows(path: Path) -> list[dict]:
    return F.read_tsv(path)


def _stage_dir(out: Path) -> Path:
    """The Welch stage directory: dea/ for a Welch primary (A-2026-10-01-46), methods_extension/welch/ for a secondary."""
    return out / "dea" if (out / "dea" / "eligibility.json").is_file() else out / "methods_extension" / "welch"


def _welch_rows(out: Path, contrast: str | None = None) -> list[dict]:
    rows = _rows(_stage_dir(out) / "welch_results.tsv")
    return [r for r in rows if contrast is None or r["contrast_id"] == contrast]


def _stable_limma(out: Path) -> list[dict]:
    return [{k: v for k, v in r.items() if k not in ("run_id", "plan_hash")} for r in _rows(out / "dea" / "zero_null.tsv")]


def _two_group_config(tmp_path: Path, *, planted=None, paper=None, paired="auto", phase: int = 5):
    values, obs = F.two_group(planted=planted if planted is not None else {0: 2.0, 1: -2.0, 2: 1.2})
    contrasts = [F.contrast("B-A", "B", "A")]
    methods = F.welch_methods(paper=paper, paired_value=paired, contrasts=None)
    return F.study_config(tmp_path / "study", values, obs, groups=["A", "B"], contrasts=contrasts, primary="welch", methods=methods, phase=phase)


@pytest.fixture(scope="module")
def two_group(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("two_group")
    cfg = _two_group_config(tmp, paper={"enabled": True, "raw_p_cutoff": 0.05})
    payload, code = _run(cfg, tmp / "run")
    assert code == 0, payload
    return tmp / "run"


# ----------------------------------------------------------------------------- V168
def test_v168_welch_statistics_equal_direct_t_test_within_1e10(two_group):
    oracle = B.r_json(T_ORACLE, str(two_group / "preprocessing" / "primary" / "matrix.tsv"), str(two_group / "preprocessing" / "primary" / "observations.tsv"), "B", "A", "two_sample")
    rows = {r["feature_id"]: r for r in _welch_rows(two_group, "B-A")}
    assert len(rows) == 200 and set(oracle) == set(rows)
    for feature, values in oracle.items():
        r = rows[feature]
        assert r["eligibility"] == "tested" and r["statistic_type"] == "Welch t" and r["engine"] == "welch"
        for key, column in (("statistic", "statistic"), ("df", "df_inference"), ("p", "p_value"), ("effect", "effect"), ("lower", "ci_lower"), ("upper", "ci_upper")):
            assert abs(float(r[column]) - values[key]) <= 1e-10, (feature, key)


def test_v168_negative_pooled_variance_and_reversed_sign_differ(two_group):
    oracle = B.r_json(T_ORACLE.replace("var.equal = FALSE", "var.equal = TRUE"), str(two_group / "preprocessing" / "primary" / "matrix.tsv"),
                      str(two_group / "preprocessing" / "primary" / "observations.tsv"), "B", "A", "two_sample")
    rows = {r["feature_id"]: r for r in _welch_rows(two_group, "B-A")}
    # the pooled df is n - 2 = 13 for every feature; the Welch df is not, because the groups have unequal variances
    assert any(abs(float(rows[f]["df_inference"]) - oracle[f]["df"]) > 1e-6 for f in oracle)
    # the effect sign follows the declared numerator (B): the planted +2 in B is positive, the planted -2 is negative
    assert float(rows["F001"]["effect"]) > 0 and float(rows["F002"]["effect"]) < 0


# ----------------------------------------------------------------------------- V170
def test_v170_bh_inside_the_declared_family_and_paper_rule_labelled(two_group):
    rows = {r["feature_id"]: r for r in _welch_rows(two_group, "B-A")}
    raw = [float(rows[f]["p_value"]) for f in sorted(rows)]
    path = two_group / "_bh_input.txt"
    path.write_text("\n".join(repr(p) for p in raw), encoding="utf-8")
    oracle_q = B.r_json(BH_ORACLE, str(path))
    for f, q in zip(sorted(rows), oracle_q):
        assert abs(float(rows[f]["q_value"]) - q) <= 1e-12
    paper = _rows(_stage_dir(two_group) / "paper_compatible.tsv")
    selected = {r["feature_id"] for r in paper}
    assert selected == {f for f in rows if float(rows[f]["p_value"]) < 0.05}                  # the raw-P set, exactly
    assert paper and all(r["claim_label"] == "not_multiplicity_controlled" for r in paper)     # every row labelled
    assert all(abs(float(r["bh_q_value"]) - float(rows[r["feature_id"]]["q_value"])) <= 1e-12 for r in paper)   # BH q in the same table
    # negative: the q column is a multiplicity adjustment, not the raw P relabelled
    assert any(abs(float(rows[f]["q_value"]) - float(rows[f]["p_value"])) > 1e-6 for f in rows)
    assert all(r["family_id"] == "welch_zero_null" for r in rows.values())


# ----------------------------------------------------------------------------- V171
def _paired_config(tmp_path: Path, *, paired="auto", planted=None):
    values, obs = F.paired(planted=planted if planted is not None else {0: 1.5, 1: -1.5})
    contrasts = [F.contrast("B-A", "B", "A")]
    design = {"blocking": {"mode": "duplicate_correlation", "subject_column": "subject_id"}}
    methods = F.welch_methods(paired_value=paired)
    return F.study_config(tmp_path / "paired", values, obs, groups=["A", "B"], contrasts=contrasts, primary="welch", methods=methods,
                          design_overrides=design)


@pytest.fixture(scope="module")
def paired_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("paired")
    payload, code = _run(_paired_config(tmp), tmp / "run")
    assert code == 0, payload
    return tmp / "run"


def test_v171_paired_variant_equals_t_test_paired_with_n_minus_one_df(paired_run):
    oracle = B.r_json(T_ORACLE, str(paired_run / "preprocessing" / "primary" / "matrix.tsv"), str(paired_run / "preprocessing" / "primary" / "observations.tsv"),
                      "B", "A", "paired")
    rows = {r["feature_id"]: r for r in _welch_rows(paired_run, "B-A")}
    for feature, values in oracle.items():
        r = rows[feature]
        assert r["statistic_type"] == "paired t" and abs(float(r["df_inference"]) - 7.0) == 0.0
        for key, column in (("statistic", "statistic"), ("df", "df_inference"), ("p", "p_value"), ("effect", "effect"), ("lower", "ci_lower"), ("upper", "ci_upper")):
            assert abs(float(r[column]) - values[key]) <= 1e-10, (feature, key)
    decision = json.loads((_stage_dir(paired_run) / "eligibility.json").read_text(encoding="utf-8"))["decisions"][0]
    assert decision["variant"] == "paired" and decision["n_pairs"] == 8


def test_v171_negative_two_sample_form_on_declared_pairs_is_refused(tmp_path):
    payload, code = _run(_paired_config(tmp_path, paired=False), tmp_path / "run")
    assert code != 0 and payload["state"] == "FAILED"
    assert payload["error"]["code"] == "E_WELCH_PAIRED_UNIT"
    assert not (tmp_path / "run" / "methods_extension" / "welch" / "welch_results.tsv").exists()


# ----------------------------------------------------------------------------- V172
@pytest.fixture(scope="module")
def covariate_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("covariate")
    base = [0.1, -0.2, 0.3, -0.1, 0.2, -0.3]
    values = {}
    sex_x = ["F", "F", "F", "F", "M", "M"]; sex_y = ["M", "M", "M", "F", "F", "F"]
    observations = [{"observation_id": f"A{i + 1}", "group": "A", "sex": sex_x[i]} for i in range(6)] + [{"observation_id": f"B{i + 1}", "group": "B", "sex": sex_y[i]} for i in range(6)]
    for f in range(40):
        row = [10 + (0.5 if o["sex"] == "M" else 0.0) + base[i % 6] + (1.0 if o["group"] == "B" and f < 5 else 0.0) for i, o in enumerate(observations)]
        values[f"F{f + 1:03d}"] = row
    design = {"categorical_covariates": [{"column": "sex", "levels": ["F", "M"], "reference": "F"}]}
    cfg = F.study_config(tmp / "cov", values, observations, groups=["A", "B"], contrasts=[F.contrast("B-A", "B", "A")], primary="welch",
                         methods=F.welch_methods(), design_overrides=design, extra_columns=("sex",))
    payload, code = _run(cfg, tmp / "run")
    assert code == 0, payload
    return tmp / "run"


def test_v172_covariate_is_run_unadjusted_and_recorded_everywhere(covariate_run):
    oracle = B.r_json(T_ORACLE, str(covariate_run / "preprocessing" / "primary" / "matrix.tsv"), str(covariate_run / "preprocessing" / "primary" / "observations.tsv"),
                      "B", "A", "two_sample")
    rows = {r["feature_id"]: r for r in _welch_rows(covariate_run, "B-A")}
    for feature, values in oracle.items():
        assert abs(float(rows[feature]["statistic"]) - values["statistic"]) <= 1e-10 and abs(float(rows[feature]["effect"]) - values["effect"]) <= 1e-10
    expected = {"analysis": "welch", "item": "covariate:sex", "requested": "adjusted", "used": "unadjusted", "reason": "welch_not_covariate_adjusting"}
    eligibility = json.loads((_stage_dir(covariate_run) / "eligibility.json").read_text(encoding="utf-8"))
    assert [{k: a[k] for k in expected} for a in eligibility["adaptations"]] == [expected]
    plan = json.loads((covariate_run / "plan.json").read_text(encoding="utf-8"))
    assert [{k: a[k] for k in expected} for a in plan["engine_eligibility"][0]["adaptations"]] == [expected]
    warnings = json.loads((covariate_run / "warnings.json").read_text(encoding="utf-8"))
    assert any(w["code"] == "W_WELCH_ADAPTED" for w in warnings)
    # negative: the sex-adjusted effect differs from the unadjusted Welch effect on this unbalanced fixture (orientation B minus A)
    adjusted = B.r_json(ADJUSTED_ORACLE, str(covariate_run / "preprocessing" / "primary" / "matrix.tsv"), str(covariate_run / "preprocessing" / "primary" / "observations.tsv"))
    assert any(abs(float(rows[f]["effect"]) + adjusted[f]) > 1e-6 for f in rows)     # the sex-adjusted group effect differs; the record is not cosmetic


# ----------------------------------------------------------------------------- V173
@pytest.fixture(scope="module")
def same_fixture_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("same_fixture")
    values, obs = F.two_group(n_num=12, n_den=12, features=200, seed=173, planted={0: 1.5}, sd_num=0.4, sd_den=1.6)
    secondary = {"id": "welch-secondary", "design_id": "joint", "engine": "welch", "role": "sensitivity", "execution_requirement": "optional",
                 "hypothesis": "zero_null", "coverage": {"policy": "available_case", "minimum_observed_per_group": 2, "minimum_fraction": 0.5, "mask_source": "original_observed"},
                 "precision_weights": {"kind": "none"}}
    families = [{"id": "primary-family", "hypothesis_type": "protein_zero_null", "model_ids": ["limma-main"], "contrast_ids": ["B-A"], "adjustment": "BH",
                 "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "primary"},
                {"id": "welch_zero_null", "hypothesis_type": "protein_zero_null", "model_ids": ["welch-secondary"], "contrast_ids": ["B-A"], "adjustment": "BH",
                 "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"}]
    cfg = F.study_config(tmp / "same", values, obs, groups=["A", "B"], contrasts=[F.contrast("B-A", "B", "A")], primary="limma",
                         methods=F.welch_methods(role="secondary"), extra_models=[secondary], families=families, phase=5)
    payload, code = _run(cfg, tmp / "run")
    assert code == 0, payload
    return tmp / "run"


def test_v173_same_effects_distinct_p_values_each_equal_its_oracle(same_fixture_run):
    limma_rows = {r["feature_id"]: r for r in _rows(same_fixture_run / "dea" / "zero_null.tsv") if r["model_id"] == "limma-main" and r["contrast_id"] == "B-A"}
    welch_rows = {r["feature_id"]: r for r in _welch_rows(same_fixture_run, "B-A")}
    limma_oracle = B.r_json(LIMMA_ORACLE, str(same_fixture_run / "preprocessing" / "primary" / "matrix.tsv"),
                            str(same_fixture_run / "preprocessing" / "primary" / "observations.tsv"))
    welch_oracle = B.r_json(T_ORACLE, str(same_fixture_run / "preprocessing" / "primary" / "matrix.tsv"),
                            str(same_fixture_run / "preprocessing" / "primary" / "observations.tsv"), "B", "A", "two_sample")
    max_gap = 0.0
    for feature, effect, p in zip(limma_oracle["feature"], limma_oracle["effect"], limma_oracle["p"]):
        assert abs(float(limma_rows[feature]["effect"]) - effect) <= 1e-10          # limma engine output equals its own oracle
        assert abs(float(limma_rows[feature]["p_value"]) - p) <= 1e-10
        assert abs(float(welch_rows[feature]["effect"]) - effect) <= 1e-10          # identical effects from both engines
        assert abs(float(welch_rows[feature]["p_value"]) - welch_oracle[feature]["p"]) <= 1e-10
        max_gap = max(max_gap, abs(float(welch_rows[feature]["p_value"]) - float(limma_rows[feature]["p_value"])))
    assert max_gap > 1e-6                                                             # negative: identical P would mean pooled variances
    assert {r["engine"] for r in welch_rows.values()} == {"welch"} and {r["engine"] for r in limma_rows.values()} == {"limma"}
    assert {r["family_id"] for r in welch_rows.values()} == {"welch_zero_null"}
    assert {r["family_id"] for r in limma_rows.values()} == {"primary-family"}


# ----------------------------------------------------------------------------- V169
def test_v169_no_methods_block_and_disabled_block_give_identical_limma_outputs(tmp_path):
    values, obs = F.three_group()
    contrasts = [F.contrast("disease-control", "disease", "control")]
    plain = F.study_config(tmp_path / "plain", values, obs, groups=["control", "disease", "treated"], contrasts=contrasts, primary="limma", phase=1)
    disabled = F.study_config(tmp_path / "disabled", values, obs, groups=["control", "disease", "treated"], contrasts=contrasts, primary="limma", phase=1,
                              methods={"welch": {"enabled": False}})
    payload_a, code_a = _run(plain, tmp_path / "run_a")
    payload_b, code_b = _run(disabled, tmp_path / "run_b")
    assert code_a == 0 and code_b == 0, (payload_a, payload_b)
    assert not (tmp_path / "run_a" / "methods_extension").exists() and not (tmp_path / "run_b" / "methods_extension").exists()
    # run_id and plan_hash are run-specific (the plan embeds the configuration); every estimate, P and q must be identical
    assert _stable_limma(tmp_path / "run_a") == _stable_limma(tmp_path / "run_b")


def test_v169_welch_secondary_leaves_limma_outputs_unchanged(same_fixture_run, tmp_path):
    values, obs = F.two_group(n_num=12, n_den=12, features=200, seed=173, planted={0: 1.5}, sd_num=0.4, sd_den=1.6)
    plain = F.study_config(tmp_path / "plain", values, obs, groups=["A", "B"], contrasts=[F.contrast("B-A", "B", "A")], primary="limma", phase=5)
    payload, code = _run(plain, tmp_path / "run")
    assert code == 0, payload
    limma_only = [r for r in _rows(same_fixture_run / "dea" / "zero_null.tsv") if r["model_id"] == "limma-main"]
    baseline = _rows(tmp_path / "run" / "dea" / "zero_null.tsv")
    by_feature = {r["feature_id"]: r for r in baseline}
    assert len(limma_only) == len(baseline)
    for r in limma_only:
        assert r["p_value"] == by_feature[r["feature_id"]]["p_value"] and r["effect"] == by_feature[r["feature_id"]]["effect"]


def test_v169_declared_welch_primary_and_second_primary_is_refused(tmp_path):
    values, obs = F.two_group()
    cfg = F.study_config(tmp_path / "two", values, obs, groups=["A", "B"], contrasts=[F.contrast("B-A", "B", "A")], primary="welch", methods=F.welch_methods(), phase=5)
    data = json.loads(cfg.read_text(encoding="utf-8"))
    second = dict(data["models"][0], id="limma-second", engine="limma", role="primary")
    data["models"].append(second)
    cfg.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ConfigurationError) as info:   # the schema's single-primary item rule refuses the second primary (E_CONFIG_SCHEMA)
        load_config(cfg)
    assert "Too many items" in info.value.message


def test_v169_welch_requested_but_not_declared_is_refused_at_validation(tmp_path):
    values, obs = F.two_group()
    cfg = F.study_config(tmp_path / "undeclared", values, obs, groups=["A", "B"], contrasts=[F.contrast("B-A", "B", "A")], primary="welch", methods=None, phase=5)
    with pytest.raises(ConfigurationError) as info:
        load_config(cfg)
    assert info.value.code == "E_WELCH_NOT_DECLARED"


def test_v169_welch_primary_run_never_emits_a_limma_result_under_the_welch_name(two_group):
    dea = _rows(two_group / "dea" / "zero_null.tsv")            # the dea/ contract is written by the welch primary (A-2026-10-01-46)
    assert {r["engine"] for r in dea} == {"welch"} and {r["engine"] for r in _welch_rows(two_group)} == {"welch"}
    assert not (two_group / "dea" / "moderation_limma-main.tsv").exists()
    assert all(r["result_type"] == "WelchResult" for r in _welch_rows(two_group))
    assert {r["statistic_type"] for r in dea if r["eligibility"] == "tested"} == {"Welch t"}
    model_status = json.loads((two_group / "dea" / "model_status.json").read_text(encoding="utf-8"))
    assert model_status["settings"]["welch-main"]["engine"] == "welch" and model_status["settings"]["welch-main"]["input_imputation"] == "none"


# ----------------------------------------------------------------------------- end to end
def test_end_to_end_three_group_welch_primary_then_verify(tmp_path, capsys):
    values, obs = F.three_group(shifts={0: {"disease": 2.0}, 1: {"treated": -2.0}})
    contrasts = [F.contrast("disease-control", "disease", "control"), F.contrast("treated-control", "treated", "control"),
                 F.contrast("treated-disease", "treated", "disease")]
    cfg = F.study_config(tmp_path / "e2e", values, obs, groups=["control", "disease", "treated"], contrasts=contrasts, primary="welch",
                         methods=F.welch_methods(paper={"enabled": True, "raw_p_cutoff": 0.05}))
    payload, code = _run(cfg, tmp_path / "run")
    assert code == 0 and payload["state"] == "COMPLETED", payload
    assert {s["stage_id"]: s["state"] for s in payload["stages"]}["welch"] == "COMPLETED"
    capsys.readouterr()
    verify_code = cli_main(["verify", "--run", str(tmp_path / "run"), "--json"])
    verified = json.loads(capsys.readouterr().out)
    assert verify_code == 0 and verified["valid"] is True and verified["run_status"]["state"] == "COMPLETED"
    oracle = B.r_json(T_ORACLE, str(tmp_path / "run" / "preprocessing" / "primary" / "matrix.tsv"), str(tmp_path / "run" / "preprocessing" / "primary" / "observations.tsv"),
                      "disease", "control", "two_sample")
    rows = {(r["contrast_id"], r["feature_id"]): r for r in _welch_rows(tmp_path / "run")}
    for feature, values_ in oracle.items():
        assert abs(float(rows[("disease-control", feature)]["p_value"]) - values_["p"]) <= 1e-10
    assert float(rows[("disease-control", "G001")]["effect"]) > 1.5                    # planted +2 in disease
    assert len(_welch_rows(tmp_path / "run")) == 3 * 60
    report = welch_report_data(_stage_dir(tmp_path / "run"))       # the report reads only the stage's own outputs
    assert report["engine"] == "welch" and report["claim_labels"] == ["welch", "not_multiplicity_controlled"]
    assert report["tested_rows"] == 3 * 60 and report["endpoint_rows"] == 3 * 60
    assert report["paper_compatible"]["claim_label"] == "not_multiplicity_controlled" and report["paper_compatible"]["selected"] >= 1
    assert {f["family_id"] for f in report["families"]} == {"welch_zero_null"}


# ----------------------------------------------------------------------------- Welch as the primary engine feeds the whole run (F1)
DEA_COLUMNS = ["model_id", "contrast_id", "feature_id", "effect", "ci_lower", "ci_upper", "statistic", "statistic_type", "df_inference", "p_value",
               "q_value", "family_id", "eligibility", "engine", "hypothesis_type"]


def test_f1_welch_primary_feeds_dea_report_figures_and_verify(tmp_path, capsys):
    values, obs = F.three_group(shifts={0: {"disease": 2.0}, 1: {"treated": -2.0}})
    contrasts = [F.contrast("disease-control", "disease", "control"), F.contrast("treated-control", "treated", "control"),
                 F.contrast("treated-disease", "treated", "disease")]
    cfg = F.study_config(tmp_path / "f1", values, obs, groups=["control", "disease", "treated"], contrasts=contrasts, primary="welch",
                         methods=F.welch_methods(paper={"enabled": True, "raw_p_cutoff": 0.05}))
    data = json.loads(cfg.read_text(encoding="utf-8"))
    data["report"]["figures"] = {"formats": ["png", "pzfx", "html"], "dpi": 300, "journal_width_mm": 120, "palette": "okabe_ito", "font": "Arial",
                                 "p_annotation": "exact"}
    cfg.write_text(json.dumps(data, indent=2), encoding="utf-8")
    payload, code = _run(cfg, tmp_path / "run")
    assert code == 0 and payload["state"] == "COMPLETED", payload
    states = {s["stage_id"]: s["state"] for s in payload["stages"]}
    assert all(states.get(k) == "COMPLETED" for k in ("welch", "report", "report_full", "figure_catalogue")), states
    out = tmp_path / "run"
    # dea/: the same differential contract as limma, with the Welch engine, statistic label, family and eligibility
    dea = _rows(out / "dea" / "zero_null.tsv")
    assert len(dea) == 3 * 60 and set(DEA_COLUMNS) <= set(dea[0].keys())
    assert {r["engine"] for r in dea} == {"welch"} and {r["statistic_type"] for r in dea if r["eligibility"] == "tested"} == {"Welch t"}
    families = _rows(out / "dea" / "families.tsv")
    assert {f["family_id"] for f in families} == {"welch_zero_null"} and families[0]["adjustment"] == "BH"
    assert all(r["family_id"] == "welch_zero_null" and r["q_value"] != "NA" for r in dea if r["eligibility"] == "tested")
    assert float(next(r for r in dea if r["contrast_id"] == "disease-control" and r["feature_id"] == "G001")["effect"]) > 1.5
    paper = _rows(out / "dea" / "paper_compatible.tsv")
    assert paper and all(r["claim_label"] == "not_multiplicity_controlled" for r in paper)
    # report: the complete-results section lists the Welch table; the methods say Welch
    report = (out / "report-full" / "index.html").read_text(encoding="utf-8")
    assert "Complete zero_null table" in report and "180 rows" in report
    methods = (out / "methods.md").read_text(encoding="utf-8")
    section = methods.split("## Differential abundance")[1].split("\n## ")[0]
    assert "Welch t-test" in section and "lmFit" not in section and "eBayes" not in section
    # R16 figure catalogue: the differential figures are produced from the Welch dea/ tables
    catalogue = {p.name for p in (out / "report" / "figure_catalogue").iterdir()}
    for stem in ("diff_volcano_labelled__disease-control", "diff_ma_plot__disease-control", "diff_protein_dotplot_brackets__disease-control",
                 "diff_upset_venn__primary"):
        assert stem + ".png" in catalogue and stem + ".pzfx" in catalogue and stem + ".html" in catalogue, stem
    # proteomics verify exit 0
    capsys.readouterr()
    verify_code = cli_main(["verify", "--run", str(out), "--json"])
    verified = json.loads(capsys.readouterr().out)
    assert verify_code == 0 and verified["valid"] is True


def test_f1_post_de_association_consumes_the_welch_dea_families(tmp_path, capsys):
    import random
    values, obs = F.three_group(shifts={0: {"disease": 2.0}, 1: {"treated": -2.0}})
    rng = random.Random(7)
    for i, o in enumerate(obs):
        o["age"] = 40 + i * 0.7
        o["score"] = 1.0 + 0.2 * i + rng.gauss(0, 0.1)
    contrasts = [F.contrast("disease-control", "disease", "control"), F.contrast("treated-control", "treated", "control"),
                 F.contrast("treated-disease", "treated", "disease")]
    cfg = F.study_config(tmp_path / "pd", values, obs, groups=["control", "disease", "treated"], contrasts=contrasts, primary="welch",
                         methods=F.welch_methods(), extra_columns=("age", "score"))
    data = json.loads(cfg.read_text(encoding="utf-8"))
    phenotype = {"column": "score", "type": "numeric", "adjust_for": [{"column": "age", "type": "continuous"}], "adjust_for_group": True, "scope": "pooled"}
    data["post_de"] = {"enabled": True, "execution_requirement": "optional", "sets": {"enabled": False}, "sensitivity": {"enabled": False},
                       "biomarker": {"enabled": False}, "networks": {"enabled": False},
                       "association": {"enabled": True, "execution_requirement": "optional", "phenotypes": [phenotype], "method": "model",
                                       "correlation": "pearson", "permutations": 99, "seed": 777}}
    cfg.write_text(json.dumps(data, indent=2), encoding="utf-8")
    payload, code = _run(cfg, tmp_path / "run")
    assert code == 0, payload
    states = {s["stage_id"]: s["state"] for s in payload["stages"]}
    assert states.get("post_de_association") == "COMPLETED" and states.get("post_de_eligibility") == "COMPLETED", states
    families = _rows(tmp_path / "run" / "dea" / "families.tsv")
    assert {f["family_id"] for f in families} == {"welch_zero_null"}                    # the association reads the Welch family definitions
    capsys.readouterr()
    assert cli_main(["verify", "--run", str(tmp_path / "run"), "--json"]) == 0


# ----------------------------------------------------------------------------- F4: the spec's "secondary" is the schema's role "sensitivity"
def _secondary_welch_config(tmp_path: Path, *, sensitivity_model: str | None = None):
    import copy as _copy
    values, obs = F.two_group(n_num=12, n_den=12, features=60, seed=173, planted={0: 1.5})
    secondary = {"id": "welch-secondary", "design_id": "joint", "engine": "welch", "role": "sensitivity", "execution_requirement": "optional",
                 "hypothesis": "zero_null", "coverage": {"policy": "available_case", "minimum_observed_per_group": 2, "minimum_fraction": 0.5, "mask_source": "original_observed"},
                 "precision_weights": {"kind": "none"}}
    families = [{"id": "primary-family", "hypothesis_type": "protein_zero_null", "model_ids": ["limma-main"], "contrast_ids": ["B-A"], "adjustment": "BH",
                 "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "primary"},
                {"id": "welch_zero_null", "hypothesis_type": "protein_zero_null", "model_ids": ["welch-secondary"], "contrast_ids": ["B-A"], "adjustment": "BH",
                 "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"}]
    cfg = F.study_config(tmp_path / "sec", values, obs, groups=["A", "B"], contrasts=[F.contrast("B-A", "B", "A")], primary="limma",
                         methods=F.welch_methods(role="secondary"), extra_models=[secondary], families=families, phase=5)
    if sensitivity_model is not None:
        data = json.loads(cfg.read_text(encoding="utf-8"))
        data.setdefault("preprocessing", {})["sensitivities"] = [{"id": "min_s", "method": "min_deterministic", "model_id": sensitivity_model}]
        cfg.write_text(json.dumps(data), encoding="utf-8")
    return cfg


def test_f4_welch_secondary_with_role_sensitivity_is_not_a_sensitivity_matrix_and_is_accepted(tmp_path):
    cfg = _secondary_welch_config(tmp_path)
    loaded = load_config(cfg)                                               # no E_SENSITIVITY_PRIMARY_MODEL and no refusal
    roles = {m["id"]: m["role"] for m in loaded["models"]}
    assert roles == {"limma-main": "primary", "welch-secondary": "sensitivity"}


def test_f4_negative_a_sensitivity_matrix_on_the_welch_model_is_refused(tmp_path):
    cfg = _secondary_welch_config(tmp_path, sensitivity_model="welch-secondary")
    with pytest.raises(ConfigurationError) as info:
        load_config(cfg)
    assert info.value.code == "E_WELCH_MATRIX_UNSUPPORTED"


def test_f4_negative_a_sensitivity_matrix_on_the_primary_model_is_still_refused(tmp_path):
    cfg = _secondary_welch_config(tmp_path, sensitivity_model="limma-main")
    with pytest.raises(ConfigurationError) as info:
        load_config(cfg)
    assert info.value.code == "E_SENSITIVITY_PRIMARY_MODEL"
