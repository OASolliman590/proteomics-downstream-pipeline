"""R03 integration through the real Python -> R stage protocol (V021, V022, V026-V028, V030)."""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from proteomics_pipeline import preprocessing_service
from proteomics_pipeline.config import load_config
from proteomics_pipeline.intake import service as intake
from proteomics_pipeline.provenance import canonical_json_sha256
from proteomics_pipeline.reporting.qc_data import qc_section

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "configs" / "examples"
FIX = ROOT / "tests" / "fixtures" / "preprocessing"


def _require_r():
    rscript = preprocessing_service.rscript_executable()
    if not shutil.which(rscript):
        pytest.skip("NOT_RUN: Rscript unavailable")
    probe = subprocess.run([rscript, "--vanilla", "-e", "quit(status=if(exists('preprocess_stage', envir=asNamespace('proteomicsCore'))) 0L else 1L)"], capture_output=True)
    if probe.returncode != 0:
        pytest.skip("NOT_RUN: installed proteomicsCore lacks the R03 handler")


def _config(tmp_path, mutate=None, name="example-independent.json"):
    raw = json.loads((EXAMPLES / name).read_text(encoding="utf-8"))
    for key in ("matrix", "observations", "features", "source_provenance"):
        raw["input"][key] = str((EXAMPLES / raw["input"][key]).resolve())
    if mutate:
        mutate(raw)
    path = tmp_path / "analysis.json"; path.write_text(json.dumps(raw), encoding="utf-8")
    return load_config(path), path


def _run(tmp_path, mutate=None, name="example-independent.json"):
    _require_r()
    config, path = _config(tmp_path, mutate, name)
    intake.run_intake(config, path.parent, tmp_path / "inputs")
    request = preprocessing_service.build_request(config, config_path=path, bundle_dir=tmp_path / "inputs", run_id="itest", output_temp_dir=tmp_path / "pre")
    result = preprocessing_service.execute(request)
    return config, result, tmp_path / "pre"


def _tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_v021_preserve_leaves_values_and_masks_unchanged(tmp_path):
    _, result, out = _run(tmp_path)
    assert result["state"] == "COMPLETED", result["message"]
    canonical = _tsv(tmp_path / "inputs" / "matrix.tsv"); primary = _tsv(out / "primary" / "matrix.tsv")
    assert [r[0] for r in canonical] == [r[0] for r in primary]
    for a, b in zip(canonical[1:], primary[1:]):
        assert [None if t == "NA" else float(t) for t in a[1:]] == [None if t == "NA" else float(t) for t in b[1:]]
    assert _tsv(out / "primary" / "observed_mask.tsv") == _tsv(tmp_path / "inputs" / "observed_mask.tsv")
    lineage = json.loads((out / "primary" / "lineage.json").read_text(encoding="utf-8"))
    assert lineage["source_sha256"] == _sha(tmp_path / "inputs" / "matrix.tsv") and lineage["display_matrix_used"] is False
    assert lineage["operations"][1]["values_changed"] is False and lineage["new_primary_imputation"] == "none"
    # display-only PCA input is a different artifact from the primary matrix
    assert _sha(out / "qc" / "pca_display_input.tsv") != _sha(out / "primary" / "matrix.tsv")
    state = json.loads((out / "qc" / "pca_state.json").read_text(encoding="utf-8"))
    assert state["primary_matrix_sha256"] == _sha(out / "primary" / "matrix.tsv")


def test_v021_v022_requested_median_is_a_distinct_artifact(tmp_path):
    _, result, out = _run(tmp_path, lambda c: c["preprocessing"].update({"normalization": "median", "normalization_sensitivities": ["quantile"]}))
    assert result["state"] == "COMPLETED", result["message"]
    lineage = json.loads((out / "primary" / "lineage.json").read_text(encoding="utf-8"))
    assert lineage["operations"][1]["values_changed"] is True and lineage["output_sha256"] != lineage["source_sha256"]
    factors = {r[0]: float(r[3]) for r in _tsv(out / "normalization_factors.tsv")[1:]}
    canonical = _tsv(tmp_path / "inputs" / "matrix.tsv")
    medians = {}
    for j, obs in enumerate(canonical[0][1:], start=1):
        vals = sorted(float(r[j]) for r in canonical[1:] if r[j] != "NA")
        n = len(vals); medians[obs] = vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2
    centre = sorted(medians.values()); centre = (centre[5] + centre[6]) / 2
    for obs, med in medians.items():
        assert abs(factors[obs] - (med - centre)) <= 1e-10
    q = json.loads((out / "sensitivity" / "quantile" / "parameters.json").read_text(encoding="utf-8"))
    assert q["role"] == "named_normalization_sensitivity_only" and sorted(q["lost_features"]) == ["P07", "P08"]


def test_v021_negative_hidden_normalization_on_linear_scale_is_refused(tmp_path):
    def mutate(c):
        c["source_scale"] = "standardized_unknown"; c["preprocessing"]["normalization"] = "median"; c["runtime"]["scope"] = "qc_only"
    _, result, _ = _run(tmp_path, mutate)
    assert result["state"] == "FAILED" and result["reason_code"] == "E_NORMALIZATION_SCALE" and result["exit_code"] == 2


def test_v027_extreme_sample_flagged_not_dropped_and_exclusion_changes_plan_hash(tmp_path):
    reason = "prespecified: synthetic sample withdrawn before analysis"
    def mutate(c):
        c["input"]["matrix"] = str(FIX / "extreme-sample-abundance.tsv")
        c["preprocessing"]["exclusions"] = [{"observation_id": "T4", "reason": reason}]
    _, result, out = _run(tmp_path, mutate)
    assert result["state"] == "COMPLETED", result["message"]
    watch = {r[0]: r for r in _tsv(out / "qc" / "watchlist.tsv")[1:]}
    assert watch["C4"][4] == "true" and watch["C4"][6] == "none_flag_only"
    retained = _tsv(out / "primary" / "matrix.tsv")[0][1:]
    assert "C4" in retained and "T4" not in retained and len(retained) == 11
    fragment = json.loads((out / "exclusions" / "fragment.json").read_text(encoding="utf-8"))
    assert fragment["exclusions"] == [{"observation_id": "T4", "reason": reason}]
    envelope = {"schema_version": "1.2.0", "kind": "synthetic_plan_envelope", "inputs": {"matrix": "fixed"}, "preprocessing_fragment": fragment}
    changed = copy.deepcopy(envelope); changed["preprocessing_fragment"]["exclusions"][0]["reason"] = "different reason"
    assert canonical_json_sha256(envelope) != canonical_json_sha256(changed)
    assert canonical_json_sha256(envelope) == canonical_json_sha256(copy.deepcopy(envelope))


def test_v027_negative_blank_exclusion_reason(tmp_path):
    _, result, _ = _run(tmp_path, lambda c: c["preprocessing"].update({"exclusions": [{"observation_id": "T4", "reason": "   "}]}))
    assert result["state"] == "FAILED" and result["reason_code"] == "E_EXCLUSION_REASON_REQUIRED"


def test_v028_stage_sensitivities_leave_primary_unchanged(tmp_path):
    # Audit 2026-10-02: sensitivities may not target the primary model, so they target a declared role=sensitivity model.
    sens = [{"id": "mindet", "method": "min_deterministic", "model_id": "limma-sens"},
            {"id": "gauss", "method": "left_shifted_gaussian", "model_id": "limma-sens", "shift_sd": 1.8, "scale_sd": 0.3}]
    def mutate(c):
        model = copy.deepcopy(c["models"][0]); model.update({"id": "limma-sens", "role": "sensitivity", "execution_requirement": "optional"})
        c["models"].append(model); c["preprocessing"].update({"sensitivities": sens})
    _, result, out = _run(tmp_path, mutate)
    assert result["state"] == "COMPLETED", result["message"]
    canonical = _tsv(tmp_path / "inputs" / "matrix.tsv")
    assert _tsv(out / "primary" / "matrix.tsv")[1:] == [[r[0]] + r[1:] for r in canonical[1:]]
    mindet = {r[0]: r for r in _tsv(out / "sensitivity" / "mindet" / "matrix.tsv")}
    p07 = [float(x) for x in canonical[7][1:] if x != "NA"]
    header = canonical[0]
    assert float(mindet["P07"][header.index("U2")]) == min(p07)
    params = json.loads((out / "sensitivity" / "gauss" / "parameters.json").read_text(encoding="utf-8"))
    assert params["seed"] == 4812026 and params["draw_order"] == "feature_id then observation_id"


def test_v030_qc_only_scope_and_zero_variable_features(tmp_path):
    def mutate(c):
        c["runtime"]["scope"] = "qc_only"; c["input"]["matrix"] = str(FIX / "constant-abundance.tsv")
    _, result, out = _run(tmp_path, mutate)
    assert result["state"] == "COMPLETED", result["message"]
    summary = json.loads((out / "preprocessing_result.json").read_text(encoding="utf-8"))
    assert summary["model_fit_performed"] is False and summary["coverage"]["state"] == "NOT_REQUESTED"
    assert summary["pca"] == {"state": "INAPPLICABLE", "reason_code": "E_QC_CONSTANT"}
    assert not (out / "qc" / "pca_scores.tsv").exists() and not (out / "qc" / "pca_variance.tsv").exists()
    section = qc_section(out)
    assert section["tables"]["pca_variance"]["state"] == "NOT_AVAILABLE"
    for name, table in section["tables"].items():
        if table["source"]:
            assert (out / table["source"].removeprefix("preprocessing/")).is_file()
    sample_n = section["tables"]["sample_n"]["rows"]
    assert [int(r["n_biological_units"]) for r in sample_n] == [4, 4, 4]


def test_v030_ordinary_example_every_qc_value_has_a_source(tmp_path):
    _, result, out = _run(tmp_path)
    section = qc_section(out)
    variance = section["tables"]["pca_variance"]["rows"]
    assert abs(sum(float(r["variance_explained"]) for r in variance) - 1) <= 1e-10
    n_by_obs = {row["observation_id"]: int(row["n"]) for row in section["tables"]["distributions"]["rows"]}
    assert n_by_obs == {o: (7 if o in ("U2", "T1") else 8) for o in n_by_obs} and len(n_by_obs) == 12
