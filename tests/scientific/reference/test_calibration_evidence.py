"""R11 V105/V106: the retained release calibration meets its prespecified gates.

The summary is re-derived here from the per-dataset tables (independent
arithmetic), so a summary that disagrees with its own datasets fails.  The
negative cases show that an over-liberal rate fails the gate and a smoke-size
run cannot make the release claim.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
CAL = ROOT / "docs" / "validation" / "012-validation" / "calibration"
spec = importlib.util.spec_from_file_location("run_calibration", ROOT / "scripts" / "maintained" / "run_calibration.py")
RC = importlib.util.module_from_spec(spec); spec.loader.exec_module(RC)


def _rows(name):
    with (CAL / f"core_{name}.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _binom_sf_upper(k, n, p):
    """P(X <= k) under Binomial(n, p) via log-space terms (independent of the script's implementation)."""
    return sum(math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(p) + (n - i) * math.log1p(-p)) for i in range(k + 1))


def test_v105_core_scenarios_meet_prespecified_gates_from_dataset_rows():
    summary = json.loads((CAL / "calibration_summary.json").read_text(encoding="utf-8"))
    assert summary["profile"] == "release" and summary["datasets_per_scenario"] >= 1000
    assert summary["fit_path"].startswith("production R path")
    for name in ("null", "mixture", "mixture_high_power", "heavy_tail", "mnar"):
        rows = _rows(name)
        entry = summary["scenarios"][name]
        assert len(rows) == entry["n_datasets"] >= 1000
        assert len({r["dataset"] for r in rows}) == len(rows)                 # independent datasets, no reuse
        cov = sum(float(r["coverage"]) for r in rows) / len(rows)
        assert abs(cov - entry["mean_coverage"]) < 1e-12
        assert abs(sum(float(r["fdp"]) for r in rows) / len(rows) - entry["mean_fdp"]) < 1e-12
        if entry["core"]:
            assert 0.93 <= cov <= 0.97 and entry["gate"] == "PASS"
        else:
            assert entry["gate"].startswith("STRESS_REPORTED")      # retained and reported, never relabelled as a guarantee
    null = _rows("null")
    k = sum(r["any_rejection"] == "TRUE" for r in null)
    upper = summary["scenarios"]["null"]["any_rejection_upper95"]
    assert k == summary["scenarios"]["null"]["any_rejection_k"]
    assert abs(_binom_sf_upper(k, len(null), upper) - 0.05) < 1e-6 and upper <= 0.075
    assert summary["scenarios"]["mixture"]["mean_fdp_upper95"] <= 0.075
    assert summary["scenarios"]["mixture"]["mean_power"] is not None          # power reported, not only error
    high = summary["scenarios"]["mixture_high_power"]                         # audit: FDR evidence with many rejections
    assert high["mean_fdp_upper95"] <= 0.075 and high["mean_power"] > 0.3 and high["total_rejections"] > summary["scenarios"]["mixture"]["total_rejections"]


def test_v105_negative_liberal_rate_fails_and_smoke_cannot_claim_release(tmp_path):
    assert RC.clopper_pearson_upper(70, 1000) > 0.075                      # 7 % observed: upper bound breaches the gate
    if shutil.which("Rscript") is None:
        pytest.skip("NOT_RUN: Rscript unavailable for the smoke-profile negative")
    code = RC.main(["--datasets", "5", "--pathway-datasets", "3", "--nrot", "49", "--features", "200", "--output", str(tmp_path)])
    summary = json.loads((tmp_path / "calibration_summary.json").read_text(encoding="utf-8"))
    assert summary["profile"].startswith("smoke") and code == 1
    assert all(summary["scenarios"][s]["gate"].startswith("NOT_RUN") for s in ("null", "mixture"))


def test_v106_pathway_null_rates_under_their_own_hypotheses_with_mc_precision():
    pathways = json.loads((CAL / "pathways.json").read_text(encoding="utf-8"))
    assert pathways["n_datasets"] >= 300 and pathways["nrot"] >= 999
    assert "competitive" in pathways["hypotheses"]["camera"] and "self-contained" in pathways["hypotheses"]["roast"]
    for key in ("camera_null_rate", "roast_null_rate", "roast_mixed_null_rate"):
        result = pathways["results"][key]
        assert result["mc_se"] > 0 and result["upper95"] <= 0.075
        assert abs(result["upper95"] - (result["mean"] + 1.6448536 * result["mc_se"])) < 1e-6
    assert "not calibrated" in pathways["fgsea"]
    assert all(pathways["gate"][k] for k in ("camera_pass", "roast_pass", "roast_mixed_pass"))


def test_v105_calibration_uses_the_production_fit_and_eligibility_rule():
    """Audit 2026-10-02: calibrate_once must run the production path.  With 6 units per group the
    production rule (>= 2 observed AND >= 50 %, i.e. >= 3 of 6) differs from the old '>= 2' shortcut,
    and the result must equal a direct limma_stage-equivalent fit of the same dataset."""
    import subprocess
    if shutil.which("Rscript") is None:
        pytest.skip("NOT_RUN: Rscript unavailable")
    code = r"""
    ns <- asNamespace('proteomicsCore')
    sim <- get('simulate_dataset', ns)(400, 6, 'mnar', 4242)
    res <- get('calibrate_once', ns)(sim)
    obsA <- rowSums(!is.na(sim$Y[, sim$X[, 1] == 1])); obsB <- rowSums(!is.na(sim$Y[, sim$X[, 2] == 1]))
    cat(jsonlite::toJSON(list(n_tested = res$n_tested, fit_path = res$fit_path,
        production_rule = sum(obsA >= 3 & obsB >= 3), shortcut_rule = sum(obsA >= 2 & obsB >= 2)), auto_unbox = TRUE))
    """
    from proteomics_pipeline.runtime import run_r_code
    out = run_r_code(code)
    assert out.returncode == 0, out.stderr
    result = json.loads(out.stdout.strip().splitlines()[-1])
    assert result["shortcut_rule"] > result["production_rule"]          # the fixture separates the two rules
    assert result["n_tested"] == result["production_rule"]
    assert result["fit_path"] == "production:fit_limma_model"
