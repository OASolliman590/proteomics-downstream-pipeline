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
    summary = json.loads((CAL / "calibration_summary.json").read_text())
    assert summary["profile"] == "release" and summary["datasets_per_scenario"] >= 1000
    for name in ("null", "mixture", "heavy_tail", "mnar"):
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


def test_v105_negative_liberal_rate_fails_and_smoke_cannot_claim_release(tmp_path):
    assert RC.clopper_pearson_upper(70, 1000) > 0.075                      # 7 % observed: upper bound breaches the gate
    if shutil.which("Rscript") is None:
        pytest.skip("NOT_RUN: Rscript unavailable for the smoke-profile negative")
    code = RC.main(["--datasets", "5", "--pathway-datasets", "3", "--nrot", "49", "--features", "200", "--output", str(tmp_path)])
    summary = json.loads((tmp_path / "calibration_summary.json").read_text())
    assert summary["profile"].startswith("smoke") and code == 1
    assert all(summary["scenarios"][s]["gate"].startswith("NOT_RUN") for s in ("null", "mixture"))


def test_v106_pathway_null_rates_under_their_own_hypotheses_with_mc_precision():
    pathways = json.loads((CAL / "pathways.json").read_text())
    assert pathways["n_datasets"] >= 300 and pathways["nrot"] >= 999
    assert "competitive" in pathways["hypotheses"]["camera"] and "self-contained" in pathways["hypotheses"]["roast"]
    for key in ("camera_null_rate", "roast_null_rate", "roast_mixed_null_rate"):
        result = pathways["results"][key]
        assert result["mc_se"] > 0 and result["upper95"] <= 0.075
        assert abs(result["upper95"] - (result["mean"] + 1.6448536 * result["mc_se"])) < 1e-6
    assert "not calibrated" in pathways["fgsea"]
    assert all(pathways["gate"][k] for k in ("camera_pass", "roast_pass", "roast_mixed_pass"))
