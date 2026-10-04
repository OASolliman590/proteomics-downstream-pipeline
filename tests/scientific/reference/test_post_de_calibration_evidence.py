"""Post-DE calibration evidence (R11 extension recorded with R14f; SM35, SM36, SM40; supports V152, V156, V165).

The retained summary is re-derived here from the per-dataset tables with independent arithmetic (log-space
binomial tails and a chi-square tail by numerical integration written in this file), so a summary that disagrees
with its own datasets, or a gate marked PASS that its numbers do not meet, fails. The negative cases show that an
over-liberal permutation P, an inflated null AUC and a non-floor planted cluster fail the gates.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
CAL = ROOT / "docs" / "validation" / "015-post-de-analysis" / "R14f" / "calibration"
spec = importlib.util.spec_from_file_location("run_post_de_calibration", ROOT / "scripts" / "maintained" / "run_post_de_calibration.py")
RC = importlib.util.module_from_spec(spec); spec.loader.exec_module(RC)


def rows(name):
    with (CAL / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def binom_upper_tail(k, n, p):
    """P(K >= k), log-space terms (independent of the script)."""
    return sum(math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(p) + (n - i) * math.log1p(-p)) for i in range(k, n + 1))


def chi2_sf(x, df, steps=200000):
    """Upper chi-square tail by Simpson integration of the density on [0, x] (independent of the script)."""
    if x <= 0:
        return 1.0
    log_c = -(df / 2) * math.log(2) - math.lgamma(df / 2)
    f = lambda t: math.exp(log_c + (df / 2 - 1) * math.log(t) - t / 2) if t > 0 else (0.0 if df > 2 else math.exp(log_c))
    h = x / steps
    total = f(0) + f(x) + sum((4 if i % 2 else 2) * f(i * h) for i in range(1, steps))
    return max(0.0, 1 - total * h / 3)


def test_post_de_calibration_summary_is_re_derived_from_its_datasets():
    summary = json.loads((CAL / "post_de_calibration_summary.json").read_text(encoding="utf-8"))
    assert summary["fit_path"].startswith("production R functions") and summary["profile"] == "post_de_release"
    s = summary["scenarios"]
    auc = [float(r["pooled_oof_auc"]) for r in rows("nested_auc_null.tsv")]
    assert len(auc) == s["nested_auc_null"]["n"] >= 400 and len({r["seed"] for r in rows("nested_auc_null.tsv")}) == len(auc)   # independent datasets
    mean = sum(auc) / len(auc); se = math.sqrt(sum((a - mean) ** 2 for a in auc) / (len(auc) - 1) / len(auc))
    assert abs(mean - s["nested_auc_null"]["mean"]) < 1e-12 and abs(se - s["nested_auc_null"]["mc_se"]) < 1e-12
    assert (s["nested_auc_null"]["gate"] == "PASS") == (abs(mean - 0.5) <= 0.03 and mean + 1.96 * se <= 0.53)
    perm = [float(r["permutation_p"]) for r in rows("permutation_p_null.tsv")]
    B = s["permutation_p_null"]["B"]
    assert len(perm) == s["permutation_p_null"]["n"] >= 400 and all(abs(p * (B + 1) - round(p * (B + 1))) < 1e-9 for p in perm)
    k = sum(p <= 0.05 + 1e-12 for p in perm)
    assert k == s["permutation_p_null"]["k_at_or_below_0_05"] and abs(binom_upper_tail(k, len(perm), 0.05) - s["permutation_p_null"]["binomial_upper_tail_p"]) < 1e-9
    counts = [sum(round(p * (B + 1)) == j for p in perm) for j in range(1, B + 2)]
    chi2 = sum((c - len(perm) / (B + 1)) ** 2 / (len(perm) / (B + 1)) for c in counts)
    assert counts == s["permutation_p_null"]["value_counts"] and abs(chi2 - s["permutation_p_null"]["chi_square"]) < 1e-9
    assert abs(chi2_sf(chi2, B) - s["permutation_p_null"]["chi_square_p"]) < 1e-5
    assert (s["permutation_p_null"]["gate"] == "PASS") == (binom_upper_tail(k, len(perm), 0.05) >= 0.01 and chi2_sf(chi2, B) >= 0.01)
    leaky = [float(r["leaky_auc"]) for r in rows("leaky_reference.tsv")]
    assert abs(sum(leaky) / len(leaky) - s["leaky_reference"]["mean"]) < 1e-12 and (s["leaky_reference"]["gate"] == "PASS") == (sum(leaky) / len(leaky) > 0.75)
    conn = [float(r["connectivity_p"]) for r in rows("connectivity_p_null.tsv")]
    kc = sum(p <= 0.05 + 1e-12 for p in conn)
    assert kc == s["connectivity_null"]["k_at_or_below_0_05"] and len(conn) == s["connectivity_null"]["n"] >= 400
    upper = s["connectivity_null"]["rate_upper95"]
    assert abs(sum(math.exp(math.lgamma(len(conn) + 1) - math.lgamma(i + 1) - math.lgamma(len(conn) - i + 1) + i * math.log(upper) + (len(conn) - i) * math.log1p(-upper)) for i in range(kc + 1)) - 0.05) < 1e-6
    assert (s["connectivity_null"]["gate"] == "PASS") == (upper <= 0.10 and abs(s["connectivity_null"]["planted_cluster_p"] - 1 / 200) < 1e-12)


def test_post_de_calibration_gates_pass():
    summary = json.loads((CAL / "post_de_calibration_summary.json").read_text(encoding="utf-8"))
    assert {k: v["gate"] for k, v in summary["scenarios"].items()} == dict.fromkeys(("nested_auc_null", "permutation_p_null", "leaky_reference", "connectivity_null",
                                                                                     "blocked_nested_auc_null", "blocked_permutation_p_null", "blocked_ungrouped_reference"), "PASS")


@pytest.mark.parametrize("scenario, mutate", [
    ("permutation_p_null", lambda perm, auc, conn, raw: [0.05] * 60 + perm[60:]),
    ("nested_auc_null", lambda perm, auc, conn, raw: [a + 0.1 for a in auc]),
    ("connectivity_null", lambda perm, auc, conn, raw: raw.update(planted_cluster_p=0.2) or conn)])
def test_post_de_calibration_negative_liberal_or_inflated_results_fail(tmp_path, scenario, mutate):
    for name in ("nested_auc_null.tsv", "permutation_p_null.tsv", "leaky_reference.tsv", "connectivity_p_null.tsv", "post_de_calibration_raw.json",
                 "blocked_nested_auc_null.tsv", "blocked_permutation_p_null.tsv", "blocked_ungrouped_reference.tsv", "post_de_calibration_blocked_raw.json"):
        (tmp_path / name).write_bytes((CAL / name).read_bytes())
    raw = json.loads((tmp_path / "post_de_calibration_raw.json").read_text(encoding="utf-8"))
    table = {"permutation_p_null": ("permutation_p_null.tsv", "permutation_p"), "nested_auc_null": ("nested_auc_null.tsv", "pooled_oof_auc"), "connectivity_null": ("connectivity_p_null.tsv", "connectivity_p")}[scenario]
    data = rows(table[0]); values = [float(r[table[1]]) for r in data]
    perm = values if scenario == "permutation_p_null" else None; auc = values if scenario == "nested_auc_null" else None; conn = values if scenario == "connectivity_null" else None
    new = mutate(perm, auc, conn, raw)
    for r, v in zip(data, new):
        r[table[1]] = repr(v)
    with (tmp_path / table[0]).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(data)
    (tmp_path / "post_de_calibration_raw.json").write_text(json.dumps(raw), encoding="utf-8")
    assert RC.summarize(tmp_path)["scenarios"][scenario]["gate"] == "FAIL"


def test_post_de_calibration_blocked_scenario_is_re_derived_and_passes():
    """Review 2026-10-05 minor 4 (D-54): subject-blocked nested CV, whole-subject permutation and the ungrouped (leaky) reference."""
    s = json.loads((CAL / "post_de_calibration_summary.json").read_text(encoding="utf-8"))["scenarios"]
    auc = [float(r["pooled_oof_auc"]) for r in rows("blocked_nested_auc_null.tsv")]
    ung = [float(r["ungrouped_auc"]) for r in rows("blocked_ungrouped_reference.tsv")]
    perm = [float(r["permutation_p"]) for r in rows("blocked_permutation_p_null.tsv")]
    assert len(auc) == len(ung) == len(perm) == s["blocked_nested_auc_null"]["n"] >= 200
    mean = sum(auc) / len(auc); se = math.sqrt(sum((a - mean) ** 2 for a in auc) / (len(auc) - 1) / len(auc))
    assert abs(mean - s["blocked_nested_auc_null"]["mean"]) < 1e-12 and (s["blocked_nested_auc_null"]["gate"] == "PASS") == (abs(mean - 0.5) <= 0.03 and mean + 1.96 * se <= 0.53)
    B = s["blocked_permutation_p_null"]["B"]; k = sum(p <= 0.05 + 1e-12 for p in perm)
    counts = [sum(round(p * (B + 1)) == j for p in perm) for j in range(1, B + 2)]
    chi2 = sum((c - len(perm) / (B + 1)) ** 2 / (len(perm) / (B + 1)) for c in counts)
    assert k == s["blocked_permutation_p_null"]["k_at_or_below_0_05"] and abs(chi2_sf(chi2, B) - s["blocked_permutation_p_null"]["chi_square_p"]) < 1e-5
    assert (s["blocked_permutation_p_null"]["gate"] == "PASS") == (binom_upper_tail(k, len(perm), 0.05) >= 0.01 and chi2_sf(chi2, B) >= 0.01)
    assert (s["blocked_ungrouped_reference"]["gate"] == "PASS") == (sum(ung) / len(ung) - mean > 0.05)
    assert {s[k]["gate"] for k in ("blocked_nested_auc_null", "blocked_permutation_p_null", "blocked_ungrouped_reference")} == {"PASS"}
