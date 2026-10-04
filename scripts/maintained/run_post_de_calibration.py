#!/usr/bin/env python3
"""Post-DE calibration (R11 extension, recorded with packet R14f; SM35, SM36, SM40).

Runs tests/scientific/calibration_post_de.R through the production R functions and writes dataset-level
tables plus a summary with prespecified gates. The gates were fixed on 2026-10-04 before any calibration
result existed; the first draft's gates were revised before the first run (see the R14f receipt):
  nested_auc_null      pooled out-of-fold AUC on pure noise: mean within 0.5 +/- 0.03 and the upper 95 % bound of the
                       mean <= 0.53 (no optimistic inflation; cross-validated AUC under the null is known to be
                       slightly pessimistic, so the interval is reported but not required to cover 0.5);
  permutation_p_null   whole-procedure permutation P (B = 19, so P takes the 20 values j/20) is uniform: the exact
                       one-sided binomial P of the count at P <= 0.05 (null rate 0.05) >= 0.01 and the chi-square
                       goodness-of-fit P over the 20 attainable values >= 0.01;
  leaky_reference      the deliberately leaky reference is inflated: mean AUC > 0.75 on the same kind of noise;
  connectivity_null    random sets drawn uniformly from the measured universe: exact upper 95 % bound of
                       P(P <= 0.05) <= 0.10 (the degree-binned null is approximate); the planted dense cluster at the floor 1/200.
Subject-blocked scenario (review 2026-10-05, D-54; gates fixed before it was run), 200 datasets of 10 vs 10 subjects x 2 observations:
  blocked_nested_auc_null      grouped nested-CV AUC: |mean - 0.5| <= 0.03 and mean + 1.96 SE <= 0.53;
  blocked_permutation_p_null   permutation of whole subjects, B = 19: binomial upper-tail P >= 0.01 and chi-square P >= 0.01;
  blocked_ungrouped_reference  observation-level folds on the same data: mean AUC exceeds the grouped mean by > 0.05.
Not calibrated: LOOCV, lasso/elastic-net selection, SVM or random-forest classifiers, tuned lambda grids.
Usage: R_LIBS_USER=$PWD/.r-lib python scripts/maintained/run_post_de_calibration.py [--auc 400] [--perm 400] [--B 19] [--ppi 400] [--cores 8] [--output DIR]
       R_LIBS_USER=$PWD/.r-lib python scripts/maintained/run_post_de_calibration.py --blocked [--blocked-n 200] [--cores 8]   (adds the blocked scenario, then re-summarises)
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "docs" / "validation" / "015-post-de-analysis" / "R14f" / "calibration"
GATES = {"nested_auc_null": "|mean - 0.5| <= 0.03 and mean + 1.96 SE <= 0.53", "permutation_p_null": "binomial P(K >= k | rate 0.05) >= 0.01 and chi-square GOF P over the B+1 values >= 0.01",
         "leaky_reference": "mean AUC > 0.75", "connectivity_null": "upper95 P(P<=0.05) <= 0.10 and planted cluster P = 1/200",
         "blocked_nested_auc_null": "|mean - 0.5| <= 0.03 and mean + 1.96 SE <= 0.53", "blocked_permutation_p_null": "binomial P(K >= k | rate 0.05) >= 0.01 and chi-square GOF P >= 0.01",
         "blocked_ungrouped_reference": "mean ungrouped AUC - mean grouped AUC > 0.05"}


def binomial_upper_tail(k: int, n: int, p: float) -> float:
    """P(K >= k) for K ~ Binomial(n, p)."""
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def chi_square_sf(x: float, df: int) -> float:
    """Upper tail of the chi-square distribution: regularized upper incomplete gamma Q(df/2, x/2)."""
    a, z = df / 2.0, x / 2.0
    if z <= 0:
        return 1.0
    if z < a + 1:   # series for P, then Q = 1 - P
        term = total = 1.0 / a; k = a
        for _ in range(10000):
            k += 1; term *= z / k; total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return max(0.0, 1.0 - total * math.exp(-z + a * math.log(z) - math.lgamma(a)))
    b = z + 1 - a; c = 1e300; d = 1 / b; h = d   # Lentz continued fraction for Q
    for i in range(1, 10000):
        an = -i * (i - a); b += 2; d = an * d + b; d = 1e-300 if abs(d) < 1e-300 else d; c = b + an / c; c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d; delta = d * c; h *= delta
        if abs(delta - 1) < 1e-15:
            break
    return math.exp(-z + a * math.log(z) - math.lgamma(a)) * h


def clopper_pearson_upper(k: int, n: int, level: float = 0.95) -> float:
    if k >= n:
        return 1.0
    def cdf(p):
        return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k + 1))
    lo, hi = k / n, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if cdf(mid) > 1 - level:
            lo = mid
        else:
            hi = mid
    return hi


def read(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize(out: Path) -> dict:
    raw = json.loads((out / "post_de_calibration_raw.json").read_text(encoding="utf-8"))
    auc = [float(r["pooled_oof_auc"]) for r in read(out / "nested_auc_null.tsv")]
    perm = [float(r["permutation_p"]) for r in read(out / "permutation_p_null.tsv")]
    leaky = [float(r["leaky_auc"]) for r in read(out / "leaky_reference.tsv")]
    conn = [float(r["connectivity_p"]) for r in read(out / "connectivity_p_null.tsv")]
    B = int(raw["B"])
    mean = sum(auc) / len(auc); se = (sum((a - mean) ** 2 for a in auc) / (len(auc) - 1)) ** 0.5 / len(auc) ** 0.5
    k_perm = sum(p <= 0.05 + 1e-12 for p in perm); k_conn = sum(p <= 0.05 + 1e-12 for p in conn)
    counts = [sum(abs(p - j / (B + 1)) < 1e-9 for p in perm) for j in range(1, B + 2)]
    if sum(counts) != len(perm):
        raise ValueError("a permutation P is not one of the B+1 attainable values j/(B+1)")
    expected = len(perm) / (B + 1)
    chi2 = sum((c - expected) ** 2 / expected for c in counts)
    s = {"nested_auc_null": {"n": len(auc), "mean": mean, "mc_se": se, "lower95": mean - 1.96 * se, "upper95": mean + 1.96 * se, "interval_covers_0_5": mean - 1.96 * se <= 0.5 <= mean + 1.96 * se},
         "permutation_p_null": {"n": len(perm), "B": B, "mean": sum(perm) / len(perm), "k_at_or_below_0_05": k_perm, "rate": k_perm / len(perm),
                                "binomial_upper_tail_p": binomial_upper_tail(k_perm, len(perm), 0.05), "value_counts": counts, "chi_square": chi2, "chi_square_df": B, "chi_square_p": chi_square_sf(chi2, B)},
         "leaky_reference": {"n": len(leaky), "mean": sum(leaky) / len(leaky)},
         "connectivity_null": {"n": len(conn), "k_at_or_below_0_05": k_conn, "rate": k_conn / len(conn), "rate_upper95": clopper_pearson_upper(k_conn, len(conn)), "planted_cluster_p": raw["planted_cluster_p"]}}
    s["nested_auc_null"]["gate"] = "PASS" if abs(mean - 0.5) <= 0.03 and s["nested_auc_null"]["upper95"] <= 0.53 else "FAIL"
    s["permutation_p_null"]["gate"] = "PASS" if s["permutation_p_null"]["binomial_upper_tail_p"] >= 0.01 and s["permutation_p_null"]["chi_square_p"] >= 0.01 else "FAIL"
    s["leaky_reference"]["gate"] = "PASS" if s["leaky_reference"]["mean"] > 0.75 else "FAIL"
    s["connectivity_null"]["gate"] = "PASS" if s["connectivity_null"]["rate_upper95"] <= 0.10 and abs(raw["planted_cluster_p"] - 1 / 200) <= 1e-12 else "FAIL"
    if (out / "blocked_nested_auc_null.tsv").is_file():   # D-54 subject-blocked scenario
        braw = json.loads((out / "post_de_calibration_blocked_raw.json").read_text(encoding="utf-8"))
        bauc = [float(r["pooled_oof_auc"]) for r in read(out / "blocked_nested_auc_null.tsv")]
        bperm = [float(r["permutation_p"]) for r in read(out / "blocked_permutation_p_null.tsv")]
        bung = [float(r["ungrouped_auc"]) for r in read(out / "blocked_ungrouped_reference.tsv")]
        bB = int(braw["B"])
        bmean = sum(bauc) / len(bauc); bse = (sum((a - bmean) ** 2 for a in bauc) / (len(bauc) - 1)) ** 0.5 / len(bauc) ** 0.5
        bk = sum(p <= 0.05 + 1e-12 for p in bperm)
        bcounts = [sum(abs(p - j / (bB + 1)) < 1e-9 for p in bperm) for j in range(1, bB + 2)]
        if sum(bcounts) != len(bperm):
            raise ValueError("a blocked permutation P is not one of the B+1 attainable values")
        bexp = len(bperm) / (bB + 1); bchi = sum((c - bexp) ** 2 / bexp for c in bcounts)
        umean = sum(bung) / len(bung)
        s["blocked_nested_auc_null"] = {"n": len(bauc), "mean": bmean, "mc_se": bse, "lower95": bmean - 1.96 * bse, "upper95": bmean + 1.96 * bse,
                                        "gate": "PASS" if abs(bmean - 0.5) <= 0.03 and bmean + 1.96 * bse <= 0.53 else "FAIL"}
        s["blocked_permutation_p_null"] = {"n": len(bperm), "B": bB, "k_at_or_below_0_05": bk, "rate": bk / len(bperm), "binomial_upper_tail_p": binomial_upper_tail(bk, len(bperm), 0.05),
                                           "value_counts": bcounts, "chi_square": bchi, "chi_square_df": bB, "chi_square_p": chi_square_sf(bchi, bB)}
        s["blocked_permutation_p_null"]["gate"] = "PASS" if s["blocked_permutation_p_null"]["binomial_upper_tail_p"] >= 0.01 and s["blocked_permutation_p_null"]["chi_square_p"] >= 0.01 else "FAIL"
        s["blocked_ungrouped_reference"] = {"n": len(bung), "mean": umean, "grouped_mean": bmean, "difference": umean - bmean, "gate": "PASS" if umean - bmean > 0.05 else "FAIL"}
        raw["blocked"] = braw
    return {"profile": "post_de_release", "gates": GATES, "seeds": raw["seeds"], "design": raw["design"], "network": raw["network"], "seconds": raw["seconds"],
            "fit_path": "production R functions: pd_bm_nested, pd_bm_permute_labels, pd_connectivity_null", "scenarios": s, "blocked": raw.get("blocked"),
            "not_calibrated": "LOOCV, lasso/elastic-net selection, SVM and random-forest classifiers, tuned lambda grids (fixed ridge lambda 0.1 and top-5 AUC selection only)"}


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--auc", type=int, default=400); parser.add_argument("--perm", type=int, default=400); parser.add_argument("--B", type=int, default=19)
    parser.add_argument("--ppi", type=int, default=400); parser.add_argument("--cores", type=int, default=8); parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--blocked", action="store_true"); parser.add_argument("--blocked-n", type=int, default=200)
    args = parser.parse_args(argv)
    rscript = shutil.which("Rscript")
    if not rscript:
        print(json.dumps({"state": "NOT_RUN", "reason": "Rscript unavailable"})); return 3
    args.output.mkdir(parents=True, exist_ok=True)
    driver = str(ROOT / "tests" / "scientific" / "calibration_post_de.R")
    if args.blocked:
        subprocess.run([rscript, "--vanilla", driver, str(args.blocked_n), "0", str(args.B), "0", str(args.cores), args.output.as_posix(), "blocked"], check=True)
    else:
        subprocess.run([rscript, "--vanilla", driver, str(args.auc), str(args.perm), str(args.B), str(args.ppi), str(args.cores), args.output.as_posix()], check=True)
    summary = summarize(args.output)
    (args.output / "post_de_calibration_summary.json").write_bytes((json.dumps(summary, indent=2) + "\n").encode("utf-8"))
    print(json.dumps({k: v["gate"] for k, v in summary["scenarios"].items()}))
    return 0 if all(v["gate"] == "PASS" for v in summary["scenarios"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
