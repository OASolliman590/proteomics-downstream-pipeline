#!/usr/bin/env python3
"""Release calibration (packet R11, V105/V106).

Runs the core null/mixture/high-power-mixture and stress scenarios with >= 1000
independent datasets each through the production R path (coverage_tables,
featurewise_estimability, fit_limma_model, adjust_family; see simulate.R), then the correlated pathway
calibration, and writes dataset-level summaries with Monte Carlo intervals.
Seeds are frozen here before any result is read.  Usage:
  R_LIBS_USER=$PWD/.r-lib python scripts/maintained/run_calibration.py [--datasets 1000] [--pathway-datasets 300] [--output DIR]
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import random
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# frozen before reading results; mixture_high_power (effect 2, 20 % alternatives) was added by the 2026-10-02 audit
# because the original mixture has low power (~0.024), so its FDR evidence rests on few rejections
SEED_BANK = {"null": 100_000, "mixture": 200_000, "mixture_high_power": 500_000, "heavy_tail": 300_000, "mnar": 400_000}
CORE = ("null", "mixture", "mixture_high_power")


def clopper_pearson_upper(k: int, n: int, level: float = 0.95) -> float:
    """One-sided exact binomial upper bound via bisection on the binomial CDF."""
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


def bootstrap_upper(values: list[float], seed: int = 1, reps: int = 2000) -> float:
    rng = random.Random(seed); n = len(values)
    means = sorted(sum(rng.choice(values) for _ in range(n)) / n for _ in range(reps))
    return means[int(0.95 * reps) - 1]


def run_scenario(rscript: str, scenario: str, n: int, features: int, per_group: int) -> list[dict]:
    code = ("a <- commandArgs(TRUE); d <- proteomicsCore:::run_core_calibration(a[1], as.integer(a[2]), as.integer(a[3]), as.integer(a[4]), as.integer(a[5])); "
            "write.table(d, stdout(), sep='\\t', quote=FALSE, row.names=FALSE)")
    from proteomics_pipeline.runtime import run_r_code
    out = run_r_code(code, [scenario, n, features, per_group, SEED_BANK[scenario]], rscript=rscript)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-2000:])
    return list(csv.DictReader(io.StringIO(out.stdout), delimiter="\t"))


def main(argv=None):
    parser = argparse.ArgumentParser(); parser.add_argument("--datasets", type=int, default=1000); parser.add_argument("--pathway-datasets", type=int, default=300)
    parser.add_argument("--features", type=int, default=2000); parser.add_argument("--per-group", type=int, default=4); parser.add_argument("--nrot", type=int, default=999)
    parser.add_argument("--output", type=Path, default=ROOT / "docs" / "validation" / "012-validation" / "calibration")
    parser.add_argument("--jobs", type=int, default=1, help="scenarios run in parallel processes (results do not depend on this)")
    args = parser.parse_args(argv)
    rscript = shutil.which("Rscript")
    if not rscript:
        print(json.dumps({"state": "NOT_RUN", "reason": "Rscript unavailable"})); return 3
    args.output.mkdir(parents=True, exist_ok=True)
    report = {"profile": "release" if args.datasets >= 1000 else "smoke (insufficient for the release claim)", "datasets_per_scenario": args.datasets,
              "fit_path": "production R path: coverage_tables -> featurewise_estimability -> fit_limma_model (trend=TRUE, robust=TRUE) -> adjust_family (BH, finite eligible)",
              "features": args.features, "units_per_group": args.per_group, "seed_bank": SEED_BANK, "scenarios": {}}
    from concurrent.futures import ThreadPoolExecutor

    def timed(scenario):
        start = time.time()
        return run_scenario(rscript, scenario, args.datasets, args.features, args.per_group), time.time() - start

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        results = dict(zip(SEED_BANK, pool.map(timed, SEED_BANK)))
    for scenario in SEED_BANK:
        rows, elapsed = results[scenario]
        anyrej = sum(r["any_rejection"] == "TRUE" for r in rows)
        fdp = [float(r["fdp"]) for r in rows]; cov = [float(r["coverage"]) for r in rows]
        power = [float(r["power"]) for r in rows if r["power"] not in ("NA", "")]
        entry = {"n_datasets": len(rows), "seconds": round(elapsed, 1), "core": scenario in CORE,
                 "mean_coverage": sum(cov) / len(cov), "coverage_mc_se": (sum((c - sum(cov) / len(cov)) ** 2 for c in cov) / (len(cov) - 1)) ** 0.5 / len(cov) ** 0.5,
                 "mean_power": (sum(power) / len(power)) if power else None, "mean_fdp": sum(fdp) / len(fdp),
                 "total_rejections": int(sum(float(r["n_rejected"]) for r in rows)), "datasets_with_rejections": sum(float(r["n_rejected"]) > 0 for r in rows)}
        if scenario == "null":
            entry.update(any_rejection_k=anyrej, any_rejection_rate=anyrej / len(rows), any_rejection_upper95=clopper_pearson_upper(anyrej, len(rows)))
            entry["gate"] = "PASS" if entry["any_rejection_upper95"] <= 0.075 and 0.93 <= entry["mean_coverage"] <= 0.97 else "FAIL"
        elif scenario in ("mixture", "mixture_high_power"):
            entry.update(mean_fdp_upper95=bootstrap_upper(fdp))
            entry["gate"] = "PASS" if entry["mean_fdp_upper95"] <= 0.075 and 0.93 <= entry["mean_coverage"] <= 0.97 else "FAIL"
        else:
            entry.update(mean_fdp_upper95=bootstrap_upper(fdp)); entry["gate"] = "STRESS_REPORTED (not a nominal guarantee)"
        report["scenarios"][scenario] = entry
        with (args.output / f"core_{scenario}.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    if args.datasets < 1000:
        for s in CORE:
            report["scenarios"][s]["gate"] = "NOT_RUN (smoke profile; release claim needs >= 1000 datasets)"
    start = time.time()
    subprocess.run([rscript, "--vanilla", str(ROOT / "tests" / "scientific" / "calibration_pathways.R"), str(args.pathway_datasets), str(args.nrot), str(args.output / "pathways.json")], check=True)
    pathways = json.loads((args.output / "pathways.json").read_text(encoding="utf-8"))
    pathways["seconds"] = round(time.time() - start, 1)
    report["pathways"] = pathways
    (args.output / "calibration_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8", newline="\n")
    print(json.dumps({s: e["gate"] for s, e in report["scenarios"].items()} | {"pathways": pathways["gate"]}))
    return 0 if all(report["scenarios"][s]["gate"] == "PASS" for s in CORE) else 1


if __name__ == "__main__":
    raise SystemExit(main())
