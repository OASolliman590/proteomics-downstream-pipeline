#!/usr/bin/env python3
"""R15c follow-up gate evaluation (D-77, D-78, D-82; FR-180, FR-181, FR-183).

Reads the outputs of tests/scientific/calibration_lpo.R and tests/scientific/calibration_rf_diagnostic.R and evaluates
every gate of the follow-up pre-registration (docs/validation/016-methods-extension-validation/R15c/followup/
gates-preregistration.json, and R15c/followup-lasso400/ for D-82), each pinned here by SHA-256: a changed or unknown gate file
fails with E_CALIBRATION_GATE_CHANGED.
A missing or short scenario file fails with E_CALIBRATION_REPLICATES. A failed gate is recorded as FAIL and never changed.
The 2026-10-13 R15c summary is not modified.

Usage: python scripts/maintained/evaluate_calibration_followup.py [--gates FILE] [--output DIR]   (writes DIR/summary.json;
DIR defaults to the calibration folder beside FILE)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FOLLOWUP = ROOT / "docs" / "validation" / "016-methods-extension-validation" / "R15c" / "followup"
GATES_PATH = FOLLOWUP / "gates-preregistration.json"
DEFAULT_OUT = FOLLOWUP / "calibration"
PREREGISTERED_SHA256 = "e98bcb2d30b37605430ef905f2c34f08542371e44d108559b06b9f4ec80d50c4"
LASSO400_GATES = FOLLOWUP.parent / "followup-lasso400" / "gates-preregistration.json"
PINNED = {GATES_PATH: PREREGISTERED_SHA256, LASSO400_GATES: "e8c3c0dac5e03f60d930fd2078e6bb852ddb7331e17aa3142a79319976b25901"}

sys.path.insert(0, str(ROOT / "scripts" / "maintained"))
import run_methods_calibration as base   # noqa: E402  (shared gate rules, TSV reader and error type)


def load_gates(path: Path = GATES_PATH, pinned: str | None = None) -> dict:
    expected = pinned or PINNED.get(Path(path).resolve(), PREREGISTERED_SHA256)
    if base.sha256(path) != expected:
        raise base.CalibrationError("E_CALIBRATION_GATE_CHANGED", "the follow-up gate file differs from its pre-registration; nothing is evaluated")
    return json.loads(path.read_text(encoding="utf-8"))


def _rows(out: Path, name: str, expected: int, family: str) -> list[dict]:
    path = out / name
    if not path.is_file():
        raise base.CalibrationError("E_CALIBRATION_REPLICATES", f"{name} is missing; the run is not accepted")
    rows = [r for r in base.read_tsv(path) if r["family"] == family]
    if len(rows) != expected:
        raise base.CalibrationError("E_CALIBRATION_REPLICATES", f"{name}: {family} has {len(rows)} replicates, {expected} were pre-registered")
    return rows


def binned_uniformity(p_values: list[float], B: int, bins: int, chi_sf) -> dict:
    """Chi-square over `bins` equal groups of the B + 1 attainable values (D-78)."""
    if (B + 1) % bins:
        raise base.CalibrationError("E_CALIBRATION_REPLICATES", f"B + 1 = {B + 1} is not divisible into {bins} equal bins")
    width = (B + 1) // bins
    counts = [0] * bins
    for p in p_values:
        counts[min(bins - 1, (round(p * (B + 1)) - 1) // width)] += 1
    expected = len(p_values) / bins
    chi2 = sum((c - expected) ** 2 / expected for c in counts)
    p = chi_sf(chi2, bins - 1)
    return {"bins": bins, "counts": counts, "expected_per_bin": expected, "chi_square": chi2, "chi_square_df": bins - 1, "chi_square_p": p,
            "gate": "PASS" if p >= 0.01 else "FAIL"}


def evaluate(out: Path, gates: dict) -> dict:
    helpers = base._r14f_helpers()
    binom, chi_sf = helpers.binomial_upper_tail, helpers.chi_square_sf
    rc = gates["replicate_counts"]
    lpo: dict[str, dict] = {}
    for fam in rc["families_lpo"]:
        null = [float(r["lpo_auc"]) for r in _rows(out, "lpo_null.tsv", rc["lpo_null"], fam)]
        lpo[fam] = {"lpo_null": base._auc_gate(null)}
        if rc.get("lpo_planted", 0):
            planted = [float(r["lpo_auc"]) for r in _rows(out, "lpo_planted.tsv", rc["lpo_planted"], fam)]
            mean_planted = sum(planted) / len(planted)
            lpo[fam]["lpo_planted"] = {"n": len(planted), "mean": mean_planted, "gate": "PASS" if mean_planted >= 0.70 else "FAIL",
                                       "note": "planted-signal rule as pre-registered (proposed, pending Maintainer confirmation)"}
    summary = {"leave_pair_out": lpo}
    if "rf_diagnostic_datasets" in rc:
        summary["random_forest_permutation_diagnostic"] = _rf(out, rc, binom, chi_sf)
    return summary


def _rf(out: Path, rc: dict, binom, chi_sf) -> dict:
    B = int(rc["rf_diagnostic_B"])
    rows = _rows(out, "rf_permutation_diagnostic.tsv", rc["rf_diagnostic_datasets"], "random_forest")
    p_values = [float(r["permutation_p"]) for r in rows]
    mid = [float(r["mid_p"]) for r in rows]
    ties = [int(r["n_tied"]) for r in rows]
    original = base._perm_gate(p_values, B, binom, chi_sf)
    binned = binned_uniformity(p_values, B, 10, chi_sf)
    rf = {"type_i": {"k_at_or_below_0_05": original["k_at_or_below_0_05"], "rate": original["rate"], "binomial_upper_tail_p": original["binomial_upper_tail_p"],
                     "gate": "PASS" if original["binomial_upper_tail_p"] >= 0.01 else "FAIL"},
          "uniformity_original": {k: original[k] for k in ("chi_square", "chi_square_df", "chi_square_p")} | {"gate": "PASS" if original["chi_square_p"] >= 0.01 else "FAIL"},
          "uniformity_binned": binned,
          "ties": {"datasets_with_ties": sum(t > 0 for t in ties), "mean_tied_permutations": sum(ties) / len(ties), "max_tied_permutations": max(ties),
                   "mean_permutation_p": sum(p_values) / len(p_values), "mean_mid_p": sum(mid) / len(mid),
                   "mid_p_at_or_below_0_05": sum(m <= 0.05 + 1e-12 for m in mid), "note": "descriptive (no gate); a uniform P has mean 0.5 + 0.5 / (B + 1)"}}
    return rf


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gates", type=Path, default=GATES_PATH)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    out = args.output or args.gates.parent / "calibration"
    try:
        gates = load_gates(args.gates)
        summary = {"gates_file": args.gates.resolve().relative_to(ROOT).as_posix(), "gates_sha256": base.sha256(args.gates), **evaluate(out, gates),
                   "note": "Earlier verdicts stay recorded; acceptance on this evidence is the Maintainer's decision (decisions.md)."}
    except base.CalibrationError as error:
        print(json.dumps({"error": error.code, "message": str(error)}), file=sys.stderr)
        return 2
    (out / "summary.json").write_bytes((json.dumps(summary, indent=1) + "\n").encode("utf-8"))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
