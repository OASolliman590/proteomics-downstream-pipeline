#!/usr/bin/env python3
"""R15c classifier calibration runner (FR-179 to FR-184; SM45; ADR 0011).

Runs tests/scientific/calibration_methods_extension.R with the replicate counts and seeds of the pre-registered
gate file (docs/validation/016-methods-extension-validation/R15c/gates-preregistration.json), then evaluates every
gate from that file and writes calibration/summary.json. The gate file is pinned by SHA-256 in this module:
a change to the gates after the pre-registration fails with E_CALIBRATION_GATE_CHANGED and nothing is run.
A reduced replicate count is refused (E_CALIBRATION_REPLICATES) unless the caller records it as an adaptation.
A failed gate is recorded as FAIL; no gate is ever changed to make a run pass.

Usage: R_LIBS_USER=$PWD/.r-lib python scripts/maintained/run_methods_calibration.py [--cores 12] [--output DIR] [--evaluate-only]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GATES_PATH = ROOT / "docs" / "validation" / "016-methods-extension-validation" / "R15c" / "gates-preregistration.json"
DEFAULT_OUT = ROOT / "docs" / "validation" / "016-methods-extension-validation" / "R15c" / "calibration"
DRIVER = ROOT / "tests" / "scientific" / "calibration_methods_extension.R"
# SHA-256 of the pre-registered gate file (commit 1c79fad). A change to the file is a change to the gates.
PREREGISTERED_SHA256 = "da7577a6cf5da61bd8b0a0d0ca0b9e3c47819cdf76be8451d6f0e06a34336a26"


class CalibrationError(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def _r14f_helpers():
    spec = importlib.util.spec_from_file_location("r14f_calibration", ROOT / "scripts" / "maintained" / "run_post_de_calibration.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gates(path: Path = GATES_PATH) -> dict:
    digest = sha256(path)
    if digest != PREREGISTERED_SHA256:
        raise CalibrationError("E_CALIBRATION_GATE_CHANGED", f"the gate file differs from the pre-registration (sha256 {digest})")
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _mean_se(values: list[float]) -> tuple[float, float]:
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    return mean, (var / n) ** 0.5


def _auc_gate(values: list[float]) -> dict:
    mean, se = _mean_se(values)
    ok = abs(mean - 0.5) <= 0.03 and mean + 1.96 * se <= 0.53
    return {"n": len(values), "mean": mean, "mc_se": se, "upper95": mean + 1.96 * se, "gate": "PASS" if ok else "FAIL"}


def _perm_gate(p_values: list[float], B: int, binom, chi_sf) -> dict:
    k = sum(p <= 0.05 + 1e-12 for p in p_values)
    counts = [sum(abs(p - j / (B + 1)) < 1e-9 for p in p_values) for j in range(1, B + 2)]
    if sum(counts) != len(p_values):
        raise CalibrationError("E_CALIBRATION_REPLICATES", "a permutation P is not one of the B + 1 attainable values")
    expected = len(p_values) / (B + 1)
    chi2 = sum((c - expected) ** 2 / expected for c in counts)
    binom_p = binom(k, len(p_values), 0.05)
    chi_p = chi_sf(chi2, B)
    ok = binom_p >= 0.01 and chi_p >= 0.01
    return {"n": len(p_values), "B": B, "k_at_or_below_0_05": k, "rate": k / len(p_values), "binomial_upper_tail_p": binom_p,
            "chi_square": chi2, "chi_square_df": B, "chi_square_p": chi_p, "gate": "PASS" if ok else "FAIL"}


def verify_replicate_counts(out: Path, gates: dict, adaptations: list[dict]) -> None:
    """A count below the pre-registered count is refused (E_CALIBRATION_REPLICATES) unless it is recorded as an adaptation."""
    expected = {"nested_auc_null": gates["replicate_counts"]["nested_auc_null"], "permutation_p_null": gates["replicate_counts"]["permutation_p_null_datasets"],
                "planted_auc": gates["replicate_counts"]["planted_auc"], "planted_permutation_p": gates["replicate_counts"]["planted_permutation_datasets"],
                "leaky_reference": gates["replicate_counts"]["leaky_reference"]}
    recorded = {(a["scenario"], a["family"]) for a in adaptations}
    for scenario, count in expected.items():
        if not (out / f"{scenario}.tsv").is_file():
            raise CalibrationError("E_CALIBRATION_REPLICATES", f"{scenario}: the scenario file is missing; the run is not accepted")
        for fam in [f["family"] for f in gates["families"]]:
            rows = [r for r in read_tsv(out / f"{scenario}.tsv") if r["family"] == fam]
            if len(rows) < count and (scenario, fam) not in recorded:
                raise CalibrationError("E_CALIBRATION_REPLICATES", f"{scenario} for {fam}: {len(rows)} of {count} pre-registered replicates and no recorded adaptation")


def evaluate(out: Path, gates: dict) -> dict:
    helpers = _r14f_helpers()
    binom, chi_sf = helpers.binomial_upper_tail, helpers.chi_square_sf
    B = int(gates["replicate_counts"]["permutation_B"])
    per_family: dict[str, dict] = {}
    families = [f["family"] for f in gates["families"]]
    for fam in families:
        entry: dict = {}
        auc = [float(r["pooled_oof_auc"]) for r in read_tsv(out / "nested_auc_null.tsv") if r["family"] == fam]
        entry["nested_auc_null"] = _auc_gate(auc)
        perm = [float(r["permutation_p"]) for r in read_tsv(out / "permutation_p_null.tsv") if r["family"] == fam]
        entry["permutation_p_null"] = _perm_gate(perm, B, binom, chi_sf)
        leaky = [float(r["leaky_auc"]) for r in read_tsv(out / "leaky_reference.tsv") if r["family"] == fam]
        lm, _ = _mean_se(leaky)
        entry["leaky_reference"] = {"n": len(leaky), "mean": lm, "gate": "PASS" if lm > 0.75 else "FAIL"}
        pl = [float(r["pooled_oof_auc"]) for r in read_tsv(out / "planted_auc.tsv") if r["family"] == fam]
        plm, _ = _mean_se(pl)
        entry["planted_signal_auc"] = {"n": len(pl), "mean": plm, "gate": "PASS" if plm >= 0.70 else "FAIL"}
        pp = [float(r["permutation_p"]) for r in read_tsv(out / "planted_permutation_p.tsv") if r["family"] == fam]
        share = sum(p <= 0.05 + 1e-12 for p in pp) / len(pp)
        entry["planted_signal_p"] = {"n": len(pp), "share_at_or_below_0_05": share, "gate": "PASS" if share >= 0.90 else "FAIL"}
        if fam in ("penalized_logistic", "lasso"):
            loo = [float(r["pooled_oof_auc"]) for r in read_tsv(out / "loocv_null.tsv") if r["family"] == fam]
            entry["loocv_null"] = _auc_gate(loo)
        if fam in ("penalized_logistic", "lasso", "elastic_net"):
            bauc = read_tsv(out / "blocked_nested_auc_null.tsv")
            rows = [r for r in bauc if r["family"] == fam]
            grouped = [float(r["pooled_oof_auc"]) for r in rows]; ungrouped = [float(r["ungrouped_auc"]) for r in rows]
            gm = _auc_gate(grouped)
            bperm = [float(r["permutation_p"]) for r in read_tsv(out / "blocked_permutation_p_null.tsv") if r["family"] == fam]
            bp = _perm_gate(bperm, B, binom, chi_sf)
            entry["blocked_null"] = {"nested": gm, "permutation": bp, "gate": "PASS" if gm["gate"] == "PASS" and bp["gate"] == "PASS" else "FAIL"}
            diff = sum(ungrouped) / len(ungrouped) - gm["mean"]
            entry["blocked_ungrouped_reference"] = {"difference": diff, "gate": "PASS" if diff > 0.05 else "FAIL"}
        per_family[fam] = entry
    return {"gates_file": GATES_PATH.relative_to(ROOT).as_posix(), "gates_sha256": PREREGISTERED_SHA256, "B": B, "families": per_family,
            "not_run": {} if _rf_present() else {"random_forest": "E_DEPENDENCY_UNAVAILABLE: randomForest is absent (NOT_RUN for this family only)"}}


def _rf_present() -> bool:
    from proteomics_pipeline.runtime import run_r_code   # the house probe (no `Rscript -e` argv)
    rscript = shutil.which("Rscript")
    if not rscript:
        return False
    proc = run_r_code("quit(status = if (requireNamespace('randomForest', quietly = TRUE)) 0L else 1L)", rscript=rscript)
    return proc.returncode == 0


def run(out: Path, cores: int, evaluate_only: bool) -> dict:
    gates = load_gates()
    counts = gates["replicate_counts"]
    if not evaluate_only:
        rscript = shutil.which("Rscript")
        if not rscript:
            raise CalibrationError("E_DEPENDENCY_UNAVAILABLE", "Rscript is unavailable; calibration is NOT_RUN")
        out.mkdir(parents=True, exist_ok=True)
        subprocess.run([rscript, "--vanilla", str(DRIVER), str(GATES_PATH), str(out), str(cores)], check=True)
    summary = evaluate(out, gates)
    raw = json.loads((out / "calibration_raw.json").read_text(encoding="utf-8"))
    summary["design"] = {"replicate_counts": counts, "seeds": gates["seeds"], "compute": {"seconds": raw["seconds"], "cores": raw["cores"]},
                         "versions": {k: raw[k] for k in ("r_version", "glmnet", "e1071", "randomForest", "proteomicsCore")}}
    summary["reduced_replicates"] = []    # every count is the pre-registered count; a reduction would be listed here as an adaptation
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cores", type=int, default=12)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--evaluate-only", action="store_true")
    args = parser.parse_args(argv)
    try:
        summary = run(args.output, args.cores, args.evaluate_only)
    except CalibrationError as error:
        print(json.dumps({"state": "REFUSED", "code": error.code, "message": str(error)}))
        return 2
    failed = [(fam, gate) for fam, entry in summary["families"].items() for gate, value in entry.items() if isinstance(value, dict) and value.get("gate") == "FAIL"]
    print(json.dumps({"families": list(summary["families"]), "failed_gates": failed}))
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
