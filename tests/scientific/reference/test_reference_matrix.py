"""R11 V104: frozen scientific golden reference matrix.

The matrix (docs/validation/reference-matrix.json) freezes package versions,
tolerance tiers and the independent oracle test for every backend case before
any candidate comparison; unsupported combinations name their typed refusal.
"""
from __future__ import annotations

import json
from pathlib import Path

from proteomics_pipeline import reproduction

ROOT = Path(__file__).resolve().parents[3]
MATRIX = ROOT / "docs" / "validation" / "reference-matrix.json"
LOCK = ROOT / "renv.lock"


def test_v104_matrix_is_frozen_consistent_with_lock_and_oracles_exist():
    result = reproduction.check_reference_matrix(MATRIX, LOCK)
    assert result["consistent"], result["differences"]
    matrix = json.loads(MATRIX.read_text())
    assert set(matrix["tolerances"]) == {"exact", "deterministic", "backend"}
    assert {c["oracle_kind"] for c in matrix["cases"]} <= {"direct_package_call", "hand_arithmetic"}
    backends = {c["backend"] for c in matrix["cases"]}
    assert {"limma", "DEqMS", "proDA", "fgsea", "vegan"} <= backends
    assert len({c["case_id"] for c in matrix["cases"]}) == len(matrix["cases"])


def test_v104_negative_version_drift_missing_oracle_and_unraised_code_are_recorded(tmp_path):
    matrix = json.loads(MATRIX.read_text())
    matrix["cases"][0]["packages"]["limma"] = "0.0.0"
    matrix["cases"][1]["oracle_test"] = "tests/scientific/test_limma.py::test_does_not_exist"
    matrix["cases"][2]["tolerance"] = "loose"
    matrix["unsupported_combinations"].append({"combination": "made-up", "reason_code": "E_NEVER_RAISED_ANYWHERE_XYZ"})
    path = tmp_path / "matrix.json"; path.write_text(json.dumps(matrix))
    result = reproduction.check_reference_matrix(path, LOCK)
    assert not result["consistent"]
    fields = {(d["case_id"], d["field"]) for d in result["differences"]}
    assert ("RM01", "packages/limma") in fields and ("RM02", "oracle_test") in fields
    assert ("RM03", "tolerance") in fields and ("made-up", "reason_code") in fields
    drift = next(d for d in result["differences"] if d["field"] == "packages/limma")
    assert drift["frozen"] == "0.0.0" and drift["observed"] == json.loads(LOCK.read_text())["Packages"]["limma"]["Version"]


def _write(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    path.write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))


def test_v104_candidate_comparison_records_actual_fieldwise_differences(tmp_path):
    base = [{"model_id": "m", "contrast_id": "c", "hypothesis_type": "zero_null", "feature_id": f"P{i}", "effect": 0.5 + i, "p_value": 0.01 * (i + 1), "run_id": "a"} for i in range(3)]
    within = [dict(r, effect=r["effect"] + 1e-12, run_id="b") for r in base]          # inside the backend tier
    outside = [dict(r, effect=r["effect"] + (0.01 if r["feature_id"] == "P1" else 0), run_id="b") for r in base]
    _write(tmp_path / "ref" / "dea" / "zero_null.tsv", base)
    _write(tmp_path / "ok" / "dea" / "zero_null.tsv", within)
    _write(tmp_path / "bad" / "dea" / "zero_null.tsv", outside)
    assert reproduction.semantic_compare(tmp_path / "ref", tmp_path / "ok")["n_differences"] == 0
    diff = reproduction.semantic_compare(tmp_path / "ref", tmp_path / "bad")
    assert diff["n_differences"] == 1
    record = diff["differences"][0]
    assert record["column"] == "effect" and record["key"][3] == "P1" and record["left"] != record["right"]
