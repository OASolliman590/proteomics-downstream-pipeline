"""R11 V101 (pinned environments) and V102 (offline semantic reproduction)."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from proteomics_pipeline import reproduction, workflow
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "configs" / "examples"
HAS_R = shutil.which("Rscript") is not None


def test_v101_python_lock_matches_the_project_environment():
    assert reproduction.check_python_lock(ROOT / "requirements.lock") == []
    assert "jsonschema" in reproduction.parse_python_lock(ROOT / "requirements.lock")


def test_v101_negative_floating_pin_and_mismatch(tmp_path):
    (tmp_path / "bad.lock").write_text("jsonschema>=4\n")
    with pytest.raises(ProteomicsError) as error:
        reproduction.parse_python_lock(tmp_path / "bad.lock")
    assert error.value.code == "E_LOCK_FLOATING"
    (tmp_path / "old.lock").write_text("jsonschema==0.1\n")
    assert reproduction.check_python_lock(tmp_path / "old.lock")[0]["installed"] != "0.1"


@pytest.mark.skipif(not HAS_R, reason="NOT_RUN: Rscript unavailable")
def test_v101_r_lock_matches_installed_library_and_api_probes(tmp_path):
    result = reproduction.check_r_lock(ROOT / "renv.lock")
    assert result["qualified"], result
    assert result["bioconductor"] == "3.23" and all(result["api_probes"].values())
    lock = json.loads((ROOT / "renv.lock").read_text())
    lock["Packages"]["limma"]["Version"] = "0.0.1"
    (tmp_path / "renv.lock").write_text(json.dumps(lock))
    bad = reproduction.check_r_lock(tmp_path / "renv.lock")
    assert not bad["qualified"] and any(m["package"] == "limma" for m in bad["mismatches"])
    lock["Packages"]["limma"]["Version"] = ""
    (tmp_path / "float.lock").write_text(json.dumps(lock))
    with pytest.raises(ProteomicsError) as error:
        reproduction.check_r_lock(tmp_path / "float.lock")
    assert error.value.code == "E_LOCK_FLOATING"


def _config(tmp_path, mutate=None):
    raw = json.loads((EXAMPLES / "example-independent.json").read_text())
    for key in ("matrix", "observations", "features", "source_provenance"):
        raw["input"][key] = str((EXAMPLES / raw["input"][key]).resolve())
    if mutate:
        mutate(raw)
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "analysis.json"; path.write_text(json.dumps(raw))
    return path


@pytest.mark.skipif(not HAS_R, reason="NOT_RUN: Rscript unavailable")
def test_v102_offline_reproduction_is_semantically_identical(tmp_path):
    path = _config(tmp_path)
    workflow.run_command(path, tmp_path / "original")
    result = reproduction.reproduce(tmp_path / "original", tmp_path / "reproduced")
    assert result["reproduced"], result
    assert "dea/zero_null.tsv" in result["semantic"]["compared_tables"] and result["same_plan_hash"]


@pytest.mark.skipif(not HAS_R, reason="NOT_RUN: Rscript unavailable")
def test_v102_negative_changed_sample_n_is_a_scientific_difference(tmp_path):
    workflow.run_command(_config(tmp_path / "a"), tmp_path / "a" / "run")
    workflow.run_command(_config(tmp_path / "b", lambda c: c["preprocessing"].update({"exclusions": [{"observation_id": "C1", "reason": "prespecified synthetic exclusion"}]})), tmp_path / "b" / "run")
    comparison = reproduction.semantic_compare(tmp_path / "a" / "run", tmp_path / "b" / "run")
    assert comparison["n_differences"] > 0
    assert any(d["table"] == "preprocessing/primary/matrix.tsv" for d in comparison["differences"]) or any(d["table"] == "dea/zero_null.tsv" for d in comparison["differences"])


def test_v102_network_attempt_is_refused():
    import socket
    with reproduction.offline():
        with pytest.raises(ProteomicsError) as error:
            socket.create_connection(("example.invalid", 80))
    assert error.value.code == "E_NETWORK_FORBIDDEN"


# --------------------------------------------------------------------------- V110 validation ledger
def test_v110_ledger_records_every_identity_with_distinct_statuses():
    ledger = reproduction.build_acceptance_ledger(ROOT)
    trace = json.loads((ROOT / "specs" / "001-downstream-proteomics" / "traceability.json").read_text())["requirements"]
    assert [r["acceptance_id"] for r in ledger] == [t["acceptance"] for t in trace]
    for record in ledger:
        reproduction.validate_evidence(record)
        if record["status"] == "PASS":
            assert record["commands"] and all(c["exit_code"] == 0 for c in record["commands"]) and record["reviewed_tree"]
        else:
            assert record["reason"]
    committed = json.loads((ROOT / "docs" / "validation" / "acceptance.json").read_text())
    assert {r["acceptance_id"]: r["status"] for r in committed["records"]} == {r["acceptance_id"]: r["status"] for r in ledger}
    assert committed["verified_commit"] is None


@pytest.mark.parametrize("mutation, message", [
    (lambda r: r.update(commands=[]), "commands"),
    (lambda r: r["commands"][0].update(exit_code=1), "exit"),
    (lambda r: r.update(status="DONE"), "status"),
    (lambda r: (r.update(status="NOT_RUN"), r.pop("reason", None)), "reason"),
    (lambda r: r.update(reviewed_tree="abc"), "tree"),
])
def test_v110_negative_incomplete_or_conflated_records_are_rejected(mutation, message):
    record = next(r for r in reproduction.build_acceptance_ledger(ROOT) if r["status"] == "PASS" and r["commands"])
    record = json.loads(json.dumps(record))
    mutation(record)
    with pytest.raises(ProteomicsError) as error:
        reproduction.validate_evidence(record)
    assert error.value.code == "E_EVIDENCE_SCHEMA", message


# --------------------------------------------------------------------------- V107 (structure only; the jobs themselves are NOT_RUN until they run on GitHub)
def test_v107_workflow_separates_smoke_ci_from_release_calibration():
    yaml = pytest.importorskip("yaml")
    workflow_file = yaml.safe_load((ROOT / ".github" / "workflows" / "maintained.yml").read_text())
    jobs = workflow_file["jobs"]
    assert set(jobs["smoke"]["strategy"]["matrix"]["os"]) == {"ubuntu-latest", "windows-latest"}
    smoke = " ".join(str(s.get("run", "")) for s in jobs["smoke"]["steps"])
    assert "pytest" in smoke and "testthat" in smoke and "proteomics_pipeline run" in smoke and "renv::restore" in smoke
    assert "run_calibration" not in smoke and "benchmark" not in smoke
    assert "workflow_dispatch" in jobs["release-calibration"]["if"]
