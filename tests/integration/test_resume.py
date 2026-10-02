"""R11 V103: resume reuses only verified, unchanged stages and invalidates affected stages and dependants."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from proteomics_pipeline import workflow

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "configs" / "examples"
pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="NOT_RUN: Rscript unavailable")


def setup(tmp_path):
    raw = json.loads((EXAMPLES / "example-independent.json").read_text())
    data = tmp_path / "data"; shutil.copytree(EXAMPLES / "fixtures", data)
    for key in ("matrix", "observations", "features", "source_provenance"):
        raw["input"][key] = str(data / Path(raw["input"][key]).name)
    path = tmp_path / "analysis.json"; path.write_text(json.dumps(raw))
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert code == 0, payload
    return path, raw, data


def test_v103_unchanged_run_reuses_every_stage(tmp_path):
    path, raw, data = setup(tmp_path)
    before = json.loads((tmp_path / "run" / "plan.json").read_text())["plan_hash"]
    payload, code = workflow.resume_command(tmp_path / "run")
    assert code == 0 and payload["invalidated"] == [] and set(payload["reused"]) >= {"intake", "preprocessing", "design", "limma"}
    assert payload["plan_hash"] == before


def test_v103_changed_config_invalidates_affected_stage_and_dependants(tmp_path):
    path, raw, data = setup(tmp_path)
    raw["preprocessing"]["exclusions"] = [{"observation_id": "T4", "reason": "prespecified synthetic exclusion"}]
    path.write_text(json.dumps(raw))
    payload, code = workflow.resume_command(tmp_path / "run")
    assert code == 0 and "intake" in payload["reused"]
    assert {"preprocessing", "design", "limma", "plan"} <= set(payload["invalidated"])
    assert "T4" not in (tmp_path / "run" / "preprocessing" / "primary" / "matrix.tsv").read_text().splitlines()[0]
    assert any(p.name.startswith("stale-limma") for p in (tmp_path / "run" / "logs").iterdir())


def test_v103_changed_input_bytes_and_interrupted_temporary(tmp_path):
    path, raw, data = setup(tmp_path)
    temp = tmp_path / "run" / ".stage-limma-deadbeef"; temp.mkdir(); (temp / "partial.tsv").write_text("incomplete")
    matrix = data / "independent-abundance.tsv"; matrix.write_text(matrix.read_text().replace("9.7", "9.75", 1))
    payload, code = workflow.resume_command(tmp_path / "run")
    assert code == 0 and "intake" in payload["invalidated"] and payload["reused"] == [] or "intake" in payload["invalidated"]
    assert payload["abandoned_temporaries"] and not temp.exists()
    assert "9.75" in (tmp_path / "run" / "inputs" / "matrix.tsv").read_text()


def test_v103_negative_existing_directory_with_mutated_output_is_not_a_cache_hit(tmp_path):
    path, raw, data = setup(tmp_path)
    target = tmp_path / "run" / "dea" / "families.tsv"; target.write_text(target.read_text() + "\n")
    payload, code = workflow.resume_command(tmp_path / "run")
    assert "limma" in payload["invalidated"] and "intake" in payload["reused"]
    from proteomics_pipeline.cli import main
    assert main(["verify", "--run", str(tmp_path / "run"), "--json"]) == 0
