"""R16e wiring (D-72): `proteomics run` with `report.figures` declared writes the figure catalogue of the run.

The fixture configuration is a copy of configs/examples/example-independent.json with ``report.figures`` declared and its synthetic
input files copied beside it (tests/fixtures/figures/run). The expected figures are the three volcano tables that the R10b report
of this run registers; every one must have the three outputs, the registered hashes must match the files, and ``proteomics verify``
must pass on the run. Negative cases: an invalid block makes only the catalogue INAPPLICABLE (the run stays COMPLETED), and a
configuration without ``report.figures`` writes no catalogue.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

from proteomics_pipeline import workflow

ROOT = Path(__file__).resolve().parents[2]
FIXTURE_CONFIG = ROOT / "tests" / "fixtures" / "figures" / "run" / "example-independent-figures.json"
EXAMPLE_CONFIG = ROOT / "configs" / "examples" / "example-independent.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("figure-run") / "run"
    payload, code = workflow.run_command(FIXTURE_CONFIG, out)
    assert code == 0, payload
    return out, payload


def _expected(out: Path) -> set[str]:
    return {p.stem for p in (out / "report-full" / "figure_sources").glob("volcano_zero_null_*.tsv")}


def test_the_run_completes_and_the_figure_stage_is_recorded_as_completed(run):
    out, payload = run
    assert payload["state"] == "COMPLETED"
    status = json.loads((out / "run_status.json").read_text(encoding="utf-8"))
    stages = {s["stage_id"]: s for s in status["stages"]}
    assert stages["figure_catalogue"]["state"] == "COMPLETED"
    assert (out / stages["figure_catalogue"]["result_path"]).is_file()


def test_every_expected_volcano_has_its_png_pzfx_and_html(run):
    out, _ = run
    expected = _expected(out)
    assert len(expected) == 3
    registry = json.loads((out / "report" / "figure_catalogue" / "registry.json").read_text(encoding="utf-8"))
    assert set(registry["figures"]) == expected
    catalogue = out / "report" / "figure_catalogue"
    for figure_id in expected:
        assert registry["figures"][figure_id]["state"] == "REGISTERED"
        for kind in ("png", "pzfx", "html"):
            path = catalogue / f"{figure_id}.{kind}"
            assert path.is_file() and path.stat().st_size > 0, path.name
        assert (catalogue / f"{figure_id}.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_the_registered_source_and_output_hashes_match_the_files(run):
    out, _ = run
    catalogue = out / "report" / "figure_catalogue"
    registry = json.loads((catalogue / "registry.json").read_text(encoding="utf-8"))
    for figure_id, record in registry["figures"].items():
        assert _sha(catalogue / record["source_path"]) == record["source_sha256"]
        assert _sha(out / "report-full" / "figure_sources" / f"{figure_id}.tsv") == record["source_sha256"]   # the producer's table
        for kind, output in record["outputs"].items():
            assert _sha(catalogue / output["path"]) == output["sha256"]
            assert output["source_sha256"] == record["source_sha256"]


def test_the_style_and_the_font_record_are_written_to_the_registry(run):
    out, _ = run
    registry = json.loads((out / "report" / "figure_catalogue" / "registry.json").read_text(encoding="utf-8"))
    assert registry["style"]["dpi"] == 300 and registry["style"]["journal_width_mm"] == 120
    assert registry["style"]["formats"] == ["png", "pzfx", "html"]
    assert "font_substitution" in registry            # null when Arial is installed, else the recorded substitution


def test_proteomics_verify_passes_and_covers_the_figure_artifacts(run):
    out, _ = run
    workflow.verify_run(out)
    completed = subprocess.run([sys.executable, "-m", "proteomics_pipeline", "verify", "--run", str(out), "--json"],
                               capture_output=True, encoding="utf-8", cwd=ROOT, check=False)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert json.loads(completed.stdout)["valid"] is True
    status = json.loads((out / "run_status.json").read_text(encoding="utf-8"))
    figure_artifacts = [a for a in status["artifacts"] if a["result_type"] == "FigureArtifact"]
    assert len(figure_artifacts) >= 3 * 3 + 3                       # the outputs and the sources of the three figures


def test_an_invalid_block_makes_only_the_catalogue_inapplicable_and_the_run_completes(tmp_path):
    config = json.loads(FIXTURE_CONFIG.read_text(encoding="utf-8"))
    config["report"]["figures"] = {"dpi": 150}
    variant = FIXTURE_CONFIG.parent / "example-independent-figures-dpi150.json"
    variant.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8", newline="\n")
    try:
        payload, code = workflow.run_command(variant, tmp_path / "run")
    finally:
        variant.unlink()
    assert payload["state"] == "COMPLETED" and code == 0
    status = json.loads((tmp_path / "run" / "run_status.json").read_text(encoding="utf-8"))
    stage = next(s for s in status["stages"] if s["stage_id"] == "figure_catalogue")
    assert stage["state"] == "INAPPLICABLE" and stage["reason_code"] == "E_FIGURE_DPI"
    assert not (tmp_path / "run" / "report" / "figure_catalogue").exists()


def test_a_configuration_without_report_figures_writes_no_catalogue(tmp_path):
    payload, code = workflow.run_command(EXAMPLE_CONFIG, tmp_path / "run")
    assert code == 0 and payload["state"] == "COMPLETED"
    assert not (tmp_path / "run" / "report" / "figure_catalogue").exists()
    status = json.loads((tmp_path / "run" / "run_status.json").read_text(encoding="utf-8"))
    assert all(s["stage_id"] != "figure_catalogue" for s in status["stages"])
