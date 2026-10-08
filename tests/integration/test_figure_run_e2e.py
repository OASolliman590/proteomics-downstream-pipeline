"""R16e wiring (D-72, C1 to C4): `proteomics run` with `report.figures` declared writes the publication figure catalogue of the run.

Two configurations (copies of the synthetic example, with their input files in tests/fixtures/figures/run):
- example-independent-figures.json (phase 1): QC (4), the differential figures of the two primary-family contrasts
  (volcano, MA, DEP heatmap, dot plot: 8) and the Venn across the two contrasts (1): 13 figures;
- example-independent-permanova.json (phase 2, PERMANOVA enabled): the same 13 and the PCoA with centroids and the PERMDISP plot
  of the PERMANOVA feature set (2): 15 figures.

The expected identifiers are written here. Every expected figure has its PNG, pzfx and HTML; the registered source and output
hashes match the files; `proteomics verify` exits 0 on the run. Negative cases: an invalid block makes only the catalogue
INAPPLICABLE (the run stays COMPLETED); a configuration without `report.figures` writes no catalogue.
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
RUN_CONFIGS = ROOT / "tests" / "fixtures" / "figures" / "run"
FIXTURE_CONFIG = RUN_CONFIGS / "example-independent-figures.json"
PERMANOVA_CONFIG = RUN_CONFIGS / "example-independent-permanova.json"
EXAMPLE_CONFIG = ROOT / "configs" / "examples" / "example-independent.json"

DIFFERENTIAL = ("diff_volcano_labelled", "diff_ma_plot", "diff_dep_heatmap_clustered", "diff_protein_dotplot_brackets")
CONTRASTS = ("disease-control", "treated-disease")
QC = ("qc_intensity_distributions", "qc_missingness_heatmap", "qc_sample_correlation_heatmap", "qc_pca_ellipses")
EXPECTED_PHASE1 = frozenset(list(QC) + [f"{f}__{c}" for f in DIFFERENTIAL for c in CONTRASTS] + ["diff_upset_venn__primary"])
EXPECTED_PERMANOVA = EXPECTED_PHASE1 | {"perm_pcoa_centroids__all", "perm_permdisp__all"}
POSTDE_NOT_RUN = frozenset({"biomarker figures", "path_enrichment_dotplot__ora_up", "path_enrichment_dotplot__ora_down", "net_module_heatmap"})


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _registry(out: Path) -> dict:
    return json.loads((out / "report" / "figure_catalogue" / "registry.json").read_text(encoding="utf-8"))


def _assert_catalogue(out: Path, expected: frozenset[str]) -> None:
    registry = _registry(out)
    catalogue = out / "report" / "figure_catalogue"
    assert set(registry["figures"]) == expected, sorted(set(registry["figures"]) ^ expected)
    for figure_id in sorted(expected):
        record = registry["figures"][figure_id]
        assert record["state"] == "REGISTERED", (figure_id, record.get("reason"))
        assert not [k for k, v in record["outputs"].items() if "code" in v], (figure_id, record["outputs"])
        for kind in ("png", "pzfx", "html"):
            path = catalogue / f"{figure_id}.{kind}"
            assert path.is_file() and path.stat().st_size > 0, path.name
            assert _sha(path) == record["outputs"][kind]["sha256"], path.name
        assert (catalogue / f"{figure_id}.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert _sha(catalogue / record["source_path"]) == record["source_sha256"]
    # only the post-DE figures of modules that did not run may be skipped; no QC, differential or PERMANOVA figure is
    skipped = {s["figure"] for s in registry["retrofit"]["run_sources"]}
    assert skipped <= POSTDE_NOT_RUN, sorted(skipped)


def _verify(out: Path) -> None:
    workflow.verify_run(out)
    completed = subprocess.run([sys.executable, "-m", "proteomics_pipeline", "verify", "--run", str(out), "--json"],
                               capture_output=True, encoding="utf-8", cwd=ROOT, check=False)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert json.loads(completed.stdout)["valid"] is True


@pytest.fixture(scope="module")
def phase1(tmp_path_factory):
    out = tmp_path_factory.mktemp("figure-run") / "run"
    payload, code = workflow.run_command(FIXTURE_CONFIG, out)
    assert code == 0, payload
    return out


@pytest.fixture(scope="module")
def permanova(tmp_path_factory):
    out = tmp_path_factory.mktemp("figure-permanova") / "run"
    payload, code = workflow.run_command(PERMANOVA_CONFIG, out)
    assert code == 0, payload
    return out


def test_the_phase1_run_completes_and_its_figure_stage_is_completed(phase1):
    status = json.loads((phase1 / "run_status.json").read_text(encoding="utf-8"))
    stages = {s["stage_id"]: s for s in status["stages"]}
    assert status["state"] == "COMPLETED"
    assert stages["figure_catalogue"]["state"] == "COMPLETED"
    assert (phase1 / stages["figure_catalogue"]["result_path"]).is_file()


def test_the_phase1_catalogue_holds_exactly_the_expected_thirteen_figures_with_their_outputs(phase1):
    _assert_catalogue(phase1, EXPECTED_PHASE1)


def test_the_phase1_run_verifies_and_covers_every_figure_file(phase1):
    _verify(phase1)
    status = json.loads((phase1 / "run_status.json").read_text(encoding="utf-8"))
    figure_artifacts = {a["relative_path"] for a in status["artifacts"] if a["result_type"] == "FigureArtifact"}
    catalogue = phase1 / "report" / "figure_catalogue"
    for figure_id in EXPECTED_PHASE1:
        for kind in ("png", "pzfx", "html"):
            assert f"report/figure_catalogue/{figure_id}.{kind}" in figure_artifacts


def test_the_permanova_run_completes_and_holds_the_expected_fifteen_figures(permanova):
    status = json.loads((permanova / "run_status.json").read_text(encoding="utf-8"))
    assert status["state"] == "COMPLETED"
    assert {s["stage_id"]: s["state"] for s in status["stages"]}["permanova"] == "COMPLETED"
    _assert_catalogue(permanova, EXPECTED_PERMANOVA)
    _verify(permanova)


def test_an_invalid_block_makes_only_the_catalogue_inapplicable_and_the_run_completes(tmp_path):
    config = json.loads(FIXTURE_CONFIG.read_text(encoding="utf-8"))
    config["report"]["figures"] = {"dpi": 150}
    variant = RUN_CONFIGS / "example-independent-figures-dpi150.json"
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
