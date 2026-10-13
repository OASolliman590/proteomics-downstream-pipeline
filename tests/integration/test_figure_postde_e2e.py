"""R16e post-DE figures of a run (D3): `proteomics run` on the synthetic post-DE study (tests/fixtures/figures/run/postde_study.py)
with report.figures declared writes exactly the expected catalogue, and `proteomics verify` exits 0 on the run.

The study has one primary contrast (so no UpSet/Venn), no down-regulated enrichment (no set with a foreground hit in that
direction) and the post-DE modules biomarker (cross-validated ROC, CV against the permutation null, single-feature forest),
enrichment (up) and networks (module heatmap, subnetwork). The expected identifiers are written here: 4 QC, 4 differential
(B-A), 6 post-DE, and 3 retrofitted R14 tables (calibration, selection stability, UpSet regions).

NOT_RUN without limma: the workflow stage needs the limma package of the pinned renv lock.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.runtime import run_r_code

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "tests" / "fixtures" / "figures" / "run" / "postde_study.py"

QC = ("qc_intensity_distributions", "qc_missingness_heatmap", "qc_sample_correlation_heatmap", "qc_pca_ellipses")
DIFFERENTIAL = ("diff_volcano_labelled", "diff_ma_plot", "diff_dep_heatmap_clustered", "diff_protein_dotplot_brackets")
POSTDE = ("bio_roc_ci_band", "bio_cv_auc_vs_permutation", "bio_panel_forest", "path_enrichment_dotplot__ora_up", "net_module_heatmap",
          "net_network_graph__primary_any")
RETROFIT = ("calibration", "selection_stability", "upset_all")
EXPECTED = frozenset(list(QC) + [f"{f}__B-A" for f in DIFFERENTIAL] + list(POSTDE) + list(RETROFIT))


def _limma_ready() -> bool:
    probe = run_r_code('cat(requireNamespace("limma", quietly = TRUE))')
    return probe.returncode == 0 and probe.stdout.strip().endswith("TRUE")


pytestmark = pytest.mark.skipif(not _limma_ready(), reason="NOT_RUN: limma (renv lock, Bioconductor) is unavailable to Rscript")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _study_config(directory: Path) -> Path:
    spec = importlib.util.spec_from_file_location("postde_study", STUDY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.study_config(directory)


@pytest.fixture(scope="module")
def postde_run(tmp_path_factory):
    base = tmp_path_factory.mktemp("postde-figures")
    out = base / "run"
    payload, code = workflow.run_command(_study_config(base / "study"), out)
    assert code == 0, payload
    return out


def test_the_postde_run_completes_and_its_catalogue_holds_exactly_the_expected_figures(postde_run):
    status = json.loads((postde_run / "run_status.json").read_text(encoding="utf-8"))
    assert status["state"] == "COMPLETED"
    assert {s["stage_id"]: s["state"] for s in status["stages"]}["figure_catalogue"] == "COMPLETED"
    catalogue = postde_run / "report" / "figure_catalogue"
    registry = json.loads((catalogue / "registry.json").read_text(encoding="utf-8"))
    assert set(registry["figures"]) == EXPECTED, sorted(set(registry["figures"]) ^ EXPECTED)
    for figure_id in sorted(EXPECTED):
        record = registry["figures"][figure_id]
        assert record["state"] == "REGISTERED", (figure_id, record.get("reason"))
        for kind in ("png", "pzfx", "html"):
            path = catalogue / f"{figure_id}.{kind}"
            assert path.is_file() and path.stat().st_size > 0, path.name
            assert _sha(path) == record["outputs"][kind]["sha256"], path.name
        assert _sha(catalogue / record["source_path"]) == record["source_sha256"]


def test_the_skipped_post_de_figures_are_recorded_with_their_reasons(postde_run):
    registry = json.loads((postde_run / "report" / "figure_catalogue" / "registry.json").read_text(encoding="utf-8"))
    skipped = {s["figure"]: s["reason"] for s in registry["retrofit"]["run_sources"]}
    assert "diff_upset_venn__primary" in skipped and "fewer than two primary contrasts" in skipped["diff_upset_venn__primary"]
    assert "path_enrichment_dotplot__ora_down" in skipped


def test_the_postde_run_verifies(postde_run):
    workflow.verify_run(postde_run)
    completed = subprocess.run([sys.executable, "-m", "proteomics_pipeline", "verify", "--run", str(postde_run), "--json"],
                               capture_output=True, encoding="utf-8", cwd=ROOT, check=False)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert json.loads(completed.stdout)["valid"] is True
