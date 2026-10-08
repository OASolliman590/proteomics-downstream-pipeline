"""Retrofit of the existing figures (D-72; R16e): the registered figure-source tables of the R10b and R14 producers are drawn
through the catalogue with the three outputs. Synthetic run trees only; the values written below are the oracle.

Each positive case reads the .pzfx tables back with xml.etree and compares them with the source rows (the producer's values, no
recomputation). Negative cases refuse with a typed code and produce no output for that figure.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from proteomics_pipeline.errors import ConfigurationError
from proteomics_pipeline.figures import run_catalogue
from proteomics_pipeline.runtime import run_r_code


def _r_ready() -> bool:
    probe = run_r_code('cat(all(vapply(c("proteomicsCore", "ggrepel", "ggprism", "ragg"), requireNamespace, logical(1), quietly = TRUE)))')
    return probe.returncode == 0 and probe.stdout.strip().endswith("TRUE")


pytestmark = pytest.mark.skipif(not _r_ready(), reason="NOT_RUN: Rscript with proteomicsCore, ggrepel, ggprism and ragg is unavailable")


def _write(path: Path, header: list[str], rows: list[list[object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)] + ["\t".join("NA" if v is None else (repr(v) if isinstance(v, float) else str(v)) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _pzfx(path: Path, index: int) -> dict:
    """Title, row titles and columns of one table of a .pzfx file, read with xml.etree (independent of the writer)."""
    table = ET.parse(path).getroot().findall("Table")[index]
    rows = None
    if table.find("RowTitlesColumn") is not None:
        rows = [d.text or "" for d in table.find("RowTitlesColumn").find("Subcolumn").findall("d")]
    columns = {}
    for child in table:
        if child.tag in ("XColumn", "YColumn"):
            columns[child.findtext("Title")] = [None if d.text in (None, "") else float(d.text) for d in child.find("Subcolumn").findall("d")]
    return {"title": table.findtext("Title"), "rows": rows, "columns": columns}


def _synthetic_run(run: Path) -> dict:
    """A run tree with the source tables of the producers, as the R10b and R14 stages write them."""
    tables = {}
    # R10b: a volcano and the response axes (the producer's own tables)
    volcano = [[f"F{i:02d}", round(0.3 * i - 2.0, 6), round(0.5 * abs(i - 6), 6), 10 ** -(0.5 * abs(i - 6)), min(1.0, 3 * 10 ** -(0.5 * abs(i - 6))), "family-1",
                "true" if i in (1, 11) else "false", 6] for i in range(1, 12)]
    tables["volcano"] = _write(run / "report-full" / "figure_sources" / "volcano_zero_null_limma-main_disease-control.tsv",
                               ["feature_id", "effect_log2", "neg_log10_p", "p_value", "q_value", "family_id", "q_at_or_below_cutoff", "n_obs_by_required_group"], volcano)
    axes = [[f"A{i}", f"F{i:02d}", round(0.2 * i, 4), round(1.5 - 0.1 * i, 4), "near_restoration" if i % 3 == 0 else "none"] for i in range(1, 8)]
    tables["response"] = _write(run / "report-full" / "figure_sources" / "response_axes.tsv", ["axis_id", "feature_id", "d", "t", "descriptive_class"], axes)
    # R14: ROC with its band, CV repeats, the permutation null, the set regions, the eigengenes, the subnetwork, the association heatmap
    roc = [[round(t, 4), round(1 - 0.09 * i, 4), round(min(1.0, 0.11 * i + 0.1), 4)] for i, t in enumerate([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])]
    tables["roc"] = _write(run / "post_de" / "biomarker" / "figure_sources" / "roc_cv.tsv", ["threshold", "fpr", "tpr"],
                           [[r[0], round(r[0] / 2.0, 4), round(min(1.0, 0.1 + 0.09 * i), 4)] for i, r in enumerate(roc)])
    band_grid = [round(0.1 * k, 4) for k in range(11)]
    tables["band"] = _write(run / "post_de" / "biomarker" / "figure_sources" / "roc_cv_band.tsv", ["fpr", "tpr_lower", "tpr_upper", "claim_label"],
                            [[g, round(max(0.0, g - 0.2), 4), round(min(1.0, g + 0.2), 4), "cross_validated_nested"] for g in band_grid])
    tables["repeats"] = _write(run / "post_de" / "biomarker" / "figure_sources" / "cv_auc_repeats.tsv", ["repeat_id", "auc", "claim_label"],
                               [[i, round(0.7 + 0.02 * i, 4), "cross_validated_nested"] for i in range(1, 6)])
    tables["null"] = _write(run / "post_de" / "biomarker" / "figure_sources" / "permutation_null.tsv", ["permutation", "pooled_oof_auc", "claim_label"],
                            [[i, round(0.5 + 0.01 * ((i * 7) % 13 - 6), 4), "cross_validated_nested"] for i in range(1, 21)])
    tables["upset"] = _write(run / "post_de" / "sets" / "figure_sources" / "upset_all.tsv", ["region_id", "defining_sets", "n_protein_groups", "n_unique_genes"],
                             [["R1", "A", 12, 10], ["R2", "A;B", 5, 5], ["R3", "B", 9, 8], ["R4", "A;B;C", 2, 2]])
    tables["venn"] = _write(run / "post_de" / "sets" / "figure_sources" / "venn_A_B.tsv", ["region_id", "defining_sets", "n_protein_groups", "n_unique_genes"],
                            [["R1", "A", 3, 3], ["R2", "B", 4, 4], ["R3", "A;B", 1, 1]])
    eigen = [[f"O{i}", "Control" if i <= 3 else "Chronic", round(0.3 * i - 0.9, 4), round(-0.2 * i + 0.8, 4), "module_level"] for i in range(1, 7)]
    tables["eigen"] = _write(run / "post_de" / "networks" / "figure_sources" / "eigengenes.tsv", ["observation_id", "group", "M1", "M2", "claim_label"], eigen)
    nodes = [["set1", f"G{i}", round(math.cos(2 * math.pi * i / 4), 6), round(math.sin(2 * math.pi * i / 4), 6), "descriptive"] for i in range(4)]
    tables["ppi_nodes"] = _write(run / "post_de" / "networks" / "figure_sources" / "ppi_subnetwork.tsv", ["set_id", "gene_id", "x", "y", "claim_label"], nodes)
    tables["ppi_edges"] = _write(run / "post_de" / "networks" / "ppi_subnetwork.tsv", ["set_id", "gene_a", "gene_b", "claim_label"],
                                 [["set1", "G0", "G1", "descriptive"], ["set1", "G1", "G2", "descriptive"], ["set1", "G0", "G2", "descriptive"]])
    heat = [[f"P{f}", f"analysis{a}", round(0.1 * f - 0.2 * a, 4), "t", f, a] for f in range(1, 4) for a in range(1, 3)]
    tables["heat"] = _write(run / "post_de" / "association" / "heatmap_source.tsv", ["feature_id", "analysis_id", "value", "statistic", "row_order", "column_order"], heat)
    tables["unknown"] = _write(run / "post_de" / "other" / "figure_sources" / "not_a_known_figure.tsv", ["x", "y"], [[1, 2]])
    return tables


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("retrofit")
    run = root / "run"
    tables = _synthetic_run(run)
    catalogue = run / "report" / "figure_catalogue"
    summary = run_catalogue.build(run, {"formats": ["png", "pzfx", "html"]}, catalogue)
    registry = json.loads((catalogue / "registry.json").read_text(encoding="utf-8"))
    return {"run": run, "catalogue": catalogue, "tables": tables, "summary": summary, "registry": registry}


EXPECTED = {
    "response_axes": "scatter",
    "upset_all": "regions",
    "venn_A_B": "regions",
    "association_heatmap": "heat_table",
}


def test_every_known_producer_table_is_retrofitted_with_its_kind_and_the_unknown_one_is_recorded(built):
    registry = built["registry"]
    assert set(registry["figures"]) == set(EXPECTED)
    for figure_id, kind in EXPECTED.items():
        assert registry["figures"][figure_id]["kind"] == kind
    skipped = {item["path"]: item["reason"] for item in registry["retrofit"]["not_retrofitted"]}
    assert set(skipped) == {"report-full/figure_sources/volcano_zero_null_limma-main_disease-control.tsv",
                            "post_de/other/figure_sources/not_a_known_figure.tsv", "post_de/biomarker/figure_sources/roc_cv_band.tsv",
                            "post_de/biomarker/figure_sources/roc_cv.tsv", "post_de/biomarker/figure_sources/cv_auc_repeats.tsv",
                            "post_de/biomarker/figure_sources/permutation_null.tsv", "post_de/networks/figure_sources/eigengenes.tsv",
                            "post_de/networks/figure_sources/ppi_subnetwork.tsv"}
    # the post-DE tables are drawn by their run-level figures (D3, figures/run_sources.py), from the same tables
    assert "bio_roc_ci_band" in skipped["post_de/biomarker/figure_sources/roc_cv.tsv"]
    assert "bio_cv_auc_vs_permutation" in skipped["post_de/biomarker/figure_sources/permutation_null.tsv"]
    assert "net_module_heatmap" in skipped["post_de/networks/figure_sources/eigengenes.tsv"]
    assert "net_network_graph" in skipped["post_de/networks/figure_sources/ppi_subnetwork.tsv"]
    # the R10b volcano is drawn by the run's own differential figure with the catalogue renderer (test_figure_run_e2e)
    assert "diff_volcano_labelled" in skipped["report-full/figure_sources/volcano_zero_null_limma-main_disease-control.tsv"]
    assert "R13" in registry["retrofit"]["producers"]


@pytest.mark.parametrize("figure_id", sorted(EXPECTED))
def test_each_figure_has_its_three_outputs_and_the_registered_hashes_match_the_files(built, figure_id):
    record = built["registry"]["figures"][figure_id]
    assert record["state"] == "REGISTERED" and not [k for k, v in record["outputs"].items() if "code" in v], record["outputs"]
    catalogue = built["catalogue"]
    for kind in ("png", "pzfx", "html"):
        path = catalogue / f"{figure_id}.{kind}"
        assert path.is_file() and path.stat().st_size > 0
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["outputs"][kind]["sha256"]
    source = catalogue / record["source_path"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == record["source_sha256"]
    assert (catalogue / f"{figure_id}.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_the_upset_regions_are_the_producer_counts_with_the_defining_sets_as_row_titles(built):
    table = _pzfx(built["catalogue"] / "upset_all.pzfx", 0)
    source = _rows(built["tables"]["upset"])
    assert table["rows"] == [r["defining_sets"] for r in source]
    assert list(table["columns"].values())[0] == [float(r["n_protein_groups"]) for r in source]


def test_the_association_cells_equal_the_producer_values(built):
    table = _pzfx(built["catalogue"] / "association_heatmap.pzfx", 0)
    source = _rows(built["tables"]["heat"])
    for r in source:
        cell = table["columns"][r["analysis_id"]][table["rows"].index(r["feature_id"])]
        assert abs(cell - float(r["value"])) <= 1e-12


def test_an_invalid_block_is_a_typed_refusal_before_anything_is_written(tmp_path):
    with pytest.raises(ConfigurationError) as caught:
        run_catalogue.build(tmp_path, {"dpi": 150}, tmp_path / "catalogue")
    assert caught.value.code == "E_FIGURE_DPI"
    assert not (tmp_path / "catalogue").exists()
