"""Post-DE figures of a run (R16e, D3): the result tables of the post-DE modules are mapped to the part-A renderers, from their own
tables. Synthetic run trees with the producers' table layout; the values written here are the oracle. No workflow run is needed.

Positive cases compare the job sources with the producer tables (no statistic is recomputed); a module that did not run is
skipped with its reason; a displayed P that differs from the producer's own (k + 1) / (B + 1) and a network edge that is not in
the interaction snapshot are refused with typed codes and no output.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from proteomics_pipeline.figures import run_catalogue, run_sources

GROUPS = {"O1": "Control", "O2": "Control", "O3": "Control", "O4": "Disease", "O5": "Disease", "O6": "Disease"}
CLAIM = "cross_validated_nested"


def _write(path: Path, header: list[str], rows: list[list[object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)] + ["\t".join("NA" if v is None else str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _tree(run: Path, *, permutation_p: str = "0.05", permutation: bool = True) -> dict:
    """A run with the post-DE tables of the biomarker, pathway, network and set modules, as the producers write them."""
    run.mkdir(parents=True, exist_ok=True)
    snap = run.parent / "snapshots" / "ppi"
    _write(snap / "edges.tsv", ["gene_a", "gene_b", "score"], [["G0", "G1", 900], ["G1", "G2", 800], ["G0", "G2", 700], ["G2", "G3", 600]])
    (snap / "manifest.json").write_text(json.dumps({"files": [{"name": "edges.tsv", "role": "ppi"}]}), encoding="utf-8", newline="\n")
    (run / "config.resolved.json").write_text(json.dumps({"resources": [{"id": "ppi", "kind": "ppi", "path": str(snap / "manifest.json")}]}),
                                              encoding="utf-8", newline="\n")
    _write(run / "preprocessing" / "primary" / "observations.tsv", ["observation_id", "group"], [[o, g] for o, g in GROUPS.items()])
    bio = run / "post_de" / "biomarker"
    _write(bio / "cv_performance.tsv", ["claim_label", "seed", "pooled_oof_auc", "n_positive_units", "n_negative_units",
                                        "threshold_bootstrap_draws", "permutation_B", "permutation_p"],
           [[CLAIM, 3, 0.9, 3, 3, 100, 19 if permutation else "NA", permutation_p if permutation else "NA"]])
    fs = bio / "figure_sources"
    _write(fs / "roc_cv.tsv", ["threshold", "fpr", "tpr", "claim_label"],
           [["Inf", 0, 0, CLAIM], ["0.8", 0, 0.67, CLAIM], ["0.5", 0.33, 1, CLAIM], ["0.1", 1, 1, CLAIM]])
    _write(fs / "roc_cv_band.tsv", ["fpr", "tpr_lower", "tpr_upper", "claim_label"],
           [[0, 0.2, 0.9, CLAIM], [0.5, 0.7, 1, CLAIM], [1, 1, 1, CLAIM]])
    _write(fs / "cv_auc_repeats.tsv", ["repeat_id", "auc", "claim_label"], [[1, 0.88, CLAIM], [2, 0.92, CLAIM]])
    _write(fs / "permutation_null.tsv", ["permutation", "pooled_oof_auc", "claim_label"],
           [[i, round(0.4 + 0.02 * ((i * 7) % 11), 4), CLAIM] for i in range(1, 20)])
    features = [f"F{i:02d}" for i in range(1, 13)]
    _write(bio / "single_feature_auc.tsv", ["feature_id", "n_positive_units", "n_negative_units", "auc", "bootstrap_lower", "bootstrap_upper", "claim_label"],
           [[f, 3, 3, round(0.5 + 0.04 * i, 4), round(0.3 + 0.04 * i, 4), round(0.7 + 0.02 * i, 4), "in_sample"] for i, f in enumerate(features)])
    _write(run / "pathways" / "ora_up.tsv", ["set_id", "N", "K", "n", "k", "p_value", "q_value", "status", "family_id", "overlap_genes"],
           [["SET_A", 40, 8, 5, 4, 0.001, 0.004, "tested", "ora_up-fam", "[]"], ["SET_B", 40, 12, 5, 0, 0.6, 0.6, "tested", "ora_up-fam", "[]"]])
    net = run / "post_de" / "networks"
    _write(net / "modules.tsv", ["feature_id", "module"], [["F01", "M1"], ["F02", "M1"], ["F03", "M2"], ["F04", "unassigned"]])
    _write(net / "eigengenes.tsv", ["observation_id", "M1", "M2", "claim_label"],
           [[o, round(0.1 * i - 0.2, 4), round(-0.1 * i + 0.3, 4), "module_level"] for i, o in enumerate(GROUPS, start=1)])
    _write(net / "ppi_subnetwork.tsv", ["set_id", "gene_a", "gene_b", "score", "claim_label"],
           [["primary_any", "G0", "G1", 900, "descriptive"], ["primary_any", "G1", "G2", 800, "descriptive"]])
    return {"fs": fs, "net": net}


def _jobs(run: Path, work: Path) -> tuple[dict[str, dict], list[dict[str, str]]]:
    work.mkdir(parents=True, exist_ok=True)
    jobs, skipped = run_sources.build_jobs(run, work)
    return {job["output_id"]: job for job in jobs}, skipped


def test_a_full_postde_tree_maps_each_module_to_its_run_level_figure(tmp_path):
    run = tmp_path / "run"
    _tree(run)
    jobs, skipped = _jobs(run, tmp_path / "work")
    assert set(jobs) == {"bio_roc_ci_band", "bio_cv_auc_vs_permutation", "bio_panel_forest", "path_enrichment_dotplot__ora_up",
                         "net_module_heatmap", "net_network_graph__primary_any"}
    assert {job["figure_id"] for job in jobs.values()} == {"bio_roc_ci_band", "bio_cv_auc_vs_permutation", "bio_panel_forest",
                                                           "path_enrichment_dotplot", "net_module_heatmap", "net_network_graph"}
    assert "path_enrichment_dotplot__ora_down" in [s["figure"] for s in skipped]     # no down table in this tree


def test_the_cv_permutation_figure_shows_the_published_p_and_the_published_null(tmp_path):
    run = tmp_path / "run"
    _tree(run)
    jobs, _ = _jobs(run, tmp_path / "work")
    rows = _rows(jobs["bio_cv_auc_vs_permutation"]["source"])
    null = [r for r in rows if r["series"] == "null"]
    repeats = [r for r in rows if r["series"] == "cv_repeat"]
    source_null = _rows(run / "post_de" / "biomarker" / "figure_sources" / "permutation_null.tsv")
    assert [float(r["value"]) for r in null] == [float(r["pooled_oof_auc"]) for r in source_null]
    assert len(null) == 19 and {r["B"] for r in rows} == {"19"}
    assert {r["P_display"] for r in rows} == {"0.05"}        # the producer's permutation P, not a recomputed one
    assert [float(r["value"]) for r in repeats] == [0.88, 0.92]


def test_the_roc_figure_is_the_producer_curve_and_band_with_the_band_draw_count(tmp_path):
    run = tmp_path / "run"
    _tree(run)
    jobs, _ = _jobs(run, tmp_path / "work")
    rows = _rows(jobs["bio_roc_ci_band"]["source"])
    curve = [r for r in rows if r["component"] == "curve"]
    band = [r for r in rows if r["component"] == "band"]
    producer = _rows(run / "post_de" / "biomarker" / "figure_sources" / "roc_cv.tsv")
    assert [float(r["fpr"]) for r in curve] == [float(r["fpr"]) for r in producer]
    assert [float(r["sensitivity"]) for r in curve] == [float(r["tpr"]) for r in producer]
    assert [float(r["band_lower"]) for r in band] == [0.2, 0.7, 1.0]
    assert {r["B"] for r in rows} == {"100"}                # min(500, threshold_bootstrap_draws)
    assert {r["n_control"] for r in rows} == {"3"} and {r["auc"] for r in rows} == {"0.9"}


def test_the_module_heatmap_holds_the_group_means_of_the_producer_eigengenes(tmp_path):
    run = tmp_path / "run"
    _tree(run)
    jobs, _ = _jobs(run, tmp_path / "work")
    rows = _rows(jobs["net_module_heatmap"]["source"])
    eigen = _rows(run / "post_de" / "networks" / "eigengenes.tsv")
    for r in rows:
        values = [float(e[r["module"]]) for e in eigen if GROUPS[e["observation_id"]] == r["group"]]
        assert abs(float(r["eigengene"]) - sum(values) / len(values)) <= 1e-12
        assert r["declared"] == "TRUE" and r["module"] in {"M1", "M2"}
    assert {r["module"] for r in rows} == {"M1", "M2"}      # the unassigned feature set is not a module


def test_the_forest_is_the_top_ten_single_feature_aucs_captioned_as_in_sample(tmp_path):
    run = tmp_path / "run"
    _tree(run)
    jobs, _ = _jobs(run, tmp_path / "work")
    rows = _rows(jobs["bio_panel_forest"]["source"])
    assert len(rows) == 10
    aucs = [float(r["estimate"]) for r in rows]
    assert aucs == sorted(aucs, reverse=True) and rows[0]["feature"] == "F12"
    assert "in-sample" in jobs["bio_panel_forest"]["options"]["caption"]


def test_the_enrichment_plot_takes_only_sets_with_a_foreground_hit(tmp_path):
    run = tmp_path / "run"
    _tree(run)
    jobs, _ = _jobs(run, tmp_path / "work")
    ids = [r["set_id"] for r in _rows(jobs["path_enrichment_dotplot__ora_up"]["source"])]
    assert ids == ["SET_A", "SET_B"]                        # the table keeps both; the builder plots only k >= 1 (SET_A)


def test_a_run_without_the_post_de_modules_skips_each_figure_with_its_reason(tmp_path):
    run = tmp_path / "run"
    (run / "preprocessing" / "primary").mkdir(parents=True)
    jobs, skipped = _jobs(run, tmp_path / "work")
    assert jobs == {}
    reasons = {s["figure"]: s["reason"] for s in skipped}
    assert "the biomarker module did not run" in reasons["biomarker figures"]
    assert "the pathway stage did not run" in reasons["path_enrichment_dotplot__ora_up"]
    assert "the network stage did not write its module tables" in reasons["net_module_heatmap"]


def test_a_disabled_permutation_test_skips_the_cv_permutation_figure(tmp_path):
    run = tmp_path / "run"
    _tree(run, permutation=False)
    jobs, skipped = _jobs(run, tmp_path / "work")
    assert "bio_cv_auc_vs_permutation" not in jobs and "bio_roc_ci_band" in jobs
    assert any(s["figure"] == "bio_cv_auc_vs_permutation" and "permutation disabled" in s["reason"] for s in skipped)


def test_a_displayed_p_that_is_not_the_rule_value_is_refused_with_no_output(tmp_path):
    run = tmp_path / "run"
    _tree(run, permutation_p="0.2")
    summary = run_catalogue.build(run, {"formats": ["pzfx"]}, run / "report" / "figure_catalogue")
    figure = next(f for f in summary["figures"] if f["figure_id"] == "bio_cv_auc_vs_permutation")
    assert figure["refusals"][0]["code"] == "E_FIGURE_PERMUTATION_P"
    assert not (run / "report" / "figure_catalogue" / "bio_cv_auc_vs_permutation.pzfx").exists()


def test_an_edge_that_is_not_in_the_interaction_snapshot_is_refused_with_no_output(tmp_path):
    run = tmp_path / "run"
    _tree(run)
    _write(run / "post_de" / "networks" / "ppi_subnetwork.tsv", ["set_id", "gene_a", "gene_b", "score", "claim_label"],
           [["primary_any", "G0", "G9", 900, "descriptive"]])
    summary = run_catalogue.build(run, {"formats": ["pzfx"]}, run / "report" / "figure_catalogue")
    figure = next(f for f in summary["figures"] if f["figure_id"] == "net_network_graph__primary_any")
    assert figure["refusals"][0]["code"] == "E_FIGURE_SOURCE_MISSING"
    assert not (run / "report" / "figure_catalogue" / "net_network_graph__primary_any.pzfx").exists()
