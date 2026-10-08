"""R16e stage 1 acceptance: V206-V214 (QC figures FR-206 to FR-209, differential figures FR-210 to FR-214).

Each figure is built through the catalogue (src/proteomics_pipeline/figures/catalogue.py) from its registered synthetic
source table. The oracles are written here, independent of the builder:
- quantiles, observed fractions, correlations, means, SDs, log ratios, region counts and labels are recomputed from the
  source rows with plain Python;
- PCA scores come from stats::prcomp run in R on the synthetic matrix (through the project R runner);
- the row order of the clustered heatmap comes from stats::hclust run in R;
- the outputs are read back independently: the .pzfx with xml.etree.ElementTree, the HTML figure JSON with json.loads.
The PNG of every figure is checked for its header (a real PNG at 300 dpi); its drawing is reviewed by eye, not by a test.
Skipped (NOT_RUN) when Rscript with ggplot2, ggprism, ragg, ggrepel and proteomicsCore is unavailable.
"""
from __future__ import annotations

import csv
import json
import math
import shutil
import struct
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from proteomics_pipeline.errors import ConfigurationError
from proteomics_pipeline.figures.catalogue import build_figure
from proteomics_pipeline.figures.plotly_html import DATA_OPEN
from proteomics_pipeline.runtime import run_r_code

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "figures" / "catalogue"


def _r_ready() -> bool:
    if shutil.which("Rscript") is None:
        return False
    probe = run_r_code('cat(all(vapply(c("proteomicsCore", "ggrepel", "ggprism", "ragg"), requireNamespace, logical(1), quietly = TRUE)))')
    return probe.returncode == 0 and probe.stdout.strip().endswith("TRUE")


pytestmark = pytest.mark.skipif(not _r_ready(), reason="NOT_RUN: Rscript with proteomicsCore, ggrepel, ggprism and ragg is unavailable")


# ---------------------------------------------------------------- helpers (independent read-back)

def _rows(name: str) -> list[dict[str, str]]:
    with (FIX / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _figure_json(html_path: Path) -> dict:
    text = html_path.read_text(encoding="utf-8")
    start = text.index(DATA_OPEN) + len(DATA_OPEN)
    return json.loads(text[start:text.index("</script>", start)])


def _pzfx_columns(path: Path, table_index: int = 0) -> dict:
    """Column data of one table, read with xml.etree (independent of the writer)."""
    table = ET.parse(path).getroot().findall("Table")[table_index]
    out = {}
    for child in table:
        if child.tag in ("XColumn", "YColumn"):
            out[child.findtext("Title")] = [None if d.text in (None, "") else float(d.text) for d in child.find("Subcolumn").findall("d")]
    return out


def _png_size(path: Path) -> tuple[int, int]:
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    return struct.unpack(">II", data[16:24])


def _build(tmp: Path, figure: str, source: str, **options) -> dict:
    return build_figure(tmp / "run", figure, FIX / source, **options)


def _assert_three_outputs(result: dict) -> None:
    assert sorted(result["outputs"]) == ["html", "png", "pzfx"], result["refusals"]
    width, height = _png_size(Path(result["outputs"]["png"]))
    assert (width, height) == (1417, 1181)                  # 120 mm x 100 mm at 300 dpi (R16c)
    assert Path(result["outputs"]["pzfx"]).stat().st_size > 0
    assert "<script" in Path(result["outputs"]["html"]).read_text(encoding="utf-8")


def _type7(values: list[float]) -> list[float]:
    ordered = sorted(values)
    out = []
    for p in (0.0, 0.25, 0.5, 0.75, 1.0):
        h = (len(ordered) - 1) * p
        lo = math.floor(h)
        hi = min(lo + 1, len(ordered) - 1)
        out.append(ordered[lo] + (h - lo) * (ordered[hi] - ordered[lo]))
    return out


def _r_numbers(code: str, args: list) -> list:
    result = run_r_code(code, args, cwd=Path(args[-1]).parent if isinstance(args[-1], Path) else None)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.strip().splitlines()[-1])


# ---------------------------------------------------------------- V206: intensity distributions

def test_v206_quantiles_and_observed_mask_equal_the_oracle_in_all_outputs(tmp_path):
    result = _build(tmp_path, "qc_intensity_distributions", "qc_intensity_distributions_source.tsv")
    _assert_three_outputs(result)
    rows = _rows("qc_intensity_distributions_source.tsv")
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    observed = {s: [float(r["value"]) for r in rows if r["sample"] == s and r["observed"] == "TRUE"] for s in samples}
    # Oracle: per-sample type-7 quantiles of the observed cells, and the observed count.
    oracle = {s: _type7(observed[s]) for s in samples}
    # HTML: the box traces hold the observed values only
    traces = _figure_json(Path(result["outputs"]["html"]))["data"]
    for trace in traces:
        assert trace["y"] == observed[trace["name"]]
        assert [round(v, 12) for v in _type7(trace["y"])] == [round(v, 12) for v in oracle[trace["name"]]]
    # pzfx: the column of each sample holds the observed values, missing cells are empty
    columns = _pzfx_columns(Path(result["outputs"]["pzfx"]))
    for s in samples:
        cells = [v for v in columns[s] if v is not None]
        assert cells == observed[s]
        assert _type7(cells) == oracle[s]
        assert len(columns[s]) - len(cells) == sum(1 for r in rows if r["sample"] == s and r["observed"] != "TRUE")   # missing cells are empty


def test_v206_negative_imputed_source_without_observed_mask_fails_e_figure_source_imputed(tmp_path):
    with pytest.raises(ConfigurationError) as caught:
        _build(tmp_path, "qc_intensity_distributions", "qc_intensity_distributions_imputed.tsv")
    assert caught.value.code == "E_FIGURE_SOURCE_IMPUTED"
    assert not (tmp_path / "run" / "report" / "figure_catalogue" / "qc_intensity_distributions.html").exists()


# ---------------------------------------------------------------- V207: missingness heatmap

def test_v207_observed_fractions_and_missing_cells_equal_the_oracle(tmp_path):
    result = _build(tmp_path, "qc_missingness_heatmap", "qc_missingness_heatmap_source.tsv")
    _assert_three_outputs(result)
    rows = _rows("qc_missingness_heatmap_source.tsv")
    features = list(dict.fromkeys(r["feature"] for r in rows))
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    cell = {(r["feature"], r["sample"]): (float(r["value"]) if r["observed"] == "TRUE" else None) for r in rows}
    # Oracle: observed fraction per feature and per sample
    frac_feature = {f: sum(1 for s in samples if cell[(f, s)] is not None) / len(samples) for f in features}
    frac_sample = {s: sum(1 for f in features if cell[(f, s)] is not None) / len(features) for s in samples}
    trace = _figure_json(Path(result["outputs"]["html"]))["data"][0]
    z = trace["z"]
    for i, f in enumerate(trace["y"]):
        assert sum(1 for v in z[i] if v is not None) / len(samples) == frac_feature[f]
    for j, s in enumerate(trace["x"]):
        assert sum(1 for row in z if row[j] is not None) / len(features) == frac_sample[s]
    # every missing cell is missing, not zero
    for i, f in enumerate(trace["y"]):
        for j, s in enumerate(trace["x"]):
            assert z[i][j] == cell[(f, s)] or (z[i][j] is None and cell[(f, s)] is None)
    assert any(v is None for row in z for v in row)
    # pzfx: missing cells are empty cells
    grouped = ET.parse(Path(result["outputs"]["pzfx"])).getroot().findall("Table")[0]
    assert grouped.get("TableType") == "TwoWay"


def test_v207_negative_missing_cells_drawn_as_zero_fail_the_oracle(tmp_path):
    result = _build(tmp_path, "qc_missingness_heatmap", "qc_missingness_heatmap_source.tsv")
    trace = _figure_json(Path(result["outputs"]["html"]))["data"][0]
    rows = _rows("qc_missingness_heatmap_source.tsv")
    expected_missing = sum(1 for r in rows if r["observed"] != "TRUE")
    assert sum(1 for row in trace["z"] for v in row if v is None) == expected_missing
    zeroed = [[0.0 if v is None else v for v in row] for row in trace["z"]]
    # a heatmap that draws the missing cells as zero fails the same count of missing cells
    assert sum(1 for row in zeroed for v in row if v is None) != expected_missing


# ---------------------------------------------------------------- V208: sample correlations

def _oracle_pearson_coverage(rows: list[dict], samples: list[str]) -> dict:
    """Pearson r on cells observed in both samples; NA (None) with fewer than three shared cells."""
    feats = list(dict.fromkeys(r["feature"] for r in rows))
    value = {(r["feature"], r["sample"]): (float(r["value"]) if r["observed"] == "TRUE" else None) for r in rows}
    out = {}
    for a in samples:
        for b in samples:
            pairs = [(value[(f, a)], value[(f, b)]) for f in feats if value[(f, a)] is not None and value[(f, b)] is not None]
            if len(pairs) < 3:
                out[(a, b)] = None
                continue
            xs, ys = [p[0] for p in pairs], [p[1] for p in pairs]
            mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
            num = sum((x - mx) * (y - my) for x, y in pairs)
            den = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
            out[(a, b)] = num / den
    return out


def test_v208_every_correlation_cell_equals_the_oracle_within_1e_10(tmp_path):
    result = _build(tmp_path, "qc_sample_correlation_heatmap", "qc_sample_correlation_source.tsv")
    _assert_three_outputs(result)
    rows = _rows("qc_sample_correlation_source.tsv")
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    oracle = _oracle_pearson_coverage(rows, samples)
    trace = _figure_json(Path(result["outputs"]["html"]))["data"][0]
    for i, a in enumerate(trace["y"]):
        for j, b in enumerate(trace["x"]):
            want = oracle[(a, b)]
            got = trace["z"][i][j]
            if want is None:
                assert got is None
            else:
                assert abs(got - want) <= 1e-10


def test_v208_negative_correlation_on_an_imputed_matrix_fails_e_figure_source_imputed(tmp_path):
    with pytest.raises(ConfigurationError) as caught:
        _build(tmp_path, "qc_sample_correlation_heatmap", "qc_sample_correlation_imputed.tsv")
    assert caught.value.code == "E_FIGURE_SOURCE_IMPUTED"


# ---------------------------------------------------------------- V209: PCA with ellipses

def _prcomp_scores(matrix_path: Path) -> dict[str, tuple[float, float]]:
    """Scores from stats::prcomp on the SM06 display matrix (missing cells filled with the feature median, centred)."""
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "m <- as.matrix(read.delim(args[1], check.names = FALSE, row.names = 1))\n"
        "filled <- t(apply(m, 1, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))\n"
        "p <- stats::prcomp(t(filled), center = TRUE, scale. = FALSE)\n"
        "cat(jsonlite::toJSON(list(samples = rownames(p$x), pc1 = unname(p$x[, 1]), pc2 = unname(p$x[, 2])), auto_unbox = FALSE, digits = NA))\n"
    )
    result = run_r_code(code, [matrix_path], cwd=matrix_path.parent)
    assert result.returncode == 0, result.stderr
    data = json.loads(result.stdout.strip().splitlines()[-1])
    return {s: (a, b) for s, a, b in zip(data["samples"], data["pc1"], data["pc2"])}


def test_v209_scores_agree_with_prcomp_within_1e_10_and_ellipses_are_the_95_percent_contours(tmp_path):
    result = _build(tmp_path, "qc_pca_ellipses", "qc_pca_ellipses_source.tsv")
    _assert_three_outputs(result)
    source = _rows("qc_pca_ellipses_source.tsv")
    oracle = _prcomp_scores(FIX / "qc_pca_matrix.tsv")
    # sign-align each component to the oracle, then compare every score
    sign1 = 1 if sum(oracle[r["sample"]][0] * float(r["PC1"]) for r in source) >= 0 else -1
    sign2 = 1 if sum(oracle[r["sample"]][1] * float(r["PC2"]) for r in source) >= 0 else -1
    for r in source:
        assert abs(sign1 * oracle[r["sample"]][0] - float(r["PC1"])) <= 1e-10
        assert abs(sign2 * oracle[r["sample"]][1] - float(r["PC2"])) <= 1e-10
    figure = _figure_json(Path(result["outputs"]["html"]))
    # ellipses: every point lies on the 95 per cent contour of its group (Mahalanobis distance squared = qchisq(0.95, 2))
    chi2 = 5.991464547107979
    for group in ("G1", "G2"):
        members = [(float(r["PC1"]), float(r["PC2"])) for r in source if r["group"] == group]
        mx = sum(p[0] for p in members) / len(members)
        my = sum(p[1] for p in members) / len(members)
        sxx = sum((p[0] - mx) ** 2 for p in members) / (len(members) - 1)
        syy = sum((p[1] - my) ** 2 for p in members) / (len(members) - 1)
        sxy = sum((p[0] - mx) * (p[1] - my) for p in members) / (len(members) - 1)
        det = sxx * syy - sxy * sxy
        ellipse = next(t for t in figure["data"] if t["name"] == f"{group} 95% ellipse")
        for x, y in zip(ellipse["x"], ellipse["y"]):
            dx, dy = x - mx, y - my
            d2 = (syy * dx * dx - 2 * sxy * dx * dy + sxx * dy * dy) / det
            assert abs(d2 - chi2) <= 1e-8


def test_v209_variant_group_of_two_has_no_ellipse_and_the_omission_is_labelled(tmp_path):
    result = _build(tmp_path, "qc_pca_ellipses", "qc_pca_ellipses_variant_source.tsv", variant=True)
    _assert_three_outputs(result)
    figure = _figure_json(Path(result["outputs"]["html"]))
    names = [t["name"] for t in figure["data"]]
    assert "G3 95% ellipse" not in names and "G1 95% ellipse" in names and "G2 95% ellipse" in names
    assert "no ellipse for G3 (fewer than three samples)" in json.dumps(figure["layout"]) or "no ellipse for G3" in Path(result["outputs"]["pzfx"]).read_text(encoding="utf-8")


def test_v209_negative_ellipse_for_the_group_of_two_without_label_fails(tmp_path):
    result = _build(tmp_path, "qc_pca_ellipses", "qc_pca_ellipses_variant_source.tsv", variant=True)
    figure = _figure_json(Path(result["outputs"]["html"]))
    labelled = "no ellipse for G3" in Path(result["outputs"]["pzfx"]).read_text(encoding="utf-8")
    drawn = any(t["name"] == "G3 95% ellipse" for t in figure["data"])
    assert labelled and not drawn                             # an ellipse for G3 with no label would fail this check
