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

def test_v207_cell_states_fractions_order_and_percentages_equal_the_oracle(tmp_path):
    result = _build(tmp_path, "qc_missingness_heatmap", "qc_missingness_heatmap_source.tsv")
    _assert_three_outputs(result)
    rows = _rows("qc_missingness_heatmap_source.tsv")
    features = list(dict.fromkeys(r["feature"] for r in rows))
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    # Oracle: the source mask (1 observed, 0 missing), the observed fraction per feature and per sample, the missing percent
    mask = {(r["feature"], r["sample"]): (1 if r["observed"] == "TRUE" else 0) for r in rows}
    frac_feature = {f: sum(mask[(f, s)] for s in samples) / len(samples) for f in features}
    missing_pct = {s: 100.0 * sum(1 - mask[(f, s)] for f in features) / len(features) for s in samples}
    trace = _figure_json(Path(result["outputs"]["html"]))["data"][0]
    z = trace["z"]
    # the heatmap rows are the features, ordered by missing fraction (most missing first)
    assert trace["y"] == sorted(features, key=lambda f: (-(1 - frac_feature[f]), f))
    column_samples = [label.split(" ")[0] for label in trace["x"]]
    assert column_samples == samples
    for i, f in enumerate(trace["y"]):
        for j, s in enumerate(column_samples):
            assert z[i][j] == mask[(f, s)]                              # the drawn state equals the source mask
        assert sum(z[i]) / len(samples) == frac_feature[f]
    for label in trace["x"]:
        s, pct = label.split(" ")
        assert pct == f"{missing_pct[s]:.0f}%"
    # pzfx: the same state table, every cell present (1 or 0)
    columns = _pzfx_columns(Path(result["outputs"]["pzfx"]))
    for s in samples:
        assert columns[s] == [float(mask[(f, s)]) for f in trace["y"]]


def test_v207_negative_missing_cells_drawn_as_observed_fail_the_oracle(tmp_path):
    result = _build(tmp_path, "qc_missingness_heatmap", "qc_missingness_heatmap_source.tsv")
    trace = _figure_json(Path(result["outputs"]["html"]))["data"][0]
    rows = _rows("qc_missingness_heatmap_source.tsv")
    mask = {(r["feature"], r["sample"]): (1 if r["observed"] == "TRUE" else 0) for r in rows}
    expected_missing = sum(1 for r in rows if r["observed"] != "TRUE")
    assert sum(1 for row in trace["z"] for v in row if v == 0) == expected_missing
    drawn_as_observed = [[1 if v == 0 else v for v in row] for row in trace["z"]]
    # a heatmap that draws the missing cells as observed no longer equals the mask
    mismatches = sum(1 for i, f in enumerate(trace["y"]) for j, s in enumerate(trace["x"]) if drawn_as_observed[i][j] != mask[(f, s.split(" ")[0])])
    assert mismatches == expected_missing and mismatches > 0


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


def _prcomp_variance_pct(matrix_path: Path) -> tuple[float, float]:
    """Variance explained by PC1 and PC2 in per cent, from the sdev of the same prcomp (SM06 display matrix)."""
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "m <- as.matrix(read.delim(args[1], check.names = FALSE, row.names = 1))\n"
        "filled <- t(apply(m, 1, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))\n"
        "p <- stats::prcomp(t(filled), center = TRUE, scale. = FALSE)\n"
        "cat(jsonlite::toJSON(as.list(100 * p$sdev^2 / sum(p$sdev^2)), auto_unbox = TRUE, digits = NA))\n"
    )
    result = run_r_code(code, [matrix_path], cwd=matrix_path.parent)
    assert result.returncode == 0, result.stderr
    values = json.loads(result.stdout.strip().splitlines()[-1])
    return values[0], values[1]


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
    # axis titles carry the variance explained by each component, from the same prcomp (one decimal)
    pct1, pct2 = _prcomp_variance_pct(FIX / "qc_pca_matrix.tsv")
    assert figure["layout"]["xaxis"]["title"]["text"] == f"PC1 ({pct1:.1f}%)"
    assert figure["layout"]["yaxis"]["title"]["text"] == f"PC2 ({pct2:.1f}%)"
    assert f"PC1 ({pct1:.1f}%)" in Path(result["outputs"]["pzfx"]).read_text(encoding="utf-8")
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


# ---------------------------------------------------------------- V210: labelled volcano

def test_v210_point_classes_and_labelled_features_equal_the_oracle(tmp_path):
    labels = [r["feature"] for r in _rows("diff_volcano_labels.tsv")]
    result = _build(tmp_path, "diff_volcano_labelled", "diff_volcano_labelled_source.tsv", labels=labels)
    _assert_three_outputs(result)
    rows = _rows("diff_volcano_labelled_source.tsv")
    cut = 0.05
    oracle = {}
    for r in rows:
        q, e = float(r["q"]), float(r["effect"])
        oracle[r["feature"]] = "up" if q < cut and e > 0 else ("down" if q < cut and e < 0 else "ns")
    figure = _figure_json(Path(result["outputs"]["html"]))
    for trace in figure["data"]:
        for customdata in trace["customdata"]:
            assert oracle[customdata[0]] == trace["name"]
    assert sorted(oracle) == sorted(c[0] for t in figure["data"] for c in t["customdata"])
    # labelled features: exactly the declared subset, all present in the source
    assert result["content"]["png"]["labels"] == labels
    assert set(labels) <= {r["feature"] for r in rows} and len(labels) <= 10


def test_v210_negative_label_for_a_feature_absent_from_the_source_fails_e_figure_source_missing(tmp_path):
    labels = [r["feature"] for r in _rows("diff_volcano_labels_absent.tsv")]
    with pytest.raises(ConfigurationError) as caught:
        _build(tmp_path, "diff_volcano_labelled", "diff_volcano_labelled_source.tsv", labels=labels)
    assert caught.value.code == "E_FIGURE_SOURCE_MISSING"
    assert "F999" in caught.value.message


# ---------------------------------------------------------------- V211: MA plot

def test_v211_a_and_m_equal_the_oracle_within_1e_10_and_the_sign_follows_the_numerator(tmp_path):
    result = _build(tmp_path, "diff_ma_plot", "diff_ma_plot_source.tsv")
    _assert_three_outputs(result)
    rows = _rows("diff_ma_plot_source.tsv")
    figure = _figure_json(Path(result["outputs"]["html"]))
    for trace in figure["data"]:
        members = [r for r in rows if r["contrast"] == trace["name"]]
        assert len(members) == len(trace["x"]) == 200
        for r, x, y in zip(members, trace["x"], trace["y"]):
            a_oracle = (float(r["mean_numerator"]) + float(r["mean_denominator"])) / 2
            m_oracle = float(r["mean_numerator"]) - float(r["mean_denominator"])
            assert abs(x - a_oracle) <= 1e-10
            assert abs(y - m_oracle) <= 1e-10
            assert r["numerator"] == trace["name"].split(" - ")[0]


def test_v211_negative_reversed_sign_of_m_fails(tmp_path):
    result = _build(tmp_path, "diff_ma_plot", "diff_ma_plot_source.tsv")
    figure = _figure_json(Path(result["outputs"]["html"]))
    rows = _rows("diff_ma_plot_source.tsv")
    trace = figure["data"][0]
    member = [r for r in rows if r["contrast"] == trace["name"]][0]
    m_oracle = float(member["mean_numerator"]) - float(member["mean_denominator"])
    assert abs(trace["y"][0] - m_oracle) <= 1e-10
    assert abs(trace["y"][0] + m_oracle) > 1e-6               # a reversed sign would not equal the oracle


# ---------------------------------------------------------------- V212: per-protein dot plots with P brackets

def _fr198_label(p: float) -> str:
    """The FR-198 exact label, written here: two significant digits (trailing zeros kept), a threshold never crossed from below."""
    if p < 0.0001:
        return "p < 0.0001"
    for digits in range(2, 18):
        exponent = int(f"{p:.{digits - 1}e}".split("e")[1])
        text = f"{p:.{max(0, digits - 1 - exponent)}f}"
        if not any(p < t <= float(text) for t in (0.05, 0.01, 0.001, 0.0001)):
            return f"p = {text}"
    return f"p = {text}"


def test_v212_group_means_sds_values_and_bracket_labels_equal_the_oracle(tmp_path):
    from proteomics_pipeline.figures.catalogue import _dotplot_brackets, read_source
    family = _rows("diff_dotplot_family.tsv")
    result = _build(tmp_path, "diff_protein_dotplot_brackets", "diff_dotplot_values.tsv", family=family)
    _assert_three_outputs(result)
    values = _rows("diff_dotplot_values.tsv")
    figure = _figure_json(Path(result["outputs"]["html"]))
    for trace in figure["data"]:
        protein = trace["name"]
        groups = {}
        for r in values:
            if r["protein"] == protein:
                groups.setdefault(r["group"], []).append(float(r["value"]))
        for g, vals in groups.items():
            got = [y for x, y in zip(trace["x"], trace["y"]) if x == g]
            assert got == vals                                  # every individual value
            mean = sum(vals) / len(vals)
            sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (len(vals) - 1))
            assert abs(sum(got) / len(got) - mean) <= 1e-12 and abs(math.sqrt(sum((v - sum(got) / len(got)) ** 2 for v in got) / (len(got) - 1)) - sd) <= 1e-12
    labels = {(b["protein"], b["group1"], b["group2"]): b["label"] for b in _dotplot_brackets(family, read_source(FIX / "diff_dotplot_values.tsv"))}
    for f in family:
        assert labels[(f["protein"], f["group_a"], f["group_b"])] == _fr198_label(float(f["P"]))


def test_v212_negative_a_bracket_with_a_p_value_not_in_the_family_table_fails(tmp_path):
    from proteomics_pipeline.figures.catalogue import _dotplot_brackets, read_source
    family = _rows("diff_dotplot_family.tsv")
    rows = read_source(FIX / "diff_dotplot_values.tsv")
    bogus = [{"protein": "P1", "group1": "Control", "group2": "Acute", "P": "0.5"}]
    with pytest.raises(ConfigurationError) as caught:
        _dotplot_brackets(family, rows, bogus)
    assert caught.value.code == "E_FIGURE_SOURCE_MISSING"


# ---------------------------------------------------------------- V213: clustered DEP heatmap

def _hclust_order(matrix_path: Path) -> list[str]:
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "m <- as.matrix(read.delim(args[1], check.names = FALSE, row.names = 1))\n"
        "scaled <- t(scale(t(m)))\n"
        "hc <- stats::hclust(stats::dist(scaled, method = 'euclidean'), method = 'complete')\n"
        "cat(jsonlite::toJSON(rownames(scaled)[hc$order], auto_unbox = FALSE))\n"
    )
    result = run_r_code(code, [matrix_path], cwd=matrix_path.parent)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.strip().splitlines()[-1])


def test_v213_row_order_and_row_scaled_values_equal_the_oracle_and_only_the_dep_set_appears(tmp_path):
    result = _build(tmp_path, "diff_dep_heatmap_clustered", "diff_dep_heatmap_source.tsv")
    _assert_three_outputs(result)
    rows = _rows("diff_dep_heatmap_source.tsv")
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    features = list(dict.fromkeys(r["feature"] for r in rows))
    matrix_path = tmp_path / "dep_matrix.tsv"
    with matrix_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("feature\t" + "\t".join(samples) + "\n")
        for f in features:
            handle.write(f + "\t" + "\t".join(repr(float(next(r for r in rows if r["feature"] == f and r["sample"] == s)["value"])) for s in samples) + "\n")
    oracle_order = _hclust_order(matrix_path)
    trace = _figure_json(Path(result["outputs"]["html"]))["data"][0]
    assert trace["y"] == oracle_order                                # row order
    assert set(trace["y"]) == set(features)                          # only the declared DEP set
    values = {f: [float(next(r for r in rows if r["feature"] == f and r["sample"] == s)["value"]) for s in samples] for f in features}
    for i, f in enumerate(trace["y"]):
        v = values[f]
        mu = sum(v) / len(v)
        sd = math.sqrt(sum((x - mu) ** 2 for x in v) / (len(v) - 1))
        for j in range(len(samples)):
            assert abs(trace["z"][i][j] - (v[j] - mu) / sd) <= 1e-10


def test_v213_negative_a_non_dep_feature_in_the_heatmap_fails_e_figure_set_mismatch(tmp_path):
    with pytest.raises(ConfigurationError) as caught:
        _build(tmp_path, "diff_dep_heatmap_clustered", "diff_dep_heatmap_with_nondep.tsv")
    assert caught.value.code == "E_FIGURE_SET_MISMATCH"


# ---------------------------------------------------------------- V214: UpSet and Venn

def _regions_oracle(rows: list[dict], sets: list[str]) -> dict[tuple, int]:
    """Exclusive region counts by brute-force enumeration of the 2^k - 1 patterns (SM31)."""
    k = len(sets)
    out = {}
    for code in range(1, 2 ** k):
        pattern = tuple((code >> (k - 1 - i)) & 1 for i in range(k))
        out[pattern] = sum(1 for r in rows if tuple(int(r[s]) for s in sets) == pattern)
    return out


def test_v214_upset_for_k3_and_k5_and_venn_for_k3_have_the_oracle_region_counts(tmp_path):
    for k, letters in ((3, ["A", "B", "C"]), (5, ["A", "B", "C", "D", "E"])):
        source = f"diff_upset_venn_k{k}.tsv"
        result = _build(tmp_path / f"upset{k}", "diff_upset_venn", source, sets=letters, venn=False)
        _assert_three_outputs(result)
        rows = _rows(source)
        oracle = _regions_oracle(rows, letters)
        bars = _figure_json(Path(result["outputs"]["html"]))["data"][0]
        assert {p: c for p, c in zip(bars["x"], bars["y"])} == {"".join(map(str, p)): c for p, c in oracle.items()}
        union = sum(1 for r in rows if any(int(r[s]) for s in letters))
        assert sum(bars["y"]) == union
    venn = _build(tmp_path / "venn3", "diff_upset_venn", "diff_upset_venn_k3.tsv", sets=["A", "B", "C"], venn=True, registry_dir=tmp_path / "venn3" / "catalogue")
    _assert_three_outputs(venn)
    rows = _rows("diff_upset_venn_k3.tsv")
    oracle = _regions_oracle(rows, ["A", "B", "C"])
    bars = _figure_json(Path(venn["outputs"]["html"]))["data"][0]
    assert sum(bars["y"]) == sum(oracle.values())


def test_v214_negative_venn_for_k5_fails_e_venn_k_and_the_upset_still_runs(tmp_path):
    with pytest.raises(ConfigurationError) as caught:
        _build(tmp_path / "venn5", "diff_upset_venn", "diff_upset_venn_k5.tsv", sets=["A", "B", "C", "D", "E"], venn=True,
               registry_dir=tmp_path / "venn5" / "catalogue")
    assert caught.value.code == "E_VENN_K"
    upset = _build(tmp_path / "upset5", "diff_upset_venn", "diff_upset_venn_k5.tsv", sets=["A", "B", "C", "D", "E"], venn=False)
    _assert_three_outputs(upset)


def _r_text_labels(source: Path, sets: list[str], venn: bool) -> list[str]:
    """The text labels drawn by the figure, read from ggplot_build of the plot that the renderer draws."""
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "suppressPackageStartupMessages(library(proteomicsCore))\n"
        "src <- utils::read.delim(args[1], stringsAsFactors = FALSE)\n"
        "sets <- strsplit(args[2], ',')[[1]]\n"
        "p <- proteomicsCore:::catalogue_upset_venn(src, sets, venn = args[3] == 'TRUE')\n"
        "cat(jsonlite::toJSON(list(labels = proteomicsCore:::catalogue_text_labels(p)), auto_unbox = FALSE))\n"
    )
    result = run_r_code(code, [source, ",".join(sets), "TRUE" if venn else "FALSE"], cwd=source.parent)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.strip().splitlines()[-1])["labels"]


def test_v214_venn_draws_every_exclusive_region_count_and_the_set_totals(tmp_path):
    from collections import Counter
    sets = ["A", "B", "C"]
    rows = _rows("diff_upset_venn_k3.tsv")
    oracle = _regions_oracle(rows, sets)
    labels = _r_text_labels(FIX / "diff_upset_venn_k3.tsv", sets, venn=True)
    drawn = Counter(labels)
    for pattern, count in oracle.items():
        assert drawn[str(count)] >= sum(1 for c in oracle.values() if c == count), (pattern, count)
    for letter in sets:
        total = sum(int(r[letter]) for r in rows)
        assert f"{letter} (n = {total})" in labels
    # the Venn itself is also written by the builder (k = 3): all three outputs
    venn = _build(tmp_path / "venn3", "diff_upset_venn", "diff_upset_venn_k3.tsv", sets=sets, venn=True, registry_dir=tmp_path / "venn3" / "catalogue")
    _assert_three_outputs(venn)


def test_v214_upset_matrix_rows_carry_the_set_names_and_the_bars_carry_their_counts(tmp_path):
    sets = ["A", "B", "C", "D", "E"]
    rows = _rows("diff_upset_venn_k5.tsv")
    oracle = _regions_oracle(rows, sets)
    labels = _r_text_labels(FIX / "diff_upset_venn_k5.tsv", sets, venn=False)
    nonzero = {p: c for p, c in oracle.items() if c}
    for letter in sets:
        assert letter in labels                                       # matrix row label (the set name)
    for count in nonzero.values():
        assert str(count) in labels                                   # intersection bar label
    assert "0" not in labels                                          # zero-count intersections are not drawn
    # set-size bars: each carries its set name and its oracle size, in the same order as the matrix rows
    expected_size = {letter: sum(int(r[letter]) for r in rows) for letter in sets}
    size_labels = [f"{letter} ({expected_size[letter]})" for letter in sets]
    for label in size_labels:
        assert label in labels
    positions = [labels.index(label) for label in size_labels]
    assert positions == sorted(positions)
    matrix_order = [labels.index(letter) for letter in sets]
    assert matrix_order == sorted(matrix_order)
