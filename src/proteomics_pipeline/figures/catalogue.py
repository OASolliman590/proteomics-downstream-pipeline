"""Figure catalogue, stage 1 (packet R16e; FR-206 to FR-214; V206-V214; ADR 0012).

Each figure is generated from its registered source table (R16a registry): the source is copied into the catalogue, its SHA-256
is recorded, and the three outputs are produced from that copy only:

- the Prism-style PNG at 300 dpi, drawn in R with the R16c theme (proteomicsCore::render_catalogue_png);
- the .pzfx file, with the figure's data tables, the D-75 statistics table where the figure has inferential statistics, and the
  Info sheet (R16b writer);
- the standalone interactive HTML (R16d writer), with the trace data built here from the same source rows.

Refusals use the typed codes of the contract: E_FIGURE_SOURCE_IMPUTED (no observed mask), E_FIGURE_SOURCE_MISSING (a label or
bracket without its source), E_FIGURE_SET_MISMATCH (a feature outside the declared set), E_VENN_K (Venn for k > 3, the UpSet
still runs). Statistics are computed here only where the figure already carries them in its source; nothing is recomputed
from another table.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import statistics as pystats
import tempfile
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from ..errors import ConfigurationError
from ..runtime import run_r_code
from .config import resolve_figure_config
from .plotly_html import write_html
from .pzfx_writer import ColumnTable, column_table, grouped_table, info_sheet, statistics_table, write_pzfx, xy_table
from .registry import FigureRegistry

ROOT = Path(__file__).resolve().parents[3]
CREATED = "2026-10-08T00:00:00Z"
COLOURS = {"up": "#D55E00", "down": "#0072B2", "ns": "#9E9E9E"}

STAGE1 = (
    "qc_intensity_distributions", "qc_missingness_heatmap", "qc_sample_correlation_heatmap", "qc_pca_ellipses",
    "diff_volcano_labelled", "diff_ma_plot", "diff_protein_dotplot_brackets", "diff_dep_heatmap_clustered", "diff_upset_venn",
)


def read_source(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _num(text: str | None) -> float | None:
    if text is None or text in ("", "NA", "NaN"):
        return None
    return float(text)


def _refuse(code: str, message: str) -> None:
    raise ConfigurationError(code, message)


def _observed_cells(rows: Sequence[Mapping[str, str]]) -> None:
    if rows and "observed" not in rows[0]:
        _refuse("E_FIGURE_SOURCE_IMPUTED", "the source has no observed mask; imputed cells cannot be plotted as observed")


def _type7_quantiles(values: Sequence[float]) -> list[float]:
    """Quantiles 0, .25, .5, .75, 1 with the type-7 rule (R's default), as used by the oracle and the figure."""
    ordered = sorted(values)
    out = []
    for p in (0.0, 0.25, 0.5, 0.75, 1.0):
        h = (len(ordered) - 1) * p
        low = math.floor(h)
        high = min(low + 1, len(ordered) - 1)
        out.append(ordered[low] + (h - low) * (ordered[high] - ordered[low]))
    return out


def _pearson(x: Sequence[float], y: Sequence[float]) -> float:
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy)


def _coverage_correlation(values: Mapping[str, Mapping[str, float | None]], samples: Sequence[str], minimum_shared: int = 3) -> dict:
    """Pearson correlation on the cells observed in both samples (the rule of R03 qc_correlations); NA below the minimum."""
    out = {}
    for a in samples:
        for b in samples:
            shared = [f for f in values if values[f].get(a) is not None and values[f].get(b) is not None]
            if len(shared) >= minimum_shared:
                out[(a, b)] = _pearson([values[f][a] for f in shared], [values[f][b] for f in shared])
            else:
                out[(a, b)] = None
    return out


def _cluster_order(matrix: Sequence[Sequence[float]], row_ids: Sequence[str], distance: str = "euclidean", linkage: str = "complete") -> list[str]:
    """Row order from stats::hclust (R) for the declared distance and linkage."""
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "m <- as.matrix(jsonlite::fromJSON(args[1]))\n"
        "hc <- stats::hclust(stats::dist(m, method = args[2]), method = args[3])\n"
        "cat(jsonlite::toJSON(hc$order, auto_unbox = FALSE))\n"
    )
    with tempfile.TemporaryDirectory(prefix="catalogue-cluster-") as tmp:
        payload = Path(tmp) / "matrix.json"
        payload.write_text(json.dumps([list(map(float, r)) for r in matrix]), encoding="utf-8", newline="\n")
        result = run_r_code(code, [payload, distance, linkage], cwd=Path(tmp))
    if result.returncode != 0:
        _refuse("E_FIGURE_OUTPUT_MISSING", f"clustering failed: {result.stderr[-200:]}")
    order = json.loads(result.stdout.strip().splitlines()[-1])
    return [row_ids[i - 1] for i in order]


def _regions(membership: Mapping[str, Sequence[int]], sets: Sequence[str]) -> list[tuple[str, int]]:
    """Exclusive region counts by brute-force enumeration of the 2^k - 1 patterns (SM31, V214 oracle)."""
    k = len(sets)
    out = []
    for code in range(1, 2 ** k):
        pattern = tuple((code >> (k - 1 - i)) & 1 for i in range(k))
        count = sum(1 for flags in membership.values() if tuple(flags) == pattern)
        out.append(("".join(str(b) for b in pattern), count))
    return out


def _info(title: str, statistics: str, test: str, source_sha: str, group_sizes: str, notes: Sequence[str]) -> dict:
    return info_sheet(title, {
        "Test": test, "Multiplicity adjustment": "BH within figure family (as in the source table)" if statistics == "table" else "none (descriptive)",
        "Effect scale": "log2 difference" if statistics == "table" else "not applicable", "Group sizes": group_sizes,
        "Pipeline version and commit": "0.1.0.dev0 (R16e stage 1)", "Source table SHA-256": source_sha, "Run ID": "catalogue-fixture",
        "Created (UTC)": CREATED, "Statistics": statistics,
    }, notes=list(notes) + ["A test run inside Prism may give a different P value from the pipeline's test."])


def _layout(title: str, xtitle: str, ytitle: str) -> dict:
    return {"title": {"text": title}, "xaxis": {"title": {"text": xtitle}, "showgrid": False, "zeroline": False},
            "yaxis": {"title": {"text": ytitle}, "showgrid": False, "zeroline": False}, "plot_bgcolor": "#FFFFFF", "paper_bgcolor": "#FFFFFF"}


def _heatmap_hover(label_x: str, label_y: str, label_z: str) -> str:
    return f"{label_y} %{{y}}<br>{label_x} %{{x}}<br>{label_z} %{{z:.3g}}<extra></extra>"


# ---------------------------------------------------------------- per-figure content

def _content_intensity(rows: list[dict], sha: str, **_: Any) -> dict:
    _observed_cells(rows)
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    observed: dict[str, list[float]] = {s: [] for s in samples}
    traces = []
    for s in samples:
        cells = [r for r in rows if r["sample"] == s]
        vals = [float(r["value"]) for r in cells if r["observed"] == "TRUE"]
        observed[s] = vals
        traces.append({"type": "box", "name": s, "y": vals, "boxpoints": "outliers", "marker": {"color": "#0072B2"},
                       "customdata": [[r["feature"], s, r["observed"], r["value"]] for r in cells if r["observed"] == "TRUE"],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>sample: %{customdata[1]}<br>observed = %{customdata[2]}<br>value %{y:.3f}<extra></extra>"})
    width = max(len(v) for v in observed.values())
    column = column_table("qc intensity distributions (observed cells)", {s: observed[s] + [None] * (width - len(observed[s])) for s in samples})
    return {"tables": [column], "statistics": None, "info": _info("qc intensity distributions", "none", "none (descriptive QC)", sha,
            "; ".join(f"{s} {len(observed[s])}" for s in samples), ["Observed cells only; missing cells are not plotted."]),
            "traces": traces, "layout": _layout("Intensity distributions", "sample", "log2 intensity (observed)"),
            "png": {}, "observed": observed}


def _content_missingness(rows: list[dict], sha: str, **_: Any) -> dict:
    _observed_cells(rows)
    features = list(dict.fromkeys(r["feature"] for r in rows))
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    cell = {(r["feature"], r["sample"]): (None if r["observed"] != "TRUE" else _num(r["value"])) for r in rows}
    fraction = {f: sum(1 for s in samples if cell[(f, s)] is not None) / len(samples) for f in features}
    z = [[cell[(f, s)] for s in samples] for f in features]
    dataset = {s: [cell[(f, s)] for f in features] for s in samples}
    table = grouped_table("qc missingness (observed log2 intensity; empty = missing)", features, dataset)
    trace = {"type": "heatmap", "x": samples, "y": features, "z": z, "colorscale": "Blues",
             "hovertemplate": _heatmap_hover("sample", "feature", "log2 intensity"), "hoverongaps": False}
    return {"tables": [table], "statistics": None, "info": _info("qc missingness", "none", "none (descriptive QC)", sha,
            "; ".join(f"{s}" for s in samples), ["Missing cells are empty in every output, never zero."]),
            "traces": [trace], "layout": _layout("Missingness", "sample", "feature"), "png": {}, "fraction": fraction, "cells": cell}


def _content_correlation(rows: list[dict], sha: str, **_: Any) -> dict:
    _observed_cells(rows)
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    features = list(dict.fromkeys(r["feature"] for r in rows))
    values = {f: {r["sample"]: (_num(r["value"]) if r["observed"] == "TRUE" else None) for r in rows if r["feature"] == f} for f in features}
    corr = _coverage_correlation(values, samples)
    z = [[corr[(a, b)] for b in samples] for a in samples]
    table = grouped_table("qc sample correlations (observed cells, minimum 3 shared)", samples, {b: [corr[(a, b)] for a in samples] for b in samples})
    trace = {"type": "heatmap", "x": samples, "y": samples, "z": z, "zmin": -1, "zmax": 1, "colorscale": "RdBu", "reversescale": True,
             "hovertemplate": _heatmap_hover("sample", "sample", "Pearson r"), "hoverongaps": False}
    return {"tables": [table], "statistics": None, "info": _info("qc sample correlations", "none", "none (descriptive QC)", sha,
            "; ".join(samples), ["Pearson r on cells observed in both samples; NA when fewer than 3 shared cells."]),
            "traces": [trace], "layout": _layout("Sample correlations", "sample", "sample"), "png": {}, "corr": corr}


def _ellipse(xs: Sequence[float], ys: Sequence[float], level_radius: float = math.sqrt(5.991464547107979), n: int = 100) -> list[tuple[float, float]]:
    """The 95 per cent ellipse of a group from its covariance (chi-square, 2 df): radius sqrt(qchisq(0.95, 2))."""
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs) / (len(xs) - 1)
    syy = sum((y - my) ** 2 for y in ys) / (len(ys) - 1)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (len(xs) - 1)
    tr, det = sxx + syy, sxx * syy - sxy * sxy
    disc = math.sqrt(max(tr * tr / 4 - det, 0.0))
    l1, l2 = tr / 2 + disc, tr / 2 - disc
    # eigenvector for l1
    if abs(sxy) > 1e-300:
        v1 = (l1 - syy, sxy)
    else:
        v1 = (1.0, 0.0) if sxx >= syy else (0.0, 1.0)
    norm = math.hypot(*v1)
    v1 = (v1[0] / norm, v1[1] / norm)
    v2 = (-v1[1], v1[0])
    a, b = level_radius * math.sqrt(max(l1, 0.0)), level_radius * math.sqrt(max(l2, 0.0))
    points = []
    for i in range(n):
        t = 2 * math.pi * i / (n - 1)
        points.append((mx + a * math.cos(t) * v1[0] + b * math.sin(t) * v2[0], my + a * math.cos(t) * v1[1] + b * math.sin(t) * v2[1]))
    return points


def _content_pca(rows: list[dict], sha: str, variant: bool = False, **_: Any) -> dict:
    groups = list(dict.fromkeys(r["group"] for r in rows))
    counts = {g: sum(1 for r in rows if r["group"] == g) for g in groups}
    small = [g for g in groups if counts[g] < 3]
    table = xy_table("qc pca scores (PC1 against PC2)", "PC1", "PC2", [float(r["PC1"]) for r in rows], [float(r["PC2"]) for r in rows])
    traces, ellipses = [], {}
    for g in groups:
        members = [r for r in rows if r["group"] == g]
        traces.append({"type": "scatter", "mode": "markers", "name": g, "x": [float(r["PC1"]) for r in members], "y": [float(r["PC2"]) for r in members],
                       "customdata": [[r["sample"], g, r["PC1"], r["PC2"]] for r in members],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>group: %{customdata[1]}<br>PC1 %{customdata[2]}<br>PC2 %{customdata[3]}<extra></extra>",
                       "marker": {"size": 8}})
        if g not in small:
            pts = _ellipse([float(r["PC1"]) for r in members], [float(r["PC2"]) for r in members])
            ellipses[g] = pts
            traces.append({"type": "scatter", "mode": "lines", "name": f"{g} 95% ellipse", "x": [p[0] for p in pts], "y": [p[1] for p in pts],
                           "customdata": [[f"{g} 95% ellipse", g, "ellipse", "95%"] for _ in pts],
                           "hovertemplate": "<b>%{customdata[0]}</b><br>group: %{customdata[1]}<br>%{customdata[2]} %{customdata[3]}<extra></extra>", "showlegend": False})
    note = "; ".join(f"no ellipse for {g} (fewer than three samples)" for g in small)
    return {"tables": [table], "statistics": None, "info": _info("qc pca", "none", "none (descriptive QC)", sha, "; ".join(f"{g} {counts[g]}" for g in groups),
            ["PCA scores of the display matrix (missing cells filled with the feature median, then centred; SM06)."] + ([note] if note else [])),
            "traces": traces, "layout": _layout("PCA with 95% ellipses", "PC1", "PC2"), "png": {"variant": variant}, "ellipses": ellipses, "small": small}


def _content_volcano(rows: list[dict], sha: str, labels: Sequence[str] | None = None, q_cutoff: float = 0.05, **_: Any) -> dict:
    if labels is not None:
        present = {r["feature"] for r in rows}
        absent = [x for x in labels if x not in present]
        if absent:
            _refuse("E_FIGURE_SOURCE_MISSING", f"label(s) {absent} are not in the source table")
        if len(labels) > 10:
            _refuse("E_FIGURE_SET_MISMATCH", "at most 10 labels per volcano (D-74)")
    classes = {}
    for r in rows:
        q, e = float(r["q"]), float(r["effect"])
        classes[r["feature"]] = "up" if q < q_cutoff and e > 0 else ("down" if q < q_cutoff and e < 0 else "ns")
    traces = []
    for cls in ("up", "down", "ns"):
        members = [r for r in rows if classes[r["feature"]] == cls]
        traces.append({"type": "scatter", "mode": "markers", "name": cls, "x": [float(r["effect"]) for r in members],
                       "y": [-math.log10(float(r["P"])) for r in members], "marker": {"color": COLOURS[cls], "size": 7, "opacity": 0.85},
                       "customdata": [[r["feature"], cls, r["P"], r["q"]] for r in members],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>class: %{customdata[1]}<br>P = %{customdata[2]}<br>q = %{customdata[3]}<extra></extra>"})
    stats = [{"row": r["feature"], "effect": float(r["effect"]), "P": float(r["P"]), "q": float(r["q"]),
              "n": {"Control": float(r["n_Control"]), "Acute": float(r["n_Acute"])}} for r in rows]
    tables = [xy_table("volcano", "log2 fold change", "-log10 P", [float(r["effect"]) for r in rows], [-math.log10(float(r["P"])) for r in rows]),
              statistics_table("volcano - pipeline statistics", stats, ["Control", "Acute"])]
    return {"tables": tables, "statistics": stats, "info": _info("volcano (labelled)", "table", "Welch t (synthetic fixture)", sha, "Control 6; Acute 6",
            [f"Family q cutoff {q_cutoff}; classes up/down/ns by q and the sign of the effect."]),
            "traces": traces, "layout": _layout("Volcano (labelled)", "log2 fold change", "-log10 P"), "png": {"labels": list(labels) if labels is not None else None, "q_cutoff": q_cutoff},
            "classes": classes}


def _content_ma(rows: list[dict], sha: str, **_: Any) -> dict:
    traces, tables, stats = [], [], []
    for contrast in dict.fromkeys(r["contrast"] for r in rows):
        members = [r for r in rows if r["contrast"] == contrast]
        a = [(float(r["mean_numerator"]) + float(r["mean_denominator"])) / 2 for r in members]
        m = [float(r["mean_numerator"]) - float(r["mean_denominator"]) for r in members]
        traces.append({"type": "scatter", "mode": "markers", "name": contrast, "x": a, "y": m,
                       "customdata": [[r["feature"], contrast, r["P"], r["q"]] for r in members],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>contrast: %{customdata[1]}<br>P = %{customdata[2]}<br>q = %{customdata[3]}<extra></extra>"})
        tables.append(xy_table(f"MA {contrast}", "A (mean of the two groups)", "M (numerator - denominator)", a, m))
        stats.extend({"row": f"{contrast}: {r['feature']}", "effect": float(r["mean_numerator"]) - float(r["mean_denominator"]), "P": float(r["P"]),
                      "q": float(r["q"]), "n": {"numerator": float(r["n_numerator"]), "denominator": float(r["n_denominator"])}} for r in members)
    tables.append(statistics_table("MA - pipeline statistics", stats, ["numerator", "denominator"]))
    return {"tables": tables, "statistics": stats, "info": _info("MA plot", "table", "Welch t (synthetic fixture)", sha, "as in the source",
            ["A and M from the group means; M follows the declared numerator."]), "traces": traces, "layout": _layout("MA plot", "A", "M"), "png": {}}


def _content_dotplot(rows: list[dict], sha: str, family: list[dict] | None = None, brackets: bool = True, **_: Any) -> dict:
    if family is None:
        _refuse("E_FIGURE_SOURCE_MISSING", "the bracket family table is missing; no P bracket is drawn")
    proteins = list(dict.fromkeys(r["protein"] for r in rows))
    groups = list(dict.fromkeys(r["group"] for r in rows))
    tables, traces = [], []
    for protein in proteins:
        cols = {g: [float(r["value"]) for r in rows if r["protein"] == protein and r["group"] == g] for g in groups}
        tables.append(column_table(protein, cols))
        traces.append({"type": "box", "name": protein, "x": [g for g in groups for _ in cols[g]], "y": [v for g in groups for v in cols[g]],
                       "boxpoints": "all", "jitter": 0.3, "pointpos": 0, "marker": {"size": 7},
                       "customdata": [[protein, r["group"], r["replicate"], r["value"]] for r in rows if r["protein"] == protein],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>group: %{customdata[1]}<br>replicate: %{customdata[2]}<br>value %{y:.3f}<extra></extra>"})
    stats = [{"row": f"{f['protein']} {f['comparison']}", "effect": float(f["effect"]), "P": float(f["P"]), "q": float(f["q"]),
              "n": {g: float(f[f"n_{g}"]) for g in groups}} for f in family]
    tables.append(statistics_table("dot plot - pipeline statistics", stats, groups))
    sizes = "; ".join(f"{g} {sum(1 for r in rows if r['group'] == g and r['protein'] == proteins[0])}" for g in groups)
    return {"tables": tables, "statistics": stats, "info": _info("per-protein dot plots with P brackets", "table", "Welch t (synthetic fixture)", sha,
            sizes, ["Brackets show the family P values with the FR-198 annotation."]),
            "traces": traces, "layout": _layout("Per-protein values", "group", "log2 abundance"), "png": {"family": family, "groups": groups}}


def _content_dep_heatmap(rows: list[dict], sha: str, **_: Any) -> dict:
    for r in rows:
        if r["declared"] != "TRUE":
            _refuse("E_FIGURE_SET_MISMATCH", f"feature {r['feature']} is outside the declared DEP set")
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    features = list(dict.fromkeys(r["feature"] for r in rows))
    values = {f: [float(r["value"]) for s in samples for r in rows if r["feature"] == f and r["sample"] == s] for f in features}
    scaled = {}
    for f in features:
        v = values[f]
        mu, sd = pystats.mean(v), pystats.stdev(v)
        scaled[f] = [(x - mu) / sd for x in v]
    order = _cluster_order([scaled[f] for f in features], features, "euclidean", "complete")
    table = grouped_table("clustered DEP heatmap (row z-scores)", order, {s: [scaled[f][i] for f in order] for i, s in enumerate(samples)})
    trace = {"type": "heatmap", "x": samples, "y": order, "z": [scaled[f] for f in order], "colorscale": "RdBu", "reversescale": True,
             "hovertemplate": _heatmap_hover("sample", "feature", "row z"), "hoverongaps": False}
    return {"tables": [table], "statistics": None, "info": _info("clustered DEP heatmap", "none", "none (descriptive)", sha, "; ".join(samples),
            ["Euclidean distance on row z-scores, complete linkage (stats::hclust)."]), "traces": [trace], "layout": _layout("Clustered DEP heatmap", "sample", "feature"),
            "png": {}, "order": order, "scaled": scaled}


def _content_upset_venn(rows: list[dict], sha: str, sets: Sequence[str] = (), **_: Any) -> dict:
    membership = {r["feature"]: [int(r[s]) for s in sets] for r in rows}
    regions = _regions(membership, sets)
    table = ColumnTable("upset regions (exclusive)", (("count", tuple(float(c) for _, c in regions)),), row_titles=tuple(p for p, _ in regions))
    trace = {"type": "bar", "x": [p for p, _ in regions], "y": [c for _, c in regions],
             "customdata": [[p, c, ",".join(s for s, b in zip(sets, p) if b == "1") or "none"] for p, c in regions],
             "hovertemplate": "<b>%{customdata[0]}</b><br>sets: %{customdata[2]}<br>count %{customdata[1]}<extra></extra>"}
    return {"tables": [table], "statistics": None, "info": _info("UpSet", "none", "none (descriptive set overlap)", sha, "; ".join(sets),
            ["Exclusive regions of the declared sets; a Venn is drawn only for k <= 3 (E_VENN_K)."]), "traces": [trace],
            "layout": _layout("UpSet (exclusive intersections)", "pattern (sets A..)", "count"), "png": {"set_names": list(sets)}, "regions": regions}


_BUILDERS: dict[str, Callable[..., dict]] = {
    "qc_intensity_distributions": _content_intensity,
    "qc_missingness_heatmap": _content_missingness,
    "qc_sample_correlation_heatmap": _content_correlation,
    "qc_pca_ellipses": _content_pca,
    "diff_volcano_labelled": _content_volcano,
    "diff_ma_plot": _content_ma,
    "diff_protein_dotplot_brackets": _content_dotplot,
    "diff_dep_heatmap_clustered": _content_dep_heatmap,
    "diff_upset_venn": _content_upset_venn,
}


def _render_png(figure_id: str, source_copy: Path, out: Path, options: Mapping[str, Any], width_mm: int = 120) -> None:
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "suppressPackageStartupMessages(library(proteomicsCore))\n"
        "opts <- jsonlite::fromJSON(args[4], simplifyVector = TRUE)\n"
        "r <- render_catalogue_png(args[1], args[2], args[3], opts, width_mm = as.numeric(args[5]))\n"
        "cat(jsonlite::toJSON(list(width_px = r$width_px, height_px = r$height_px), auto_unbox = TRUE))\n"
    )
    with tempfile.TemporaryDirectory(prefix="catalogue-png-") as tmp:
        options_file = Path(tmp) / "options.json"
        options_file.write_text(json.dumps(dict(options)), encoding="utf-8", newline="\n")
        result = run_r_code(code, [figure_id, source_copy, out, options_file, width_mm], cwd=Path(tmp))
    if result.returncode != 0:
        match = re.search(r"(E_[A-Z_]+):\s*([^\n]*)", result.stderr)
        if match:
            raise ConfigurationError(match.group(1), match.group(2).strip())
        raise ConfigurationError("E_FIGURE_OUTPUT_MISSING", f"the PNG of {figure_id} was not written: {result.stderr.strip()[-200:]}")


def _dotplot_brackets(family: Sequence[Mapping[str, Any]], rows: Sequence[Mapping[str, str]], brackets: Sequence[Mapping[str, Any]] | None = None) -> list[dict]:
    """FR-198 labels of the family P values, one bracket per family comparison, stacked above each protein's data.

    Explicit ``brackets`` (protein, group1, group2, P) must each match a family row; a bracket whose P is not in the family
    table is refused (E_FIGURE_SOURCE_MISSING) and no bracket is drawn for it."""
    from .style import p_exact
    if brackets is not None:
        for b in brackets:
            match = [f for f in family if f["protein"] == b["protein"] and f["group_a"] == b["group1"] and f["group_b"] == b["group2"]
                     and math.isclose(float(f["P"]), float(b["P"]), rel_tol=0, abs_tol=0)]
            if not match:
                _refuse("E_FIGURE_SOURCE_MISSING", f"bracket {b['protein']} {b['group1']}-{b['group2']} with P {b['P']} is not in the family table")
    out = []
    for protein in dict.fromkeys(f["protein"] for f in family):
        top = max(float(r["value"]) for r in rows if r["protein"] == protein)
        step = 0.12 * (top - min(float(r["value"]) for r in rows if r["protein"] == protein))
        members = [f for f in family if f["protein"] == protein]
        for index, f in enumerate(members):
            out.append({"protein": protein, "group1": f["group_a"], "group2": f["group_b"], "label": p_exact(float(f["P"])),
                        "y.position": top + step * (index + 1)})
    return out


def build_figure(run_dir: str | Path, figure_id: str, source: str | Path, *, registry_dir: str | Path | None = None, **options: Any) -> dict:
    """Produce the three outputs of a stage-1 figure from its registered source copy. Returns the outputs and the refusals.

    Options: ``labels`` and ``q_cutoff`` (volcano), ``variant`` (PCA), ``family`` (dot plot brackets), ``sets`` and ``venn``
    (UpSet and Venn; ``venn=True`` with more than three sets raises E_VENN_K before any Venn output is written).
    """
    if figure_id not in _BUILDERS:
        _refuse("E_FIGURE_SOURCE_MISSING", f"{figure_id} is not a stage-1 catalogue figure")
    if figure_id == "diff_upset_venn" and options.get("venn") and len(options.get("sets", ())) > 3:
        _refuse("E_VENN_K", f"a circle Venn needs k <= 3 sets; {len(options['sets'])} were given (UpSet is still produced)")
    style = resolve_figure_config(None)
    catalogue_dir = Path(registry_dir) if registry_dir is not None else Path(run_dir) / "report" / "figure_catalogue"
    registry = FigureRegistry(catalogue_dir, style)
    registry.register(figure_id, source)
    copy = catalogue_dir / "sources" / f"{figure_id}.tsv"
    rows = read_source(copy)
    sha = hashlib.sha256(copy.read_bytes()).hexdigest()
    builder_options = {k: v for k, v in options.items() if k not in ("sets", "venn")}
    if figure_id == "diff_upset_venn":
        builder_options["sets"] = list(options.get("sets", ()))
    if figure_id == "diff_protein_dotplot_brackets":
        family_rows = [dict(f) for f in options.get("family", [])]
        builder_options["family"] = family_rows
    content = _BUILDERS[figure_id](rows, sha, **builder_options)
    png_options: dict[str, Any] = {}
    if figure_id == "diff_volcano_labelled":
        png_options = {"labels": content["png"]["labels"], "q_cutoff": content["png"]["q_cutoff"]}
    elif figure_id == "qc_pca_ellipses":
        png_options = {"variant": bool(options.get("variant"))}
    elif figure_id == "diff_protein_dotplot_brackets":
        png_options = {"brackets": _dotplot_brackets(content["png"]["family"], rows, options.get("brackets"))}
    elif figure_id == "diff_upset_venn":
        png_options = {"set_names": list(options.get("sets", ())), "venn": bool(options.get("venn"))}
    html_title = figure_id.replace("_", " ")

    def render_png(source_path: Path, out_path: Path) -> None:
        _render_png(figure_id, source_path, out_path, png_options)

    def render_pzfx(source_path: Path, out_path: Path) -> None:
        write_pzfx(out_path, content["tables"], content["info"])

    def render_html(source_path: Path, out_path: Path) -> None:
        write_html(out_path, html_title, content["traces"], content["layout"])

    outputs: dict[str, str] = {}
    refusals: list[dict[str, str]] = []
    for kind, render in (("png", render_png), ("pzfx", render_pzfx), ("html", render_html)):
        try:
            registry.produce(figure_id, kind, render)
            outputs[kind] = str(catalogue_dir / f"{figure_id}.{kind}")
        except ConfigurationError as error:
            registry.record_missing_output(figure_id, kind, error.message)
            refusals.append({"kind": kind, "code": error.code, "reason": error.message})
    registry.write()
    return {"figure_id": figure_id, "outputs": outputs, "refusals": refusals, "content": content, "registry": registry}
