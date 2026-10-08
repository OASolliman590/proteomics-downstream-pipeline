"""Figure catalogue, stages 1 and 2 (packet R16e; FR-206 to FR-222; V206-V222; ADR 0012).

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
STAGE2 = (
    "bio_roc_ci_band", "bio_panel_forest", "bio_cv_auc_vs_permutation", "path_enrichment_dotplot",
    "net_module_heatmap", "net_network_graph", "perm_pcoa_centroids", "perm_permdisp",
)
# Claim labels of the biomarker results (post_de_biomarker.R). The cross-validated label is the only one a ROC or a
# cross-validated AUC may carry; a panel plot may carry any frozen label (E_BIOMARKER_CLAIM otherwise).
CV_CLAIM = "cross_validated_nested"
PANEL_CLAIMS = ("in_sample", "cross_validated_nested", "fixed_panel_cv", "independently_validated")
CLAIM_NOTE = "claim label {claim}: performance estimate from the pipeline, not a diagnostic or clinical-utility claim"


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


def _info(title: str, statistics: str, test: str, source_sha: str, group_sizes: str, notes: Sequence[str], *,
          multiplicity: str | None = None, effect_scale: str | None = None, stage: str = "R16e stage 1") -> dict:
    return info_sheet(title, {
        "Test": test,
        "Multiplicity adjustment": multiplicity or ("BH within figure family (as in the source table)" if statistics == "table" else "none (descriptive)"),
        "Effect scale": effect_scale or ("log2 difference" if statistics == "table" else "not applicable"), "Group sizes": group_sizes,
        "Pipeline version and commit": f"0.1.0.dev0 ({stage})", "Source table SHA-256": source_sha, "Run ID": "catalogue-fixture",
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
    """A binary observed / missing map: state 1 = observed, 0 = missing (never an intensity). Features are ordered by their
    missing fraction (most missing first); the per-sample missing percentage is written in the column label."""
    _observed_cells(rows)
    features = list(dict.fromkeys(r["feature"] for r in rows))
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    state = {(r["feature"], r["sample"]): (1 if r["observed"] == "TRUE" else 0) for r in rows}
    fraction = {f: sum(state[(f, s)] for s in samples) / len(samples) for f in features}
    missing_pct = {s: 100.0 * sum(1 - state[(f, s)] for f in features) / len(features) for s in samples}
    ordered = sorted(features, key=lambda f: (-(1 - fraction[f]), f))
    z = [[state[(f, s)] for s in samples] for f in ordered]
    table = grouped_table("qc missingness (1 = observed, 0 = missing)", ordered, {s: [state[(f, s)] for f in ordered] for s in samples})
    trace = {"type": "heatmap", "x": [f"{s} {missing_pct[s]:.0f}%" for s in samples], "y": ordered, "z": z,
             "zmin": 0, "zmax": 1, "showscale": False, "colorscale": [[0.0, "#000000"], [0.5, "#000000"], [0.5, "#D9D9D9"], [1.0, "#D9D9D9"]],
             "hovertemplate": _heatmap_hover("sample", "feature", "state (1 observed, 0 missing)"),
             "hoverongaps": False}
    return {"tables": [table], "statistics": None, "info": _info("qc missingness", "none", "none (descriptive QC)", sha,
            "; ".join(samples), ["Binary observed / missing state of every cell; features ordered by missing fraction."]),
            "traces": [trace], "layout": _layout("Missingness", "sample (missing %)", "feature (ordered by missing fraction)"),
            "png": {}, "fraction": fraction, "state": state, "ordered": ordered, "missing_pct": missing_pct, "samples": samples}


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
    pct1, pct2 = float(rows[0]["PC1_var_pct"]), float(rows[0]["PC2_var_pct"])
    x_title, y_title = f"PC1 ({pct1:.1f}%)", f"PC2 ({pct2:.1f}%)"
    table = xy_table("qc pca scores (PC1 against PC2)", x_title, y_title, [float(r["PC1"]) for r in rows], [float(r["PC2"]) for r in rows])
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
            "traces": traces, "layout": _layout("PCA with 95% ellipses", x_title, y_title), "png": {"variant": variant}, "ellipses": ellipses, "small": small,
            "variance": (pct1, pct2)}


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


# ---------------------------------------------------------------- stage 2 (V215-V222)

def _claim_check(rows: Sequence[Mapping[str, str]], allowed: Sequence[str]) -> str:
    """The one claim label of the figure, which must be frozen (allowed). E_BIOMARKER_CLAIM otherwise."""
    labels = {r.get("claim_label", "") for r in rows}
    if len(labels) != 1 or next(iter(labels)) not in allowed:
        _refuse("E_BIOMARKER_CLAIM", f"the claim label is {sorted(labels)!r}; the figure needs one of {list(allowed)}")
    return next(iter(labels))


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values)


def _content_roc(rows: list[dict], sha: str, **_: Any) -> dict:
    """V215: ROC curve and stratified bootstrap band (pROC, seed recorded). Only out-of-fold scores with the cross-validated
    label may be drawn (a training-fold ROC is refused as E_BIOMARKER_CLAIM)."""
    claim = _claim_check(rows, (CV_CLAIM,))
    if any(r["score_source"] != "out_of_fold" for r in rows):
        _refuse("E_BIOMARKER_CLAIM", "the ROC is drawn from training-fold scores; only out-of-fold scores carry the cross-validated claim")
    head = rows[0]
    curve = [r for r in rows if r["component"] == "curve"]
    band = [r for r in rows if r["component"] == "band"]
    auc, lo, hi = float(head["auc"]), float(head["auc_lower"]), float(head["auc_upper"])
    seed, boots = int(head["seed"]), int(head["B"])
    n_control, n_case = int(head["n_control"]), int(head["n_case"])
    cfpr = [float(r["fpr"]) for r in curve]
    csens = [float(r["sensitivity"]) for r in curve]
    bfpr = [float(r["fpr"]) for r in band]
    blo = [float(r["band_lower"]) for r in band]
    bhi = [float(r["band_upper"]) for r in band]
    hover = "<b>%{customdata[0]}</b><br>claim: %{customdata[1]}<br>%{customdata[2]}<extra></extra>"
    traces = [
        {"type": "scatter", "mode": "lines", "name": "95% CI band, lower", "x": bfpr, "y": blo, "line": {"width": 0}, "showlegend": False,
         "customdata": [["95% CI band (stratified bootstrap)", claim, f"lower sensitivity {v:.3f}; seed {seed}, B {boots}"] for v in blo],
         "hovertemplate": hover},
        {"type": "scatter", "mode": "lines", "name": "95% CI band", "x": bfpr, "y": bhi, "line": {"width": 0}, "fill": "tonexty",
         "fillcolor": "rgba(213,94,0,0.18)",
         "customdata": [["95% CI band (stratified bootstrap)", claim, f"upper sensitivity {v:.3f}; seed {seed}, B {boots}"] for v in bhi],
         "hovertemplate": hover},
        {"type": "scatter", "mode": "lines", "name": "chance", "x": [0.0, 1.0], "y": [0.0, 1.0], "line": {"dash": "dash", "color": "#9E9E9E", "width": 1},
         "customdata": [["chance diagonal", claim, "no discrimination"]] * 2, "hovertemplate": hover},
        {"type": "scatter", "mode": "lines", "name": "ROC (out-of-fold)", "x": cfpr, "y": csens, "line": {"shape": "hv", "color": "#D55E00", "width": 2.5},
         "customdata": [["ROC (out-of-fold)", claim, f"AUC {auc:.3f} (95% CI {lo:.3f}-{hi:.3f})"] for _ in cfpr],
         "hovertemplate": "<b>%{customdata[0]}</b><br>claim: %{customdata[1]}<br>1 - specificity %{x:.3f}"
                          "<br>sensitivity %{y:.3f}<br>%{customdata[2]}<extra></extra>"},
    ]
    tables = [xy_table("ROC curve (out-of-fold)", "1 - specificity", "sensitivity", cfpr, csens),
              xy_table("ROC CI band, lower", "1 - specificity", "sensitivity (lower)", bfpr, blo),
              xy_table("ROC CI band, upper", "1 - specificity", "sensitivity (upper)", bfpr, bhi)]
    notes = [f"AUC {auc:.3f} (95% CI {lo:.3f}-{hi:.3f}, DeLong); band from the stratified bootstrap with seed {seed} and B = {boots}.",
             CLAIM_NOTE.format(claim=claim)]
    return {"tables": tables, "statistics": None, "info": _info("ROC with CI band", "none", "pROC ROC; DeLong AUC CI; stratified bootstrap band",
            sha, f"Control {n_control}; Case {n_case}", notes, multiplicity="none (performance estimate; no multiplicity family)",
            effect_scale="not applicable", stage="R16e stage 2"),
            "traces": traces, "layout": _layout(f"ROC with 95% CI band ({claim})", "1 - specificity", "sensitivity"), "png": {}}


def _content_forest(rows: list[dict], sha: str, **_: Any) -> dict:
    """V216: panel effects with their 95 per cent intervals (Welch, per feature). The claim label must be frozen."""
    claim = _claim_check(rows, PANEL_CLAIMS)
    features = [r["feature"] for r in rows]
    estimate = [float(r["estimate"]) for r in rows]
    lo = [float(r["ci_lower"]) for r in rows]
    hi = [float(r["ci_upper"]) for r in rows]
    p_values = [float(r["P"]) for r in rows]
    q_values = [float(r["q"]) for r in rows]
    n_control, n_case = int(rows[0]["n_Control"]), int(rows[0]["n_Case"])
    table = ColumnTable("panel forest (effect and 95% CI)", (("estimate", tuple(estimate)), ("ci_lower", tuple(lo)), ("ci_upper", tuple(hi))),
                        row_titles=tuple(features))
    stats = [{"row": f"{f}: Case - Control", "effect": e, "P": p, "q": q, "n": {"Control": n_control, "Case": n_case}}
             for f, e, p, q in zip(features, estimate, p_values, q_values)]
    tables = [table, statistics_table("panel forest - pipeline statistics", stats, ["Control", "Case"])]
    traces = [
        {"type": "scatter", "mode": "markers", "name": "panel effect", "x": estimate, "y": features,
         "error_x": {"type": "data", "symmetric": False, "array": [h - e for h, e in zip(hi, estimate)],
                     "arrayminus": [e - l for l, e in zip(lo, estimate)], "color": "#000000"},
         "marker": {"color": "#0072B2", "size": 8},
         "customdata": [[f, claim, f"P = {p:.3g}", f"q = {q:.3g}"] for f, p, q in zip(features, p_values, q_values)],
         "hovertemplate": "<b>%{customdata[0]}</b><br>claim: %{customdata[1]}<br>effect %{x:.3f}<br>%{customdata[2]}"
                          "<br>%{customdata[3]}<extra></extra>"},
        {"type": "scatter", "mode": "lines", "name": "no difference", "x": [0.0, 0.0], "y": [features[0], features[-1]],
         "line": {"dash": "dash", "color": "#9E9E9E", "width": 1}, "showlegend": False,
         "customdata": [["no difference", claim, "effect 0", "zero line"]] * 2,
         "hovertemplate": "<b>%{customdata[0]}</b><br>claim: %{customdata[1]}<br>%{customdata[2]}<extra></extra>"},
    ]
    notes = ["Effects are Case minus Control on the log2 scale (Welch t, 95% CI); q is BH across the panel.", CLAIM_NOTE.format(claim=claim)]
    return {"tables": tables, "statistics": stats, "info": _info("panel forest", "table", "Welch t per panel feature (95% CI from t.test)", sha,
            f"Control {n_control}; Case {n_case}", notes, multiplicity="BH within the panel family (q)", effect_scale="log2 difference", stage="R16e stage 2"),
            "traces": traces, "layout": _layout(f"Panel forest plot ({claim})", "log2 difference (Case - Control), 95% CI", "panel feature"), "png": {}}


def _content_cv_permutation(rows: list[dict], sha: str, **_: Any) -> dict:
    """V217: the observed cross-validated AUC (mean of the nested-CV repeats) against the permutation null; the displayed P
    must equal (k + 1) / (B + 1), otherwise E_FIGURE_PERMUTATION_P."""
    claim = _claim_check(rows, (CV_CLAIM,))
    seed, boots = int(rows[0]["seed"]), int(rows[0]["B"])
    cv = [float(r["value"]) for r in rows if r["series"] == "cv_repeat"]
    null = [float(r["value"]) for r in rows if r["series"] == "null"]
    if len(null) != boots:
        _refuse("E_FIGURE_SOURCE_MISSING", f"the permutation null has {len(null)} values for B = {boots}")
    observed = _mean(cv)
    k = sum(1 for v in null if v >= observed)
    p_value = (k + 1) / (boots + 1)
    shown = {float(r["P_display"]) for r in rows}
    if len(shown) != 1 or not math.isclose(next(iter(shown)), p_value, rel_tol=0, abs_tol=1e-12):
        _refuse("E_FIGURE_PERMUTATION_P", f"the displayed permutation P {sorted(shown)} differs from (k + 1) / (B + 1) = {p_value:.6g} "
                f"with k = {k} and B = {boots}")
    table = column_table("CV AUC against permutation null", {"cv_repeat": cv + [None] * (boots - len(cv)), "null": null})
    stats = [{"row": "observed CV AUC (mean of repeats) against permutation null", "effect": observed - _mean(null), "P": p_value, "q": None,
              "n": {"cv_repeat": len(cv), "null": boots}}]
    tables = [table, statistics_table("CV AUC - pipeline statistics", stats, ["cv_repeat", "null"])]
    hover = "<b>%{customdata[0]}</b><br>claim: %{customdata[1]}<br>AUC %{customdata[2]:.3f}<br>B = %{customdata[3]}<extra></extra>"
    traces = [
        {"type": "scatter", "mode": "markers", "name": "permutation null", "x": ["permutation null"] * len(null), "y": null,
         "marker": {"color": "#9E9E9E", "size": 7, "opacity": 0.85},
         "customdata": [["permutation null", claim, v, boots] for v in null], "hovertemplate": hover},
        {"type": "scatter", "mode": "markers", "name": "nested-CV repeat", "x": ["observed CV"] * len(cv), "y": cv,
         "marker": {"color": "#0072B2", "size": 8},
         "customdata": [["nested-CV repeat", claim, v, boots] for v in cv], "hovertemplate": hover},
        {"type": "scatter", "mode": "markers", "name": "observed (mean of repeats)", "x": ["observed CV"], "y": [observed],
         "marker": {"color": "#000000", "size": 12, "symbol": "diamond"},
         "customdata": [["observed (mean of repeats)", claim, observed, boots]], "hovertemplate": hover},
    ]
    notes = [f"P = (k + 1) / (B + 1) = ({k} + 1) / ({boots} + 1) = {p_value:.4g}, one-sided (null AUC at or above the observed mean); seed {seed}.",
             CLAIM_NOTE.format(claim=claim)]
    return {"tables": tables, "statistics": stats, "info": _info("CV AUC against permutation null", "table",
            "permutation test of the observed cross-validated AUC", sha, f"cv_repeat {len(cv)}; null {boots}", notes,
            multiplicity="none (one observed statistic per figure)", effect_scale="AUC difference from the null mean", stage="R16e stage 2"),
            "traces": traces, "layout": _layout(f"CV AUC against permutation null (permutation P = {p_value:.3f}, B = {boots})", "", "AUC"),
            "png": {}}


def _content_enrichment(rows: list[dict], sha: str, show_sets: Sequence[str] | None = None, **_: Any) -> dict:
    """V218: enrichment dot plot. Gene ratio k / n and -log10 q of the sets with at least one foreground hit. A set with no
    foreground hit is never plotted as enriched (E_FIGURE_SET_MISMATCH when it is requested)."""
    by_id = {r["set_id"]: r for r in rows}
    if show_sets is None:
        shown = [r["set_id"] for r in rows if int(r["k"]) >= 1]
    else:
        unknown = [s for s in show_sets if s not in by_id]
        if unknown:
            _refuse("E_FIGURE_SOURCE_MISSING", f"set(s) {unknown} are not in the ORA source")
        empty = [s for s in show_sets if int(by_id[s]["k"]) < 1]
        if empty:
            _refuse("E_FIGURE_SET_MISMATCH", f"set(s) {empty} have no foreground hit; they are not enriched and are not plotted")
        shown = list(show_sets)
    plotted = [by_id[s] for s in shown]
    gene_ratio = [int(r["k"]) / int(r["n"]) for r in plotted]
    neglog_q = [-math.log10(float(r["q_value"])) for r in plotted]
    hits = [int(r["k"]) for r in plotted]
    n_fg, universe = int(rows[0]["n"]), int(rows[0]["N"])
    table = ColumnTable("enrichment (gene ratio and -log10 q)", (("gene_ratio", tuple(gene_ratio)), ("neglog10_q", tuple(neglog_q))),
                        row_titles=tuple(shown))
    stats = [{"row": r["set_id"], "effect": g, "P": float(r["p_value"]), "q": float(r["q_value"]),
              "n": {"foreground": int(r["n"]), "set": int(r["K"])}} for r, g in zip(plotted, gene_ratio)]
    traces = [{"type": "scatter", "mode": "markers", "name": "enriched set (foreground hit)", "x": gene_ratio, "y": shown,
               "marker": {"size": [4 + 2 * h for h in hits], "color": neglog_q, "colorscale": "Reds", "showscale": False,
                          "line": {"color": "#000000", "width": 0.5}},
               "customdata": [[r["set_id"], f"{int(r['k'])} of {int(r['n'])} foreground genes", float(r["q_value"]), int(r["K"])] for r in plotted],
               "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}<br>gene ratio %{x:.3f}<br>q = %{customdata[2]:.3g}"
                                "<br>set size %{customdata[3]}<extra></extra>"}]
    notes = [f"Gene ratio = k / n (foreground hits over foreground size {n_fg}); universe N = {universe}; q is BH across the ORA family.",
             "Sets without a foreground hit are not plotted as enriched."]
    return {"tables": [table, statistics_table("enrichment - pipeline statistics", stats, ["foreground", "set"])], "statistics": stats,
            "info": _info("enrichment dot plot", "table", "one-sided hypergeometric over-representation (ORA)", sha,
                          f"foreground {n_fg}; universe {universe}", notes, multiplicity="BH within the ORA family (q)",
                          effect_scale="gene ratio (k / n)", stage="R16e stage 2"),
            "traces": traces, "layout": _layout("Enrichment dot plot (family q)", "gene ratio (k / n)", "gene set"),
            "png": {"show_sets": shown}}


def _content_module_heatmap(rows: list[dict], sha: str, modules: Sequence[str] = (), **_: Any) -> dict:
    """V219: module eigengene group means for the declared modules only (E_FIGURE_SET_MISMATCH otherwise)."""
    declared = list(modules)
    if not declared:
        _refuse("E_FIGURE_SET_MISMATCH", "no declared module set was given; no module is plotted")
    for r in rows:
        if r["declared"] != "TRUE" or r["module"] not in declared:
            _refuse("E_FIGURE_SET_MISMATCH", f"module {r['module']} is outside the declared set {declared}")
    missing = [m for m in declared if not any(r["module"] == m for r in rows)]
    if missing:
        _refuse("E_FIGURE_SOURCE_MISSING", f"declared module(s) {missing} are not in the source")
    groups = list(dict.fromkeys(r["group"] for r in rows))
    cell = {(r["module"], r["group"]): float(r["eigengene"]) for r in rows}
    table = column_table("module eigengenes (group means)", {g: [cell[(m, g)] for m in declared] for g in groups})
    table = ColumnTable(table.title, table.columns, row_titles=tuple(declared))
    z = [[cell[(m, g)] for g in groups] for m in declared]
    n_per_group = {r["group"]: int(r["n"]) for r in rows}
    trace = {"type": "heatmap", "x": groups, "y": declared, "z": z, "zmid": 0, "colorscale": "RdBu", "reversescale": True,
             "colorbar": {"title": {"text": "eigengene"}}, "hovertemplate": _heatmap_hover("group", "module", "eigengene"), "hoverongaps": False}
    notes = ["Cell = mean over the samples of the group of the module eigengene (PC1 of the scaled module, sign-aligned).",
             f"Declared modules only: {', '.join(declared)}."]
    return {"tables": [table], "statistics": None, "info": _info("module heatmap", "none", "PC1 of each scaled module, group means", sha,
            "; ".join(f"{g} {n_per_group[g]}" for g in groups), notes, multiplicity="none (descriptive)", effect_scale="eigengene units",
            stage="R16e stage 2"), "traces": [trace], "layout": _layout("Module eigengenes", "group", "module"), "png": {"modules": declared}}


def _snapshot_pairs(path: str | Path) -> set[frozenset[str]]:
    """The undirected edges of the interaction snapshot (node_a, node_b columns)."""
    return {frozenset((r["node_a"], r["node_b"])) for r in read_source(path)}


def _network_layout(nodes: Sequence[Mapping[str, str]]) -> dict[str, tuple[float, float]]:
    """Deterministic layout coordinates computed in R (proteomicsCore:::catalogue_network_layout). Presentation only."""
    payload = [{"id": n["id"], "module": n["module"]} for n in nodes]
    with tempfile.TemporaryDirectory(prefix="catalogue-layout-") as tmp:
        nodes_file = Path(tmp) / "nodes.json"
        nodes_file.write_text(json.dumps(payload), encoding="utf-8", newline="\n")
        code = (
            "args <- commandArgs(trailingOnly = TRUE)\n"
            "suppressPackageStartupMessages(library(proteomicsCore))\n"
            "nodes <- jsonlite::fromJSON(args[1], simplifyVector = TRUE)\n"
            "xy <- proteomicsCore:::catalogue_network_layout(nodes)\n"
            "cat(jsonlite::toJSON(xy, digits = NA))\n"
        )
        result = run_r_code(code, [nodes_file], cwd=Path(tmp))
    if result.returncode != 0:
        _refuse("E_FIGURE_OUTPUT_MISSING", f"network layout failed: {result.stderr[-200:]}")
    coords = json.loads(result.stdout.strip().splitlines()[-1])
    return {c["id"]: (float(c["x"]), float(c["y"])) for c in coords}


def _content_network(rows: list[dict], sha: str, snapshot: set[frozenset[str]] | None = None, **_: Any) -> dict:
    """V220: nodes and edges of the declared subnetwork, exactly as in the interaction snapshot. An edge that is not in the
    snapshot is refused (E_FIGURE_SOURCE_MISSING). Layout coordinates are not data and are not compared."""
    if snapshot is None:
        _refuse("E_FIGURE_SOURCE_MISSING", "no interaction snapshot was given; the edges cannot be checked against it")
    nodes = [r for r in rows if r["kind"] == "node"]
    edges = [r for r in rows if r["kind"] == "edge"]
    node_ids = [n["id"] for n in nodes]
    for e in edges:
        if e["from"] not in node_ids or e["to"] not in node_ids:
            _refuse("E_FIGURE_SOURCE_MISSING", f"edge {e['from']}-{e['to']} has an endpoint that is not a node of the figure")
        if frozenset((e["from"], e["to"])) not in snapshot:
            _refuse("E_FIGURE_SOURCE_MISSING", f"edge {e['from']}-{e['to']} is not in the interaction snapshot")
    coords = _network_layout(nodes)
    degree = {n: 0 for n in node_ids}
    for e in edges:
        degree[e["from"]] += 1
        degree[e["to"]] += 1
    edge_x, edge_y, edge_custom = [], [], []
    for e in edges:
        (x1, y1), (x2, y2) = coords[e["from"]], coords[e["to"]]
        label = [f"{e['from']}-{e['to']}", "interaction (snapshot edge)", "present in the interaction snapshot"]
        edge_x += [x1, x2, None]
        edge_y += [y1, y2, None]
        edge_custom += [label, label, label]
    traces = [{"type": "scatter", "mode": "lines", "name": "edges", "x": edge_x, "y": edge_y, "line": {"color": "#9E9E9E", "width": 1.2},
               "customdata": edge_custom, "hoverinfo": "text",
               "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}<br>%{customdata[2]}<extra></extra>"}]
    for module in dict.fromkeys(n["module"] for n in nodes):
        members = [n for n in nodes if n["module"] == module]
        traces.append({"type": "scatter", "mode": "markers+text", "name": f"module {module}", "x": [coords[m["id"]][0] for m in members],
                       "y": [coords[m["id"]][1] for m in members], "text": [m["id"] for m in members], "textposition": "top center",
                       "textfont": {"size": 9}, "marker": {"size": 9},
                       "customdata": [[m["id"], module, degree[m["id"]]] for m in members],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>module: %{customdata[1]}<br>degree: %{customdata[2]}<extra></extra>"})
    position = {n: i + 1 for i, n in enumerate(node_ids)}
    tables = [ColumnTable("node degree in the declared subnetwork", (("degree", tuple(float(degree[n]) for n in node_ids)),), row_titles=tuple(node_ids)),
              xy_table("edge list (node indices)", "from (node index)", "to (node index)", [position[e["from"]] for e in edges],
                       [position[e["to"]] for e in edges])]
    layout = {"title": {"text": "Declared subnetwork (interaction snapshot)"}, "xaxis": {"visible": False},
              "yaxis": {"visible": False, "scaleanchor": "x"}, "plot_bgcolor": "#FFFFFF", "paper_bgcolor": "#FFFFFF", "showlegend": False}
    notes = [f"{len(node_ids)} nodes and {len(edges)} edges, exactly as in the interaction snapshot (induced on the declared subnetwork).",
             "Layout coordinates are presentation only (deterministic, computed in R) and are not part of the data."]
    return {"tables": tables, "statistics": None, "info": _info("network graph", "none", "none (descriptive; interaction snapshot)", sha,
            f"nodes {len(node_ids)}; edges {len(edges)}", notes, multiplicity="none (descriptive)", effect_scale="not applicable",
            stage="R16e stage 2"), "traces": traces, "layout": layout, "png": {}}


def _content_pcoa(rows: list[dict], sha: str, **_: Any) -> dict:
    """V221: principal coordinates with the group centroids (means of each group's members; never the full sample)."""
    groups = list(dict.fromkeys(r["group"] for r in rows))
    pct1, pct2 = float(rows[0]["percent_axis1"]), float(rows[0]["percent_axis2"])
    x_title, y_title = f"PCoA 1 ({pct1:.1f}%)", f"PCoA 2 ({pct2:.1f}%)"
    axis1 = [float(r["axis1"]) for r in rows]
    axis2 = [float(r["axis2"]) for r in rows]
    centroid = {g: (_mean([float(r["axis1"]) for r in rows if r["group"] == g]), _mean([float(r["axis2"]) for r in rows if r["group"] == g]))
                for g in groups}
    tables = [xy_table("PCoA scores", x_title, y_title, axis1, axis2),
              xy_table("PCoA group centroids (means of members)", x_title, y_title, [centroid[g][0] for g in groups], [centroid[g][1] for g in groups])]
    traces = []
    for g in groups:
        members = [r for r in rows if r["group"] == g]
        cx, cy = centroid[g]
        seg_x, seg_y, seg_custom = [], [], []
        for r in members:
            seg_x += [float(r["axis1"]), cx, None]
            seg_y += [float(r["axis2"]), cy, None]
            seg_custom += [[r["observation_id"], g, "segment to the group centroid"]] * 3
        traces.append({"type": "scatter", "mode": "lines", "name": f"{g} to centroid", "x": seg_x, "y": seg_y,
                       "line": {"width": 0.8, "color": "#BDBDBD"}, "showlegend": False, "customdata": seg_custom,
                       "hovertemplate": "<b>%{customdata[0]}</b><br>group: %{customdata[1]}<br>%{customdata[2]}<extra></extra>"})
        traces.append({"type": "scatter", "mode": "markers", "name": g, "x": [float(r["axis1"]) for r in members],
                       "y": [float(r["axis2"]) for r in members], "marker": {"size": 7, "opacity": 0.85},
                       "customdata": [[r["observation_id"], g, float(r["axis1"]), float(r["axis2"])] for r in members],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>group: %{customdata[1]}<br>PCoA 1 %{customdata[2]:.3f}"
                                        "<br>PCoA 2 %{customdata[3]:.3f}<extra></extra>"})
        traces.append({"type": "scatter", "mode": "markers", "name": f"{g} centroid", "x": [cx], "y": [cy],
                       "marker": {"size": 13, "symbol": "diamond", "line": {"color": "#000000", "width": 1}},
                       "customdata": [[f"centroid {g}", g, f"mean of {len(members)} members"]],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>group: %{customdata[1]}<br>%{customdata[2]}<extra></extra>"})
    sizes = "; ".join(f"{g} {sum(1 for r in rows if r['group'] == g)}" for g in groups)
    notes = ["Principal coordinates of the Euclidean distance (cmdscale); axis labels give the percent of positive eigenvalues.",
             "Each centroid is the mean of its own group's members; the full-sample mean is not drawn."]
    return {"tables": tables, "statistics": None, "info": _info("PCoA with group centroids", "none", "principal coordinates (cmdscale) of the Euclidean distance",
            sha, sizes, notes, multiplicity="none (descriptive)", effect_scale="not applicable", stage="R16e stage 2"),
            "traces": traces, "layout": _layout("PCoA with group centroids", x_title, y_title), "png": {}}


def _content_permdisp(rows: list[dict], sha: str, **_: Any) -> dict:
    """V222: distances to group centroids with the PERMDISP P, the number of permutations and the seed. A P without B or
    seed is refused (E_FIGURE_SOURCE_MISSING)."""
    if any(not r["permutations"] or not r["seed"] for r in rows):
        _refuse("E_FIGURE_SOURCE_MISSING", "the PERMDISP P is shown without its number of permutations B and its seed")
    head = rows[0]
    boots, seed = int(head["permutations"]), int(head["seed"])
    p_value = float(head["permdisp_p"])
    groups = list(dict.fromkeys(r["group"] for r in rows))
    distances = {g: [float(r["distance_to_centroid"]) for r in rows if r["group"] == g] for g in groups}
    tables = [column_table("PERMDISP distances to group centroid", distances)]
    stats = [{"row": f"{g} - {groups[0]}", "effect": _mean(distances[g]) - _mean(distances[groups[0]]), "P": p_value, "q": None,
              "n": {gg: len(distances[gg]) for gg in groups}} for g in groups[1:]]
    if stats:
        tables.append(statistics_table("PERMDISP - pipeline statistics", stats, groups))
    traces = []
    for g in groups:
        members = [r for r in rows if r["group"] == g]
        traces.append({"type": "box", "name": g, "x": [g] * len(members), "y": [float(r["distance_to_centroid"]) for r in members],
                       "boxpoints": "all", "jitter": 0.3, "pointpos": 0, "marker": {"color": "#0072B2", "size": 6}, "line": {"color": "#000000"},
                       "customdata": [[r["observation_id"], g, f"PERMDISP P = {p_value:.3f}; B = {boots}; seed {seed}"] for r in members],
                       "hovertemplate": "<b>%{customdata[0]}</b><br>group: %{customdata[1]}<br>distance %{y:.3f}"
                                        "<br>%{customdata[2]}<extra></extra>"})
    sizes = "; ".join(f"{g} {len(distances[g])}" for g in groups)
    notes = [f"PERMDISP P = {p_value:.4g} from permutest with B = {boots} permutations and seed {seed}; Euclidean distance to the group centroid.",
             "Statistics table: the dispersion difference of each group from the first group; P is the PERMDISP family P."]
    return {"tables": tables, "statistics": stats or None, "info": _info("PERMDISP distances to centroid", "table" if stats else "none",
            f"PERMDISP (vegan betadisper and permutest), B = {boots}, seed {seed}", sha, sizes, notes,
            multiplicity="none (one PERMDISP family per figure)", effect_scale="distance to centroid difference", stage="R16e stage 2"),
            "traces": traces, "layout": _layout(f"PERMDISP: distances to group centroids (P = {p_value:.3f}, B = {boots}, seed {seed})",
                                                "group", "distance to group centroid"), "png": {}}


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
    "bio_roc_ci_band": _content_roc,
    "bio_panel_forest": _content_forest,
    "bio_cv_auc_vs_permutation": _content_cv_permutation,
    "path_enrichment_dotplot": _content_enrichment,
    "net_module_heatmap": _content_module_heatmap,
    "net_network_graph": _content_network,
    "perm_pcoa_centroids": _content_pcoa,
    "perm_permdisp": _content_permdisp,
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
    if figure_id == "net_network_graph":
        if not options.get("snapshot"):
            _refuse("E_FIGURE_SOURCE_MISSING", "net_network_graph needs its interaction snapshot (option snapshot)")
        builder_options["snapshot"] = _snapshot_pairs(options["snapshot"])
    if figure_id == "net_module_heatmap":
        builder_options["modules"] = list(options.get("modules") or ())
    if figure_id == "path_enrichment_dotplot" and options.get("show_sets") is not None:
        builder_options["show_sets"] = list(options["show_sets"])
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
    elif figure_id == "path_enrichment_dotplot":
        png_options = {"show_sets": content["png"]["show_sets"]}
    elif figure_id == "net_module_heatmap":
        png_options = {"modules": content["png"]["modules"]}
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
