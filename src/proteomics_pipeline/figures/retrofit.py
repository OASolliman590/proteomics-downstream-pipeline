"""Retrofit of the existing figures to the three outputs (D-72, ADR 0012; FR-195 to FR-205 apply).

The R10b and R14 producers already write a registered figure-source table for each figure (``figure_sources/<stem>.tsv``, and
``heatmap_source.tsv`` for the association heatmap). A rule, chosen from the stem of that table, gives the figure kind and the
column names. The figure is then drawn from the registered copy only, and the three outputs are written next to the catalogue
figures: the .pzfx file (data table and Info sheet), the Prism-style PNG at 300 dpi (R16c theme) and the standalone HTML.

Transformations are limited to those named by a rule: the absolute value of the treatment effect of the response axes, group
means of the eigengenes, equal-width bins of a histogram in the HTML, and the subnetwork edges read from the companion table of
the same stage. No figure recomputes an analysis. The R13 PERMANOVA producer writes no figure-source table, so it has no rule here.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..errors import ConfigurationError
from ..runtime import run_r_code
from .catalogue import _heatmap_hover, _info, _layout, _refuse, read_source
from .plotly_html import write_html
from .pzfx_writer import ColumnTable, column_table, write_pzfx, xy_table
from .registry import FigureRegistry

STAGE_LABEL = "R16e retrofit (D-72)"
COLUMN_OPTIONS = ("x", "y", "label", "value", "fill", "flag_col", "band_x", "band_lo", "band_hi", "edge_a", "edge_b")
PNG_COLOUR = "#D55E00"


@dataclass(frozen=True)
class Rule:
    pattern: str                      # full match against the stem of the registered table
    kind: str                         # figure kind (the R renderer and the Python content builder)
    options: Mapping[str, Any]        # column names and labels; the R renderer reads them
    skip: bool = False                # the table is a companion of another figure, or covered by a run-level figure
    companion: str | None = None      # stem of the companion table; "parent:<stem>" looks in the parent directory
    note: str = ""                    # why a skipped table is not drawn on its own


RULES: tuple[Rule, ...] = (
    Rule(r"volcano_(?:zero_null|treat)_.+", "diff_volcano_labelled", {}, skip=True,
         note="drawn by the run's differential figure diff_volcano_labelled__<contrast> from the same table (catalogue renderer)"),
    Rule(r"response_axes", "scatter", {"x": "d", "y": "t", "abs_y": True, "flag_col": "descriptive_class",
         "flag_values": ["near_restoration", "partial_return"], "label": "feature_id", "xlab": "Disease effect d (log2)",
         "ylab": "|treatment effect t| (log2)", "title": "Descriptive response axes", "caption": "red: near restoration or partial return; no P values"}),
    Rule(r"roc_cv", "curve_band", {"x": "fpr", "y": "tpr", "band_x": "fpr", "band_lo": "tpr_lower", "band_hi": "tpr_upper",
         "xlab": "1 - specificity", "ylab": "sensitivity", "title": "ROC with bootstrap band"}, companion="roc_cv_band", skip=True,
         note="drawn by the run's post-DE figure bio_roc_ci_band from the same table (D3)"),
    Rule(r"roc_cv_band", "curve_band", {}, skip=True, note="the band of roc_cv, drawn with it"),
    Rule(r"cv_auc_repeats", "strip", {"x": "auc", "xlab": "nested cross-validation repeat AUC", "title": "Nested cross-validation: repeat AUCs"},
         skip=True, note="drawn by the run's post-DE figure bio_cv_auc_vs_permutation from the same tables (D3)"),
    Rule(r"permutation_null", "histogram", {"x": "pooled_oof_auc", "xlab": "pooled out-of-fold AUC under permutation",
         "title": "Permutation null of pooled AUC"}, skip=True,
         note="drawn by the run's post-DE figure bio_cv_auc_vs_permutation from the same tables (D3)"),
    Rule(r"calibration", "scatter", {"x": "mean_predicted", "y": "observed_fraction", "diagonal": True, "xlab": "mean predicted probability",
         "ylab": "observed fraction", "title": "Calibration (cross-validated)"}),
    Rule(r"selection_stability", "bar", {"label": "feature_id", "value": "selection_frequency", "xlab": "selection frequency across outer folds",
         "title": "Feature selection stability"}),
    Rule(r"upset_all|venn_.+", "regions", {"label": "defining_sets", "value": "n_protein_groups", "xlab": "protein groups in the region",
         "title": "Exclusive set regions"}),
    Rule(r"concordance_.+", "scatter", {"x": "left_effect", "y": "right_effect", "flag_col": "significant_either", "flag_values": ["true"],
         "xlab": "left effect", "ylab": "right effect", "title": "Effect concordance (descriptive)", "caption": "red: significant in either set"}),
    Rule(r"comparison_.+", "scatter", {"x": "primary_effect", "y": "sensitivity_effect", "diagonal": True, "xlab": "primary effect",
         "ylab": "sensitivity effect", "title": "Sensitivity versus primary (descriptive)"}),
    Rule(r"matched_n_discoveries", "histogram", {"x": "n_discoveries", "xlab": "discoveries per resampled draw",
         "title": "Matched-n resampling: discoveries per draw"}),
    Rule(r"influence_max_effect_change", "scatter", {"x": "rank", "y": "max_abs_effect_change", "xlab": "influence rank",
         "ylab": "max absolute effect change", "title": "Influence of single units", "caption": "descriptive"}),
    Rule(r"eigengenes", "eigen_means", {"xlab": "group", "ylab": "module", "fill_title": "eigengene (group mean)",
         "title": "Module eigengenes by group (group means)"}, skip=True,
         note="drawn by the run's post-DE figure net_module_heatmap (declared modules, group means; D3)"),
    Rule(r"ppi_subnetwork", "network", {"x": "x", "y": "y", "label": "gene_id", "edge_a": "gene_a", "edge_b": "gene_b",
         "title": "Declared subnetwork (descriptive)"}, companion="parent:ppi_subnetwork", skip=True,
         note="drawn by the run's post-DE figure net_network_graph from the subnetwork edges, checked against the snapshot (D3)"),
    Rule(r"association_heatmap", "heat_table", {"x": "analysis_id", "y": "feature_id", "fill": "value", "fill_title": "statistic",
         "xlab": "analysis", "ylab": "feature", "title": "Association statistics"}),
)


def rule_for(stem: str) -> Rule | None:
    for rule in RULES:
        if re.fullmatch(rule.pattern, stem):
            return rule
    return None


COMPANION_OPTIONS = ("band_x", "band_lo", "band_hi", "edge_a", "edge_b")


def required_columns(rule: Rule) -> list[str]:
    """The columns of the registered table that the rule needs (the companion's columns are checked against the companion)."""
    columns = [str(rule.options[key]) for key in COLUMN_OPTIONS if key in rule.options and key not in COMPANION_OPTIONS]
    if rule.kind == "eigen_means":
        columns.append("group")
    return columns


def companion_columns(rule: Rule) -> list[str]:
    return [str(rule.options[key]) for key in COMPANION_OPTIONS if key in rule.options]


def _f(text: str | None) -> float | None:
    if text is None or text in ("", "NA", "NaN"):
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _claims(rows: Sequence[Mapping[str, str]]) -> list[str]:
    return sorted({r["claim_label"] for r in rows if r.get("claim_label")}) if rows and "claim_label" in rows[0] else []


def _equal_bins(values: Sequence[float], count: int = 30) -> tuple[list[float], list[int]]:
    """Centres and counts of equal-width bins from the minimum to the maximum (the HTML view of a histogram)."""
    lo, hi = min(values), max(values)
    width = (hi - lo) / count if hi > lo else 1.0
    counts = [0] * count
    for v in values:
        index = min(int((v - lo) / width), count - 1)
        counts[index] += 1
    return [lo + (i + 0.5) * width for i in range(count)], counts


def _info_for(title: str, sha: str, rows: Sequence[Mapping[str, str]], test: str, notes: list[str], claims: list[str]) -> dict:
    notes = list(notes) + ([f"claim label(s) of the source: {', '.join(claims)}"] if claims else [])
    return _info(title, "none", test, sha, f"rows {len(rows)}", notes, multiplicity="none (descriptive; retrofit of a registered table)",
                 effect_scale="not applicable", stage=STAGE_LABEL)


def _content(rule: Rule, rows: list[dict], sha: str, companion: list[dict] | None, options: Mapping[str, Any], title: str) -> dict:
    o = rule.options
    claims = sorted(set(_claims(rows)) | set(_claims(companion or [])))
    notes = ["Drawn from the registered figure-source table of the producer; the values are those of the table."]
    if rule.kind == "scatter":
        xs = [(_f(r[o["x"]]), _f(r[o["y"]]), r) for r in rows]
        if o.get("abs_y"):
            xs = [(x, (abs(y) if y is not None else None), r) for x, y, r in xs]
        pts = [(x, y, r) for x, y, r in xs if x is not None and y is not None]
        flag_values = set(o.get("flag_values") or [])
        is_flag = [bool(flag_values) and str(r.get(o.get("flag_col", ""), "")) in flag_values for _, _, r in pts]
        flagged = [p for p, f in zip(pts, is_flag) if f]
        plain = [p for p, f in zip(pts, is_flag) if not f]
        label = (lambda r: str(r.get(o.get("label", ""), ""))) if o.get("label") else (lambda r: "")
        tables = [xy_table(title, o["xlab"], o["ylab"], [x for x, _, _ in pts], [y for _, y, _ in pts])]
        traces = []
        for name, members, colour in (("flagged", flagged, PNG_COLOUR), ("other", plain, "#9E9E9E")):
            if not members and name == "flagged":
                continue
            traces.append({"type": "scatter", "mode": "markers", "name": name, "x": [x for x, _, _ in members], "y": [y for _, y, _ in members],
                           "marker": {"color": colour, "size": 6, "opacity": 0.85},
                           "customdata": [[label(r) or "unit", name, f"{o['x']} {x:.4g}; {o['ylab']} {y:.4g}"] for x, y, r in members],
                           "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}<br>%{customdata[2]}<extra></extra>"})
        if o.get("diagonal"):
            traces.append({"type": "scatter", "mode": "lines", "name": "diagonal", "x": [0.0, 1.0], "y": [0.0, 1.0], "line": {"dash": "dash", "color": "#9E9E9E"},
                           "customdata": [["diagonal", "reference line", "y = x"]] * 2,
                           "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}<br>%{customdata[2]}<extra></extra>"})
        layout = _layout(title, o["xlab"], o["ylab"])
        info = _info_for(title, sha, rows, "retrofit of the producer's table (no inference recomputed)", notes, claims)
        return {"tables": tables, "traces": traces, "layout": layout, "info": info, "statistics": None}
    if rule.kind == "curve_band":
        band = companion or []
        cx = [_f(r[o["x"]]) for r in rows]
        cy = [_f(r[o["y"]]) for r in rows]
        curve = [(x, y) for x, y in zip(cx, cy) if x is not None and y is not None]
        bx = [_f(r[o["band_x"]]) for r in band]
        blo = [_f(r[o["band_lo"]]) for r in band]
        bhi = [_f(r[o["band_hi"]]) for r in band]
        bpts = [(x, lo, hi) for x, lo, hi in zip(bx, blo, bhi) if x is not None and lo is not None and hi is not None]
        tables = [xy_table("ROC curve", o["xlab"], o["ylab"], [x for x, _ in curve], [y for _, y in curve]),
                  xy_table("ROC band, lower", o["xlab"], "lower bound", [x for x, _, _ in bpts], [lo for _, lo, _ in bpts]),
                  xy_table("ROC band, upper", o["xlab"], "upper bound", [x for x, _, _ in bpts], [hi for _, _, hi in bpts])]
        traces = [
            {"type": "scatter", "mode": "lines", "name": "band, lower", "x": [x for x, _, _ in bpts], "y": [lo for _, lo, _ in bpts], "line": {"width": 0},
             "showlegend": False, "customdata": [["95% band (bootstrap)", "lower bound", f"{lo:.3f}"] for _, lo, _ in bpts],
             "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}: %{customdata[2]}<extra></extra>"},
            {"type": "scatter", "mode": "lines", "name": "band", "x": [x for x, _, _ in bpts], "y": [hi for _, _, hi in bpts], "line": {"width": 0},
             "fill": "tonexty", "fillcolor": "rgba(213,94,0,0.18)", "customdata": [["95% band (bootstrap)", "upper bound", f"{hi:.3f}"] for _, _, hi in bpts],
             "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}: %{customdata[2]}<extra></extra>"},
            {"type": "scatter", "mode": "lines", "name": "ROC", "x": [x for x, _ in curve], "y": [y for _, y in curve],
             "line": {"shape": "hv", "color": PNG_COLOUR, "width": 2.5}, "customdata": [["ROC", "curve point", f"{y:.3f}"] for _, y in curve],
             "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}: sensitivity %{customdata[2]}<extra></extra>"},
        ]
        info = _info_for(title, sha, rows, "retrofit of the producer's ROC table and band (no inference recomputed)", notes, claims)
        return {"tables": tables, "traces": traces, "layout": _layout(title, o["xlab"], o["ylab"]), "info": info, "statistics": None}
    if rule.kind == "histogram":
        values = [v for v in (_f(r[o["x"]]) for r in rows) if v is not None]
        if not values:
            _refuse("E_FIGURE_SOURCE_MISSING", f"no finite values in {o['x']}")
        centres, counts = _equal_bins(values)
        tables = [column_table(title, {o["x"]: values})]
        traces = [{"type": "bar", "name": "count", "x": centres, "y": counts, "marker": {"color": "#BDBDBD"},
                   "customdata": [[f"bin around {c:.4g}", "count", n] for c, n in zip(centres, counts)],
                   "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}: %{customdata[2]}<extra></extra>"}]
        if o.get("ref") is not None:
            traces.append({"type": "scatter", "mode": "lines", "name": "reference", "x": [float(o["ref"])] * 2, "y": [0, max(counts)],
                           "line": {"color": PNG_COLOUR}, "customdata": [["reference", "value", f"{float(o['ref']):.4g}"]] * 2,
                           "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}: %{customdata[2]}<extra></extra>"})
        info = _info_for(title, sha, rows, "retrofit of the producer's table; the HTML bins are equal-width", notes, claims)
        return {"tables": tables, "traces": traces, "layout": _layout(title, o["xlab"], "count"), "info": info, "statistics": None}
    if rule.kind == "strip":
        values = [v for v in (_f(r[o["x"]]) for r in rows) if v is not None]
        tables = [column_table(title, {o["x"]: values})]
        traces = [{"type": "scatter", "mode": "markers", "name": "repeats", "x": ["repeats"] * len(values), "y": values,
                   "marker": {"color": "#0072B2", "size": 8}, "customdata": [[f"repeat {i + 1}", "AUC", f"{v:.4f}"] for i, v in enumerate(values)],
                   "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}: %{customdata[2]}<extra></extra>"}]
        info = _info_for(title, sha, rows, "retrofit of the producer's table", notes, claims)
        return {"tables": tables, "traces": traces, "layout": _layout(title, "", o["xlab"]), "info": info, "statistics": None}
    if rule.kind in ("bar", "regions"):
        labels = [str(r[o["label"]]) for r in rows]
        values = [_f(r[o["value"]]) for r in rows]
        keep = [(l, v) for l, v in zip(labels, values) if v is not None]
        tables = [ColumnTable(title, (("value", tuple(v for _, v in keep)),), row_titles=tuple(l for l, _ in keep))]
        traces = [{"type": "bar", "orientation": "h", "name": title, "x": [v for _, v in keep], "y": [l for l, _ in keep],
                   "marker": {"color": "#4c78a8" if rule.kind == "bar" else "#000000"},
                   "customdata": [[l, o["value"], f"{v:.4g}"] for l, v in keep],
                   "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}: %{customdata[2]}<extra></extra>"}]
        info = _info_for(title, sha, rows, "retrofit of the producer's table (no inference recomputed)", notes, claims)
        return {"tables": tables, "traces": traces, "layout": _layout(title, o["xlab"], ""), "info": info, "statistics": None}
    if rule.kind in ("heat_table", "eigen_means"):
        if rule.kind == "heat_table":
            cells = [(str(r[o["x"]]), str(r[o["y"]]), _f(r[o["fill"]])) for r in rows]
        else:
            modules = [c for c in rows[0] if c not in ("observation_id", "group", "claim_label")]
            groups = list(dict.fromkeys(r["group"] for r in rows))
            cells = [(g, m, sum(_f(r[m]) for r in rows if r["group"] == g and _f(r[m]) is not None) /
                      max(1, sum(1 for r in rows if r["group"] == g and _f(r[m]) is not None)))
                     for m in modules for g in groups]
        xs = list(dict.fromkeys(c[0] for c in cells))
        ys = list(dict.fromkeys(c[1] for c in cells))
        table_cells = {(c[0], c[1]): c[2] for c in cells}
        z = [[table_cells.get((x, y)) for x in xs] for y in ys]
        tables = [ColumnTable(title, tuple((x, tuple(table_cells.get((x, y)) for y in ys)) for x in xs), row_titles=tuple(ys))]
        trace = {"type": "heatmap", "x": xs, "y": ys, "z": z, "zmid": 0, "colorscale": "RdBu", "reversescale": True,
                 "colorbar": {"title": {"text": o.get("fill_title", "value")}}, "hovertemplate": _heatmap_hover(o.get("xlab", "x"), o.get("ylab", "y"), o.get("fill_title", "value")),
                 "hoverongaps": False}
        info = _info_for(title, sha, rows, "retrofit of the producer's table (no inference recomputed)", notes, claims)
        return {"tables": tables, "traces": [trace], "layout": _layout(title, o.get("xlab", ""), o.get("ylab", "")), "info": info, "statistics": None}
    if rule.kind == "network":
        if companion is None:
            _refuse("E_FIGURE_SOURCE_MISSING", "the subnetwork edges (companion table) are missing")
        ids = [str(r[o["label"]]) for r in rows]
        coords = {str(r[o["label"]]): (_f(r[o["x"]]), _f(r[o["y"]])) for r in rows}
        edges = [(str(e[o["edge_a"]]), str(e[o["edge_b"]])) for e in companion]
        for a, b in edges:
            if a not in coords or b not in coords:
                _refuse("E_FIGURE_SOURCE_MISSING", f"edge {a}-{b} has an endpoint that is not a node of the figure")
        degree = {n: 0 for n in ids}
        for a, b in edges:
            degree[a] += 1
            degree[b] += 1
        ex, ey, ec = [], [], []
        for a, b in edges:
            label = [f"{a}-{b}", "subnetwork edge", "present in the producer's table"]
            (x1, y1), (x2, y2) = coords[a], coords[b]
            ex += [x1, x2, None]
            ey += [y1, y2, None]
            ec += [label, label, label]
        traces = [{"type": "scatter", "mode": "lines", "name": "edges", "x": ex, "y": ey, "line": {"color": "#9E9E9E", "width": 1.2},
                   "customdata": ec, "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}<br>%{customdata[2]}<extra></extra>"},
                  {"type": "scatter", "mode": "markers+text", "name": "nodes", "x": [coords[n][0] for n in ids], "y": [coords[n][1] for n in ids],
                   "text": ids, "textposition": "top center", "textfont": {"size": 9}, "marker": {"size": 9, "color": "#0072B2"},
                   "customdata": [[n, "node", degree[n]] for n in ids],
                   "hovertemplate": "<b>%{customdata[0]}</b><br>%{customdata[1]}; degree %{customdata[2]}<extra></extra>"}]
        node_index = {n: i + 1 for i, n in enumerate(ids)}
        tables = [ColumnTable(title + " (degree)", (("degree", tuple(float(degree[n]) for n in ids)),), row_titles=tuple(ids)),
                  xy_table("edge list (node indices)", "from (node index)", "to (node index)", [node_index[a] for a, _ in edges],
                           [node_index[b] for _, b in edges])]
        layout = {"title": {"text": title}, "xaxis": {"visible": False}, "yaxis": {"visible": False, "scaleanchor": "x"},
                  "plot_bgcolor": "#FFFFFF", "paper_bgcolor": "#FFFFFF", "showlegend": False}
        info = _info_for(title, sha, rows, "none (descriptive; subnetwork as in the producer's table)", notes, claims)
        return {"tables": tables, "traces": traces, "layout": layout, "info": info, "statistics": None}
    raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"no retrofit content for kind {rule.kind}")


def _render_png(kind: str, source_copy: Path, out: Path, options: Mapping[str, Any], width_mm: int) -> None:
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "suppressPackageStartupMessages(library(proteomicsCore))\n"
        "opts <- jsonlite::fromJSON(args[4], simplifyVector = TRUE)\n"
        "r <- proteomicsCore:::render_retrofit_png(args[1], args[2], args[3], opts, width_mm = as.numeric(args[5]))\n"
        "cat(jsonlite::toJSON(list(width_px = r$width_px, height_px = r$height_px), auto_unbox = TRUE))\n"
    )
    with tempfile.TemporaryDirectory(prefix="retrofit-png-") as tmp:
        options_file = Path(tmp) / "options.json"
        options_file.write_text(json.dumps(dict(options)), encoding="utf-8", newline="\n")
        result = run_r_code(code, [kind, source_copy, out, options_file, width_mm], cwd=Path(tmp))
    if result.returncode != 0:
        match = re.search(r"(E_[A-Z_]+):\s*([^\n]*)", result.stderr)
        if match:
            raise ConfigurationError(match.group(1), match.group(2).strip())
        raise ConfigurationError("E_FIGURE_OUTPUT_MISSING", f"the PNG of the retrofitted figure was not written: {result.stderr.strip()[-200:]}")


def build_retrofit(registry: FigureRegistry, figure_id: str, source: Path, rule: Rule, formats: Sequence[str], width_mm: int,
                   companion: Path | None = None) -> dict:
    """Produce the requested outputs of one retrofitted figure from its registered source copy (and its companion copy)."""
    registry.register(figure_id, source, retrofit_kind=rule.kind)
    catalogue_dir = registry.directory
    copy = catalogue_dir / "sources" / f"{figure_id}.tsv"
    rows = read_source(copy)
    missing = [c for c in required_columns(rule) if rows and c not in rows[0]]
    if missing:
        raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"{source.name} lacks the column(s) {missing} that the rule of {rule.kind} needs")
    companion_rows = None
    companion_copy = None
    if rule.companion is not None:
        if companion is None:
            raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"{source.name} needs its companion table {rule.companion}")
        companion_copy = catalogue_dir / "sources" / f"{figure_id}.{companion.stem}.tsv"
        companion_copy.write_bytes(companion.read_bytes())
        companion_rows = read_source(companion_copy)
        absent = [c for c in companion_columns(rule) if companion_rows and c not in companion_rows[0]]
        if absent:
            raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"{companion.name} lacks the column(s) {absent} that the rule of {rule.kind} needs")
        registry.figures[figure_id]["companions"][companion.name] = {"path": companion_copy.relative_to(catalogue_dir).as_posix(),
                                                                     "sha256": hashlib.sha256(companion_copy.read_bytes()).hexdigest()}
    sha = registry.figures[figure_id]["source_sha256"]
    title = figure_id.replace("_", " ")
    content = _content(rule, rows, sha, companion_rows, rule.options, title)
    png_options = dict(rule.options)
    claims = sorted(set(_claims(rows)) | set(_claims(companion_rows or [])))
    notes_caption = str(rule.options.get("caption", ""))
    claim_caption = f"claim label: {', '.join(claims)}" if claims else ""
    png_options["caption"] = "; ".join(part for part in (notes_caption, claim_caption) if part) or "retrofit of the registered producer table"
    if companion_copy is not None:
        png_options["band_path" if rule.kind == "curve_band" else "edges_path"] = str(companion_copy)
    if rule.kind == "histogram" and rule.options.get("ref") is not None:
        png_options["ref"] = rule.options["ref"]

    def render_png(source_path: Path, out_path: Path) -> None:
        _render_png(rule.kind, source_path, out_path, png_options, width_mm)

    def render_pzfx(source_path: Path, out_path: Path) -> None:
        write_pzfx(out_path, content["tables"], content["info"])

    def render_html(source_path: Path, out_path: Path) -> None:
        write_html(out_path, figure_id, content["traces"], content["layout"])

    renderers = {"png": render_png, "pzfx": render_pzfx, "html": render_html}
    outputs: dict[str, str] = {}
    refusals: list[dict[str, str]] = []
    for kind in ("png", "pzfx", "html"):
        if kind not in formats:
            continue
        try:
            registry.produce(figure_id, kind, renderers[kind])
            outputs[kind] = str(catalogue_dir / f"{figure_id}.{kind}")
        except ConfigurationError as error:
            registry.record_missing_output(figure_id, kind, error.message)
            refusals.append({"kind": kind, "code": error.code, "reason": error.message})
    return {"figure_id": figure_id, "kind": rule.kind, "outputs": outputs, "refusals": refusals}
