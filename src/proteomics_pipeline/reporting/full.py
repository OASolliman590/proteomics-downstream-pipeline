"""Full offline report (packet R10b: FR-095 to FR-099).

Reads the same verified ReportData as the R10a stub (the stub is reused, not
changed) and adds: every planned endpoint row of every differential table,
figures drawn only from saved figure-source tables (SVG in Python; PDF/PNG
through the R plotting helpers when requested), pathway/response/multivariate
views and generated methods/limitations.  Failed, partial or cancelled work
is shown as such; nothing is substituted with zeros.
"""
from __future__ import annotations

import csv
import html
import json
import math
import os
import shutil
import subprocess
from pathlib import Path
from string import Template

from ..provenance import sha256_file
from .assemble import assemble, validate_report_data
from .methods import render_full_methods
from .stub import Raw, e, link, render as render_stub, state_span, table

CAPABILITY = "report_full"
PALETTE = {"sig": "#b2182b", "ns": "#7f7f7f"}


def capabilities():
    return [{"id": CAPABILITY, "dependencies": [], "required_r_packages": []}]


def _table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(["NA" if v is None else v for v in row])


def _num(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def scatter_svg(points: list[dict], *, x_label: str, y_label: str, title: str, hline: float | None = None) -> str:
    """Plain SVG scatter from exact source coordinates (no rounding of plotted data)."""
    W, H, L, R, T, Bm = 520, 380, 64, 16, 40, 52
    xs = [p["x"] for p in points if p["x"] is not None and p["y"] is not None]
    ys = [p["y"] for p in points if p["x"] is not None and p["y"] is not None]
    if not xs:
        return f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="80" role="img" aria-label="{html.escape(title)}"><text x="10" y="40">{html.escape(title)}: no finite points (a real empty result)</text></svg>'
    x0, x1 = min(xs + [0.0]), max(xs + [0.0]); y0, y1 = 0.0, max(ys + [hline or 0.0]) or 1.0
    if x0 == x1:
        x0, x1 = x0 - 1, x1 + 1
    sx = lambda x: L + (x - x0) / (x1 - x0) * (W - L - R)
    sy = lambda y: H - Bm - (y - y0) / (y1 - y0) * (H - T - Bm)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="{html.escape(title)}">',
             f'<rect width="{W}" height="{H}" fill="white"/>', f'<text x="{W / 2}" y="22" text-anchor="middle" font-size="13" font-family="sans-serif">{html.escape(title)}</text>',
             f'<line x1="{L}" y1="{H - Bm}" x2="{W - R}" y2="{H - Bm}" stroke="black"/>', f'<line x1="{L}" y1="{T}" x2="{L}" y2="{H - Bm}" stroke="black"/>',
             f'<text x="{(L + W - R) / 2}" y="{H - 14}" text-anchor="middle" font-size="12" font-family="sans-serif">{html.escape(x_label)}</text>',
             f'<text transform="translate(16,{(T + H - Bm) / 2}) rotate(-90)" text-anchor="middle" font-size="12" font-family="sans-serif">{html.escape(y_label)}</text>']
    for tick in (x0, (x0 + x1) / 2, x1):
        parts.append(f'<text x="{sx(tick):.1f}" y="{H - Bm + 16}" text-anchor="middle" font-size="10" font-family="sans-serif">{tick:.3g}</text>')
    for tick in (y0, (y0 + y1) / 2, y1):
        parts.append(f'<text x="{L - 6}" y="{sy(tick) + 3:.1f}" text-anchor="end" font-size="10" font-family="sans-serif">{tick:.3g}</text>')
    if hline is not None:
        parts.append(f'<line x1="{L}" y1="{sy(hline):.1f}" x2="{W - R}" y2="{sy(hline):.1f}" stroke="#555" stroke-dasharray="4 3"/>')
    for p in points:
        if p["x"] is None or p["y"] is None:
            continue
        colour = PALETTE["sig"] if p.get("highlight") else PALETTE["ns"]
        parts.append(f'<circle cx="{sx(p["x"]):.2f}" cy="{sy(p["y"]):.2f}" r="3" fill="{colour}" fill-opacity="0.75"><title>{html.escape(str(p["label"]))}</title></circle>')
    parts.append("</svg>")
    return "".join(parts)


def _r_figure(source: Path, kind: str, target: Path, fmt: str) -> bool:
    rscript = os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript")
    if not rscript:
        return False
    code = "a <- commandArgs(TRUE); proteomicsCore:::plot_figure_source(a[1], a[2], a[3], a[4])"
    from ..runtime import run_r_code
    result = run_r_code(code, [source, kind, target, fmt], rscript=rscript)
    return result.returncode == 0 and target.is_file()


def build_figures(root: Path, out: Path, data: dict, formats: list[str]) -> list[dict]:
    figures = []
    families = {f["family_id"]: float(f["q_cutoff"]) for f in (data["dea"].get("families") or [])} if data["dea"].get("state") == "COMPLETED" else {}
    if data["dea"].get("state") == "COMPLETED":
        for name in ("zero_null", "treat"):
            source = root / "dea" / f"{name}.tsv"
            if not source.is_file():
                continue
            rows = _table(source)
            for key in sorted({(r["model_id"], r["contrast_id"]) for r in rows}):
                sel = [r for r in rows if (r["model_id"], r["contrast_id"]) == key]
                tested = [r for r in sel if r["eligibility"] == "tested"]
                stem = f"volcano_{name}_{key[0]}_{key[1]}"
                fam = sel[0]["family_id"]; cutoff = families.get(fam)
                src_rows = []
                for r in tested:
                    p = _num(r["p_value"]); q = _num(r["q_value"])
                    src_rows.append([r["feature_id"], r["effect"], repr(-math.log10(p)) if p and p > 0 else "NA", r["p_value"], r["q_value"], fam,
                                     "true" if (q is not None and cutoff is not None and q <= cutoff) else "false", r["n_obs_by_required_group"]])
                src = out / "figure_sources" / f"{stem}.tsv"
                _write_tsv(src, ["feature_id", "effect_log2", "neg_log10_p", "p_value", "q_value", "family_id", "q_at_or_below_cutoff", "n_obs_by_required_group"], src_rows)
                points = [{"x": _num(r[1]), "y": _num(r[2]), "label": r[0], "highlight": r[6] == "true"} for r in src_rows]
                svg = scatter_svg(points, x_label="Effect (log2 difference)", y_label="-log10 P", title=f"{key[1]} ({key[0]}, {name}); {len(tested)} tested of {len(sel)} planned")
                (out / "figures").mkdir(parents=True, exist_ok=True)
                (out / "figures" / f"{stem}.svg").write_text(svg, encoding="utf-8")
                record = {"stem": stem, "source": f"figure_sources/{stem}.tsv", "files": [f"figures/{stem}.svg"],
                          "caption": f"Effect versus -log10 P for {key[1]} ({name}); red points have family q at or below {cutoff} in family {fam}. Unadjusted P is not a discovery criterion.",
                          "not_run": []}
                for fmt in formats:
                    if fmt == "svg":
                        continue
                    target = out / "figures" / f"{stem}.{fmt}"
                    if _r_figure(src, "volcano", target, fmt):
                        record["files"].append(f"figures/{stem}.{fmt}")
                    else:
                        record["not_run"].append(fmt)
                figures.append(record)
    response = root / "response" / "descriptive.tsv"
    if response.is_file():
        rows = [r for r in _table(response) if r["eligibility"] == "tested"]
        src = out / "figure_sources" / "response_axes.tsv"
        _write_tsv(src, ["axis_id", "feature_id", "d", "t", "descriptive_class"], [[r["axis_id"], r["feature_id"], r["d"], r["t"], r["descriptive_class"]] for r in rows])
        pts = [{"x": _num(r["d"]), "y": abs(_num(r["t"]) or 0) if _num(r["t"]) is not None else None, "label": f"{r['feature_id']} ({r['descriptive_class']})",
                "highlight": r["descriptive_class"] in ("near_restoration", "partial_return")} for r in rows]
        (out / "figures").mkdir(parents=True, exist_ok=True)
        (out / "figures" / "response_axes.svg").write_text(scatter_svg(pts, x_label="Disease effect d (log2)", y_label="|treatment effect t| (log2)", title="Descriptive response axes (no P values)"), encoding="utf-8")
        figures.append({"stem": "response_axes", "source": "figure_sources/response_axes.tsv", "files": ["figures/response_axes.svg"], "not_run": [],
                        "caption": "Descriptive disease and treatment effects; highlighted classes are descriptive, not equivalence or rescue."})
    return figures


def _complete_table_html(root: Path, name: str, rows: list[dict]) -> str:
    columns = ["model_id", "contrast_id", "hypothesis_type", "role", "feature_id", "eligibility", "reason_code", "estimable", "effect", "ci_lower", "ci_upper", "p_value", "q_value", "family_id"]
    body = [[r.get(c) if r.get(c) not in ("NA",) else None for c in columns] for r in rows]
    return table(f"Complete {name} table: every planned endpoint ({len(rows)} rows, all contrasts and roles)", columns, body, f"dea/{name}.tsv")


def render(data: dict, root: Path, figures: list[dict]) -> str:
    stub_page = render_stub(data)
    main = stub_page.split("<main>", 1)[1].split("</main>", 1)[0]
    nav = stub_page.split('<nav aria-label="Report sections"><ul>', 1)[1].split("</ul></nav>", 1)[0]
    extra = []
    if data["dea"].get("state") == "COMPLETED":
        content = ""
        for name, info in data["dea"]["tables"].items():
            content += _complete_table_html(root, name, _table(root / info["source"]))
        if (root / "assay_engines" / "assay_results.tsv").is_file():
            content += _complete_table_html(root, "assay_engines", _table(root / "assay_engines" / "assay_results.tsv")).replace('href="../dea/assay_engines.tsv"', 'href="../assay_engines/assay_results.tsv"').replace("dea/assay_engines.tsv", "assay_engines/assay_results.tsv")
        extra.append(("complete-dea", "Complete differential results", content))
    else:
        extra.append(("complete-dea", "Complete differential results", f"<p>{state_span(data['dea'].get('state'))}: no differential rows exist for this run.</p>"))
    fig_html = ""
    for fig in figures:
        svg = next(f for f in fig["files"] if f.endswith(".svg"))
        others = ", ".join(f'<a href="{html.escape(f)}">{html.escape(Path(f).suffix[1:].upper())}</a>' for f in fig["files"])
        missing = f' Not produced: {", ".join(fig["not_run"])} (renderer unavailable).' if fig["not_run"] else ""
        fig_html += (f'<figure><img src="{html.escape(svg)}" alt="{html.escape(fig["caption"])}"><figcaption>{html.escape(fig["caption"])} '
                     f'Files: {others}. Source data: <a href="{html.escape(fig["source"])}">{html.escape(fig["source"])}</a>.{missing}</figcaption></figure>')
    extra.append(("figures", "Figures and their source data", fig_html or "<p>No figures: no completed differential or response results.</p>"))
    extra.append(("methods-full", "Methods and limitations", '<p>Full methods and limitations generated from the executed run: <a href="methods_full.md">methods_full.md</a>.</p>'))
    nav += "".join(f'<li><a href="#{a}">{html.escape(t)}</a></li>' for a, t, _ in extra)
    main += "".join(f'<section id="{a}" aria-labelledby="{a}-h"><h2 id="{a}-h">{html.escape(t)}</h2>{c}</section>' for a, t, c in extra)
    template = Template((Path(__file__).parent / "full_templates" / "report.html").read_text(encoding="utf-8"))
    return template.substitute(title=html.escape(f"Proteomics full report — {data['run']['run_id']}"), phase_label=html.escape("Full offline report generated from verified artifacts."), nav=nav, body=main)


def build_report(run_root: str | Path, snapshot: dict, out: Path, formats: list[str]) -> list[str]:
    root = Path(run_root)
    data = assemble(root, snapshot)
    validate_report_data(data)
    plan = json.loads((root / "plan.json").read_text(encoding="utf-8")) if (root / "plan.json").is_file() else None
    figures = build_figures(root, out, data, formats)
    (out / "report_data.json").write_text(json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    (out / "methods_full.md").write_text(render_full_methods(data, plan), encoding="utf-8")
    from .stub_methods import render_methods
    (out / "methods.md").write_text(render_methods(data), encoding="utf-8")
    (out / "figures.json").write_text(json.dumps(figures, indent=2), encoding="utf-8")
    (out / "index.html").write_text(render(data, root, figures), encoding="utf-8")
    return sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file())


def execute(request: dict) -> dict:
    from ..errors import ProteomicsError
    from ..runtime import stage_result, utc_now, validate_stage_request
    validate_stage_request(request)
    out = Path(request["output_temp_dir"]); out.mkdir(parents=True, exist_ok=False)
    started = utc_now()
    try:
        files = build_report(request["parameters"]["run_root"], request["parameters"]["run_snapshot"], out, list(request["parameters"].get("figure_formats", [])))
    except ProteomicsError as error:
        result = stage_result(request["run_id"], request["stage_id"], CAPABILITY, "FAILED", plan_hash=request["plan_hash"], exit_code=error.exit_code, reason_code=error.code, message=error.message)
        result["started_at"] = started; return result
    except Exception as error:   # a renderer defect fails the report stage only
        result = stage_result(request["run_id"], request["stage_id"], CAPABILITY, "FAILED", plan_hash=request["plan_hash"], exit_code=4, reason_code="E_REPORT_RENDER", message=f"{type(error).__name__}: {error}")
        result["started_at"] = started; return result
    result = stage_result(request["run_id"], request["stage_id"], CAPABILITY, "COMPLETED", plan_hash=request["plan_hash"], exit_code=0, message="full offline report rendered")
    result["started_at"] = started
    result["outputs"] = [{"artifact_id": "full_" + f.replace("/", "_").replace(".", "_")[:100], "relative_path": f, "sha256": sha256_file(out / f), "result_type": "ReportArtifact"} for f in files]
    return result
