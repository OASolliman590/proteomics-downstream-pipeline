"""R11 V109: automated numerical checks of report figures and tables.

An independent re-derivation (written here, not imported from the renderer)
maps every figure-source row to SVG pixel coordinates and compares it with the
rendered marks; tick labels, axis labels, source links and raster/vector files
are checked.  Human visual inspection is recorded separately in the R11
receipt (the independent-reviewer inspection is NOT_RUN).
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import re
import struct
from pathlib import Path

import pytest

from proteomics_pipeline import workflow

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore unavailable")

W, H, L, R, T, BM = 520, 380, 64, 16, 40, 52     # frozen layout of the report scatter
CIRCLE = re.compile(r'<circle cx="([-0-9.]+)" cy="([-0-9.]+)" r="3" fill="([^"]+)"[^>]*><title>([^<]*)</title>')


def _tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _f(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def check_figure(svg_text: str, source_rows: list[dict], xcol: str, ycol: str, flag: str | None) -> list[str]:
    """Return a list of findings (empty = consistent) for one scatter figure."""
    findings = []
    pts = [(r["feature_id"], _f(r[xcol]), _f(r[ycol]), r.get(flag) == "true" if flag else False) for r in source_rows]
    pts = [p for p in pts if p[1] is not None and p[2] is not None]
    marks = {m.group(4): (float(m.group(1)), float(m.group(2)), m.group(3)) for m in CIRCLE.finditer(svg_text)}
    if len(marks) != len(pts):
        findings.append(f"{len(marks)} marks for {len(pts)} finite source rows")
    if not pts:
        return findings
    xs, ys = [p[1] for p in pts], [p[2] for p in pts]
    x0, x1 = min(xs + [0.0]), max(xs + [0.0])
    y0, y1 = 0.0, max(ys) or 1.0
    if x0 == x1:
        x0, x1 = x0 - 1, x1 + 1
    colours = set()
    for label, x, y, hi in pts:
        if label not in marks:
            findings.append(f"{label}: no mark"); continue
        cx, cy, colour = marks[label]
        ex = L + (x - x0) / (x1 - x0) * (W - L - R)
        ey = H - BM - (y - y0) / (y1 - y0) * (H - T - BM)
        if abs(cx - ex) > 0.006 or abs(cy - ey) > 0.006:
            findings.append(f"{label}: mark at ({cx},{cy}) but source gives ({ex:.3f},{ey:.3f})")
        colours.add((hi, colour))
    if len({c for h, c in colours if h} & {c for h, c in colours if not h}):
        findings.append("highlighted and non-highlighted points share a colour")
    ticks = re.findall(r'font-size="10" font-family="sans-serif">([^<]+)</text>', svg_text)
    for value in (x0, x1, y1):
        if f"{value:.3g}" not in ticks:
            findings.append(f"tick {value:.3g} missing")
    return findings


def _png_size(path: Path) -> tuple[int, int]:
    data = path.read_bytes()[:24]
    assert data[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    return struct.unpack(">II", data[16:24])


def _run(tmp_path):
    raw = json.loads((B.EXAMPLES / "example-independent.json").read_text())
    for key in ("matrix", "observations", "features", "source_provenance"):
        raw["input"][key] = str((B.EXAMPLES / raw["input"][key]).resolve())
    raw["report"]["figure_formats"] = ["png", "pdf", "svg"]
    path = tmp_path / "analysis.json"; path.write_text(json.dumps(raw))
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert code == 0, payload
    return tmp_path / "run"


def test_v109_figures_labels_coordinates_links_and_files(tmp_path):
    run = _run(tmp_path)
    full = run / "report-full"
    figures = json.loads((full / "figures.json").read_text())
    assert figures, "no figures produced"
    page = (full / "index.html").read_text(encoding="utf-8")
    for fig in figures:
        svg = (full / next(f for f in fig["files"] if f.endswith(".svg"))).read_text(encoding="utf-8")
        rows = _tsv(full / fig["source"])
        if fig["stem"].startswith("volcano_"):
            assert check_figure(svg, rows, "effect_log2", "neg_log10_p", "q_at_or_below_cutoff") == []
            # -log10 P in the source equals the DEA table exactly
            name = "treat" if fig["stem"].startswith("volcano_treat_") else "zero_null"
            dea = {r["feature_id"]: r for r in _tsv(run / "dea" / f"{name}.tsv")
                   if r["eligibility"] == "tested" and f"volcano_{name}_{r['model_id']}_{r['contrast_id']}" == fig["stem"]}
            assert len(dea) == len(rows)
            for r in rows:
                assert r["p_value"] == dea[r["feature_id"]]["p_value"] and r["effect_log2"] == dea[r["feature_id"]]["effect"]
                assert abs(float(r["neg_log10_p"]) + math.log10(float(r["p_value"]))) < 1e-12
            assert "Effect (log2 difference)" in svg and "-log10 P" in svg
        assert f'href="{fig["source"]}"' in page
        for f in fig["files"]:
            assert (full / f).is_file() and (full / f).stat().st_size > 0
            if f.endswith(".png"):
                w, h = _png_size(full / f)
                assert w >= 600 and h >= 400, (f, w, h)
            if f.endswith(".pdf"):
                assert (full / f).read_bytes()[:5] == b"%PDF-"
        assert fig["not_run"] == []
        # readability: every text element fits inside the canvas at ~0.6 em per character
        for size, text in re.findall(r'font-size="(\d+)" font-family="sans-serif">([^<]+)<', svg):
            assert int(size) >= 10 and len(text) * int(size) * 0.6 <= W, text
    for url in re.findall(r'(?:href|src)="([^"#:]+)"', page):
        assert (full / url).exists(), f"broken link {url}"


def test_v109_negative_tampered_source_or_mark_is_detected(tmp_path):
    run = _run(tmp_path)
    full = run / "report-full"
    fig = next(f for f in json.loads((full / "figures.json").read_text()) if f["stem"].startswith("volcano_zero_null"))
    svg = (full / fig["source"].replace("figure_sources", "figures").replace(".tsv", ".svg")).read_text(encoding="utf-8")
    rows = _tsv(full / fig["source"])
    tampered = [dict(r) for r in rows]
    finite = [r for r in tampered if _f(r["effect_log2"]) is not None and _f(r["neg_log10_p"]) is not None]
    mid = sorted(finite, key=lambda r: float(r["effect_log2"]))[len(finite) // 2]
    mid["effect_log2"] = repr(float(mid["effect_log2"]) + 0.05)
    assert any(mid["feature_id"] in f for f in check_figure(svg, tampered, "effect_log2", "neg_log10_p", "q_at_or_below_cutoff"))
    dropped = re.sub(r'<circle [^>]*><title>[^<]*</title></circle>', "", svg, count=1)
    assert check_figure(dropped, rows, "effect_log2", "neg_log10_p", "q_at_or_below_cutoff")
