"""R16d acceptance: V204 (embedded trace data and hover fields) and V205 (offline HTML and pinned plotly.js), FR-204, FR-205.

Oracles, written here and independent of the writer:
- the synthetic volcano source tests/fixtures/figures/plotly/volcano_source.tsv, read with the csv module;
- the trace arrays parsed from the HTML figure JSON with json.loads (not with the writer's own helpers);
- the SHA-256 of vendor/plotly/plotly.min.js computed with hashlib, compared with vendor/plotly/PINNED.json;
- a static scan of the HTML outside the library block, written here with plain string checks.
The HTML has no network reference. The library block is compared with the pinned file by hash, because the minified
library contains URL strings of its own (see plotly_html.py).
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
from pathlib import Path

import pytest

from proteomics_pipeline.errors import ConfigurationError
from proteomics_pipeline.figures.plotly_html import (
    DATA_OPEN, VENDOR_DIR, check_hover, library_block, outside_library, render_html, scan_html, verified_library,
    volcano_traces, write_html,
)

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "tests" / "fixtures" / "figures" / "plotly" / "volcano_source.tsv"
GROUPS = ["Control", "Acute"]
COLOURS = {"Control": "#0072B2", "Acute": "#E69F00"}          # Okabe-Ito, first two (contracts/figures.md)
NAMESPACES = ("http://www.w3.org/2000/svg", "http://www.w3.org/1999/xhtml", "http://www.w3.org/1999/xlink")


def _source() -> list[dict]:
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _embedded_figure(document: str) -> dict:
    start = document.index(DATA_OPEN) + len(DATA_OPEN)
    end = document.index("</script>", start)
    return json.loads(document[start:end])


def _figure_html(tmp_path: Path) -> tuple[Path, str]:
    rows = _source()
    traces = volcano_traces(rows, COLOURS, GROUPS)
    path = write_html(tmp_path / "volcano.html", "Volcano: synthetic", traces, {"xaxis": {"title": "log2 fold change"}})
    return path, path.read_text(encoding="utf-8")


# ---------------------------------------------------------------- V204: embedded data and hover fields

def test_v204_embedded_trace_arrays_equal_the_source_rows(tmp_path):
    path, document = _figure_html(tmp_path)
    figure = _embedded_figure(document)
    rows = _source()
    assert [t["name"] for t in figure["data"]] == GROUPS
    for trace in figure["data"]:
        members = [r for r in rows if r["group"] == trace["name"]]
        # Oracle: x is log2FC, y is -log10 P, customdata is (feature, group, P, q), all from the source rows.
        assert trace["x"] == [float(r["log2FC"]) for r in members]
        assert trace["y"] == [-math.log10(float(r["P"])) for r in members]
        assert trace["customdata"] == [[r["feature"], r["group"], float(r["P"]), float(r["q"])] for r in members]
        assert len(trace["x"]) == len(members) == 20                        # one point per source row


def test_v204_hover_template_names_feature_group_and_p_or_q(tmp_path):
    _, document = _figure_html(tmp_path)
    for trace in _embedded_figure(document)["data"]:
        template = trace["hovertemplate"]
        assert "customdata[0]" in template                                  # feature identifier
        assert "customdata[1]" in template                                  # group
        assert "customdata[2]" in template and "customdata[3]" in template  # P and q
    assert (tmp_path / "volcano.html").is_file()


def test_v204_negative_a_hover_template_without_the_feature_identifier_fails_e_plotly_hover():
    traces = volcano_traces(_source(), COLOURS, GROUPS)
    broken = [dict(traces[0], hovertemplate="group %{customdata[1]} P = %{customdata[2]}<extra></extra>")]
    with pytest.raises(ConfigurationError) as caught:
        check_hover(broken)
    assert caught.value.code == "E_PLOTLY_HOVER"
    with pytest.raises(ConfigurationError) as written:
        render_html("broken", broken, {})
    assert written.value.code == "E_PLOTLY_HOVER"


# ---------------------------------------------------------------- V205: offline HTML and pinned plotly.js

def test_v205_vendored_hash_equals_the_pin_and_the_licence_is_present():
    pin = json.loads((VENDOR_DIR / "PINNED.json").read_text(encoding="utf-8"))
    digest = hashlib.sha256((VENDOR_DIR / "plotly.min.js").read_bytes()).hexdigest()
    assert digest == pin["sha256"]
    assert pin["licence"] == "MIT"
    licence = (VENDOR_DIR / "LICENSE").read_text(encoding="utf-8")
    assert "Permission is hereby granted, free of charge" in licence
    assert hashlib.sha256((VENDOR_DIR / "LICENSE").read_bytes()).hexdigest() == pin["licence_sha256"]


def test_v205_static_scan_of_the_html_outside_the_library_block(tmp_path):
    _, document = _figure_html(tmp_path)
    outside = outside_library(document)
    # Oracle: plain string and regex checks written here.
    assert not re.search(r"<script\b[^>]*\bsrc\s*=", outside, re.I)
    assert not re.search(r"<link\b[^>]*\bhref\s*=", outside, re.I)
    assert "fetch(" not in outside and "XMLHttpRequest" not in outside
    assert not re.search(r"\bimport\s*\(", outside)
    for url in re.findall(r"https?://[^\s\"'<>)]*", outside):
        assert url.startswith(NAMESPACES), url
    assert "Permission is hereby granted, free of charge" in document          # the MIT notice is in the HTML head


def test_v205_library_block_equals_the_pinned_file_by_hash(tmp_path):
    _, document = _figure_html(tmp_path)
    pin = json.loads((VENDOR_DIR / "PINNED.json").read_text(encoding="utf-8"))
    assert hashlib.sha256(library_block(document).encode("utf-8")).hexdigest() == pin["sha256"]
    summary = scan_html(document)
    assert summary["library_sha256_matches_pin"] is True


def test_v205_negative_a_cdn_script_tag_fails_e_plotly_network(tmp_path):
    _, document = _figure_html(tmp_path)
    injected = document.replace("</body>", '<script src="https://cdn.plot.ly/plotly-2.25.2.min.js"></script>\n</body>')
    with pytest.raises(ConfigurationError) as caught:
        scan_html(injected)
    assert caught.value.code == "E_PLOTLY_NETWORK"


def test_v205_negative_an_external_fetch_outside_the_library_fails_e_plotly_network(tmp_path):
    _, document = _figure_html(tmp_path)
    injected = document.replace("</body>", "<script>fetch('https://example.invalid/x');</script>\n</body>")
    with pytest.raises(ConfigurationError) as caught:
        scan_html(injected)
    assert caught.value.code == "E_PLOTLY_NETWORK"


def test_v205_negative_a_modified_vendored_file_fails_e_plotly_hash(tmp_path):
    vendor = tmp_path / "vendor"
    shutil.copytree(VENDOR_DIR, vendor)
    (vendor / "plotly.min.js").write_bytes((vendor / "plotly.min.js").read_bytes() + b"\n;")
    with pytest.raises(ConfigurationError) as caught:
        verified_library(vendor)
    assert caught.value.code == "E_PLOTLY_HASH"
    out = tmp_path / "out.html"
    with pytest.raises(ConfigurationError):
        write_html(out, "modified", volcano_traces(_source(), COLOURS, GROUPS), {}, vendor_dir=vendor)
    assert not out.exists()
