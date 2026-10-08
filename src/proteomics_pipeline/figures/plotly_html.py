"""Standalone interactive HTML per figure (packet R16d, FR-204, FR-205, V204, V205; ADR 0012).

The HTML inlines the vendored plotly.js (vendor/plotly/plotly.min.js) after its SHA-256 is compared with the pin in
vendor/plotly/PINNED.json (E_PLOTLY_HASH otherwise). The figure data are written as JSON built here, without the Python
plotly package. The MIT notice of plotly.js is included in the HTML head.

The static scan (contracts/figures.md, "Vendored plotly.js") covers the document outside the library block. The library
block itself is checked by hash equality with the pinned file, because the minified library contains URL strings of its
own (map tiles, documentation links). Outside the block the document may contain only the XML namespace strings
http://www.w3.org/2000/svg, http://www.w3.org/1999/xhtml and http://www.w3.org/1999/xlink.
"""
from __future__ import annotations

import hashlib
import html as html_lib
import json
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from ..errors import ConfigurationError

VENDOR_DIR = Path(__file__).resolve().parents[3] / "vendor" / "plotly"
LIBRARY_FILE = "plotly.min.js"
PIN_FILE = "PINNED.json"
LICENCE_FILE = "LICENSE"
LIBRARY_OPEN = '<script id="vendored-plotly-js" type="text/javascript">\n'
LIBRARY_CLOSE = "\n</script>\n"
DATA_OPEN = '<script id="figure-data" type="application/json">'
NAMESPACE_ALLOWLIST = ("http://www.w3.org/2000/svg", "http://www.w3.org/1999/xhtml", "http://www.w3.org/1999/xlink")
_URL = re.compile(r"https?://[^\s\"'<>)]*")
_FORBIDDEN = [
    (re.compile(r"<script\b[^>]*\bsrc\s*=", re.I), "a <script src> reference"),
    (re.compile(r"<link\b[^>]*\bhref\s*=", re.I), "a <link href> reference"),
    (re.compile(r"\bfetch\s*\("), "a fetch() call"),
    (re.compile(r"\bXMLHttpRequest\b"), "an XMLHttpRequest"),
    (re.compile(r"\bimport\s*\("), "a dynamic import()"),
]


def read_pin(vendor_dir: Path = VENDOR_DIR) -> dict[str, Any]:
    pin = vendor_dir / PIN_FILE
    if not pin.is_file():
        raise ConfigurationError("E_PLOTLY_HASH", "the pinned plotly.js record PINNED.json is missing", "/vendor/plotly/PINNED.json")
    return json.loads(pin.read_text(encoding="utf-8"))


def verified_library(vendor_dir: Path = VENDOR_DIR) -> bytes:
    """The vendored plotly.js bytes, after the hash check against the pin (FR-205). Raises E_PLOTLY_HASH on a mismatch."""
    library = vendor_dir / LIBRARY_FILE
    if not library.is_file():
        raise ConfigurationError("E_PLOTLY_HASH", "the vendored plotly.js file is missing", "/vendor/plotly/plotly.min.js")
    data = library.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    pinned = read_pin(vendor_dir)["sha256"]
    if digest != pinned:
        raise ConfigurationError("E_PLOTLY_HASH", f"vendored plotly.js sha256 {digest} differs from the pin {pinned}", "/vendor/plotly/plotly.min.js")
    return data


def licence_text(vendor_dir: Path = VENDOR_DIR) -> str:
    return (vendor_dir / LICENCE_FILE).read_text(encoding="utf-8")


def volcano_traces(rows: Sequence[Mapping[str, Any]], group_colours: Mapping[str, str], groups: Sequence[str]) -> list[dict[str, Any]]:
    """One scatter trace per group, in declared group order. x is the effect, y is -log10 P, and customdata holds
    feature, group, P and q for every point. Hover text names the feature, the group, P and q (FR-204)."""
    traces = []
    for group in groups:
        members = [r for r in rows if r["group"] == group]
        if not members:
            continue
        traces.append({
            "type": "scatter", "mode": "markers", "name": group,
            "x": [float(r["log2FC"]) for r in members],
            "y": [-_log10(float(r["P"])) for r in members],
            "customdata": [[str(r["feature"]), str(r["group"]), float(r["P"]), float(r["q"])] for r in members],
            "hovertemplate": ("<b>%{customdata[0]}</b><br>group: %{customdata[1]}<br>P = %{customdata[2]:.3g}"
                              "<br>q = %{customdata[3]:.3g}<extra></extra>"),
            "marker": {"color": group_colours[group], "size": 7, "opacity": 0.85},
        })
    return traces


def _log10(value: float) -> float:
    if not value > 0:
        raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"P value {value!r} cannot be plotted on a -log10 scale", "/P")
    return math.log10(value)


def check_hover(traces: Sequence[Mapping[str, Any]]) -> None:
    """Every trace's hover template names the feature identifier (customdata[0]), the group (customdata[1]) and P or q
    (customdata[2] or customdata[3]). Otherwise E_PLOTLY_HOVER."""
    for index, trace in enumerate(traces):
        template = str(trace.get("hovertemplate", ""))
        if trace.get("type") == "heatmap":
            # a heatmap names its feature and sample through the axis positions of the cell (R16e amendment A-2026-10-01-37)
            if "%{x}" not in template or "%{y}" not in template:
                raise ConfigurationError("E_PLOTLY_HOVER", f"trace {index} hover text lacks the row or column identifier", f"/data/{index}/hovertemplate")
            continue
        if "customdata[0]" not in template or "customdata[1]" not in template or ("customdata[2]" not in template and "customdata[3]" not in template):
            raise ConfigurationError("E_PLOTLY_HOVER", f"trace {index} hover text lacks the feature, group or P/q field", f"/data/{index}/hovertemplate")


def figure_json(traces: Sequence[Mapping[str, Any]], layout: Mapping[str, Any]) -> str:
    """The figure as JSON. '<', '>' and '&' are escaped so the data cannot end the script element."""
    text = json.dumps({"data": list(traces), "layout": dict(layout)}, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    return text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def render_html(title: str, traces: Sequence[Mapping[str, Any]], layout: Mapping[str, Any], vendor_dir: Path = VENDOR_DIR) -> str:
    library = verified_library(vendor_dir)
    check_hover(traces)
    pin = read_pin(vendor_dir)
    notice = licence_text(vendor_dir)
    document = (
        "<!doctype html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{html_lib.escape(title)}</title>\n"
        f"<!--\nplotly.js {pin['version']} (vendor/plotly/plotly.min.js), MIT licence notice (vendor/plotly/LICENSE, unmodified):\n"
        f"{notice}-->\n"
        "<style>html,body{margin:0;background:#fff;font-family:Arial,Helvetica,sans-serif}#figure{width:100vw;height:100vh}</style>\n"
        "</head>\n<body>\n<div id=\"figure\"></div>\n"
        + LIBRARY_OPEN + library.decode("utf-8") + LIBRARY_CLOSE
        + DATA_OPEN + figure_json(traces, layout) + "</script>\n"
        + '<script type="text/javascript">\n'
        + '(function () { var f = JSON.parse(document.getElementById("figure-data").textContent); '
          'Plotly.newPlot("figure", f.data, f.layout, {responsive: true, displaylogo: false}); })();\n'
        + "</script>\n</body>\n</html>\n"
    )
    return document


def library_block(document: str) -> str:
    """The text of the vendored library block, between its markers."""
    start = document.find(LIBRARY_OPEN)
    if start < 0:
        raise ConfigurationError("E_PLOTLY_NETWORK", "the library block marker is missing", "/head")
    body_start = start + len(LIBRARY_OPEN)
    end = document.find(LIBRARY_CLOSE, body_start)
    if end < 0:
        raise ConfigurationError("E_PLOTLY_NETWORK", "the library block is not closed", "/head")
    return document[body_start:end]


def outside_library(document: str) -> str:
    start = document.find(LIBRARY_OPEN)
    end = document.find(LIBRARY_CLOSE, start + len(LIBRARY_OPEN))
    return document[:start] + document[end + len(LIBRARY_CLOSE):]


def scan_html(document: str, vendor_dir: Path = VENDOR_DIR) -> dict[str, Any]:
    """The static scan. The library block must equal the pinned file (E_PLOTLY_HASH). Outside it, there must be no
    external script, stylesheet, fetch, XMLHttpRequest or dynamic import (E_PLOTLY_NETWORK), and every http(s) string
    must be one of the XML namespace strings. Returns what was checked."""
    block = library_block(document)
    library = verified_library(vendor_dir)
    if block.encode("utf-8") != library:
        raise ConfigurationError("E_PLOTLY_HASH", "the library block in the HTML differs from the vendored plotly.js", "/head")
    rest = outside_library(document)
    for pattern, what in _FORBIDDEN:
        if pattern.search(rest):
            raise ConfigurationError("E_PLOTLY_NETWORK", f"the HTML contains {what}", "/body")
    for match in _URL.finditer(rest):
        url = match.group(0)
        if not any(url.startswith(ns) for ns in NAMESPACE_ALLOWLIST):
            raise ConfigurationError("E_PLOTLY_NETWORK", f"the HTML contains an external URL outside the library block: {url[:60]}", "/body")
    return {"library_bytes": len(library), "library_sha256_matches_pin": True, "outside_library_chars": len(rest)}


def write_html(path: str | Path, title: str, traces: Sequence[Mapping[str, Any]], layout: Mapping[str, Any], vendor_dir: Path = VENDOR_DIR) -> Path:
    """Write the standalone HTML (UTF-8, LF). The static scan runs over the complete document before it is written."""
    document = render_html(title, traces, layout, vendor_dir)
    scan_html(document, vendor_dir)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(document.encode("utf-8"))
    return destination
