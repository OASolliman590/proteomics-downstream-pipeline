"""Figure configuration: the optional ``report.figures`` block (packet R16a, FR-195, V195).

The schema (analysis.schema.json) checks the shape and the keys of the block. The values are checked here with the
typed errors of the contract: an unknown format is E_FIGURE_FORMAT, a dpi other than 300 is E_FIGURE_DPI and an invalid
palette is E_FIGURE_PALETTE. Unknown keys are refused as E_CONFIG_SCHEMA. The defaults and the allowed values are those
of specs/017-visualization/contracts/figures.md (width 120 mm by D-74). ``report.figure_formats`` (R10b) is a different key
and is never read here.
"""
from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any

from ..errors import ConfigurationError

FIGURE_FORMATS = ("png", "pzfx", "html")
DPI = 300
WIDTHS_MM = (85, 120, 180)
OKABE_ITO = ("#0072B2", "#E69F00", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#000000")
P_ANNOTATIONS = ("stars", "exact")
FONT = "Arial"
KEYS = ("formats", "dpi", "journal_width_mm", "palette", "group_colours", "font", "p_annotation")
DEFAULTS: dict[str, Any] = {
    "formats": ["png", "pzfx", "html"],
    "dpi": DPI,
    "journal_width_mm": 120,
    "palette": "okabe_ito",
    "group_colours": None,
    "font": FONT,
    "p_annotation": "exact",
}
_HEX = re.compile(r"^#[0-9A-Fa-f]{6}$")


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _resolve_formats(value: Any) -> list[str]:
    if not isinstance(value, (list, tuple)) or not value:
        raise ConfigurationError("E_FIGURE_FORMAT", "report.figures.formats must be a non-empty list of output formats", "/report/figures/formats")
    for item in value:
        if item not in FIGURE_FORMATS:
            raise ConfigurationError("E_FIGURE_FORMAT", f"unknown figure format {item!r}; allowed {list(FIGURE_FORMATS)}", "/report/figures/formats")
    if len(set(value)) != len(value):
        raise ConfigurationError("E_FIGURE_FORMAT", "report.figures.formats lists a format twice", "/report/figures/formats")
    return list(value)


def _resolve_palette(value: Any) -> str | list[str]:
    if value == "okabe_ito":
        return value
    if not isinstance(value, (list, tuple)) or not value:
        raise ConfigurationError("E_FIGURE_PALETTE", "report.figures.palette must be 'okabe_ito' or a non-empty list of hex colours", "/report/figures/palette")
    for colour in value:
        if not isinstance(colour, str) or not _HEX.match(colour):
            raise ConfigurationError("E_FIGURE_PALETTE", f"palette entry {colour!r} is not a #RRGGBB hex colour", "/report/figures/palette")
    if len(set(value)) != len(value):
        raise ConfigurationError("E_FIGURE_PALETTE", "palette lists a colour twice", "/report/figures/palette")
    return list(value)


def _resolve_group_colours(value: Any) -> dict[str, str] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ConfigurationError("E_FIGURE_PALETTE", "report.figures.group_colours must map group identifiers to hex colours", "/report/figures/group_colours")
    for group, colour in value.items():
        if not isinstance(colour, str) or not _HEX.match(colour):
            raise ConfigurationError("E_FIGURE_PALETTE", f"group_colours[{group!r}] is not a #RRGGBB hex colour", f"/report/figures/group_colours/{group}")
    return dict(value)


def resolve_figure_config(block: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return the resolved style. ``None`` (no ``report.figures`` block) resolves to the contract defaults.

    Raises ConfigurationError with a typed code; nothing is written and no figure is produced on refusal.
    """
    if block is None:
        return {key: (list(value) if isinstance(value, list) else value) for key, value in DEFAULTS.items()}
    if not isinstance(block, Mapping):
        raise ConfigurationError("E_CONFIG_SCHEMA", "report.figures must be an object", "/report/figures")
    unknown = sorted(str(key) for key in block if key not in KEYS)
    if unknown:
        raise ConfigurationError("E_CONFIG_SCHEMA", f"unknown key(s) in report.figures: {unknown}; allowed {list(KEYS)}", f"/report/figures/{unknown[0]}")
    values = {key: block.get(key, DEFAULTS[key]) for key in KEYS}
    if not _is_int(values["dpi"]) or values["dpi"] != DPI:
        raise ConfigurationError("E_FIGURE_DPI", f"dpi {values['dpi']!r} is refused; the only allowed value is {DPI}", "/report/figures/dpi")
    if not _is_int(values["journal_width_mm"]) or values["journal_width_mm"] not in WIDTHS_MM:
        raise ConfigurationError("E_CONFIG_SCHEMA", f"journal_width_mm must be one of {list(WIDTHS_MM)}", "/report/figures/journal_width_mm")
    if values["font"] != FONT:
        raise ConfigurationError("E_CONFIG_SCHEMA", f"font must be {FONT!r} (a substitution is recorded in the registry, never declared)", "/report/figures/font")
    if values["p_annotation"] not in P_ANNOTATIONS:
        raise ConfigurationError("E_CONFIG_SCHEMA", f"p_annotation must be one of {list(P_ANNOTATIONS)}", "/report/figures/p_annotation")
    return {
        "formats": _resolve_formats(values["formats"]),
        "dpi": DPI,
        "journal_width_mm": int(values["journal_width_mm"]),
        "palette": _resolve_palette(values["palette"]),
        "group_colours": _resolve_group_colours(values["group_colours"]),
        "font": FONT,
        "p_annotation": values["p_annotation"],
    }
