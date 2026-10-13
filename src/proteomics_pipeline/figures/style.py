"""Style system: stable group colours (FR-196, V196) and P-value annotation (FR-198, V198).

Group colours are built once per run from the declared group order and an explicit ``group_colours`` mapping, and every
figure of the run reads that same map. The contract (figures.md) fixes the Okabe-Ito order and the P-value rules.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from ..errors import ConfigurationError
from .config import OKABE_ITO, P_ANNOTATIONS


def group_colours(groups: Sequence[str], style: Mapping[str, Any]) -> dict[str, str]:
    """Colour of each group, in declared group order: Okabe-Ito (or the declared palette) by position, then overridden
    by ``group_colours``. Duplicate groups, unknown override groups and duplicate colours are refused (E_FIGURE_PALETTE)."""
    groups = list(groups)
    if len(set(groups)) != len(groups):
        raise ConfigurationError("E_FIGURE_PALETTE", "declared group order lists a group twice", "/groups")
    palette = list(OKABE_ITO) if style["palette"] == "okabe_ito" else list(style["palette"])
    if len(groups) > len(palette):
        raise ConfigurationError("E_FIGURE_PALETTE", f"the palette has {len(palette)} colours for {len(groups)} groups", "/report/figures/palette")
    overrides = dict(style.get("group_colours") or {})
    unknown = sorted(set(overrides) - set(groups))
    if unknown:
        raise ConfigurationError("E_FIGURE_PALETTE", f"group_colours names group(s) that are not declared: {unknown}", "/report/figures/group_colours")
    colours = {group: overrides.get(group, palette[index]) for index, group in enumerate(groups)}
    if len(set(colours.values())) != len(colours):
        raise ConfigurationError("E_FIGURE_PALETTE", "two groups would share one colour", "/report/figures/group_colours")
    return colours


class RunStyle:
    """The resolved style and group colours of one run. Built once; every figure reads the same colour map."""

    def __init__(self, style: Mapping[str, Any], groups: Sequence[str]):
        self.style = dict(style)
        self.groups = tuple(groups)
        self._colours = group_colours(self.groups, self.style)

    def colours_for_figure(self, figure_id: str) -> dict[str, str]:
        """The run's colour map. The figure identifier does not change the map (FR-196)."""
        if not isinstance(figure_id, str) or not figure_id:
            raise ConfigurationError("E_FIGURE_SOURCE_MISSING", "a figure needs an identifier before colours are assigned", "/figure_id")
        return dict(self._colours)


def assert_colours_stable(maps: Mapping[str, Mapping[str, str]]) -> None:
    """Every named colour map must be identical (V196). ``maps`` is keyed by a label such as 'order A / volcano'."""
    reference_label = next(iter(maps), None)
    if reference_label is None:
        return
    reference = dict(maps[reference_label])
    differing = sorted(label for label, colours in maps.items() if dict(colours) != reference)
    if differing:
        raise ConfigurationError("E_FIGURE_COLOUR_UNSTABLE", f"group colours differ from {reference_label!r} in {differing}", "/colours")


def _p_value(p: Any) -> float:
    if p is None or isinstance(p, bool) or not isinstance(p, (int, float)) or math.isnan(float(p)):
        raise ConfigurationError("E_FIGURE_SOURCE_MISSING", "no P value in the source table for this annotation; none is drawn", "/p")
    value = float(p)
    if not 0.0 <= value <= 1.0:
        raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"P value {p!r} in the source table is outside [0, 1]", "/p")
    return value


def p_stars(p: Any) -> str:
    """Stars annotation of contract FR-198: ns, *, **, *** or ****."""
    value = _p_value(p)
    if value >= 0.05:
        return "ns"
    if value >= 0.01:
        return "*"
    if value >= 0.001:
        return "**"
    if value >= 0.0001:
        return "***"
    return "****"


_THRESHOLDS = (0.05, 0.01, 0.001, 0.0001)


def _significant(value: float, digits: int) -> str:
    """``value`` rounded to ``digits`` significant digits, in plain decimal, trailing zeros kept."""
    exponent = int(f"{value:.{digits - 1}e}".split("e")[1])
    return f"{value:.{max(0, digits - 1 - exponent)}f}"


def p_exact(p: Any) -> str:
    """Exact annotation (FR-198, figures.md note): two significant digits in plain decimal, trailing zeros kept
    (0.2 -> 0.20, 0.01 -> 0.010, 0.001 -> 0.0010, 0.05 -> 0.050), floor at 0.0001.

    A value below a significance threshold (0.05, 0.01, 0.001, 0.0001) is never shown as reaching it: if the two-digit
    rounding would land on or above the threshold, more significant digits are used until it stays below
    (0.04999 -> 0.04999, 0.0499 -> 0.0499). A value at or above a threshold is not rounded below it.
    """
    value = _p_value(p)
    if value < 0.0001:
        return "p < 0.0001"
    for digits in range(2, 18):
        text = _significant(value, digits)
        rounded = float(text)
        if not any(value < threshold <= rounded for threshold in _THRESHOLDS):
            break
    return f"p = {text}"


def p_annotation(p: Any, mode: str = "exact") -> str:
    """The annotation of one P value in the declared mode. Stars and exact values never mix."""
    if mode == "stars":
        return p_stars(p)
    if mode == "exact":
        return p_exact(p)
    raise ConfigurationError("E_CONFIG_SCHEMA", f"p_annotation must be one of {list(P_ANNOTATIONS)}", "/report/figures/p_annotation")


def annotate_rows(rows: Iterable[Mapping[str, Any]], mode: str = "exact", p_field: str = "p") -> list[str]:
    """Annotations of a source table's rows, in row order. A row without its P value is refused, never annotated."""
    return [p_annotation(row.get(p_field), mode) for row in rows]
