"""Visualization outputs (ADR 0012, specs/017-visualization): figure configuration, style, registry and writers.

R16a provides the configuration (FR-195), stable group colours (FR-196), the figure registry (FR-197) and P-value
annotation (FR-198). The Prism writer (R16b) lives in ``pzfx_writer``.
"""
from .config import DEFAULTS, FIGURE_FORMATS, OKABE_ITO, resolve_figure_config
from .registry import CATALOGUE, FigureRegistry
from .style import RunStyle, annotate_rows, assert_colours_stable, group_colours, p_annotation

__all__ = [
    "CATALOGUE", "DEFAULTS", "FIGURE_FORMATS", "FigureRegistry", "OKABE_ITO", "RunStyle",
    "annotate_rows", "assert_colours_stable", "group_colours", "p_annotation", "resolve_figure_config",
]
