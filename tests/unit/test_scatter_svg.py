"""Open item 10 (V109 layout findings): the SVG scatter/volcano marks overplotted identical points, pads the axis
ranges and draws a zero reference line; every plotted coordinate is still an exact source value."""
from __future__ import annotations

import re

from proteomics_pipeline.reporting.full import scatter_svg


def _svg(points):
    return scatter_svg(points, x_label="effect", y_label="-log10 P", title="t", hline=1.3)


def test_identical_points_are_marked_with_their_count():
    svg = _svg([{"x": 0.5, "y": 2.0, "label": "F1"}, {"x": 0.5, "y": 2.0, "label": "F2"}, {"x": -1.0, "y": 0.4, "label": "F3"}])
    circles = re.findall(r"<circle [^>]*>.*?</circle>", svg)
    assert len(circles) == 2                                                      # one mark per distinct coordinate
    dup = next(c for c in circles if 'data-n="2"' in c)
    assert "2 identical points: F1, F2" in dup and 'stroke="black"' in dup


def test_ranges_are_padded_and_zero_is_marked():
    svg = _svg([{"x": -1.0, "y": 0.0, "label": "a"}, {"x": 2.0, "y": 3.0, "label": "b"}])
    xs = [float(v) for v in re.findall(r'<circle cx="([0-9.]+)"', svg)]
    assert min(xs) > 64 and max(xs) < 520 - 16                                    # no point sits on the plot frame
    assert re.search(r'<line [^>]*class="zero"', svg)                             # x = 0 reference line
