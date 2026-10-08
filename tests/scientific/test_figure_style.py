"""R16a acceptance: V196 (stable Okabe-Ito group colours) and V198 (P-value annotation), FR-196 and FR-198.

Oracles are written here: the Okabe-Ito hex values and the declared group order come from contracts/figures.md, and the
P-value labels are the rows of the annotation table of the same contract. The P-value fixture is synthetic.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from proteomics_pipeline.errors import ConfigurationError
from proteomics_pipeline.figures import RunStyle, annotate_rows, resolve_figure_config
from proteomics_pipeline.figures.style import assert_colours_stable, p_annotation

ROOT = Path(__file__).resolve().parents[2]
P_FIXTURE = ROOT / "tests" / "fixtures" / "figures" / "style" / "p_values.tsv"

# Oracle: Okabe-Ito palette in contract order (figures.md, "Style defaults").
OKABE_ITO_ORACLE = ["#0072B2", "#E69F00", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]
# Oracle: declared group order of the synthetic fixture.
GROUPS = ["Control", "Acute", "Chronic"]
FIGURE_TYPES = ["diff_volcano_labelled", "diff_protein_dotplot_brackets", "bio_roc_ci_band"]


def _maps(run: RunStyle, order: list[str]) -> dict[str, dict[str, str]]:
    return {f"order {'-'.join(order)} / {fid}": run.colours_for_figure(fid) for fid in order}


def test_v196_colours_equal_the_oracle_and_are_identical_across_figure_orders():
    run = RunStyle(resolve_figure_config(None), GROUPS)
    expected = {group: OKABE_ITO_ORACLE[index] for index, group in enumerate(GROUPS)}
    forward = _maps(run, FIGURE_TYPES)
    backward = _maps(run, list(reversed(FIGURE_TYPES)))
    for label, colours in {**forward, **backward}.items():
        assert colours == expected, label
    assert_colours_stable({**forward, **backward})
    # the first three groups are blue, orange and vermillion, as previewed (D-74)
    assert [expected[g] for g in GROUPS] == ["#0072B2", "#E69F00", "#D55E00"]


def _order_dependent_colours(order: list[str]) -> dict[str, dict[str, str]]:
    """Deliberately wrong assigner for the negative case: a figure's colours depend on its position in the run."""
    return {fid: {group: OKABE_ITO_ORACLE[(index + position) % len(OKABE_ITO_ORACLE)] for index, group in enumerate(GROUPS)}
            for position, fid in enumerate(order)}


def test_v196_negative_order_dependent_colours_fail_e_figure_colour_unstable():
    wrong = {**{f"forward / {k}": v for k, v in _order_dependent_colours(FIGURE_TYPES).items()},
             **{f"backward / {k}": v for k, v in _order_dependent_colours(list(reversed(FIGURE_TYPES))).items()}}
    with pytest.raises(ConfigurationError) as caught:
        assert_colours_stable(wrong)
    assert caught.value.code == "E_FIGURE_COLOUR_UNSTABLE"


def test_v196_explicit_group_colours_override_one_group_and_stay_stable():
    style = resolve_figure_config({"group_colours": {"Chronic": "#009E73"}})
    run = RunStyle(style, GROUPS)
    expected = {"Control": "#0072B2", "Acute": "#E69F00", "Chronic": "#009E73"}
    assert run.colours_for_figure(FIGURE_TYPES[0]) == expected
    assert run.colours_for_figure(FIGURE_TYPES[2]) == expected


def test_v196_negative_more_groups_than_palette_colours_fail_e_figure_palette():
    too_many = [f"G{i}" for i in range(len(OKABE_ITO_ORACLE) + 1)]
    with pytest.raises(ConfigurationError) as caught:
        RunStyle(resolve_figure_config(None), too_many)
    assert caught.value.code == "E_FIGURE_PALETTE"


def _p_fixture() -> list[dict]:
    with P_FIXTURE.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row["p"] = float(row["p"])
    return rows


def test_v198_fixture_annotations_equal_the_contract_table():
    rows = _p_fixture()
    assert [row["p"] for row in rows] == [0.2, 0.04, 0.0004, 1e-7]
    # Oracle: stars ns, *, ***, **** and exact p = 0.20, p = 0.04, p = 0.0004, p < 0.0001 (contract FR-198 table).
    assert annotate_rows(rows, "stars") == ["ns", "*", "***", "****"]
    assert annotate_rows(rows, "exact") == ["p = 0.20", "p = 0.04", "p = 0.0004", "p < 0.0001"]


@pytest.mark.parametrize(("p", "stars"), [
    (0.05, "ns"), (0.2, "ns"), (0.0499, "*"), (0.01, "*"), (0.0099, "**"), (0.001, "**"),
    (0.0009, "***"), (0.0001, "***"), (0.00009, "****"), (1e-12, "****"),
])
def test_v198_star_thresholds_follow_the_contract(p, stars):
    assert p_annotation(p, "stars") == stars


@pytest.mark.parametrize(("p", "exact"), [
    (0.05, "p = 0.05"), (0.0001, "p = 0.0001"), (0.00009, "p < 0.0001"), (1.0, "p = 1.00"), (0.0999, "p = 0.10"),
])
def test_v198_exact_rounding_follows_the_contract(p, exact):
    assert p_annotation(p, "exact") == exact


def test_v198_negative_stars_never_appear_in_exact_mode():
    exact = annotate_rows(_p_fixture(), "exact")
    assert all("*" not in label and label != "ns" for label in exact)
    assert exact != annotate_rows(_p_fixture(), "stars")


@pytest.mark.parametrize("missing", [None, float("nan"), True, -0.1, 1.5])
def test_v198_negative_no_annotation_without_a_p_value_in_the_source(missing):
    with pytest.raises(ConfigurationError) as caught:
        p_annotation(missing, "exact")
    assert caught.value.code == "E_FIGURE_SOURCE_MISSING"
    with pytest.raises(ConfigurationError):
        annotate_rows([{"label": "no p column"}], "stars")
