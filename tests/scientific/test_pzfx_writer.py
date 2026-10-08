"""R16b acceptance: V199 (Column tables), V200 (XY tables) and V201 (Grouped tables), FR-199 to FR-201.

Oracles, written in this file and independent of the writer:
- the source values, read from the synthetic TSVs in tests/fixtures/figures/pzfx with the csv module;
- the .pzfx XML, parsed here with xml.etree.ElementTree (not with the writer's read-back).
The optional R read-back with the CRAN package pzfx is in r/proteomicsCore/tests/testthat/test-pzfx-readback.R.
Opening the files in GraphPad Prism is a manual Maintainer check and is NOT_RUN in automated tests.
"""
from __future__ import annotations

import csv
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from proteomics_pipeline.figures.pzfx_writer import column_table, grouped_table, pzfx_bytes, write_pzfx, xy_table

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "figures" / "pzfx"

# Oracle: Prism's TableType codes for the three kinds of table (Column = OneWay, XY = XY, Grouped = TwoWay).
PRISM_CODE = {"Column": "OneWay", "XY": "XY", "Grouped": "TwoWay"}


def _tsv(name: str) -> list[dict[str, str]]:
    with (FIXTURES / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _oracle_tables(path: Path) -> list[dict]:
    """Independent parse of a .pzfx file: every table with its type, its title, its data sets and its row titles."""
    tables = []
    for node in ET.parse(path).getroot().findall("Table"):
        columns: dict[str, list] = {}
        order: list[str] = []
        for child in node:
            if child.tag in ("XColumn", "YColumn"):
                title = child.find("Title").text
                cells = [d.text for d in child.find("Subcolumn").findall("d")]
                columns[title] = [None if c is None else float(c) for c in cells]
                order.append(title)
        row_node = node.find("RowTitlesColumn")
        rows = [d.text for d in row_node.find("Subcolumn").findall("d")] if row_node is not None else None
        tables.append({"id": node.get("ID"), "type": node.get("TableType"), "xformat": node.get("XFormat"),
                       "title": node.find("Title").text, "columns": columns, "order": order, "rows": rows})
    return tables


def _mismatches(expected: dict[str, list], parsed: dict[str, list]) -> list[str]:
    """Every difference between the source (expected) and the file (parsed), exact comparison."""
    problems = []
    for name in sorted(set(expected) - set(parsed)):
        problems.append(f"data set {name!r} is missing from the file")
    for name in sorted(set(parsed) - set(expected)):
        problems.append(f"data set {name!r} is not in the source")
    for name in sorted(set(expected) & set(parsed)):
        if len(expected[name]) != len(parsed[name]):
            problems.append(f"{name!r}: {len(parsed[name])} values, source has {len(expected[name])}")
            continue
        for index, (want, got) in enumerate(zip(expected[name], parsed[name])):
            if want is None and got is None:
                continue
            if want is None or got is None or want != got:
                problems.append(f"{name!r}[{index}]: file {got!r} != source {want!r}")
    return problems


def _fixture_bytes(name: str) -> bytes:
    return (FIXTURES / name).read_bytes()


# ---------------------------------------------------------------- V199: column tables

def _column_source() -> tuple[list[str], list[str], dict[str, dict[str, list[float]]]]:
    rows = _tsv("column_source.tsv")
    proteins = list(dict.fromkeys(r["protein"] for r in rows))
    groups = list(dict.fromkeys(r["group"] for r in rows))
    values = {p: {g: [float(r["value"]) for r in rows if r["protein"] == p and r["group"] == g] for g in groups} for p in proteins}
    return proteins, groups, values


def test_v199_column_tables_hold_the_source_values_exactly(tmp_path):
    proteins, groups, source = _column_source()
    tables = [column_table(p, {g: source[p][g] for g in groups}) for p in proteins]
    path = write_pzfx(tmp_path / "column_groups.pzfx", tables)
    parsed = _oracle_tables(path)
    assert [t["type"] for t in parsed] == [PRISM_CODE["Column"]] * len(proteins)
    for protein, table in zip(proteins, parsed):
        assert table["title"] == protein and table["order"] == groups
        assert _mismatches({g: source[protein][g] for g in groups}, table["columns"]) == []
        assert len(table["columns"]["Control"]) == 6          # six replicate rows per group
    assert path.read_bytes() == _fixture_bytes("column_groups.pzfx"), "the committed fixture is the writer output"


def test_v199_round_trip_precision_is_full_float(tmp_path):
    # repr() of a float reads back as the same float; a 17-digit value survives unchanged
    values = [0.1 + 0.2, 1e-300, 123456789.123456789, -2.5e-7]
    path = write_pzfx(tmp_path / "precision.pzfx", [column_table("precision", {"G": values})])
    assert _oracle_tables(path)[0]["columns"]["G"] == values


def test_v199_missing_values_are_empty_cells(tmp_path):
    path = write_pzfx(tmp_path / "missing.pzfx", [column_table("missing", {"A": [1.0, float("nan"), 3.0], "B": [None, 2.0, 2.5]})])
    parsed = _oracle_tables(path)[0]["columns"]
    assert parsed == {"A": [1.0, None, 3.0], "B": [None, 2.0, 2.5]}
    assert _mismatches({"A": [1.0, None, 3.0], "B": [None, 2.0, 2.5]}, parsed) == []


def test_v199_negative_single_changed_value_fails_the_comparison(tmp_path):
    proteins, groups, source = _column_source()
    path = write_pzfx(tmp_path / "column.pzfx", [column_table(proteins[0], {g: source[proteins[0]][g] for g in groups})])
    parsed = _oracle_tables(path)[0]["columns"]
    assert _mismatches(source[proteins[0]], parsed) == []
    changed = {k: list(v) for k, v in parsed.items()}
    changed["Acute"][2] = changed["Acute"][2] + 1e-9
    assert _mismatches(source[proteins[0]], changed) != []


def test_v199_negative_rounding_to_three_decimals_fails_the_comparison(tmp_path):
    proteins, groups, source = _column_source()
    rounded = {g: [round(v, 3) for v in source[proteins[0]][g]] for g in groups}
    path = write_pzfx(tmp_path / "rounded.pzfx", [column_table(proteins[0], rounded)])
    assert _mismatches(source[proteins[0]], _oracle_tables(path)[0]["columns"]) != []


# ---------------------------------------------------------------- V200: XY tables

def _volcano_and_roc() -> tuple[dict, dict]:
    volcano = _tsv("volcano_source.tsv")
    roc = _tsv("roc_source.tsv")
    return ({"log2FC": [float(r["log2FC"]) for r in volcano], "neg_log10P": [float(r["neg_log10P"]) for r in volcano]},
            {"one_minus_specificity": [float(r["one_minus_specificity"]) for r in roc], "sensitivity": [float(r["sensitivity"]) for r in roc]})


def _xy_tables(volcano: dict, roc: dict) -> list:
    return [xy_table("volcano", "log2FC", "neg_log10P", volcano["log2FC"], volcano["neg_log10P"]),
            xy_table("roc", "one_minus_specificity", "sensitivity", roc["one_minus_specificity"], roc["sensitivity"])]


def test_v200_xy_tables_hold_every_source_row_exactly(tmp_path):
    volcano, roc = _volcano_and_roc()
    assert len(volcano["log2FC"]) == 200 and len(roc["sensitivity"]) == 30
    path = write_pzfx(tmp_path / "xy.pzfx", _xy_tables(volcano, roc))
    parsed = _oracle_tables(path)
    assert [t["type"] for t in parsed] == [PRISM_CODE["XY"], PRISM_CODE["XY"]]
    assert parsed[0]["order"] == ["log2FC", "neg_log10P"] and parsed[1]["order"] == ["one_minus_specificity", "sensitivity"]
    assert _mismatches(volcano, parsed[0]["columns"]) == []
    assert _mismatches(roc, parsed[1]["columns"]) == []
    assert len(parsed[0]["columns"]["log2FC"]) == 200                 # one point per source row
    assert path.read_bytes() == _fixture_bytes("xy_volcano_roc.pzfx"), "the committed fixture is the writer output"


def test_v200_negative_values_rounded_to_three_decimals_fail(tmp_path):
    volcano, roc = _volcano_and_roc()
    rounded_volcano = {k: [round(v, 3) for v in vals] for k, vals in volcano.items()}
    path = write_pzfx(tmp_path / "rounded_xy.pzfx", _xy_tables(rounded_volcano, roc))
    parsed = _oracle_tables(path)
    assert _mismatches(volcano, parsed[0]["columns"]) != []
    assert _mismatches(roc, parsed[1]["columns"]) == []


# ---------------------------------------------------------------- V201: grouped tables

def _grouped_source() -> tuple[list[str], list[str], list[str], dict[str, list[float]]]:
    rows = _tsv("grouped_source.tsv")
    replicates = sorted({r["replicate"] for r in rows}, key=int)
    proteins = list(dict.fromkeys(r["protein"] for r in rows))
    groups = list(dict.fromkeys(r["group"] for r in rows))
    datasets = {f"{p} {g}": [float(r["value"]) for rep in replicates for r in rows
                             if r["protein"] == p and r["group"] == g and r["replicate"] == rep] for p in proteins for g in groups}
    return replicates, proteins, groups, datasets


def test_v201_grouped_table_rows_are_replicates_and_datasets_are_protein_by_group(tmp_path):
    replicates, proteins, groups, datasets = _grouped_source()
    assert len(proteins) == 4 and len(groups) == 2 and len(replicates) == 6
    table = grouped_table("panel", [f"rep{r}" for r in replicates], datasets)
    path = write_pzfx(tmp_path / "grouped.pzfx", [table])
    parsed = _oracle_tables(path)[0]
    assert parsed["type"] == PRISM_CODE["Grouped"]
    assert parsed["rows"] == [f"rep{r}" for r in replicates]
    assert len(parsed["order"]) == len(proteins) * len(groups)          # one data set per protein and group
    assert _mismatches(datasets, parsed["columns"]) == []
    assert path.read_bytes() == _fixture_bytes("grouped_panel.pzfx"), "the committed fixture is the writer output"


def test_v201_negative_missing_protein_data_set_fails_the_comparison(tmp_path):
    replicates, proteins, groups, datasets = _grouped_source()
    without_p2 = {name: values for name, values in datasets.items() if not name.startswith("P2 ")}
    path = write_pzfx(tmp_path / "missing_p2.pzfx", [grouped_table("panel", [f"rep{r}" for r in replicates], without_p2)])
    problems = _mismatches(datasets, _oracle_tables(path)[0]["columns"])
    assert problems and all("P2" in p for p in problems)


def test_v201_writer_refuses_a_row_count_that_does_not_match():
    with pytest.raises(ValueError, match="1 values for 2 rows"):
        grouped_table("panel", ["rep1", "rep2"], {"P1 Control": [1.0]})


def test_v199_to_v201_writer_is_deterministic():
    _, _, source = _column_source()
    first = pzfx_bytes([column_table("P1", {"Control": source["P1"]["Control"]})])
    second = pzfx_bytes([column_table("P1", {"Control": source["P1"]["Control"]})])
    assert first == second and b"\r" not in first and first.startswith(b'<?xml version="1.0" encoding="UTF-8"?>\n<GraphPadPrismFile')
    assert not any(math.isnan(v) for v in source["P1"]["Control"])
