"""R16b acceptance: V199 (Column tables with the D-75 statistics table and Info sheet), V200 (XY tables) and V201 (Grouped
tables), FR-199 to FR-201, and the Info-sheet rules of contracts/figures.md (D-75).

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

from proteomics_pipeline.errors import ConfigurationError
from proteomics_pipeline.figures.pzfx_writer import (
    INFO_CONSTANTS, column_table, grouped_table, info_sheet, pzfx_bytes, statistics_table, write_pzfx, xy_table,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "figures" / "pzfx"
GROUPS = ["Control", "Acute"]

# Oracle: Prism's TableType codes for the three kinds of table (Column = OneWay, XY = XY, Grouped = TwoWay).
PRISM_CODE = {"Column": "OneWay", "XY": "XY", "Grouped": "TwoWay"}


def _tsv(name: str) -> list[dict[str, str]]:
    with (FIXTURES / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _sha(name: str) -> str:
    import hashlib
    return hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest()


def _info(statistics: str = "table", test: str = "Welch t (fixture)", source_sha: str = "0" * 64) -> dict:
    return info_sheet("fixture figure", {"Test": test, "Multiplicity adjustment": "BH within fixture-family (synthetic)",
                                         "Effect scale": "log2 difference", "Group sizes": "Control 6; Acute 6",
                                         "Pipeline version and commit": "0.1.0.dev0 (fixture)", "Source table SHA-256": source_sha,
                                         "Run ID": "fixture-run", "Created (UTC)": "2026-10-08T00:00:00Z", "Statistics": statistics},
                      notes=["Synthetic fixture for the pzfx writer.", "A test run inside Prism may give a different P value from the pipeline's test."])


def _oracle_tables(path: Path) -> list[dict]:
    """Independent parse of a .pzfx file: every table with its type, title, data sets and row titles."""
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


def _oracle_info(path: Path) -> dict:
    """Independent parse of the Info sheet: title, notes and the name/value constants."""
    info = ET.parse(path).getroot().find("Info")
    return {"title": info.find("Title").text, "notes": [f.text for f in info.find("Notes").findall("Font")],
            "constants": {c.find("Name").text: c.find("Value").text for c in info.findall("Constant")}}


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


def _stats_comparisons(rows: list[dict[str, str]], groups: list[str]) -> list[dict]:
    """The statistics source of the fixture as comparison mappings (values parsed from the TSV text)."""
    return [{"row": r["comparison"], "effect": float(r["effect"]), "P": float(r["P"]), "q": float(r["q"]),
             "n": {g: float(r[f"n_{g}"]) for g in groups}} for r in rows]


def _expected_stats(rows: list[dict[str, str]]) -> dict[str, list]:
    columns = {"effect": [float(r["effect"]) for r in rows], "P": [float(r["P"]) for r in rows], "q": [float(r["q"]) for r in rows]}
    for g in GROUPS:
        columns[f"n_{g}"] = [float(r[f"n_{g}"]) for r in rows]
    return columns


# ---------------------------------------------------------------- V199: column tables, statistics table, Info sheet

def _column_source() -> tuple[list[str], dict[str, dict[str, list[float]]]]:
    rows = _tsv("column_source.tsv")
    proteins = list(dict.fromkeys(r["protein"] for r in rows))
    groups = list(dict.fromkeys(r["group"] for r in rows))
    values = {p: {g: [float(r["value"]) for r in rows if r["protein"] == p and r["group"] == g] for g in groups} for p in proteins}
    return proteins, values


def test_v199_column_tables_hold_the_source_values_exactly_and_the_statistics_table_and_info_sheet_follow(tmp_path):
    proteins, source = _column_source()
    stats = _tsv("stats_proteins.tsv")
    tables = [column_table(p, {g: source[p][g] for g in GROUPS}) for p in proteins]
    tables.append(statistics_table("fixture figure - pipeline statistics", _stats_comparisons(stats, GROUPS), GROUPS))
    path = write_pzfx(tmp_path / "column_groups.pzfx", tables, _info("table", source_sha=_sha("stats_proteins.tsv")))
    parsed = _oracle_tables(path)
    assert [t["type"] for t in parsed] == [PRISM_CODE["Column"]] * (len(proteins) + 1)
    for protein, table in zip(proteins, parsed):
        assert table["title"] == protein and table["order"] == GROUPS
        assert _mismatches({g: source[protein][g] for g in GROUPS}, table["columns"]) == []
        assert len(table["columns"]["Control"]) == 6                      # six replicate rows per group
    # the statistics table: comparison row titles; effect, P, q and n per group equal the source statistics exactly
    stats_table = parsed[-1]
    assert stats_table["rows"] == [r["comparison"] for r in stats]
    assert _mismatches(_expected_stats(stats), stats_table["columns"]) == []
    # the Info sheet names every required constant
    info = _oracle_info(path)
    assert set(INFO_CONSTANTS) <= set(info["constants"])
    assert info["constants"]["Statistics"] == "table"
    assert info["constants"]["Source table SHA-256"] == _sha("stats_proteins.tsv")
    assert len(info["notes"]) >= 2
    assert path.read_bytes() == _fixture_bytes("column_groups.pzfx"), "the committed fixture is the writer output"


def test_v199_round_trip_precision_is_full_float(tmp_path):
    values = [0.1 + 0.2, 1e-300, 123456789.123456789, -2.5e-7]
    path = write_pzfx(tmp_path / "precision.pzfx", [column_table("precision", {"G": values})], _info("none"))
    assert _oracle_tables(path)[0]["columns"]["G"] == values


def test_v199_missing_values_are_empty_cells(tmp_path):
    path = write_pzfx(tmp_path / "missing.pzfx", [column_table("missing", {"A": [1.0, float("nan"), 3.0], "B": [None, 2.0, 2.5]})], _info("none"))
    parsed = _oracle_tables(path)[0]["columns"]
    assert parsed == {"A": [1.0, None, 3.0], "B": [None, 2.0, 2.5]}
    assert _mismatches({"A": [1.0, None, 3.0], "B": [None, 2.0, 2.5]}, parsed) == []


def test_v199_negative_single_changed_value_fails_the_comparison(tmp_path):
    proteins, source = _column_source()
    path = write_pzfx(tmp_path / "column.pzfx", [column_table(proteins[0], {g: source[proteins[0]][g] for g in GROUPS})], _info("none"))
    parsed = _oracle_tables(path)[0]["columns"]
    assert _mismatches(source[proteins[0]], parsed) == []
    changed = {k: list(v) for k, v in parsed.items()}
    changed["Acute"][2] = changed["Acute"][2] + 1e-9
    assert _mismatches(source[proteins[0]], changed) != []


def test_v199_negative_rounding_to_three_decimals_fails_the_comparison(tmp_path):
    proteins, source = _column_source()
    rounded = {g: [round(v, 3) for v in source[proteins[0]][g]] for g in GROUPS}
    path = write_pzfx(tmp_path / "rounded.pzfx", [column_table(proteins[0], rounded)], _info("none"))
    assert _mismatches(source[proteins[0]], _oracle_tables(path)[0]["columns"]) != []


def test_v199_negative_a_statistics_value_mismatch_fails_the_comparison(tmp_path):
    stats = _tsv("stats_proteins.tsv")
    expected = _expected_stats(stats)
    tables = [statistics_table("stats", _stats_comparisons(stats, GROUPS), GROUPS)]
    path = write_pzfx(tmp_path / "stats.pzfx", tables, _info("table"))
    parsed = _oracle_tables(path)[0]["columns"]
    assert _mismatches(expected, parsed) == []
    tampered = {k: list(v) for k, v in parsed.items()}
    tampered["P"][1] = tampered["P"][1] * (1 + 1e-12)
    assert _mismatches(expected, tampered) != []


def test_v199_negative_missing_info_constant_is_refused_and_no_file_is_written(tmp_path):
    info = _info("none")
    del info["constants"]["Source table SHA-256"]
    with pytest.raises(ConfigurationError) as caught:
        pzfx_bytes([column_table("x", {"A": [1.0]})], info)
    assert caught.value.code == "E_PZFX_INFO_MISSING"
    assert "Source table SHA-256" in str(caught.value)
    target = tmp_path / "no_info.pzfx"
    with pytest.raises(ConfigurationError):
        write_pzfx(target, [column_table("x", {"A": [1.0]})], info)
    assert not target.exists()


def test_v199_statistics_none_is_recorded_in_the_info_sheet(tmp_path):
    path = write_pzfx(tmp_path / "none.pzfx", [column_table("x", {"A": [1.0, 2.0]})], _info("none"))
    assert _oracle_info(path)["constants"]["Statistics"] == "none"
    assert len(_oracle_tables(path)) == 1                              # no statistics table is added


# ---------------------------------------------------------------- V200: XY tables

def _volcano() -> tuple[dict, dict]:
    rows = _tsv("volcano_source.tsv")
    roc = _tsv("roc_source.tsv")
    return ({"log2FC": [float(r["log2FC"]) for r in rows], "neg_log10P": [float(r["neg_log10P"]) for r in rows]},
            {"one_minus_specificity": [float(r["one_minus_specificity"]) for r in roc], "sensitivity": [float(r["sensitivity"]) for r in roc]})


def test_v200_xy_tables_hold_every_source_row_exactly(tmp_path):
    volcano, roc = _volcano()
    assert len(volcano["log2FC"]) == 200 and len(roc["sensitivity"]) == 30
    path = write_pzfx(tmp_path / "xy_volcano.pzfx", [xy_table("volcano", "log2FC", "neg_log10P", volcano["log2FC"], volcano["neg_log10P"])], _info("none"))
    parsed = _oracle_tables(path)
    assert parsed[0]["type"] == PRISM_CODE["XY"] and parsed[0]["order"] == ["log2FC", "neg_log10P"]
    assert _mismatches(volcano, parsed[0]["columns"]) == []
    assert len(parsed[0]["columns"]["log2FC"]) == 200                    # one point per source row
    roc_path = write_pzfx(tmp_path / "xy_roc.pzfx", [xy_table("roc", "one_minus_specificity", "sensitivity", roc["one_minus_specificity"], roc["sensitivity"])], _info("none"))
    roc_parsed = _oracle_tables(roc_path)
    assert roc_parsed[0]["type"] == PRISM_CODE["XY"]
    assert _mismatches(roc, roc_parsed[0]["columns"]) == []


def test_v200_volcano_statistics_table_has_one_row_per_feature(tmp_path):
    rows = _tsv("stats_features.tsv")
    assert len(rows) == 200
    table = statistics_table("volcano - pipeline statistics", _stats_comparisons(rows, GROUPS), GROUPS)
    parsed = _oracle_tables(write_pzfx(tmp_path / "v.pzfx", [table], _info("table")))[0]
    assert parsed["rows"] == [r["comparison"] for r in rows]
    assert _mismatches(_expected_stats(rows), parsed["columns"]) == []


def test_v200_negative_values_rounded_to_three_decimals_fail(tmp_path):
    volcano, _ = _volcano()
    rounded = {k: [round(v, 3) for v in vals] for k, vals in volcano.items()}
    path = write_pzfx(tmp_path / "rounded_xy.pzfx", [xy_table("volcano", "log2FC", "neg_log10P", rounded["log2FC"], rounded["neg_log10P"])], _info("none"))
    assert _mismatches(volcano, _oracle_tables(path)[0]["columns"]) != []


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
    path = write_pzfx(tmp_path / "grouped.pzfx", [table], _info("none", test="none (multi-protein replicate panel)", source_sha=_sha("grouped_source.tsv")))
    parsed = _oracle_tables(path)[0]
    assert parsed["type"] == PRISM_CODE["Grouped"]
    assert parsed["rows"] == [f"rep{r}" for r in replicates]
    assert len(parsed["order"]) == len(proteins) * len(groups)           # one data set per protein and group
    assert _mismatches(datasets, parsed["columns"]) == []
    assert path.read_bytes() == _fixture_bytes("grouped_panel.pzfx"), "the committed fixture is the writer output"


def test_v201_negative_missing_protein_data_set_fails_the_comparison(tmp_path):
    replicates, proteins, groups, datasets = _grouped_source()
    without_p2 = {name: values for name, values in datasets.items() if not name.startswith("P2 ")}
    path = write_pzfx(tmp_path / "missing_p2.pzfx", [grouped_table("panel", [f"rep{r}" for r in replicates], without_p2)], _info("none"))
    problems = _mismatches(datasets, _oracle_tables(path)[0]["columns"])
    assert problems and all("P2" in p for p in problems)


def test_v201_writer_refuses_a_row_count_that_does_not_match():
    with pytest.raises(ValueError, match="1 values for 2 rows"):
        grouped_table("panel", ["rep1", "rep2"], {"P1 Control": [1.0]})


def test_v199_to_v201_writer_is_deterministic():
    _, source = _column_source()
    first = pzfx_bytes([column_table("P1", {"Control": source["P1"]["Control"]})], _info("none"))
    second = pzfx_bytes([column_table("P1", {"Control": source["P1"]["Control"]})], _info("none"))
    assert first == second and b"\r" not in first and first.startswith(b'<?xml version="1.0" encoding="UTF-8"?>\n<GraphPadPrismFile')
    assert not any(math.isnan(v) for v in source["P1"]["Control"])
