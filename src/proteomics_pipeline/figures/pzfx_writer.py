"""GraphPad Prism .pzfx writer (packet R16b, FR-199 to FR-201, V199 to V201).

The XML is written here with the standard library only (xml.etree.ElementTree); the R package pzfx is not used to write.
Three table kinds follow the contract (figures.md, "Prism .pzfx table rules"):

- ``Column``: one column per group, holding the individual replicate values (FR-199);
- ``XY``: one point per source row, X and Y exactly as in the source (FR-200);
- ``Grouped``: one row per replicate and one data set per protein and group (FR-201).

Numbers are written with ``repr`` (the shortest string that reads back as the same float, D-12). A missing value is an
empty cell. Before a file is returned, the writer parses the bytes it produced and compares every value with the source
(E_PZFX_VALUE_MISMATCH). Prism's own table codes are used in the XML: Column = OneWay, XY = XY, Grouped = TwoWay.
"""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union

from ..errors import ConfigurationError

PRISM_TABLE_TYPE = {"Column": "OneWay", "XY": "XY", "Grouped": "TwoWay"}
XML_DECLARATION = b'<?xml version="1.0" encoding="UTF-8"?>\n'
Cell = Union[float, None]


@dataclass(frozen=True)
class ColumnTable:
    title: str
    columns: tuple[tuple[str, tuple[Cell, ...]], ...]


@dataclass(frozen=True)
class XYTable:
    title: str
    x_title: str
    y_title: str
    x: tuple[Cell, ...]
    y: tuple[Cell, ...]


@dataclass(frozen=True)
class GroupedTable:
    title: str
    row_titles: tuple[str, ...]
    datasets: tuple[tuple[str, tuple[Cell, ...]], ...]


PzfxTable = Union[ColumnTable, XYTable, GroupedTable]


def _cell(value: Any, where: str) -> Cell:
    """A source value as written: a finite float, or None for a missing value (None or NaN)."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError("E_PZFX_VALUE_MISMATCH", f"{where}: value {value!r} is not numeric", where)
    number = float(value)
    if math.isnan(number):
        return None
    if math.isinf(number):
        raise ConfigurationError("E_PZFX_VALUE_MISMATCH", f"{where}: infinite values cannot be stored in a Prism table", where)
    return number


def _cells(values: Sequence[Any], where: str) -> tuple[Cell, ...]:
    return tuple(_cell(v, f"{where}[{i}]") for i, v in enumerate(values))


def column_table(title: str, groups: Mapping[str, Sequence[Any]]) -> ColumnTable:
    """Group comparison: one column per group, in the declared group order."""
    if not groups:
        raise ValueError("a column table needs at least one group")
    return ColumnTable(title, tuple((str(g), _cells(v, f"{title}/{g}")) for g, v in groups.items()))


def xy_table(title: str, x_title: str, y_title: str, x: Sequence[Any], y: Sequence[Any]) -> XYTable:
    """XY table: one point per source row."""
    if len(x) != len(y):
        raise ValueError("an XY table needs one Y value per X value")
    return XYTable(title, x_title, y_title, _cells(x, f"{title}/x"), _cells(y, f"{title}/y"))


def grouped_table(title: str, row_titles: Sequence[str], datasets: Mapping[str, Sequence[Any]]) -> GroupedTable:
    """Grouped table: rows are replicates, one data set per protein and group (every data set has one value per row)."""
    rows = tuple(str(r) for r in row_titles)
    if not rows or any(not r for r in rows) or len(set(rows)) != len(rows):
        raise ValueError("grouped rows need unique, non-empty titles")
    for name, values in datasets.items():
        if len(values) != len(rows):
            raise ValueError(f"data set {name!r} has {len(values)} values for {len(rows)} rows")
    return GroupedTable(title, rows, tuple((str(n), _cells(v, f"{title}/{n}")) for n, v in datasets.items()))


def _expected(table: PzfxTable) -> dict[str, Any]:
    """The table as the source defines it; the self-check compares the written XML with this."""
    if isinstance(table, ColumnTable):
        return {"kind": "Column", "title": table.title, "rows": None, "columns": [(name, list(values)) for name, values in table.columns]}
    if isinstance(table, XYTable):
        return {"kind": "XY", "title": table.title, "rows": None,
                "columns": [(table.x_title, list(table.x)), (table.y_title, list(table.y))]}
    if isinstance(table, GroupedTable):
        return {"kind": "Grouped", "title": table.title, "rows": list(table.row_titles),
                "columns": [(name, list(values)) for name, values in table.datasets]}
    raise TypeError(f"not a Prism table: {type(table).__name__}")


def _write_cells(parent: ET.Element, values: Sequence[Cell]) -> None:
    for value in values:
        cell = ET.SubElement(parent, "d")
        if value is not None:
            cell.text = repr(float(value))


def _build(tables: Sequence[PzfxTable]) -> ET.Element:
    root = ET.Element("GraphPadPrismFile", {"PrismXMLVersion": "5.00"})
    created = ET.SubElement(root, "Created")
    ET.SubElement(created, "OriginalVersion", {"CreatedByProgram": "proteomics_pipeline", "CreatedByVersion": "R16b"})
    sequence = ET.SubElement(root, "TableSequence")
    for index in range(len(tables)):
        ET.SubElement(sequence, "Ref", {"ID": f"Table{index}", **({"Selected": "1"} if index == 0 else {})})
    for index, table in enumerate(tables):
        kind = _expected(table)["kind"]
        is_xy = kind == "XY"
        grouped = kind == "Grouped"
        attributes = {"ID": f"Table{index}", "XFormat": "numbers" if is_xy else "none", "YFormat": "replicates",
                      "Replicates": "1", "TableType": PRISM_TABLE_TYPE[kind], "EVFormat": "AsteriskAfterNumber"}
        node = ET.SubElement(root, "Table", attributes)
        ET.SubElement(node, "Title").text = table.title
        if grouped:
            row_column = ET.SubElement(node, "RowTitlesColumn")
            _write_cells_text(row_column, table.row_titles)
        if isinstance(table, XYTable):
            x_column = ET.SubElement(node, "XColumn", {"Subcolumns": "1"})
            ET.SubElement(x_column, "Title").text = table.x_title
            _write_cells(ET.SubElement(x_column, "Subcolumn"), table.x)
            y_column = ET.SubElement(node, "YColumn", {"Subcolumns": "1"})
            ET.SubElement(y_column, "Title").text = table.y_title
            _write_cells(ET.SubElement(y_column, "Subcolumn"), table.y)
        else:
            columns = table.columns if isinstance(table, ColumnTable) else table.datasets
            for name, values in columns:
                y_column = ET.SubElement(node, "YColumn", {"Subcolumns": "1"})
                ET.SubElement(y_column, "Title").text = name
                _write_cells(ET.SubElement(y_column, "Subcolumn"), values)
    return root


def _write_cells_text(parent: ET.Element, titles: Sequence[str]) -> None:
    subcolumn = ET.SubElement(parent, "Subcolumn")
    for title in titles:
        ET.SubElement(subcolumn, "d").text = title


def _read_cell(node: ET.Element) -> Cell:
    text = node.text
    return None if text is None or text == "" else float(text)


def _parsed(data: bytes) -> list[dict[str, Any]]:
    """Read back the bytes the writer produced (the self-check reads the XML, not the Python objects)."""
    root = ET.fromstring(data)
    parsed = []
    for node in root.findall("Table"):
        kind = {v: k for k, v in PRISM_TABLE_TYPE.items()}[node.get("TableType")]
        rows = None
        if node.find("RowTitlesColumn") is not None:
            rows = [d.text or "" for d in node.find("RowTitlesColumn").find("Subcolumn").findall("d")]
        columns = []
        for child in node:
            if child.tag in ("XColumn", "YColumn"):
                values = [_read_cell(d) for d in child.find("Subcolumn").findall("d")]
                columns.append((child.findtext("Title"), values))
        parsed.append({"kind": kind, "title": node.findtext("Title"), "rows": rows, "columns": columns})
    return parsed


def _self_check(data: bytes, expected: Sequence[dict[str, Any]]) -> None:
    parsed = _parsed(data)
    if len(parsed) != len(expected):
        raise ConfigurationError("E_PZFX_VALUE_MISMATCH", f"the file holds {len(parsed)} tables, the source defines {len(expected)}")
    for index, (found, wanted) in enumerate(zip(parsed, expected)):
        if found["kind"] != wanted["kind"] or found["title"] != wanted["title"]:
            raise ConfigurationError("E_PZFX_VALUE_MISMATCH", f"table {index} kind or title differs from the source", f"/Table{index}")
        if found["rows"] != wanted["rows"]:
            raise ConfigurationError("E_PZFX_VALUE_MISMATCH", f"table {index} row titles differ from the source", f"/Table{index}")
        if [name for name, _ in found["columns"]] != [name for name, _ in wanted["columns"]]:
            raise ConfigurationError("E_PZFX_VALUE_MISMATCH", f"table {index} data set titles differ from the source", f"/Table{index}")
        for (name, got), (_, want) in zip(found["columns"], wanted["columns"]):
            if got != [None if v is None or (isinstance(v, float) and math.isnan(v)) else float(v) for v in want]:
                raise ConfigurationError("E_PZFX_VALUE_MISMATCH", f"table {index}, data set {name!r}: a value differs from the source", f"/Table{index}")


def pzfx_bytes(tables: Sequence[PzfxTable]) -> bytes:
    """The complete .pzfx document for the tables, after the self-check. Deterministic for the same input."""
    tables = list(tables)
    if not tables:
        raise ValueError("a .pzfx file needs at least one table")
    expected = [_expected(t) for t in tables]
    root = _build(tables)
    ET.indent(root, space="  ")
    data = XML_DECLARATION + ET.tostring(root, encoding="utf-8")
    _self_check(data, expected)
    return data


def write_pzfx(path: str | Path, tables: Sequence[PzfxTable]) -> Path:
    """Write the .pzfx file (LF line endings; the bytes are written as produced)."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(pzfx_bytes(tables))
    return destination
