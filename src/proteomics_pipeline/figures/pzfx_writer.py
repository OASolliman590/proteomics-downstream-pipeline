"""GraphPad Prism .pzfx writer (packets R16b and R16e; FR-199 to FR-201, D-75).

The XML is written here with the standard library only (xml.etree.ElementTree); the R package pzfx is not used to write.
Table kinds (figures.md, "Prism .pzfx table rules"):

- ``Column``: one column per group holding the replicate values (FR-199), or the statistics table of D-75 (rows are the
  comparisons, columns effect, P, q and n_<group>);
- ``XY``: one point per source row, X and Y exactly as in the source (FR-200);
- ``Grouped``: one row per replicate and one data set per protein and group (FR-201).

Every file also carries an Info sheet (D-75): the Prism ``InfoSequence`` and ``Info`` elements with ``Constant`` name/value
pairs and ``Notes``. A missing required constant is refused (E_PZFX_INFO_MISSING).

Numbers are written with ``repr`` (the shortest string that reads back as the same float, D-12). A missing value is an empty
cell. Before a file is returned, the writer parses the bytes it produced and compares every value, title, row title and
Info constant with the source (E_PZFX_VALUE_MISMATCH). Prism's table codes: Column = OneWay, XY = XY, Grouped = TwoWay.
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
# Required Info constants (figures.md, "Statistics table and Info sheet (D-75)"). "Statistics" is "none" when the figure has
# no inferential statistics.
INFO_CONSTANTS = ("Test", "Multiplicity adjustment", "Effect scale", "Group sizes", "Pipeline version and commit",
                  "Source table SHA-256", "Run ID", "Created (UTC)", "Statistics")
STATISTICS_COLUMNS = ("effect", "P", "q")
Cell = Union[float, None]


@dataclass(frozen=True)
class ColumnTable:
    title: str
    columns: tuple[tuple[str, tuple[Cell, ...]], ...]
    row_titles: tuple[str, ...] | None = None


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


def statistics_table(title: str, comparisons: Sequence[Mapping[str, Any]], groups: Sequence[str]) -> ColumnTable:
    """The D-75 statistics table. Each comparison is a mapping with ``row`` (its title), ``effect``, ``P``, ``q`` and
    ``n`` (a mapping from group to size). Columns are effect, P, q and n_<group> for each group, in declared order."""
    if not comparisons:
        raise ValueError("a statistics table needs at least one comparison")
    rows = tuple(str(c["row"]) for c in comparisons)
    if len(set(rows)) != len(rows):
        raise ValueError("statistics comparisons must have unique row titles")
    columns: list[tuple[str, tuple[Cell, ...]]] = []
    for name in STATISTICS_COLUMNS:
        columns.append((name, _cells([c[name] for c in comparisons], f"{title}/{name}")))
    for group in groups:
        columns.append((f"n_{group}", _cells([c["n"][group] for c in comparisons], f"{title}/n_{group}")))
    return ColumnTable(title, tuple(columns), row_titles=rows)


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


def info_sheet(title: str, constants: Mapping[str, Any], notes: Sequence[str] = ()) -> dict[str, Any]:
    """The Info sheet of D-75. ``constants`` must hold every name in INFO_CONSTANTS; ``notes`` repeat the caption."""
    return {"title": str(title), "constants": {str(k): str(v) for k, v in constants.items()}, "notes": [str(n) for n in notes]}


def _expected(table: PzfxTable) -> dict[str, Any]:
    """The table as the source defines it; the self-check compares the written XML with this."""
    if isinstance(table, ColumnTable):
        return {"kind": "Column", "title": table.title, "rows": list(table.row_titles) if table.row_titles else None,
                "columns": [(name, list(values)) for name, values in table.columns]}
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


def _write_titles(parent: ET.Element, titles: Sequence[str]) -> None:
    subcolumn = ET.SubElement(parent, "Subcolumn")
    for title in titles:
        ET.SubElement(subcolumn, "d").text = title


def _info_elements(info: Mapping[str, Any]) -> list[ET.Element]:
    refs = ET.Element("InfoSequence")
    ET.SubElement(refs, "Ref", {"ID": "Info0", "Selected": "1"})
    node = ET.Element("Info", {"ID": "Info0"})
    ET.SubElement(node, "Title").text = info["title"]
    notes = ET.SubElement(node, "Notes")
    for line in info["notes"]:
        font = ET.SubElement(notes, "Font", {"Color": "#000000", "Face": "Helvetica"})
        font.text = line
    for name, value in info["constants"].items():
        constant = ET.SubElement(node, "Constant")
        ET.SubElement(constant, "Name").text = name
        ET.SubElement(constant, "Value").text = value
    return [refs, node]


def _build(tables: Sequence[PzfxTable], info: Mapping[str, Any]) -> ET.Element:
    root = ET.Element("GraphPadPrismFile", {"PrismXMLVersion": "5.00"})
    created = ET.SubElement(root, "Created")
    ET.SubElement(created, "OriginalVersion", {"CreatedByProgram": "proteomics_pipeline", "CreatedByVersion": "R16e"})
    for element in _info_elements(info):
        root.append(element)
    sequence = ET.SubElement(root, "TableSequence")
    for index in range(len(tables)):
        ET.SubElement(sequence, "Ref", {"ID": f"Table{index}", **({"Selected": "1"} if index == 0 else {})})
    for index, table in enumerate(tables):
        kind = _expected(table)["kind"]
        attributes = {"ID": f"Table{index}", "XFormat": "numbers" if kind == "XY" else "none", "YFormat": "replicates",
                      "Replicates": "1", "TableType": PRISM_TABLE_TYPE[kind], "EVFormat": "AsteriskAfterNumber"}
        node = ET.SubElement(root, "Table", attributes)
        ET.SubElement(node, "Title").text = table.title
        has_rows = isinstance(table, GroupedTable) or (isinstance(table, ColumnTable) and table.row_titles)
        if has_rows:
            _write_titles(ET.SubElement(node, "RowTitlesColumn"), table.row_titles)
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


def _read_cell(node: ET.Element) -> Cell:
    text = node.text
    return None if text is None or text == "" else float(text)


def _parsed(data: bytes) -> tuple[list[dict[str, Any]], dict[str, Any]]:
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
    info_node = root.find("Info")
    info_read: dict[str, Any] = {"title": None, "constants": {}, "notes": []}
    if info_node is not None:
        info_read["title"] = info_node.findtext("Title")
        info_read["notes"] = [font.text or "" for font in info_node.find("Notes").findall("Font")]
        for constant in info_node.findall("Constant"):
            info_read["constants"][constant.findtext("Name")] = constant.findtext("Value")
    return parsed, info_read


def _self_check(data: bytes, expected: Sequence[dict[str, Any]], info: Mapping[str, Any]) -> None:
    parsed, info_read = _parsed(data)
    if len(parsed) != len(expected):
        raise ConfigurationError("E_PZFX_VALUE_MISMATCH", f"the file holds {len(parsed)} tables, the source defines {len(expected)}")
    if info_read != {"title": info["title"], "constants": info["constants"], "notes": info["notes"]}:
        raise ConfigurationError("E_PZFX_VALUE_MISMATCH", "the Info sheet differs from the source", "/Info0")
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


def _check_info(info: Mapping[str, Any]) -> None:
    if not isinstance(info, Mapping) or not all(key in info for key in ("title", "constants", "notes")):
        raise ConfigurationError("E_PZFX_INFO_MISSING", "the Info sheet needs a title, constants and notes", "/Info0")
    missing = [name for name in INFO_CONSTANTS if name not in info["constants"]]
    if missing:
        raise ConfigurationError("E_PZFX_INFO_MISSING", f"the Info sheet lacks the constant(s) {missing}", "/Info0")


def pzfx_bytes(tables: Sequence[PzfxTable], info: Mapping[str, Any]) -> bytes:
    """The complete .pzfx document for the tables and the Info sheet, after the self-check. Deterministic for the same input."""
    tables = list(tables)
    if not tables:
        raise ValueError("a .pzfx file needs at least one table")
    _check_info(info)
    expected = [_expected(t) for t in tables]
    root = _build(tables, info)
    ET.indent(root, space="  ")
    data = XML_DECLARATION + ET.tostring(root, encoding="utf-8")
    _self_check(data, expected, info)
    return data


def write_pzfx(path: str | Path, tables: Sequence[PzfxTable], info: Mapping[str, Any]) -> Path:
    """Write the .pzfx file (LF line endings; the bytes are written as produced)."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(pzfx_bytes(tables, info))
    return destination
