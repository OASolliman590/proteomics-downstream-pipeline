"""Strict TSV parsing for canonical wide/long protein-level matrices (FR-011, FR-020).

Values are parsed exactly once; every rejection carries the file, 1-based
row number and field name so a user can find the defect.  Nothing here
transforms scale or aggregates observations.
"""
from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from ..errors import ProteomicsError

# Column names that identify a sub-protein grain.  A table containing any of
# them is not a protein-level table and is rejected (SM01, V016/V020).
SUBPROTEIN_COLUMNS = frozenset(
    name.lower()
    for name in (
        "peptide", "peptide_sequence", "sequence", "stripped_sequence", "stripped.sequence",
        "modified_sequence", "modified.sequence", "modified sequence", "precursor", "precursor_id",
        "precursor.id", "eg.precursorid", "pep.strippedsequence", "psm", "psm_id", "spectrum",
        "site", "phosphosite", "site_position", "modification_site", "charge",
    )
)


class InputError(ProteomicsError):
    """Typed input rejection with a precise file/row/field location."""

    def __init__(self, code: str, message: str, *, file: str | None = None, row: int | None = None,
                 field: str | None = None, exit_code: int = 2):
        location = []
        if file is not None:
            location.append(f"file={Path(file).name}")
        if row is not None:
            location.append(f"row={row}")
        if field is not None:
            location.append(f"field={field}")
        full = message if not location else f"{message} ({', '.join(location)})"
        super().__init__(code, full, exit_code=exit_code)
        self.file, self.row, self.field = (str(file) if file is not None else None), row, field

    def as_dict(self) -> dict:
        value = super().as_dict()
        for key in ("file", "row", "field"):
            item = getattr(self, key)
            if item is not None:
                value[key] = Path(item).name if key == "file" else item
        return value


def read_text_utf8(path: str | Path) -> str:
    """Decode a file as UTF-8 (BOM tolerated); any other encoding is rejected."""
    path = Path(path)
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise InputError("E_INPUT_READ", str(exc), file=str(path)) from exc
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        line = raw[: exc.start].count(b"\n") + 1
        raise InputError("E_INPUT_ENCODING", f"file is not valid UTF-8 at byte {exc.start}", file=str(path), row=line) from exc


def read_tsv(path: str | Path) -> tuple[list[str], list[list[str]]]:
    """Read a TAB-delimited table with standard double-quote escaping.

    Returns the header and the data rows.  Row numbers in errors count the
    header as row 1.  Ragged rows, empty or duplicated header names are
    rejected instead of padded or renamed.
    """
    text = read_text_utf8(path)
    if not text.strip():
        raise InputError("E_INPUT_EMPTY", "table is empty", file=str(path))
    reader = csv.reader(io.StringIO(text, newline=""), delimiter="\t", quotechar='"', doublequote=True, strict=True)
    try:
        rows = list(reader)
    except csv.Error as exc:
        raise InputError("E_INPUT_QUOTING", f"malformed quoting: {exc}", file=str(path), row=reader.line_num) from exc
    while rows and rows[-1] == []:
        rows.pop()
    header, body = rows[0], rows[1:]
    for index, name in enumerate(header):
        if name == "":
            raise InputError("E_INPUT_HEADER", "header contains an empty column name", file=str(path), row=1, field=f"column {index + 1}")
    seen: set[str] = set()
    for name in header:
        if name in seen:
            raise InputError("E_ID_DUPLICATE", f"duplicate header {name!r}", file=str(path), row=1, field=name)
        seen.add(name)
    for number, row in enumerate(body, start=2):
        if row == []:
            raise InputError("E_INPUT_RAGGED", "blank line inside table", file=str(path), row=number)
        if len(row) != len(header):
            raise InputError("E_INPUT_RAGGED", f"row has {len(row)} fields but header has {len(header)}", file=str(path), row=number)
    return header, body


def parse_number(token: str, missing: Iterable[str], *, file: str, row: int, field: str) -> float | None:
    """Parse one numeric cell; declared missing tokens become None.

    Non-finite values (Inf, NaN) are rejected rather than treated as missing
    or clipped (V020 negative case).
    """
    if token in missing:
        return None
    try:
        value = float(token)
    except ValueError:
        raise InputError("E_INPUT_NONNUMERIC", f"value {token!r} is not numeric and is not a declared missing token", file=file, row=row, field=field) from None
    if not math.isfinite(value):
        raise InputError("E_INPUT_NONFINITE", f"non-finite value {token!r} is not allowed", file=file, row=row, field=field)
    return value


def reject_subprotein_grain(header: Iterable[str], file: str) -> None:
    for name in header:
        if name.strip().lower() in SUBPROTEIN_COLUMNS:
            raise InputError("E_UNSUPPORTED_SCOPE", f"column {name!r} indicates a peptide/precursor/site-level table; only protein-level input is supported", file=file, row=1, field=name)


@dataclass
class RawMatrix:
    """Parsed protein x observation values on the source scale."""

    feature_ids: list[str]
    observation_ids: list[str]
    values: dict[tuple[str, str], float | None] = field(default_factory=dict)
    source_rows: dict[str, int] = field(default_factory=dict)

    def get(self, feature_id: str, observation_id: str) -> float | None:
        return self.values[(feature_id, observation_id)]


def parse_wide(path: str | Path, missing: Iterable[str], id_column: str = "feature_id") -> RawMatrix:
    path = str(path)
    missing = set(missing)
    header, body = read_tsv(path)
    reject_subprotein_grain(header, path)
    if header[0] != id_column:
        raise InputError("E_INPUT_HEADER", f"first column must be {id_column!r}", file=path, row=1, field=header[0])
    observation_ids = header[1:]
    if not observation_ids:
        raise InputError("E_INPUT_HEADER", "wide matrix has no observation columns", file=path, row=1)
    matrix = RawMatrix(feature_ids=[], observation_ids=list(observation_ids))
    for number, row in enumerate(body, start=2):
        feature_id = row[0]
        if feature_id == "":
            raise InputError("E_ID_MISSING", "feature_id is empty", file=path, row=number, field=id_column)
        if feature_id in matrix.source_rows:
            raise InputError("E_ID_DUPLICATE", f"duplicate feature_id {feature_id!r} (first at row {matrix.source_rows[feature_id]})", file=path, row=number, field=id_column)
        matrix.source_rows[feature_id] = number
        matrix.feature_ids.append(feature_id)
        for observation_id, token in zip(observation_ids, row[1:]):
            matrix.values[(feature_id, observation_id)] = parse_number(token, missing, file=path, row=number, field=observation_id)
    if not matrix.feature_ids:
        raise InputError("E_INPUT_EMPTY", "wide matrix has no feature rows", file=path)
    return matrix


LONG_COLUMNS = ("feature_id", "observation_id", "abundance")


def parse_long(path: str | Path, missing: Iterable[str]) -> RawMatrix:
    """Parse a long table whose feature/observation pairs must be unique and complete."""
    path = str(path)
    missing = set(missing)
    header, body = read_tsv(path)
    reject_subprotein_grain(header, path)
    if tuple(header) != LONG_COLUMNS:
        raise InputError("E_INPUT_HEADER", f"long table header must be exactly {list(LONG_COLUMNS)}", file=path, row=1)
    matrix = RawMatrix(feature_ids=[], observation_ids=[])
    seen_features: set[str] = set()
    seen_observations: set[str] = set()
    pair_rows: dict[tuple[str, str], int] = {}
    for number, (feature_id, observation_id, token) in enumerate(body, start=2):
        if feature_id == "" or observation_id == "":
            raise InputError("E_ID_MISSING", "feature_id and observation_id are required", file=path, row=number)
        key = (feature_id, observation_id)
        if key in pair_rows:
            raise InputError("E_ID_DUPLICATE", f"duplicate feature-observation pair {key!r} (first at row {pair_rows[key]})", file=path, row=number, field="feature_id,observation_id")
        pair_rows[key] = number
        if feature_id not in seen_features:
            seen_features.add(feature_id)
            matrix.feature_ids.append(feature_id)
            matrix.source_rows[feature_id] = number
        if observation_id not in seen_observations:
            seen_observations.add(observation_id)
            matrix.observation_ids.append(observation_id)
        matrix.values[key] = parse_number(token, missing, file=path, row=number, field="abundance")
    for feature_id in matrix.feature_ids:
        for observation_id in matrix.observation_ids:
            if (feature_id, observation_id) not in matrix.values:
                raise InputError("E_ID_MISSING", f"missing feature-observation pair ({feature_id!r}, {observation_id!r}); absent pairs must be explicit missing tokens", file=path, field="feature_id,observation_id")
    if not matrix.feature_ids:
        raise InputError("E_INPUT_EMPTY", "long table has no rows", file=path)
    return matrix


def parse_bool_mask(path: str | Path, feature_ids: list[str], observation_ids: list[str], name: str) -> dict[tuple[str, str], bool]:
    """Read a wide true/false mask whose keys must equal the matrix keys exactly."""
    path = str(path)
    header, body = read_tsv(path)
    if header[0] != "feature_id":
        raise InputError("E_INPUT_HEADER", f"{name} first column must be 'feature_id'", file=path, row=1, field=header[0])
    if set(header[1:]) != set(observation_ids) or len(header) - 1 != len(observation_ids):
        raise InputError("E_ID_ALIGNMENT", f"{name} observation columns do not match the matrix observations", file=path, row=1)
    rows: dict[str, list[str]] = {}
    for number, row in enumerate(body, start=2):
        if row[0] in rows:
            raise InputError("E_ID_DUPLICATE", f"duplicate feature_id {row[0]!r} in {name}", file=path, row=number, field="feature_id")
        rows[row[0]] = row
    if set(rows) != set(feature_ids):
        raise InputError("E_ID_ALIGNMENT", f"{name} feature IDs do not match the matrix features", file=path)
    mask: dict[tuple[str, str], bool] = {}
    for feature_id in feature_ids:
        row = rows[feature_id]
        number = list(rows).index(feature_id) + 2
        for observation_id, token in zip(header[1:], row[1:]):
            if token not in ("true", "false"):
                raise InputError("E_MASK_VALUE", f"{name} values must be 'true' or 'false', got {token!r}", file=path, row=number, field=observation_id)
            mask[(feature_id, observation_id)] = token == "true"
    return mask


def format_number(value: float | None) -> str:
    """Serialize a double exactly (shortest round-trip repr); missing is NA."""
    if value is None:
        return "NA"
    if not math.isfinite(value):
        raise ValueError("non-finite values cannot be serialized")
    return "0.0" if value == 0 else repr(float(value))


def write_tsv(path: str | Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    """Write UTF-8, TAB-delimited, LF-terminated TSV with minimal quoting."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(["NA" if value is None else value for value in row])
