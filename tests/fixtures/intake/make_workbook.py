"""Write a tiny synthetic four-sheet XLSX in the historical workbook layout (V017).

Values are invented: abundance(feature_index, workbook_label, replicate) =
10 + feature_index + label_offset + replicate/10, with one cell (accession
A3, sample 'PDM 2') left empty in every sheet that contains it.
"""
from __future__ import annotations

import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

SHEETS = [("DM_vs_CTL", "CTL", "DM"), ("PDM_vs_CTL", "CTL", "PDM"), ("DM-Treated_vs_DM", "DM", "DM-Treated"), ("PDM-Treated_vs_PDM", "PDM", "PDM-Treated")]
LABEL_OFFSET = {"CTL": 0.0, "DM": 1.0, "PDM": 0.5, "DM-Treated": 0.25, "PDM-Treated": 0.125}
ACCESSIONS = [("A1", "Gn1"), ("A2", "Gn2"), ("A3", "")]
MISSING = ("A3", "PDM 2")


def abundance(feature_index: int, label: str, replicate: int) -> float:
    return 10 + feature_index + LABEL_OFFSET[label] + replicate / 10


def _col(index: int) -> str:
    name = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        name = chr(65 + rem) + name
    return name


def write_workbook(path: str | Path) -> dict:
    strings: list[str] = []

    def s(text: str) -> int:
        if text not in strings:
            strings.append(text)
        return strings.index(text)

    sheets_xml = []
    expected: dict[tuple[str, str], float | None] = {}
    for sheet, first, second in SHEETS:
        samples = [f"{first} {r}" for r in range(1, 5)] + [f"{second} {r}" for r in range(1, 5)]
        header = ["Accession", "Gene", "Protein"] + samples + ["log2fc", "pvalue", "p.adj", "nlog10p"]
        rows = [[("s", h) for h in header]]
        for fi, (accession, gene) in enumerate(ACCESSIONS):
            row = [("s", accession), ("s", gene), ("s", f"synthetic protein {accession}")]
            for sample in samples:
                label, rep = sample.rsplit(" ", 1)
                value = None if (accession, sample) == MISSING else abundance(fi, label, int(rep))
                row.append(("n", value))
                expected[(accession, sample)] = value
            row += [("n", LABEL_OFFSET[second] - LABEL_OFFSET[first]), ("n", 0.5), ("n", 0.9), ("n", 0.30102999566398120)]
            rows.append(row)
        cells = []
        for ri, row in enumerate(rows, start=1):
            parts = []
            for ci, (kind, value) in enumerate(row):
                ref = f"{_col(ci)}{ri}"
                if kind == "s":
                    parts.append(f'<c r="{ref}" t="s"><v>{s(value)}</v></c>')
                elif value is None:
                    parts.append(f'<c r="{ref}"/>')
                else:
                    parts.append(f'<c r="{ref}"><v>{value!r}</v></c>')
            cells.append(f'<row r="{ri}">{"".join(parts)}</row>')
        sheets_xml.append((sheet, '<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + "".join(cells) + "</sheetData></worksheet>"))
    shared = '<?xml version="1.0" encoding="UTF-8"?><sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' + "".join(f"<si><t>{escape(t)}</t></si>" for t in strings) + "</sst>"
    workbook = ('<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
                + "".join(f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>' for i, (name, _) in enumerate(sheets_xml, start=1)) + "</sheets></workbook>")
    rels = ('<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + "".join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1, len(sheets_xml) + 1))
            + "</Relationships>")
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", rels)
        archive.writestr("xl/sharedStrings.xml", shared)
        for i, (_, xml) in enumerate(sheets_xml, start=1):
            archive.writestr(f"xl/worksheets/sheet{i}.xml", xml)
    return expected
