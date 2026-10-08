"""Synthetic fixtures for the R16b .pzfx writer (V199-V201). Run from the repository root:

    .venv/bin/python tests/fixtures/figures/pzfx/make_fixtures.py

Values come from fixed trigonometric formulas with planted group and protein offsets; nothing is study data. Each
source table is written with repr() so that the round trip is tested at full float precision. The .pzfx files are the
writer's output for these tables and are checked byte for byte by tests/scientific/test_pzfx_writer.py.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))

from proteomics_pipeline.figures.pzfx_writer import column_table, grouped_table, write_pzfx, xy_table  # noqa: E402

PROTEINS = ["P1", "P2", "P3"]                 # V199 column tables: three proteins
PANEL_PROTEINS = ["P1", "P2", "P3", "P4"]  # V201 grouped panel: four proteins
GROUPS = ["Control", "Acute"]
REPLICATES = 6


def _tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(v if isinstance(v, str) else repr(float(v)) if isinstance(v, float) else str(v) for v in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def column_value(k: int, g: int, r: int) -> float:
    return 20.0 + 2.5 * k + 1.5 * g + 0.731 * math.sin(1.7 * r + 0.9 * k + 0.4 * g)


def grouped_value(k: int, g: int, r: int) -> float:
    return 18.0 + 1.1 * k + 0.9 * g + 0.6 * math.cos(1.3 * r + k)


def main() -> None:
    column_rows = [[PROTEINS[k], GROUPS[g], r + 1, column_value(k, g, r)] for k in range(len(PROTEINS)) for g in range(len(GROUPS)) for r in range(REPLICATES)]
    _tsv(HERE / "column_source.tsv", ["protein", "group", "replicate", "value"], column_rows)

    volcano = [[f"F{i:03d}", 2.0 * math.sin(0.61 * i) + 0.05 * math.cos(i), -math.log10(10 ** -(1.0 + 3.0 * abs(math.sin(0.23 * i)))) ]
               for i in range(1, 201)]
    _tsv(HERE / "volcano_source.tsv", ["feature", "log2FC", "neg_log10P"], volcano)

    roc = [[i / 29, min(1.0, math.sqrt(i / 29) + 0.03 * math.sin(i))] for i in range(30)]
    _tsv(HERE / "roc_source.tsv", ["one_minus_specificity", "sensitivity"], roc)

    grouped_rows = [[r + 1, PANEL_PROTEINS[k], GROUPS[g], grouped_value(k, g, r)] for r in range(REPLICATES) for k in range(len(PANEL_PROTEINS)) for g in range(len(GROUPS))]
    _tsv(HERE / "grouped_source.tsv", ["replicate", "protein", "group", "value"], [[r, p, g, v] for r, p, g, v in grouped_rows])

    column_tables = []
    for k, protein in enumerate(PROTEINS):
        column_tables.append(column_table(protein, {GROUPS[g]: [column_value(k, g, r) for r in range(REPLICATES)] for g in range(len(GROUPS))}))
    write_pzfx(HERE / "column_groups.pzfx", column_tables)

    xy = [xy_table("volcano", "log2FC", "neg_log10P", [row[1] for row in volcano], [row[2] for row in volcano]),
          xy_table("roc", "one_minus_specificity", "sensitivity", [row[0] for row in roc], [row[1] for row in roc])]
    write_pzfx(HERE / "xy_volcano_roc.pzfx", xy)

    datasets = {f"{protein} {group}": [grouped_value(k, g, r) for r in range(REPLICATES)]
                for k, protein in enumerate(PANEL_PROTEINS) for g, group in enumerate(GROUPS)}
    write_pzfx(HERE / "grouped_panel.pzfx", [grouped_table("panel", [f"rep{r + 1}" for r in range(REPLICATES)], datasets)])


if __name__ == "__main__":
    main()
