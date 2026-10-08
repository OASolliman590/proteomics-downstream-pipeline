"""Synthetic fixtures for the .pzfx writer (R16b, V199-V201, D-75). Run from the repository root:

    .venv/bin/python tests/fixtures/figures/pzfx/make_fixtures.py

Values come from fixed trigonometric formulas with planted group and protein offsets; nothing is study data. Each
source table is written with repr() so that the round trip is tested at full float precision. The .pzfx files are the
writer's output for these tables, with the statistics table and the Info sheet of D-75, and are checked byte for byte by
tests/scientific/test_pzfx_writer.py.
"""
from __future__ import annotations

import hashlib
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))

from proteomics_pipeline.figures.pzfx_writer import (  # noqa: E402
    column_table, grouped_table, info_sheet, statistics_table, write_pzfx, xy_table,
)

PROTEINS = ["P1", "P2", "P3"]                 # V199 column tables: three proteins
PANEL_PROTEINS = ["P1", "P2", "P3", "P4"]  # V201 grouped panel: four proteins
GROUPS = ["Control", "Acute"]
REPLICATES = 6
CREATED = "2026-10-08T00:00:00Z"


def _tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(v if isinstance(v, str) else repr(float(v)) if isinstance(v, float) else str(v) for v in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def column_value(k: int, g: int, r: int) -> float:
    return 20.0 + 2.5 * k + 1.5 * g + 0.731 * math.sin(1.7 * r + 0.9 * k + 0.4 * g)


def grouped_value(k: int, g: int, r: int) -> float:
    return 18.0 + 1.1 * k + 0.9 * g + 0.6 * math.cos(1.3 * r + k)


def info_for(source: Path, statistics: str, test: str) -> dict:
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    return info_sheet(
        "fixture figure",
        {"Test": test, "Multiplicity adjustment": "BH within fixture-family (synthetic)", "Effect scale": "log2 difference",
         "Group sizes": f"Control {REPLICATES}; Acute {REPLICATES}", "Pipeline version and commit": "0.1.0.dev0 (fixture)",
         "Source table SHA-256": digest, "Run ID": "fixture-run", "Created (UTC)": CREATED, "Statistics": statistics},
        notes=["Synthetic fixture for the pzfx writer.", "A test run inside Prism may give a different P value from the pipeline's test."])


def _comparison(row: list) -> dict:
    return {"row": row[0], "effect": row[1], "P": row[2], "q": row[3], "n": {"Control": row[4], "Acute": row[5]}}


def main() -> None:
    column_rows = [[PROTEINS[k], GROUPS[g], r + 1, column_value(k, g, r)] for k in range(len(PROTEINS)) for g in range(len(GROUPS)) for r in range(REPLICATES)]
    _tsv(HERE / "column_source.tsv", ["protein", "group", "replicate", "value"], column_rows)

    # statistics of the three protein comparisons (effect, P, q and group sizes); synthetic values
    protein_stats = []
    for k, protein in enumerate(PROTEINS):
        p_value = 10 ** -(1.0 + 1.5 * k)
        protein_stats.append([f"{protein} - Control", 1.5 + 0.2 * k, p_value, min(1.0, p_value * 3), REPLICATES, REPLICATES])
    stats_source = HERE / "stats_proteins.tsv"
    _tsv(stats_source, ["comparison", "effect", "P", "q", "n_Control", "n_Acute"], protein_stats)
    column_tables = [column_table(protein, {GROUPS[g]: [column_value(k, g, r) for r in range(REPLICATES)] for g in range(len(GROUPS))})
                     for k, protein in enumerate(PROTEINS)]
    column_tables.append(statistics_table("fixture figure - pipeline statistics", [_comparison(r) for r in protein_stats], GROUPS))
    write_pzfx(HERE / "column_groups.pzfx", column_tables, info_for(stats_source, "table", "Welch t (fixture)"))

    volcano = [[f"F{i:03d}", 2.0 * math.sin(0.61 * i) + 0.05 * math.cos(i), -math.log10(10 ** -(1.0 + 3.0 * abs(math.sin(0.23 * i))))]
               for i in range(1, 201)]
    _tsv(HERE / "volcano_source.tsv", ["feature", "log2FC", "neg_log10P"], volcano)
    feature_stats = [[row[0], row[1], 10 ** -row[2], min(1.0, 10 ** -row[2] * 40), REPLICATES, REPLICATES] for row in volcano]
    feature_source = HERE / "stats_features.tsv"
    _tsv(feature_source, ["comparison", "effect", "P", "q", "n_Control", "n_Acute"], feature_stats)
    volcano_tables = [
        xy_table("volcano", "log2FC", "neg_log10P", [row[1] for row in volcano], [row[2] for row in volcano]),
        statistics_table("volcano - pipeline statistics", [_comparison(r) for r in feature_stats], GROUPS),
    ]
    write_pzfx(HERE / "xy_volcano.pzfx", volcano_tables, info_for(feature_source, "table", "Welch t (fixture)"))

    roc = [[i / 29, min(1.0, math.sqrt(i / 29) + 0.03 * math.sin(i))] for i in range(30)]
    roc_source = HERE / "roc_source.tsv"
    _tsv(roc_source, ["one_minus_specificity", "sensitivity"], roc)
    write_pzfx(HERE / "xy_roc.pzfx", [xy_table("roc", "one_minus_specificity", "sensitivity", [row[0] for row in roc], [row[1] for row in roc])],
               info_for(roc_source, "none", "none (ROC curve, no inferential statistics)"))

    grouped_rows = [[r + 1, PANEL_PROTEINS[k], GROUPS[g], grouped_value(k, g, r)] for r in range(REPLICATES) for k in range(len(PANEL_PROTEINS)) for g in range(len(GROUPS))]
    grouped_source = HERE / "grouped_source.tsv"
    _tsv(grouped_source, ["replicate", "protein", "group", "value"], [[r, p, g, v] for r, p, g, v in grouped_rows])
    datasets = {f"{protein} {group}": [grouped_value(k, g, r) for r in range(REPLICATES)]
                for k, protein in enumerate(PANEL_PROTEINS) for g, group in enumerate(GROUPS)}
    write_pzfx(HERE / "grouped_panel.pzfx", [grouped_table("panel", [f"rep{r + 1}" for r in range(REPLICATES)], datasets)],
               info_for(grouped_source, "none", "none (multi-protein replicate panel)"))
    (HERE / "xy_volcano_roc.pzfx").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
