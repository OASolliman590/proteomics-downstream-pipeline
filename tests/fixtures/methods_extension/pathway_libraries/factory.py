"""Synthetic GMT libraries and Welch-shaped result tables for R15b (V174-V178).

Every GMT row carries the marker SYNTHETIC_FIXTURE in its description field; the repository may contain only
such marked fixtures (FR-175). Gene names are synthetic identifiers (G0001 ...). Nothing here is a real pathway
library, and no KEGG or MSigDB content is reproduced.
"""
from __future__ import annotations

import csv
from pathlib import Path

MARKER = "SYNTHETIC_FIXTURE"


def write_gmt(path: Path, sets: dict[str, list[str]], *, note: str = "curated test set", marked: bool = True) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for name, genes in sets.items():
        description = f"{MARKER}: {note}" if marked else note
        lines.append("\t".join([name, description, *genes]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def write_welch_table(path: Path, rows: list[dict]) -> Path:
    """A Welch-shaped result table: only the columns the R15b stage reads are required; the rest mirror R15a."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["feature_id", "gene_id", "model_id", "contrast_id", "eligibility", "effect", "statistic", "p_value", "q_value", "engine"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            base = {"model_id": "welch-main", "eligibility": "tested", "engine": "welch", "p_value": "NA", "q_value": "NA", "statistic": "NA"}
            writer.writerow({**base, **{k: ("NA" if v is None else v) for k, v in row.items()}})
    return path


def ora_universe_table(path: Path, *, genes: list[str], foreground: set[str], contrast: str = "B-A") -> Path:
    """One Welch row per measured gene; foreground genes are significant (q 0.01) and up; the rest are not significant."""
    rows = []
    for gene in genes:
        hit = gene in foreground
        rows.append({"feature_id": f"P_{gene}", "gene_id": gene, "contrast_id": contrast, "effect": 1.0 if hit else -0.1,
                     "statistic": 2.0 if hit else -0.2, "q_value": 0.01 if hit else 0.9, "p_value": 0.001 if hit else 0.5})
    return write_welch_table(path, rows)


def gene_names(n: int, prefix: str = "G") -> list[str]:
    return [f"{prefix}{i:04d}" for i in range(1, n + 1)]
