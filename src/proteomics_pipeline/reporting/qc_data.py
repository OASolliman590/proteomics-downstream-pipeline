"""QC section of ReportData built from verified R03 preprocessing outputs.

Every number exposed here is read from a published QC table, and each
entry carries the relative path of that source table (SM06, V030/V094).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def qc_section(stage_dir: str | Path, *, prefix: str = "preprocessing") -> dict:
    root = Path(stage_dir)
    result = json.loads((root / "preprocessing_result.json").read_text(encoding="utf-8"))
    section = {"state": "COMPLETED", "summary": {k: result[k] for k in ("scope", "assay", "scale", "normalization", "primary_values_changed", "original_mask_state", "n_features", "n_observations")},
               "excluded_observations": result.get("excluded_observations", []), "reference_only_channels": result.get("reference_only_channels", []),
               "pca": result.get("pca"), "detection": result.get("detection"), "coverage": result.get("coverage"), "tables": {}}
    for key, relative in (("sample_n", "qc/sample_n.tsv"), ("missingness_by_group", "missingness/by_group.tsv"), ("normalization_factors", "normalization_factors.tsv"),
                          ("exclusions", "exclusions/exclusions.tsv"), ("watchlist", "qc/watchlist.tsv"), ("pca_variance", "qc/pca_variance.tsv"),
                          ("distributions", "qc/distributions.tsv")):
        path = root / relative
        if path.is_file():
            section["tables"][key] = {"source": f"{prefix}/{relative}", "rows": _table(path)}
        else:
            section["tables"][key] = {"source": None, "rows": None, "state": "NOT_AVAILABLE"}
    statement = root / "missingness" / "statement.json"
    section["missingness_statement"] = json.loads(statement.read_text(encoding="utf-8")) if statement.is_file() else None
    return section
