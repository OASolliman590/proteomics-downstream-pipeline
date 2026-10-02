"""Shared builders for R02 tests.  Oracles live in the test functions, not here."""
from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXAMPLES = ROOT / "configs" / "examples"
FIXTURES = ROOT / "tests" / "fixtures" / "intake"
MAPPINGS = ROOT / "configs" / "mappings"


def example_config(name: str = "example-independent.json") -> dict:
    config = json.loads((EXAMPLES / name).read_text(encoding="utf-8"))
    for key in ("matrix", "observations", "features", "source_provenance"):
        config["input"][key] = str((EXAMPLES / config["input"][key]).resolve())
    return config


def write_config(tmp_path: Path, config: dict) -> Path:
    path = tmp_path / "analysis.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    return path


def write_tsv(path: Path, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow(row)
    return path


def read_tsv(path: Path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance(tmp_path: Path, **fields) -> str:
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "provenance.json"
    path.write_text(json.dumps({"kind": "synthetic_test_fixture", **fields}), encoding="utf-8")
    return str(path)


def minimal_config(tmp_path: Path, *, matrix: Path, observations: Path, features: Path | str, source_scale="log2", transform="none",
                   zero_as_missing=False, technical=None, fmt="wide_tsv", profile="canonical_protein", mapping=None,
                   missing=("NA", ""), group_levels=None, prior_imputation="none_documented", original_mask=None,
                   prior_mask=None, provenance_fields=None) -> dict:
    config = copy.deepcopy(example_config())
    config["source_scale"] = source_scale
    config["preprocessing"]["transform"] = transform
    config["preprocessing"]["technical_replicates"] = technical or {"mode": "none", "biological_unit_column": "biological_unit_id"}
    config["input"].update({"matrix": str(matrix), "observations": str(observations), "features": str(features), "format": fmt,
                            "profile": profile, "missing_values": list(missing), "linear_zero_as_missing": zero_as_missing,
                            "prior_imputation": prior_imputation, "original_observed_mask": None if original_mask is None else str(original_mask),
                            "prior_imputed_mask": None if prior_mask is None else str(prior_mask),
                            "source_provenance": provenance(tmp_path, **(provenance_fields or {}))})
    if mapping is not None:
        config["input"]["mapping"] = str(mapping)
    if group_levels is not None:
        config["design"]["group_levels"] = list(group_levels)
    return config


def features_rows(ids):
    rows = [["feature_id", "accessions", "gene_ids", "gene_symbols", "is_decoy", "is_contaminant", "protein_group_ambiguous"]]
    rows += [[i, json.dumps([i]), "[]", "[]", "unknown", "unknown", "false"] for i in ids]
    return rows


def observation_rows(entries):
    """entries: (observation_id, unit, subject, technical, group)."""
    return [["observation_id", "biological_unit_id", "subject_id", "technical_replicate_id", "group"]] + [list(e) for e in entries]
