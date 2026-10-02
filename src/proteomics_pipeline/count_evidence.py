"""Count-evidence table checks for DEqMS (packet R06, SM10, V051).

The table is genuine peptide/PSM evidence keyed by feature_id (optionally
observation_id/plex_id), with recorded source, aggregation and pseudocount
policy.  Observed-sample counts or abundance-derived proxies are rejected;
missing counts are never treated as zero.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from .errors import ProteomicsError

REQUIRED = ("feature_id", "count_type", "count_value", "count_source", "count_aggregation", "pseudocount_policy")
PROXIES = {"observed_sample_count", "abundance_derived", "intensity_derived", "imputed"}


def _schema() -> dict:
    return json.loads((Path(__file__).parent / "schemas" / "count-evidence.schema.json").read_text(encoding="utf-8"))


def check_count_table(path: str | Path, *, model_id: str, aggregation: str) -> list[dict]:
    from jsonschema import Draft202012Validator
    validator = Draft202012Validator(_schema())
    pointer = f"/models[{model_id}]/count_evidence"
    if not Path(path).is_file():
        raise ProteomicsError("E_DEQMS_COUNT_EVIDENCE", f"count evidence file {Path(path).name!r} does not exist", pointer)
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(c not in rows[0] for c in REQUIRED):
        raise ProteomicsError("E_DEQMS_COUNT_EVIDENCE", f"count evidence must have columns {list(REQUIRED)}", pointer)
    for number, row in enumerate(rows, start=2):
        errors = list(validator.iter_errors(row))
        if errors:
            raise ProteomicsError("E_DEQMS_COUNT_EVIDENCE", f"row {number}: {errors[0].message}", pointer)
        if row["count_source"] in PROXIES:
            raise ProteomicsError("E_DEQMS_COUNT_EVIDENCE", f"row {number}: {row['count_source']} is a proxy, not peptide/PSM count evidence", pointer)
        if row["count_aggregation"] != aggregation:
            raise ProteomicsError("E_DEQMS_COUNT_EVIDENCE", f"row {number}: recorded aggregation {row['count_aggregation']} differs from the model's {aggregation}", pointer)
    return rows
