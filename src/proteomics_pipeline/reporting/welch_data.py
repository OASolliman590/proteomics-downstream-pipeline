"""Report data of the Welch engine stage (packet R15a, ADR 0011).

Reads only the stage's own outputs (welch_results.tsv, eligibility.json, families.tsv and, when the
paper-compatible rule is declared, paper_compatible.tsv). Every value shown downstream keeps the
frozen claim label: a Welch result is labelled ``welch`` and a raw-P selection is labelled
``not_multiplicity_controlled`` with BH q alongside (SM43).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def welch_report_data(stage_dir: str | Path) -> dict | None:
    stage = Path(stage_dir)
    results = stage / "welch_results.tsv"
    if not results.is_file():
        return None
    rows = _rows(results)
    tested = [r for r in rows if r["eligibility"] == "tested"]
    eligibility = json.loads((stage / "eligibility.json").read_text(encoding="utf-8"))
    data = {
        "engine": "welch",
        "claim_labels": list(eligibility.get("claim_labels", [])),
        "endpoint_rows": len(rows),
        "tested_rows": len(tested),
        "families": _rows(stage / "families.tsv") if (stage / "families.tsv").is_file() else [],
        "adaptations": eligibility.get("adaptations", []),
        "decisions": eligibility.get("decisions", []),
    }
    paper = stage / "paper_compatible.tsv"
    if paper.is_file():
        selected = _rows(paper)
        data["paper_compatible"] = {"selected": len(selected), "claim_label": "not_multiplicity_controlled",
                                    "bh_q_column": "bh_q_value", "raw_p_cutoff": eligibility.get("paper_compatible", {}).get("raw_p_cutoff")}
    return data
