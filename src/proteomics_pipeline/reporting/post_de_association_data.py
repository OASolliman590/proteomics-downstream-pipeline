"""Post-DE protein-phenotype association section of ReportData (packet R14c, FR-148, SM25).

Built only from verified post_de/association outputs: one row per analysis with its family, complete-case n,
adjustment set and method; refused analyses (for example an aliased pooled analysis) and the Simpson warning are shown.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _table(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def post_de_association_section(directory: str | Path, *, prefix: str = "post_de/association") -> dict:
    root = Path(directory)
    eligibility = json.loads((root / "eligibility.json").read_text(encoding="utf-8"))
    summary = _table(root / "association_summary.tsv")
    section = {"state": "COMPLETED", "title": "Post-DE protein-phenotype association", "claim_label": None,
               "summary": ["Proteins are the response and each phenotype is a design term; " + eligibility.get("inference", ""),
                           "Phenotype values are never imputed: units with a missing phenotype are excluded for that phenotype (complete case)."],
               "tables": [{"caption": "Association analyses (one secondary family per phenotype; not part of the group-contrast families)",
                           "columns": ["analysis", "phenotype", "scope", "group", "family", "method", "n units", "missing phenotype", "adjusted for", "tested", "q < 0.05", "Simpson flags"],
                           "rows": [[s["analysis_id"], s["phenotype"], s["scope"], s["group"] if s["group"] != "NA" else "", s["family_id"], s["method"], s["n_units"], s["n_missing_phenotype"],
                                     s["adjusted_for"], s["n_tested"], s["n_q_below_0_05"], s["n_simpson_flag"]] for s in summary],
                           "source": f"{prefix}/association_summary.tsv", "numeric": ["n units", "missing phenotype", "tested", "q < 0.05", "Simpson flags"]}],
               "figures": [], "notes": [], "methods": [eligibility.get("rule", "")], "values": {"n_analyses": len(summary)}}
    if eligibility.get("simpson_warning"):
        section["notes"].append(eligibility["simpson_warning"])
    refusals = _table(root / "refusals.tsv")
    if refusals:
        section["tables"].append({"caption": "Requested association analyses that were refused", "columns": ["analysis", "item", "reason code", "reason"],
                                  "rows": [[r["analysis"], r["item"], r["reason_code"], r["reason"]] for r in refusals], "source": f"{prefix}/refusals.tsv"})
    for record in eligibility.get("figures", []):
        if record.get("state") == "COMPLETED" and record["file"].endswith(".png"):
            section["figures"].append({"src": f"{prefix}/{record['file']}", "alt": "association heatmap",
                                       "caption": "Feature x analysis association statistics, rows clustered; every cell is in the source table.", "source": f"{prefix}/heatmap_source.tsv"})
    return section
