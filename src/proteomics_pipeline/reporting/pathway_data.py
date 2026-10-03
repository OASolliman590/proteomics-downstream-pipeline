"""Pathway section of ReportData (packet R10b, FR-096).

Each method keeps its own null type, family and universe; fgsea and ORA are
labelled exploratory; nothing claims cross-null confirmation, sample-level
replication from gene-set randomization, or a causal mechanism.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

NULL_LABELS = {"competitive_enrichment": "CAMERA competitive null (estimated inter-gene correlation)",
               "self_contained_directional": "ROAST rotation, self-contained directional null",
               "self_contained_mixed": "ROAST rotation, self-contained mixed null",
               "preranked_gene_set": "fgsea preranked gene-set randomization (exploratory; not sample-level inference)",
               "ora_up": "ORA hypergeometric over the measured background (exploratory, up)",
               "ora_down": "ORA hypergeometric over the measured background (exploratory, down)"}


def _table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _f(v):
    try:
        return f"{float(v):.3g}"
    except (TypeError, ValueError):
        return v


def pathway_section(directory: str | Path, *, prefix: str = "pathways") -> dict:
    root = Path(directory)
    result = json.loads((root / "pathway_result.json").read_text(encoding="utf-8"))
    tables = []
    for htype, label in NULL_LABELS.items():
        path = root / f"{htype}.tsv"
        if not path.is_file():
            continue
        rows = _table(path)
        cols = ["contrast_id", "set_id", "p_value", "q_value", "native_q_value", "family_id"]
        extra = {"competitive_enrichment": ["correlation", "direction"], "self_contained_directional": ["direction", "nrot"], "self_contained_mixed": ["nrot"],
                 "preranked_gene_set": ["nes", "size", "leading_edge"], "ora_up": ["N", "K", "n", "k"], "ora_down": ["N", "K", "n", "k"]}[htype]
        tables.append({"caption": f"{label} — all eligible sets, including non-significant", "source": f"{prefix}/{htype}.tsv", "columns": cols + extra,
                       "rows": [[_f(r.get(c)) for c in cols + extra] for r in rows]})
    if (root / "leading_edge_overlap.tsv").is_file():
        ov = _table(root / "leading_edge_overlap.tsv")
        tables.append({"caption": "Leading-edge overlap (redundancy display only; tests unchanged)", "source": f"{prefix}/leading_edge_overlap.tsv",
                       "columns": ["set_a", "set_b", "intersection", "jaccard"], "rows": [[r["set_a"], r["set_b"], r["intersection"], _f(r["jaccard"])] for r in ov]})
    return {"state": "COMPLETED", "title": "Pathways and gene-set enrichment", "tables": tables, "figures": [],
            "summary": [result["dispatch"], f"Gene-level model: {result['gene_model']['model_id']} ({result['gene_model']['role']})."],
            "notes": ["Different null types answer different questions and are never pooled or used as confirmation of each other.",
                      "A pathway name is not evidence of a biological mechanism; overlapping sets are not independent discoveries."],
            "methods": [f"Methods: {', '.join(result['methods'])}; limma {result['versions']['limma']}" + (f", fgsea {result['versions']['fgsea']}" if result["versions"].get("fgsea") else "") + "."],
            "values": {"methods": result["methods"]}}
