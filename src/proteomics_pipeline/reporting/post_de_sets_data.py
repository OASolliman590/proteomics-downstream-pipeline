"""Post-DE result-structure section of ReportData (packet R14a, FR-131-FR-136, SM25).

Built only from the verified post_de/sets outputs.  Every table names its source file and claim label;
refused sub-analyses (for example a Venn for k > 3) are listed, never shown as empty results.
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


def _f(value, digits=3):
    try:
        return f"{float(value):.{digits}g}"
    except (TypeError, ValueError):
        return value


def post_de_sets_section(directory: str | Path, *, prefix: str = "post_de/sets") -> dict:
    root = Path(directory)
    eligibility = json.loads((root / "eligibility.json").read_text(encoding="utf-8"))
    sets = _table(root / "sets.tsv")
    regions = _table(root / "regions.tsv")
    section = {"state": "COMPLETED", "title": "Post-DE result structure: declared sets, regions and concordance",
               "claim_label": eligibility.get("claim_label"),
               "summary": [f"{len(sets)} declared set(s) built from completed families with the finite set grammar; region counts are exact and sum to the size of the union.",
                           "Concordance and direction-aware overlap are descriptive; a hypergeometric overlap P is reported only for contrasts that share no biological units."],
               "tables": [{"caption": "Declared sets (rule, criterion, families, null type, claim label)",
                           "columns": ["set", "rule", "criteria", "threshold", "families", "null type", "matrix", "members", "unique genes", "claim"],
                           "rows": [[s["set_id"], s["rule"], s["criteria"], s["thresholds"], s["family_ids"], s["null_type"], s["input_matrix"], s["n_members"], s["n_unique_genes"], s["claim_label"]] for s in sets],
                           "source": f"{prefix}/sets.tsv", "numeric": ["members", "unique genes"]},
                          {"caption": "Exact exclusive regions (protein groups and unique genes are counted separately)",
                           "columns": ["region", "in sets", "protein groups", "unique genes", "without gene", "claim"],
                           "rows": [[r["region_id"], r["defining_sets"], r["n_protein_groups"], r["n_unique_genes"], r["n_without_gene"], r["claim_label"]] for r in regions],
                           "source": f"{prefix}/regions.tsv", "numeric": ["protein groups", "unique genes", "without gene"]}],
               "figures": [], "notes": [], "methods": [], "values": {"n_sets": len(sets), "n_regions": len(regions)}}
    direction = _table(root / "direction_overlap.tsv")
    if direction:
        section["tables"].append({"caption": "Direction-aware overlap of single-endpoint sets (descriptive)",
                                  "columns": ["set A", "set B", "concordant up", "concordant down", "discordant", "only A", "only B"],
                                  "rows": [[d["set_a"], d["set_b"], d["n_concordant_up"], d["n_concordant_down"], d["n_discordant"], d["n_only_a"], d["n_only_b"]] for d in direction],
                                  "source": f"{prefix}/direction_overlap.tsv"})
    concordance = _table(root / "concordance.tsv")
    if concordance:
        section["tables"].append({"caption": "Effect concordance on common tested features (descriptive; no P values)",
                                  "columns": ["pair", "left", "right", "common tested", "Pearson", "Spearman", "origin slope", "sign agreement", "claim"],
                                  "rows": [[c["pair_id"], c["left"], c["right"], c["n_common_tested"], _f(c["pearson"]), _f(c["spearman"]), _f(c["origin_slope"]), _f(c["sign_agreement"]), c["claim_label"]] for c in concordance],
                                  "source": f"{prefix}/concordance.tsv"})
    eligibility_rows = _table(root / "overlap_eligibility.tsv")
    if eligibility_rows:
        section["tables"].append({"caption": "Overlap inference eligibility (shared biological units make overlap descriptive only)",
                                  "columns": ["set A", "set B", "shared units", "overlap P eligible", "reason"],
                                  "rows": [[e["set_a"], e["set_b"], e["n_shared_units"], e["overlap_p_eligible"], e["reason"]] for e in eligibility_rows],
                                  "source": f"{prefix}/overlap_eligibility.tsv"})
    tests = _table(root / "overlap_test.tsv")
    if tests:
        section["tables"].append({"caption": "Hypergeometric overlap (unit-disjoint contrasts only; own descriptive-overlap family, BH)",
                                  "columns": ["set A", "set B", "universe", "n A", "n B", "overlap", "expected", "P", "q"],
                                  "rows": [[t["set_a"], t["set_b"], t["universe_n"], t["n_a"], t["n_b"], t["overlap"], _f(t["expected"]), _f(t["p_value"]), _f(t["q_value"])] for t in tests],
                                  "source": f"{prefix}/overlap_test.tsv"})
    refusals = _table(root / "refusals.tsv")
    if refusals:
        section["tables"].append({"caption": "Requested sub-analyses that were refused", "columns": ["analysis", "item", "reason code", "reason"],
                                  "rows": [[r["analysis"], r["item"], r["reason_code"], r["reason"]] for r in refusals], "source": f"{prefix}/refusals.tsv"})
    for record in eligibility.get("figures", []):
        if record.get("state") != "COMPLETED" or not record["file"].endswith(".png"):
            continue
        stem = Path(record["file"]).stem
        section["figures"].append({"src": f"{prefix}/{record['file']}", "alt": stem.replace("_", " "), "caption": stem.replace("_", " ").capitalize() + " (descriptive).",
                                   "source": f"{prefix}/figure_sources/{stem}.tsv"})
    section["methods"] = [eligibility.get("rule", "")]
    return section
