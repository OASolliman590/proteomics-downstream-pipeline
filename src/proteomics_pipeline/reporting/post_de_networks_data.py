"""Post-DE co-abundance and network section of ReportData (packet R14e, FR-161-FR-165, SM25).

Built only from verified post_de/networks outputs.  Module-trait results are labelled module-level; hub and degree
results are descriptive; the snapshot release, species, score type and hash and the measured-universe null are shown.
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


def post_de_networks_section(directory: str | Path, *, prefix: str = "post_de/networks") -> dict:
    root = Path(directory)
    eligibility = json.loads((root / "eligibility.json").read_text(encoding="utf-8"))
    summary = eligibility.get("summary") or {}
    section = {"state": "COMPLETED", "title": "Post-DE co-abundance modules and interaction networks", "claim_label": eligibility.get("claim_label"),
               "summary": [], "tables": [], "figures": [], "notes": list(eligibility.get("limitations", [])), "methods": [eligibility.get("rule", "")], "values": {}}
    co = summary.get("coabundance")
    if co:
        section["summary"].append(f"Co-abundance ({co['rule']}, soft threshold {co['soft_threshold']}): {co['n_modules']} module(s) from {co['n_features']} complete features of {co['n_units']} units.")
        stability = _table(root / "module_stability.tsv")
        if stability:
            section["tables"].append({"caption": "Modules and unit-bootstrap stability (descriptive)", "columns": ["module", "features", "variance explained", "bootstrap", "seed", "mean Jaccard"],
                                      "rows": [[s["module"], s["n_features"], _f(s["variance_explained"]), s["bootstrap"], s["seed"], _f(s["mean_jaccard"])] for s in stability],
                                      "source": f"{prefix}/module_stability.tsv"})
        trait = _table(root / "module_trait.tsv")
        if trait:
            section["tables"].append({"caption": "Module-trait association on eigengenes [module_level]; not inference on individual proteins",
                                      "columns": ["trait", "contrast", "module", "effect", "t", "P", "q", "family", "claim"],
                                      "rows": [[t["trait"], t["contrast_id"], t["module"], _f(t["effect"]), _f(t["statistic"]), _f(t["p_value"]), _f(t["q_value"]), t["family_id"], t["claim_label"]] for t in trait],
                                      "source": f"{prefix}/module_trait.tsv"})
    ppi = summary.get("ppi")
    if ppi:
        section["summary"].append(f"Interaction snapshot {ppi['snapshot_id']} ({ppi.get('source')}, release {ppi['release']}, taxonomy {ppi['species']}, score {ppi['score_type']}, SHA-256 {ppi['snapshot_sha256'][:12]}...): "
                                  f"measured mapped universe of {ppi['universe_size']} genes, minimum score {ppi['min_score']}.")
        conn = _table(root / "connectivity.tsv")
        if conn:
            section["tables"].append({"caption": "Connectivity of declared sets against degree-preserving random sets from the measured mapped universe",
                                      "columns": ["set", "set genes", "universe", "observed edges", "null mean", "draws", "seed", "P"],
                                      "rows": [[c["set_id"], c["n_set_genes_in_universe"], c["universe_size"], c["observed_edges"], _f(c["null_mean_edges"]), c["null_draws"], c["seed"], _f(c["p_value"])] for c in conn],
                                      "source": f"{prefix}/connectivity.tsv"})
    refusals = _table(root / "refusals.tsv")
    if refusals:
        section["tables"].append({"caption": "Refused network analyses", "columns": ["analysis", "item", "reason code", "reason"],
                                  "rows": [[r["analysis"], r["item"], r["reason_code"], r["reason"]] for r in refusals], "source": f"{prefix}/refusals.tsv"})
    for record in eligibility.get("figures", []):
        if record.get("state") == "COMPLETED" and record["file"].endswith(".png"):
            stem = Path(record["file"]).stem
            section["figures"].append({"src": f"{prefix}/{record['file']}", "alt": stem.replace("_", " "), "caption": stem.replace("_", " ").capitalize() + ".",
                                       "source": f"{prefix}/figure_sources/{stem}.tsv"})
    return section
