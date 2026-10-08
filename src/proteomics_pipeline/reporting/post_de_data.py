"""Integrated post-DE eligibility and dependency section of ReportData (packet R14f, FR-166, SM41, SM25).

Every requested post-DE module is listed with its plan eligibility, actual state, typed reason and claim label.
NOT_RUN and INAPPLICABLE are rendered as different states with their meanings; a refused module never appears as a
zero result.  The dependency summary links each primary discovery to its post-DE outcomes.
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


def post_de_section(directory: str | Path, *, prefix: str = "post_de/eligibility") -> dict:
    root = Path(directory)
    data = json.loads((root / "eligibility.json").read_text(encoding="utf-8"))
    modules = data["modules"]
    deps = _table(root / "dependency_summary.tsv")
    rows = []
    for m in modules:
        reasons = "; ".join(f"{r['reason_code']}: {r['reason']}" for r in m.get("plan_reasons") or []) or (m.get("reason") or "")
        if m["state"] != "COMPLETED":
            claim = "none (no values)"
        else:   # a completed module without a post-DE claim label (association) shows its declared family-level inference instead
            claim = m.get("claim_label") or (f"no post-DE label; {m['inference']}" if m.get("inference") else "no post-DE label")
        hashes = m.get("input_hashes") or []
        inputs = f"{len(hashes)} verified input(s)" if hashes else ("none (not run)" if m["state"] != "COMPLETED" else "none recorded")
        rows.append([m["module"], "required" if m["required"] else "optional", m["plan_eligibility"], m["state"], m.get("reason_code") or "", reasons[:400], claim, inputs, m["state_meaning"]])
    section = {"state": "COMPLETED", "title": "Post-DE eligibility and what each conclusion depends on", "claim_label": "descriptive",
               "summary": [f"{len(modules)} post-DE module(s) requested; states: " + ", ".join(f"{m['module']} {m['state']}" for m in modules) + ".",
                           "NOT_RUN means software or an implementation was unavailable (not a scientific refusal); INAPPLICABLE means the design or data are scientifically ineligible (typed reason). Neither shows values."],
               "tables": [{"caption": "Post-DE modules: eligibility decided before computation and the executed state",
                           "columns": ["module", "requirement", "plan eligibility", "state", "reason code", "reason", "claim label", "inputs", "meaning"], "rows": rows, "source": f"{prefix}/eligibility.json"}],
               "figures": [], "notes": [], "methods": ["Eligibility rules are published per module (SM41); the planner records ELIGIBLE or INAPPLICABLE with a typed reason before any computation."],
               "values": {"n_modules": len(modules), "n_primary_discoveries": len(deps)}}
    sub = [[m["module"], s["analysis"], s["item"], s["state"], s["reason_code"], s["reason"][:300]] for m in modules for s in m.get("subanalyses") or []]
    if sub:
        section["tables"].append({"caption": "Refused or inapplicable sub-analyses (decided by the planner)", "columns": ["module", "analysis", "item", "state", "reason code", "reason"],
                                  "rows": sub, "source": f"{prefix}/eligibility.json"})
    adapted = [[m["module"], a.get("analysis", ""), str(a.get("item", "")), json.dumps(a.get("requested"), ensure_ascii=False).strip('"'), json.dumps(a.get("used"), ensure_ascii=False).strip('"'),
                str(a.get("reason", ""))[:300]] for m in modules for a in m.get("adaptations") or []]
    if adapted:   # D-59: adaptations are never silent
        section["tables"].append({"caption": "Adaptations to the data (requested versus used, with the reason)", "columns": ["module", "analysis", "item", "requested", "used", "reason"],
                                  "rows": adapted, "source": f"{prefix}/eligibility.json"})
        section["values"]["n_adaptations"] = len(adapted)
    if deps:
        section["tables"].append({"caption": "What each primary discovery depends on (descriptive; module states are shown where a module did not run)",
                                  "columns": ["contrast", "feature", "q", "sets", "robustness fraction", "sign stability", "influence flag", "phenotype associations", "biomarker selection", "module"],
                                  "rows": [[d["contrast_id"], d["feature_id"], d["q_value"], d["sets"], d["robustness_fraction"], d["sign_stability"], d["influence_flag"], d["phenotype_associations"],
                                            d["biomarker_selection_frequency"], d["coabundance_module"]] for d in deps],
                                  "source": f"{prefix}/dependency_summary.tsv"})
    else:
        section["notes"].append("The primary families have no discoveries at their cutoffs; the dependency summary is empty (a real completed result, not evidence of no effect).")
    return section
