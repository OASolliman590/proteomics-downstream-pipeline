"""PERMANOVA section of ReportData (packet R13, FR-129).

Built only from the verified PERMANOVA stage outputs; every table and
figure carries its source path.  The interpretation rule and the
selection-circularity flags are always shown.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _f(value, digits=3):
    try:
        return f"{float(value):.{digits}g}"
    except (TypeError, ValueError):
        return value


def permanova_section(directory: str | Path, *, prefix: str = "permanova") -> dict:
    root = Path(directory)
    result = json.loads((root / "permanova_result.json").read_text(encoding="utf-8"))
    tests = _table(root / "tests.tsv")
    rows = [[t["feature_set_id"], t["analysis"], t["comparison"], t["term"], t["n"], t["n_features"], t["df"], _f(t["r2"]), _f(t["pseudo_f"]),
             t["p_display"], t["p_adjusted_display"] if t["p_adjusted_display"] != "NA" else "", t["permutation_scheme"],
             t["permdisp_p_display"] if t["permdisp_p_display"] != "NA" else "not applicable", t["interpretation"]] for t in tests]
    section = {"state": "COMPLETED", "title": "Multivariate separation (PERMANOVA and PERMDISP)",
               "summary": [f"Distance: {result['settings']['metric']} on {result['settings']['scaling']}-scaled genuinely observed log2 abundance; "
                           f"{result['settings']['permutations_requested']} permutations, seed {result['settings']['seed']}; pairwise adjustment {result['settings']['adjustment']}.",
                           result["interpretation_rule"]],
               "tables": [{"caption": "PERMANOVA terms with their PERMDISP results", "columns": ["feature set", "analysis", "comparison", "term", "n", "features", "Df", "R2", "pseudo-F",
                                                                                                    "P", "adjusted P", "permutation scheme", "PERMDISP P", "interpretation"],
                           "rows": rows, "source": f"{prefix}/tests.tsv", "numeric": ["n", "features", "Df", "R2", "pseudo-F"]}],
               "figures": [], "notes": list(result["limitations"]), "methods": [], "values": {"n_tests": len(tests)}}
    identity = result.get("identity", [])
    if identity:
        section["tables"].append({"caption": "Mean per-feature R2 versus multivariate R2 (identity under z-scored Euclidean distance)",
                                  "columns": ["feature set", "multivariate R2", "mean per-feature R2", "state"],
                                  "rows": [[i["feature_set_id"], _f(i["multivariate_r2"], 6), _f(i["mean_univariate_r2"], 6), i["state"]] for i in identity],
                                  "source": f"{prefix}/identity_check.json"})
    context_path = root / "selection_context.tsv"
    if context_path.is_file():
        context = _table(context_path)
        section["tables"].append({"caption": "Selection-circularity context: in-sample R2 is optimistic for sets chosen on these data",
                                  "columns": ["feature set", "k", "set R2", "null median", "null 95th pct", "fraction of null ≥ set", "best possible R2", "flag"],
                                  "rows": [[c["feature_set_id"], c["k"], _f(c["set_r2"]), _f(c["null_median"]), _f(c["null_q95"]), _f(c["fraction_null_at_or_above_set"]),
                                            _f(c["best_possible_r2"]), "in-sample R2 is optimistic" if c["in_sample_optimistic"] == "true" else ""] for c in context],
                                  "source": f"{prefix}/selection_context.tsv"})
    refusals_path = root / "refusals.tsv"   # audit 2026-10-02: untested comparisons/terms are shown, never silently dropped
    refusals = _table(refusals_path) if refusals_path.is_file() else []
    if refusals:
        section["tables"].append({"caption": "Planned comparisons and terms that were not tested", "columns": ["feature set", "analysis", "comparison", "reason code", "reason"],
                                  "rows": [[r["feature_set_id"], r["analysis"], r["comparison"], r["reason_code"], r["reason"]] for r in refusals],
                                  "source": f"{prefix}/refusals.tsv"})
    family = {(t["adjustment_family_size"], t["adjustment_family_planned"]) for t in tests if t["analysis"] == "pairwise" and t.get("adjustment_family_size") not in (None, "", "NA")}
    for size, planned in sorted(family):
        section["notes"].append(f"Pairwise {result['settings']['adjustment']} adjustment covered {size} tested pair(s) of {planned} planned.")
    states = result.get("feature_sets", {})
    refused = [f"{k}: {v.get('reason')}" for k, v in states.items() if v.get("state") != "COMPLETED"]
    if refused:
        section["notes"].append("Feature sets not tested: " + "; ".join(refused) + ".")
    sources = {"ordination": "ordination.tsv", "permdisp": "ordination.tsv", "feature_r2": "feature_r2.tsv", "random_null": "random_set_null.tsv", "null": "permutation_null.tsv", "r2_summary": "tests.tsv"}
    for record in result.get("figures", []):
        if record["state"] != "COMPLETED" or not record["file"].endswith(".png"):
            continue
        stem = Path(record["file"]).stem
        if stem.startswith("null_") and "_global_" not in stem:
            continue  # every permutation null is a published figure; the report embeds the global ones
        kind = next((k for k in sources if stem.startswith(k)), "r2_summary")
        section["figures"].append({"src": f"{prefix}/{record['file']}", "alt": stem.replace("_", " "), "caption": stem.replace("_", " ").capitalize() + ".", "source": f"{prefix}/{sources[kind]}"})
    section["methods"] = [f"PERMANOVA (vegan {result['versions']['vegan']}, permute {result['versions']['permute']}) with marginal sums of squares; covariate and interaction terms permuted within the primary group; "
                          f"group terms from unrestricted (or subject-blocked) permutations; PERMDISP (betadisper/permutest) with the same scheme.",
                          f"Seed {result['settings']['seed']} ({result['settings']['rng_kind']}); P floored at 1/(nperm+1)."]
    return section
