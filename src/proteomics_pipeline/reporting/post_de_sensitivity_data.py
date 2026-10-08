"""Post-DE robustness and sensitivity section of ReportData (packet R14b, FR-137-FR-143, SM25).

Built only from verified post_de/sensitivity outputs.  The non-interaction warning, the estimability
classes, refused analyses and the descriptive labels are always shown.
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


def post_de_sensitivity_section(directory: str | Path, *, prefix: str = "post_de/sensitivity") -> dict:
    root = Path(directory)
    eligibility = json.loads((root / "eligibility.json").read_text(encoding="utf-8"))
    engine = eligibility.get("primary_engine", {})
    section = {"state": "COMPLETED", "title": "Post-DE robustness and sensitivity", "claim_label": "descriptive",
               "summary": [f"Sensitivity models are named secondary models fitted with the primary engine settings (limma, trend={engine.get('trend')}, robust={engine.get('robust')}), each in its own family; primary families are unchanged.",
                           eligibility.get("warning", "")],
               "tables": [], "figures": [], "notes": [], "methods": [eligibility.get("rule", "")], "values": {}}
    classes = _table(root / "estimability_classes.tsv")
    if classes:
        section["tables"].append({"caption": "Estimability of each requested adjustment (decided before any refit)", "columns": ["model", "kind", "class", "reason code", "detail"],
                                  "rows": [[c["model_id"], c["kind"], c["class"], c["reason_code"], c["detail"]] for c in classes], "source": f"{prefix}/estimability_classes.tsv"})
    imbalance = _table(root / "imbalance.tsv")
    if imbalance:
        section["tables"].append({"caption": "Covariate x group cross-tabulation (biological units)", "columns": ["covariate", "level", "group", "units", "aliased with group"],
                                  "rows": [[r["covariate"], r["level"], r["group"], r["n_units"], r.get("aliased_with_group", "")] for r in imbalance], "source": f"{prefix}/imbalance.tsv"})
    comparison = _table(root / "comparison.tsv")
    if comparison:
        section["tables"].append({"caption": "Primary versus sensitivity models (descriptive)", "columns": ["model", "contrast", "criterion", "retained", "lost", "gained", "Jaccard", "effect r", "attenuation slope", "interaction family"],
                                  "rows": [[c["model_id"], c["contrast_id"], c["criterion"], c["n_retained"], c["n_lost"], c["n_gained"], _f(c["jaccard"]), _f(c["effect_correlation"]), _f(c["attenuation_slope"]),
                                            c["interaction_family"] if c["interaction_family"] != "NA" else "not estimable / not requested"] for c in comparison],
                                  "source": f"{prefix}/comparison.tsv"})
    settings = _table(root / "matched_n_settings.tsv")
    if settings:
        s = settings[0]
        section["notes"].append(f"Matched-n resampling: {s['draws']} draws, seed {s['seed']} ({s['rng_kind']}), units = {s['unit_level']}, targets {s['targets']} (descriptive).")
    influence = _table(root / "influence.tsv")
    if influence:
        top = sorted((r for r in influence if r["state"] == "refit"), key=lambda r: int(r["rank"]))[:10]
        section["tables"].append({"caption": "Most influential biological units (leave-one-unit-out; descriptive)", "columns": ["unit", "contrast", "omitted observations", "lost", "gained", "max |effect change|", "rank"],
                                  "rows": [[r["unit_id"], r["contrast_id"], r["omitted_observation_ids"], r["n_lost"], r["n_gained"], _f(r["max_abs_effect_change"]), r["rank"]] for r in top],
                                  "source": f"{prefix}/influence.tsv"})
    robust = _table(root / "robustness_summary.tsv")
    if robust:
        members = [r for r in robust if r["primary_member"] == "true"]
        section["tables"].append({"caption": "Robustness of primary discoveries (fraction of analyses meeting the criterion; descriptive, never re-adjusted)",
                                  "columns": ["contrast", "feature", "robustness fraction", "sensitivity models", "influence retained", "matched-n frequency", "sign stability", "influence flag"],
                                  "rows": [[r["contrast_id"], r["feature_id"], _f(r["robustness_fraction"]), _f(r["sensitivity_model_fraction"]), _f(r["influence_retained_fraction"]),
                                            _f(r["matched_n_selection_frequency"]), _f(r["sign_stability"]), r["influence_flag"]] for r in members],
                                  "source": f"{prefix}/robustness_summary.tsv"})
    refusals = _table(root / "refusals.tsv")
    if refusals:
        section["tables"].append({"caption": "Requested sensitivity analyses that were refused", "columns": ["analysis", "item", "reason code", "reason"],
                                  "rows": [[r["analysis"], r["item"], r["reason_code"], r["reason"]] for r in refusals], "source": f"{prefix}/refusals.tsv"})
    for record in eligibility.get("figures", []):
        if record.get("state") == "COMPLETED" and record["file"].endswith(".png"):
            stem = Path(record["file"]).stem
            section["figures"].append({"src": f"{prefix}/{record['file']}", "alt": stem.replace("_", " "), "caption": stem.replace("_", " ").capitalize() + " (descriptive).",
                                       "source": f"{prefix}/figure_sources/{stem}.tsv"})
    return section
