"""Post-DE biomarker discrimination section of ReportData (packet R14d, FR-149, FR-160, SM34, SM25).

Built only from verified post_de/biomarker outputs.  Every table and the section carry their claim labels; n per class,
the CV scheme, seeds and the permutation B are always shown, and "not externally validated" is stated when no
independent cohort was declared.  No diagnostic or clinical-use claim is made.
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


def post_de_biomarker_section(directory: str | Path, *, prefix: str = "post_de/biomarker") -> dict:
    root = Path(directory)
    result = json.loads((root / "biomarker_result.json").read_text(encoding="utf-8"))
    eligibility = json.loads((root / "eligibility.json").read_text(encoding="utf-8"))
    cv = _table(root / "cv_performance.tsv")[0]
    perm = result.get("permutation") or {}
    section = {"state": "COMPLETED", "title": "Post-DE biomarker discrimination (performance estimates)", "claim_label": eligibility["claim_label"],
               "summary": [f"Binary contrast {result['contrast']}: {cv['n_positive_units']} {result['positive_class']} vs {cv['n_negative_units']} {result['negative_class']} biological units "
                           f"(grouped by {cv['grouped_by']}); classifier {result['classifier']}, selection {result['selection']['method']}.",
                           f"Nested CV ({cv['scheme']}, outer k = {cv['outer_k']}, inner k = {cv['inner_k']}, repeats = {cv['repeats']}; seed {result['seed']}): pooled out-of-fold AUC "
                           f"{_f(cv['pooled_oof_auc'])}, per-repeat mean {_f(cv['repeat_auc_mean'])} (SD {_f(cv['repeat_auc_sd'])}) [cross_validated_nested].",
                           (f"Whole-procedure permutation test (selection and tuning rerun; B = {cv['permutation_B']}, seed {result['seed']}, scheme {cv['permutation_scheme']}): "
                            f"P = {_f(cv['permutation_p'])} [cross_validated_nested]." if perm.get("enabled") else "No permutation test was requested."),
                           f"External validation: {result['external_validation']}."],
               "tables": [{"caption": "Nested cross-validated performance [cross_validated_nested]",
                           "columns": ["pooled AUC", "repeat mean", "repeat SD", "Brier", "threshold rule", "sensitivity", "specificity", "permutation P", "claim"],
                           "rows": [[_f(cv["pooled_oof_auc"]), _f(cv["repeat_auc_mean"]), _f(cv["repeat_auc_sd"]), _f(cv["brier"]), cv["threshold_rule"],
                                     f"{_f(cv['sensitivity'])} ({_f(cv['sensitivity_lower'])}-{_f(cv['sensitivity_upper'])})",
                                     f"{_f(cv['specificity'])} ({_f(cv['specificity_lower'])}-{_f(cv['specificity_upper'])})", _f(cv["permutation_p"]), cv["claim_label"]]],
                           "source": f"{prefix}/cv_performance.tsv"}],
               "figures": [], "notes": list(result.get("limitations", [])), "methods": [eligibility.get("rule", ""), result.get("rng", "")],
               "values": {"pooled_oof_auc": cv["pooled_oof_auc"]}}
    single = _table(root / "single_feature_auc.tsv")
    if single:
        top = sorted(single, key=lambda r: -abs(float(r["auc"]) - 0.5))[:15]
        section["tables"].append({"caption": "Single-feature AUC on biological units [in_sample]" + ("; auto-direction values are direction-optimistic" if single[0]["direction_optimistic"] == "true" else ""),
                                  "columns": ["feature", "AUC", "direction", "DeLong CI", "bootstrap CI", "claim"],
                                  "rows": [[r["feature_id"], _f(r["auc"]), r["direction"], f"{_f(r['delong_lower'])}-{_f(r['delong_upper'])}", f"{_f(r['bootstrap_lower'])}-{_f(r['bootstrap_upper'])}", r["claim_label"]] for r in top],
                                  "source": f"{prefix}/single_feature_auc.tsv"})
    panels = _table(root / "fixed_panels.tsv")
    if panels:
        section["tables"].append({"caption": "Fixed panels [fixed_panel_cv]; same-data or unknown panels are selection-optimistic and shown beside a nested reselection of the same size",
                                  "columns": ["panel", "provenance", "features", "fixed-panel AUC", "optimistic", "nested reselection AUC", "claim"],
                                  "rows": [[p["panel_id"], p["provenance"], p["n_features"], _f(p["pooled_oof_auc"]), p["selection_optimistic"], _f(p["nested_reselection_auc"]), p["claim_label"]] for p in panels],
                                  "source": f"{prefix}/fixed_panels.tsv"})
    validation = _table(root / "validation.tsv")
    if validation:
        v = validation[0]
        section["tables"].append({"caption": f"External validation [{v['claim_label']}]", "columns": ["state", "AUC", "locked model SHA-256", "claim"],
                                  "rows": [[v["state"], _f(v.get("auc")), v.get("locked_model_sha256"), v["claim_label"]]], "source": f"{prefix}/validation.tsv"})
    refusals = _table(root / "refusals.tsv")
    if refusals:
        section["tables"].append({"caption": "Refused or not applicable sub-analyses", "columns": ["analysis", "item", "reason code", "reason"],
                                  "rows": [[r["analysis"], r["item"], r["reason_code"], r["reason"]] for r in refusals], "source": f"{prefix}/refusals.tsv"})
    for record in eligibility.get("figures", []):
        if record.get("state") == "COMPLETED" and record["file"].endswith(".png"):
            stem = Path(record["file"]).stem
            section["figures"].append({"src": f"{prefix}/{record['file']}", "alt": stem.replace("_", " "), "caption": stem.replace("_", " ").capitalize() + " [cross_validated_nested].",
                                       "source": f"{prefix}/figure_sources/{stem}.tsv"})
    return section
