"""Small factual methods/status summary for the R10a report (FR-094 scope).

Generated only from ReportData (resolved configuration plus executed stage
metadata).  It names no engine that did not run and asserts nothing that the
input declarations leave unknown.  Full publication methods are R10b.
"""
from __future__ import annotations


def render_methods(data: dict) -> str:
    run = data["run"]
    lines = [f"# Methods and status summary — {run['run_id']}", "", f"Run state: **{run['state']}** (exit code {run['exit_code']}). {run['phase_label']}", ""]
    inputs = data["inputs"].get("values")
    if inputs:
        scale = inputs["scale"]
        tech = inputs["technical_replicates"]
        lines += ["## Input", "",
                  f"- Declared {inputs['grain']}-level {inputs['assay']} input ({inputs['input_format']}); {inputs['n_features']} features.",
                  f"- Source scale `{scale['source_scale']}`, transform `{scale['applied_transform']}`, analysis scale `{scale['output_scale']}`.",
                  f"- {inputs['n_source_observations']} source observations; {inputs['n_canonical_observations']} canonical biological units; technical handling `{tech['mode']}`"
                  + (f" ({tech['method']})." if tech.get("method") else "."), ""]
    qc = data["qc"]
    if qc.get("state") == "COMPLETED":
        summary = qc["summary"]
        lines += ["## Preprocessing", "", f"- Normalization `{summary['normalization']}`; primary values changed: {str(summary['primary_values_changed']).lower()}; no new primary imputation.",
                  f"- Original-observed mask: {summary['original_mask_state']}. Declared exclusions: {len(qc['excluded_observations'])}.", ""]
    design = data["design"]
    if design.get("state") == "COMPLETED":
        lines += ["## Design", ""]
        for d in design["designs"]:
            lines.append(f"- Design `{d['design_id']}`: {d['n_coefficients']} coefficients, rank {d['rank']}, blocking `{d['blocking_mode']}`.")
        lines.append("")
    dea = data["dea"]
    lines += ["## Differential abundance", ""]
    if dea.get("state") == "COMPLETED":
        for model_id, settings in dea["settings"].items():
            lines.append(f"- Model `{model_id}`: limma {settings['limma_version']}, eBayes trend={str(settings['trend_applied']).lower()}, robust={str(settings['robust_applied']).lower()} "
                         f"(robust empirical-Bayes hyperparameters), exact contrast-as-coefficient refits, genuinely observed cells only.")
        for family in dea["families"]:
            lines.append(f"- Family `{family['family_id']}` ({family['hypothesis_type']}, {family['role']}): {family['adjustment']} over {family['n_finite']} finite of {family['n_planned']} planned tests; "
                         f"{family['completeness']}; {family['rejection_count']} at q ≤ {family['q_cutoff']}.")
    else:
        lines.append(f"- Not available: {dea.get('state')}" + (f" ({dea.get('reason_code')})." if dea.get("reason_code") else "."))
    lines.append("")
    for name, section in data.get("sections", {}).items():
        lines += [f"## {section.get('title', name)}", ""]
        for item in section.get("methods", []) if section.get("state") == "COMPLETED" else [f"Not available: {section.get('state')} {section.get('reason_code') or ''}".strip()]:
            lines.append(f"- {item}")
        lines.append("")
    lines += ["## Stage status", ""]
    for stage in data["stages"]:
        lines.append(f"- `{stage['stage_id']}`: {stage['display_state']}" + (f" ({stage['reason_code']})" if stage.get("reason_code") else ""))
    lines += ["", "## Limitations", ""] + [f"- {item}" for item in data["limitations"]] + [""]
    return "\n".join(lines)
