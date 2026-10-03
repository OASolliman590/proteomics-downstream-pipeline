"""Full methods and limitations (packet R10b, FR-098).

Generated only from the frozen plan (resolved configuration, families,
eligibility), the executed stage states and the recorded package inventory.
Engines that did not run are named only with their actual state; unknown
provenance stays explicit; no canned biological or rescue claim is made.
"""
from __future__ import annotations


def _versions(plan: dict | None) -> str:
    if not plan:
        return "unknown (no plan was frozen)"
    env = plan.get("environment", {})
    packages = ", ".join(f"{k} {v}" for k, v in sorted((env.get("r_packages") or {}).items()) if v)
    r = (env.get("rscript") or {}).get("version", "unknown R")
    return f"Python {env.get('python', {}).get('version', '?')}; {r}; {packages}"


def render_full_methods(data: dict, plan: dict | None) -> str:
    run = data["run"]
    executed = {s["stage_id"]: s for s in data["stages"]}
    lines = [f"# Methods and limitations — {run['run_id']}", "", f"Run state **{run['state']}** (exit {run['exit_code']}); plan `{run['plan_hash']}`.", "",
             "## Software actually executed", "", f"- {_versions(plan)}.", ""]
    inputs = data["inputs"].get("values")
    lines += ["## Input and design", ""]
    if inputs:
        s = inputs["scale"]
        lines.append(f"- {inputs['grain'].capitalize()}-level {inputs['assay']} data ({inputs['input_format']}), source scale {s['source_scale']}, transform {s['applied_transform']}, analysis scale {s['output_scale']}; "
                     f"{inputs['n_features']} features, {inputs['n_canonical_observations']} biological units from {inputs['n_source_observations']} source observations.")
    else:
        lines.append("- Input stage did not complete; no input description is available.")
    if plan:
        for model in plan["models"]:
            state = "fitted" if executed.get("limma", {}).get("state") == "COMPLETED" and model["engine"] == "limma" else executed.get(f"model.{model['model_id']}", executed.get("assay_engines", {})).get("state", "not executed")
            lines.append(f"- Model `{model['model_id']}` ({model['engine']}, {model['role']}, {model['execution_requirement']}): {state}.")
        for family in plan["families"]:
            members = ", ".join(f"{m['model_id']}/{m['contrast_id']}" for m in family["members"])
            lines.append(f"- Family `{family['family_id']}` ({family['hypothesis_type']}, {family['role']}, {family['adjustment']}; {family['dependence_assumption']}): {members}.")
    lines.append("")
    dea = data["dea"]
    lines += ["## Differential abundance", ""]
    if dea.get("state") == "COMPLETED":
        for model_id, settings in dea["settings"].items():
            lines.append(f"- `{model_id}`: limma {settings['limma_version']} lmFit/eBayes(trend={settings['trend_applied']}, robust={settings['robust_applied']}) on genuinely observed cells; "
                         "each contrast estimated by an exact contrast-as-coefficient refit; 95% intervals use the moderated-t quantile.")
        for family in dea["families"]:
            lines.append(f"- `{family['family_id']}`: {family['n_finite']} of {family['n_planned']} planned tests finite; {family['rejection_count']} at q ≤ {family['q_cutoff']}; {family['completeness']}.")
    else:
        lines.append(f"- Not available ({dea.get('state')}).")
    for name, section in data.get("sections", {}).items():
        lines += ["", f"## {section.get('title', name)}", ""]
        if section.get("state") == "COMPLETED":
            lines += [f"- {m}" for m in section.get("methods", [])]
        else:
            lines.append(f"- Not executed: {section.get('state')} ({section.get('reason_code') or 'no reason recorded'}).")
    lines += ["", "## Limitations", ""] + [f"- {item}" for item in data["limitations"]]
    lines += ["- Results are statistical associations within the declared design; they do not establish mechanism or causal drug effects.",
              "- Calibration of error rates on simulated data is a separate validation gate and is not implied by this run.", ""]
    return "\n".join(lines)
