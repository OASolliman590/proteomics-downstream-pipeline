"""Thin offline HTML report stage (packet R10a: FR-091 to FR-094).

Renders ReportData into one self-contained HTML page: embedded CSS, plain
accessible tables, relative links to the exact source tables, no scripts
and no network resources.  Absolute local paths are redacted in the
presentation only.
"""
from __future__ import annotations

import html
import json
import re
from pathlib import Path
from string import Template

from ..provenance import sha256_file
from .assemble import assemble, validate_report_data
from .stub_methods import render_methods

CAPABILITY = "report_stub"
_ABSOLUTE = re.compile(r"(?:[A-Za-z]:\\|/)(?:[^\s\"'<>|]+[/\\])+([^\s\"'<>|/\\]+)")


class Raw(str):
    """Already-escaped HTML fragment produced by this module."""


def capabilities():
    return [{"id": CAPABILITY, "dependencies": [], "required_r_packages": []}]


def redact(text) -> str:
    """Replace absolute local paths by their final component for presentation."""
    return _ABSOLUTE.sub(lambda m: f"…/{m.group(1)}", str(text))


def e(value) -> str:
    if value is None:
        return '<span class="muted">unknown</span>'
    if isinstance(value, bool):
        return "true" if value else "false"
    return html.escape(redact(value))


def link(source: str | None, label: str | None = None) -> str:
    if not source:
        return Raw('<span class="muted">no source table</span>')
    return Raw(f'<a href="../{html.escape(source)}">{html.escape(label or source)}</a>')


def table(caption: str, columns: list[str], rows: list[list], source: str | None = None, numeric: tuple[str, ...] = ()) -> str:
    head = "".join(f'<th scope="col">{html.escape(c)}</th>' for c in columns)
    body = []
    for row in rows:
        cells = []
        for column, value in zip(columns, row):
            cls = ' class="num"' if column in numeric else ""
            cells.append(f"<td{cls}>{value if isinstance(value, Raw) else e(value)}</td>")
        body.append(f"<tr>{''.join(cells)}</tr>")
    if not body:
        body.append(f'<tr><td colspan="{len(columns)}" class="muted">No rows (a real empty result, not missing data).</td></tr>')
    src = f'<p class="muted">Source: {link(source)}</p>' if source else ""
    return f'<div class="table-wrap"><table><caption>{html.escape(caption)}</caption><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>{src}'


def state_span(state: str | None) -> str:
    state = state or "NOT_RUN"
    return Raw(f'<span class="state-{html.escape(state)}">{html.escape(state)}</span>')


def _fmt(value, digits=4):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return value
    return f"{number:.{digits}g}"


def render(data: dict) -> str:
    run = data["run"]
    sections: list[tuple[str, str, str]] = []
    status = (f'<div class="banner" role="status"><p>Run state: {state_span(run["state"])} — exit code {e(run["exit_code"])}'
              + (f" — reason {e(run['reason_code'])}" if run["reason_code"] else "") + "</p>"
              + f'<p class="muted">Run {e(run["run_id"])}; plan hash <code>{e(run["plan_hash"])}</code>; requested phase {e(run["requested_phase"])}; scope {e(run["scope"])}.</p></div>')
    stage_rows = [[s["stage_id"], s["capability"], "required" if s["required"] else "optional", state_span(s["display_state"]), s["reason_code"] or "none",
                   "verified" if s["verified"] else (s.get("verification_problem") or "not applicable")] for s in data["stages"]]
    sections.append(("status", "Status", status + table("Executed stages and their actual states", ["stage", "capability", "requiredness", "state", "reason", "manifest"], stage_rows)))

    inputs = data["inputs"]
    if inputs.get("values"):
        v = inputs["values"]
        scale = v["scale"]
        body = (f"<p>Declared {e(v['grain'])}-level {e(v['assay'])} input ({e(v['input_format'])}); source scale <code>{e(scale['source_scale'])}</code>, transform "
                f"<code>{e(scale['applied_transform'])}</code>, analysis scale <code>{e(scale['output_scale'])}</code>. {e(v['n_features'])} features; "
                f"{e(v['n_source_observations'])} source observations (injections) became {e(v['n_canonical_observations'])} biological units from {e(v['n_subjects'])} subjects.</p>")
        body += table("Biological versus technical n by group", ["group", "injections", "biological units", "subjects", "small n"],
                      [[g["group"], g["n_injections"], g["n_biological_units"], g["n_subjects"], g["small_n"]] for g in v["groups"]], inputs["sources"]["sample_n"], ("injections", "biological units", "subjects"))
        m = v["missingness"]
        body += table("Missing values", ["total cells", "numerically missing", "originally observed", "prior-imputed", "original mask"],
                      [[m["total_cells"], m["numeric_missing_cells"], m["originally_observed_cells"], m["prior_imputed_cells"], m["original_observed_mask"]["state"]]], inputs["sources"]["missingness"])
        gaps = v["provenance_gaps"]
        body += table("Unknown provenance (not asserted)", ["field", "state", "consequence"], [[g["field"], g["state"], g["consequence"]] for g in gaps], inputs["source"])
        body += f'<p class="muted">Observation lineage: {link(inputs["lineage_source"])} ({inputs["n_lineage_rows"]} source observations mapped).</p>'
    else:
        body = f"<p>Input stage {state_span(inputs.get('state'))}: no input values are shown.</p>"
    sections.append(("inputs", "Input and biological identity", body))

    qc = data["qc"]
    if qc.get("state") == "COMPLETED":
        s = qc["summary"]
        t = qc["tables"]
        body = (f"<p>Normalization <code>{e(s['normalization'])}</code> (primary values changed: {e(s['primary_values_changed'])}); original-observed mask {e(s['original_mask_state'])}; "
                f"{e(s['n_features'])} features × {e(s['n_observations'])} analysed observations. No new primary imputation.</p>")
        rows = t["sample_n"]["rows"] or []
        body += table("Analysed biological units by group (after declared exclusions)", ["group", "biological units", "injections", "subjects"],
                      [[r["group"], r["n_biological_units"], r["n_injections"], r["n_subjects"]] for r in rows], t["sample_n"]["source"], ("biological units", "injections", "subjects"))
        excl = t["exclusions"]["rows"] or []
        body += table("Prespecified exclusions and their recorded reasons", ["observation", "reason"], [[r["observation_id"], r["reason"]] for r in excl], t["exclusions"]["source"])
        watch = [r for r in (t["watchlist"]["rows"] or []) if r["flagged"] == "true"]
        body += table("QC watchlist flags (flag only; no sample removed)", ["observation", "PCA robust z", "reasons", "action"],
                      [[r["observation_id"], _fmt(r["pca_robust_z"]), r["reasons"], r["action"]] for r in watch], t["watchlist"]["source"])
        miss = t["missingness_by_group"]["rows"] or []
        body += table("Missingness by group", ["group", "units", "cells", "numeric", "originally observed", "prior-imputed"],
                      [[r["group"], r["n_units"], r["n_cells"], r["n_numeric"], r["n_original_observed"], r["n_prior_imputed"]] for r in miss], t["missingness_by_group"]["source"])
        norm = t["normalization_factors"]["rows"] or []
        body += table("Normalization factors (subtracted, log2)", ["observation", "method", "factor"], [[r["observation_id"], r["method"], _fmt(r["factor"])] for r in norm], t["normalization_factors"]["source"])
        pca = qc["pca"]
        if t["pca_variance"]["rows"]:
            body += table("PCA (display-only median fill, centred; never an inference input)", ["component", "variance explained"],
                          [[r["component"], _fmt(r["variance_explained"])] for r in t["pca_variance"]["rows"]], t["pca_variance"]["source"])
        else:
            body += f"<p>PCA {state_span(pca.get('state'))}: {e(pca.get('reason_code'))}. No variance values are shown.</p>"
        if qc.get("missingness_statement"):
            body += f'<p class="muted">{e(qc["missingness_statement"]["note"])}</p>'
    else:
        body = f"<p>Preprocessing/QC {state_span(qc.get('state'))}: no QC values are shown.</p>"
    sections.append(("qc", "Quality control and inclusion", body))

    design = data["design"]
    if design.get("state") == "COMPLETED":
        body = table("Designs", ["design", "observations", "coefficients", "rank", "full rank", "blocking"],
                     [[d["design_id"], d["n_observations"], d["n_coefficients"], d["rank"], d["full_rank"], d["blocking_mode"]] for d in design["designs"]], design["source"])
        body += table("Planned contrasts (all roles)", ["contrast", "role", "required groups", "estimable", "biological units", "subjects"],
                      [[c["contrast_id"], c["role"], ", ".join(c["required_groups"]), c["estimable"], json.dumps(c["units_by_required_group"]), json.dumps(c["subjects_by_required_group"])] for c in design["contrasts"]])
    else:
        body = f"<p>Design {state_span(design.get('state'))}.</p>"
    sections.append(("design", "Design and contrasts", body))

    dea = data["dea"]
    if dea.get("state") == "COMPLETED":
        body = table("Multiplicity families (prespecified)", ["family", "hypothesis", "role", "adjustment", "planned", "eligible", "finite", "numerical failures", "completeness", "q cutoff", "rejections"],
                     [[f["family_id"], f["hypothesis_type"], f["role"], f["adjustment"], f["n_planned"], f["n_eligible"], f["n_finite"], f["n_numerical_failure"], f["completeness"], f["q_cutoff"], f["rejection_count"]] for f in dea["families"]],
                     dea["families_source"], ("planned", "eligible", "finite", "numerical failures", "rejections"))
        body += table("Every planned endpoint group, including secondary contrasts", ["model", "contrast", "hypothesis", "role", "family", "planned", "tested", "excluded", "nonestimable", "numerical failure", "q ≤ cutoff"],
                      [[c["model_id"], c["contrast_id"], c["hypothesis_type"], c["role"], c["family_id"], c["n_planned"], c["n_tested"], c["n_excluded"], c["n_nonestimable"], c["n_numerical_failure"], c["n_q_at_or_below_cutoff"]] for c in dea["contrast_summaries"]],
                      None, ("planned", "tested", "excluded", "nonestimable", "numerical failure", "q ≤ cutoff"))
        body += "<p>Complete result tables: " + ", ".join(link(t["source"], f"{name} ({t['n_rows']} rows)") for name, t in dea["tables"].items()) + ".</p>"
        if dea["influence_sources"]:
            body += "<p>Whole-subject influence (descriptive, not held-out validation): " + ", ".join(link(s) for s in dea["influence_sources"]) + ".</p>"
        body += f'<p class="muted">{e(dea["note"])}</p>'
    else:
        body = f"<p>Differential analysis {state_span(dea.get('state'))}" + (f" ({e(dea.get('reason_code'))})" if dea.get("reason_code") else "") + f". {e(dea.get('note'))}</p>"
    sections.append(("dea", "Differential abundance", body))

    for name, section in data["sections"].items():
        if section.get("state") != "COMPLETED":
            body = f"<p>{state_span(section.get('state'))}" + (f" ({e(section.get('reason_code'))})" if section.get("reason_code") else "") + ". No values are shown.</p>"
        else:
            body = "".join(f"<p>{e(p)}</p>" for p in section.get("summary", []))
            for t in section.get("tables", []):
                body += table(t["caption"], t["columns"], t["rows"], t.get("source"), tuple(t.get("numeric", ())))
            for fig in section.get("figures", []):
                body += f'<figure><img src="../{html.escape(fig["src"])}" alt="{html.escape(fig["alt"])}"><figcaption>{e(fig["caption"])} Source data: {link(fig.get("source"))}</figcaption></figure>'
            body += "".join(f'<p class="muted">{e(n)}</p>' for n in section.get("notes", []))
        sections.append((f"section-{name}", section.get("title", name), body))

    warnings = data["warnings"]
    body = table("Warnings", ["code", "stage", "severity", "message"], [[w.get("code"), w.get("stage_id"), w.get("severity"), w.get("message")] for w in warnings])
    body += "<ul>" + "".join(f"<li>{e(item)}</li>" for item in data["limitations"]) + "</ul>"
    body += '<p>Methods and status summary: <a href="methods.md">methods.md</a>; machine-readable report data: <a href="report_data.json">report_data.json</a>.</p>'
    sections.append(("limitations", "Warnings and limitations", body))

    nav = "".join(f'<li><a href="#{anchor}">{html.escape(title)}</a></li>' for anchor, title, _ in sections)
    body = "".join(f'<section id="{anchor}" aria-labelledby="{anchor}-h"><h2 id="{anchor}-h">{html.escape(title)}</h2>{content}</section>' for anchor, title, content in sections)
    template = Template((Path(__file__).parent / "stub_template.html").read_text(encoding="utf-8"))
    return template.substitute(title=html.escape(f"Proteomics run report — {run['run_id']}"), phase_label=html.escape(run["phase_label"]), nav=nav, body=body)


def execute(request: dict) -> dict:
    from ..errors import ProteomicsError
    from ..runtime import stage_result, validate_stage_request, utc_now
    validate_stage_request(request)
    out = Path(request["output_temp_dir"])
    out.mkdir(parents=True, exist_ok=False)
    started = utc_now()
    try:
        data = assemble(request["parameters"]["run_root"], request["parameters"]["run_snapshot"])
        validate_report_data(data)
        (out / "report_data.json").write_text(json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True), encoding="utf-8")
        (out / "methods.md").write_text(render_methods(data), encoding="utf-8")
        (out / "index.html").write_text(render(data), encoding="utf-8")
    except ProteomicsError as error:
        result = stage_result(request["run_id"], request["stage_id"], CAPABILITY, "FAILED", plan_hash=request["plan_hash"], exit_code=error.exit_code, reason_code=error.code, message=error.message)
        result["started_at"] = started
        return result
    except Exception as error:  # renderer bug: report fails, analysis state is not upgraded
        result = stage_result(request["run_id"], request["stage_id"], CAPABILITY, "FAILED", plan_hash=request["plan_hash"], exit_code=4, reason_code="E_REPORT_RENDER", message=f"{type(error).__name__}: {error}")
        result["started_at"] = started
        return result
    result = stage_result(request["run_id"], request["stage_id"], CAPABILITY, "COMPLETED", plan_hash=request["plan_hash"], exit_code=0, message="thin offline report rendered")
    result["started_at"] = started
    result["outputs"] = [{"artifact_id": a, "relative_path": p, "sha256": sha256_file(out / p), "result_type": t}
                         for a, p, t in (("report_data", "report_data.json", "ReportData"), ("methods", "methods.md", "MethodsSummary"), ("report_html", "index.html", "ReportHTML"))]
    return result
