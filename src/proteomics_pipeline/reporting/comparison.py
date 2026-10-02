"""Run comparison (packet R10b, FR-100: ``proteomics compare``).

Both runs are verified first (manifests and plan hashes).  The comparison is
keyed by stable identifiers (model, contrast, hypothesis, feature); it reports
configuration differences, stage states, family counts and per-endpoint effect
and P differences, and never edits either run.
"""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path

from ..errors import CollisionError, ProteomicsError
from ..runtime import validate_run_status


def _flatten(value, prefix=""):
    if isinstance(value, dict):
        out = {}
        for k, v in value.items():
            out.update(_flatten(v, f"{prefix}/{k}"))
        return out
    if isinstance(value, list):
        out = {}
        for i, v in enumerate(value):
            out.update(_flatten(v, f"{prefix}/{i}"))
        return out
    return {prefix: value}


def _rows(path: Path) -> dict:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        return {(r["model_id"], r["contrast_id"], r["hypothesis_type"], r["feature_id"]): r for r in csv.DictReader(handle, delimiter="\t")}


def _load_run(root: Path) -> dict:
    from ..workflow import verify_run
    status_path = root / "run_status.json"
    if not status_path.is_file():
        raise ProteomicsError("E_INTEGRITY", f"{root.name} has no run_status.json", exit_code=5)
    status = json.loads(status_path.read_text(encoding="utf-8"))
    try:
        validate_run_status(status, root)
    except ValueError as exc:
        raise ProteomicsError("E_INTEGRITY", f"{root.name}: {exc}", exit_code=5) from exc
    verify_run(root)
    plan = json.loads((root / "plan.json").read_text(encoding="utf-8")) if (root / "plan.json").is_file() else None
    return {"status": status, "plan": plan}


def compare(left: str | Path, right: str | Path, output: str | Path) -> dict:
    left, right, output = Path(left).resolve(), Path(right).resolve(), Path(output).resolve()
    if output.exists():
        raise CollisionError(f"comparison output already exists: {output.name}")
    a, b = _load_run(left), _load_run(right)
    ca = _flatten(a["plan"]["config"]) if a["plan"] else {}
    cb = _flatten(b["plan"]["config"]) if b["plan"] else {}
    config_diff = [{"field": k, "left": ca.get(k), "right": cb.get(k)} for k in sorted(set(ca) | set(cb)) if ca.get(k) != cb.get(k)]
    stages = {s["stage_id"]: s["state"] for s in a["status"]["stages"]}, {s["stage_id"]: s["state"] for s in b["status"]["stages"]}
    stage_diff = [{"stage_id": k, "left": stages[0].get(k, "absent"), "right": stages[1].get(k, "absent")} for k in sorted(set(stages[0]) | set(stages[1]))]
    ra, rb = _rows(left / "dea" / "zero_null.tsv"), _rows(right / "dea" / "zero_null.tsv")
    keys = sorted(set(ra) | set(rb))
    output.mkdir(parents=True)
    with (output / "dea_differences.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["model_id", "contrast_id", "hypothesis_type", "feature_id", "left_eligibility", "right_eligibility", "left_effect", "right_effect", "effect_difference", "left_p_value", "right_p_value", "left_q_value", "right_q_value"])
        for key in keys:
            l, r = ra.get(key, {}), rb.get(key, {})
            try:
                diff = repr(float(r["effect"]) - float(l["effect"]))
            except (KeyError, ValueError):
                diff = "NA"
            writer.writerow([*key, l.get("eligibility", "absent"), r.get("eligibility", "absent"), l.get("effect", "NA"), r.get("effect", "NA"), diff,
                             l.get("p_value", "NA"), r.get("p_value", "NA"), l.get("q_value", "NA"), r.get("q_value", "NA")])
    summary = {"result_type": "RunComparison", "left": {"run_id": a["status"]["run_id"], "plan_hash": a["status"]["plan_hash"], "state": a["status"]["state"]},
               "right": {"run_id": b["status"]["run_id"], "plan_hash": b["status"]["plan_hash"], "state": b["status"]["state"]},
               "same_plan": a["status"]["plan_hash"] == b["status"]["plan_hash"], "config_differences": config_diff, "stage_states": stage_diff,
               "n_endpoints_left": len(ra), "n_endpoints_right": len(rb), "n_shared_endpoints": len(set(ra) & set(rb)),
               "note": "Keyed by model, contrast, hypothesis and feature; differences are descriptive and do not rank either run."}
    (output / "comparison.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    rows = "".join(f"<tr><td>{html.escape(d['field'])}</td><td>{html.escape(str(d['left']))}</td><td>{html.escape(str(d['right']))}</td></tr>" for d in config_diff) or '<tr><td colspan="3">No configuration differences.</td></tr>'
    srows = "".join(f"<tr><td>{html.escape(s['stage_id'])}</td><td>{s['left']}</td><td>{s['right']}</td></tr>" for s in stage_diff)
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Run comparison</title><style>body{{font:14px sans-serif;margin:24px;}}table{{border-collapse:collapse;margin:8px 0}}'
            f'td,th{{border:1px solid #ccc;padding:3px 6px}}</style></head><body><h1>Run comparison</h1><p>Left {html.escape(summary["left"]["run_id"])} ({summary["left"]["state"]}); right '
            f'{html.escape(summary["right"]["run_id"])} ({summary["right"]["state"]}); same plan: {summary["same_plan"]}.</p>'
            f'<table><caption>Configuration differences</caption><thead><tr><th scope="col">field</th><th scope="col">left</th><th scope="col">right</th></tr></thead><tbody>{rows}</tbody></table>'
            f'<table><caption>Stage states</caption><thead><tr><th scope="col">stage</th><th scope="col">left</th><th scope="col">right</th></tr></thead><tbody>{srows}</tbody></table>'
            f'<p>Endpoint-level differences: <a href="dea_differences.tsv">dea_differences.tsv</a> ({len(keys)} keyed endpoints).</p></body></html>')
    (output / "index.html").write_text(page, encoding="utf-8")
    return summary
