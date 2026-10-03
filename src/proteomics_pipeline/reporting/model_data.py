"""Alternative-engine section of ReportData (packet R06).

Shows each declared DEqMS/proDA model with its actual state, separate family
and the stable-ID comparison with the frozen primary model; nothing here can
promote an alternative engine to primary.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def model_section(directory: str | Path, *, prefix: str = "assay_engines") -> dict:
    root = Path(directory)
    status = json.loads((root / "model_status.json").read_text(encoding="utf-8"))
    tables = [{"caption": "Alternative-engine models and their actual states", "source": f"{prefix}/model_status.json", "columns": ["model", "engine", "state", "reason"],
               "rows": [[m["model_id"], m["engine"], m["state"], m.get("reason_code") or "none"] for m in status["models"]]}]
    if (root / "families.tsv").is_file():
        fam = _table(root / "families.tsv")
        tables.append({"caption": "Alternative-engine families (separate from the primary family)", "source": f"{prefix}/families.tsv",
                       "columns": ["family", "planned", "finite", "completeness", "rejections"], "rows": [[f["family_id"], f["n_planned"], f["n_finite"], f["completeness"], f["rejection_count"]] for f in fam]})
    if (root / "engine_comparison.tsv").is_file():
        comp = _table(root / "engine_comparison.tsv")
        matched = sum(1 for r in comp if r["in_alternative_universe"] == "true" and r["in_primary_universe"] == "true")
        tables.append({"caption": "Stable-ID comparison with the primary model (sensitivity, not primary)", "source": f"{prefix}/engine_comparison.tsv",
                       "columns": ["matched features", "alternative only", "primary only"],
                       "rows": [[matched, sum(1 for r in comp if r["in_alternative_universe"] == "true" and r["in_primary_universe"] != "true"),
                                 sum(1 for r in comp if r["in_primary_universe"] == "true" and r["in_alternative_universe"] != "true")]]})
    return {"state": "COMPLETED", "title": "Alternative engines (DEqMS / proDA)", "tables": tables, "figures": [], "summary": [status["note"]], "notes": [],
            "methods": [f"{k}: {v.get('engine')} {v.get('version')}" for k, v in status.get("diagnostics", {}).items()], "values": {"n_models": len(status["models"])}}
