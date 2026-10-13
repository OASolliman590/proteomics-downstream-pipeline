"""Design-specific acceptance checks of a study run (packet R15d: FR-186, FR-189 to FR-193; ADR 0011).

Each check reads only the run's own outputs and reports its own state. A design with no tested endpoint is reported
as FAIL for the check, never as a healthy zero (SM25)."""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def acceptance(run_root: str | Path, design: str) -> dict:
    root = Path(run_root)
    table = root / "dea" / "zero_null.tsv"
    if not table.is_file():
        return {"check": f"{design}_endpoints", "state": "FAIL", "reason": "no differential table was produced", "tested": 0}
    rows = _rows(table)
    tested = sum(1 for r in rows if r["eligibility"] == "tested")
    state = "PASS" if tested > 0 else "FAIL"
    detail = {"check": f"{design}_endpoints", "state": state, "tested": tested, "rows": len(rows)}
    if design in ("paired", "repeated"):
        plan = json.loads((root / "plan.json").read_text(encoding="utf-8"))
        blocked = any((b or {}).get("mode", "none") != "none" for b in plan.get("blocking", []))
        detail["subject_structure_preserved"] = blocked
        detail["state"] = "PASS" if state == "PASS" and blocked else "FAIL"
    return detail
