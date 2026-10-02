"""Treatment-response section of ReportData (packet R09, SM19-SM24).

Descriptive movement, equivalence, formal conjunction and independent score
tests are kept in separate tables; descriptive classes never carry P values
or equivalence language, and no causal rescue wording is generated.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _f(v):
    try:
        return f"{float(v):.3g}"
    except (TypeError, ValueError):
        return v


def response_section(directory: str | Path, *, prefix: str = "response") -> dict:
    root = Path(directory)
    result = json.loads((root / "response_result.json").read_text(encoding="utf-8"))
    tables = []
    desc = _table(root / "descriptive.tsv") if (root / "descriptive.tsv").is_file() else []
    if desc:
        tables.append({"caption": "Descriptive response axes (no P values; classes are descriptive)", "source": f"{prefix}/descriptive.tsv",
                       "columns": ["axis", "feature", "d", "t", "r", "Cov(d,t)", "RI", "class", "crossed control", "ratio interval"],
                       "rows": [[r["axis_id"], r["feature_id"], _f(r["d"]), _f(r["t"]), _f(r["r"]), _f(r["cov_dt"]), _f(r["reversal_index"]), r["descriptive_class"],
                                 r["crossed_control"], r.get("ratio_interval_kind", "")] for r in desc]})
    for name, caption, cols in (("equivalence", "Residual equivalence (TOST; separate family)", ["feature_id", "residual_estimate", "margin", "p_lower", "p_upper", "p_value", "q_value"]),
                                ("formal_rescue", "Formal rescue conjunction (intersection-union; independent directions)", ["feature_id", "direction", "p_disease", "p_treatment", "p_value", "q_value"]),
                                ("independent_score_tests", "Independent score randomization tests", ["score_id", "method", "statistic", "k", "N", "B", "p_value", "q_value"]),
                                ("descriptive_scores", "Descriptive scores (no inference: selection overlap or unverified independence)", ["score_id", "observation_id", "group", "value", "score_eligibility_reason"])):
        path = root / f"{name}.tsv"
        if path.is_file():
            rows = _table(path)
            tables.append({"caption": caption, "source": f"{prefix}/{name}.tsv", "columns": cols, "rows": [[_f(r.get(c)) for c in cols] for r in rows]})
    return {"state": "COMPLETED", "title": "Treatment response (descriptive axes, equivalence, conjunction, scores)", "tables": tables, "figures": [],
            "summary": [f"Response mode: {result['response_mode']}; dmin = {result['dmin']}; ratio uncertainty: {result['ratio_uncertainty']}.", result["vocabulary"]],
            "notes": ["Opposition between disease and treatment contrasts that share a control group is descriptive; the estimates are negatively correlated, not independent confirmation.",
                      "A near_restoration class is not statistical equivalence; only a passing EquivalenceResult supports equivalence within the declared margin."],
            "methods": ["Same-model axes d, t, r = d + t with exact per-feature covariance; TOST on r with model-correct SE/df; conjunction as the maximum of one-sided component P values; "
                        "fixed independent scores tested by exact or Monte Carlo randomization only when training participants are disjoint from tested units."],
            "values": {"axes": list(result.get("axes", {}))}}
