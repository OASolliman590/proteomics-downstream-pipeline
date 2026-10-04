"""R12 V116: spec-task-evidence reconciliation (FR-116, SM25).

The oracle is independent of traceability.json: identities are parsed from every slice's spec.md (FR definitions and
V anchors) and tasks.md (T lines), then compared by set, cardinality and key with the traceability rows, the
committed acceptance ledger and the evidence on disk. The negative cases apply the same checks to mutated copies.
"""
from __future__ import annotations

import copy
import json
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPECS = ROOT / "specs"
TRACE = SPECS / "001-downstream-proteomics" / "traceability.json"
LEDGER = ROOT / "docs" / "validation" / "acceptance.json"
STATUSES = {"PASS", "FAIL", "NOT_RUN", "SKIPPED", "INAPPLICABLE"}
CORE = range(1, 121)                       # FR-001-FR-120 / T001-T120 / V001-V120 (frozen kit)
AMENDMENTS = {"014-multivariate-permanova": range(121, 131), "015-post-de-analysis": range(131, 168)}   # ADR 0006, ADR 0009


def spec_identities() -> dict[str, dict[str, list[str]]]:
    out = {}
    for spec in sorted(SPECS.glob("0[0-9][0-9]-*/spec.md")):
        if spec.parent.name == "001-downstream-proteomics":
            continue
        text = spec.read_text(encoding="utf-8")
        tasks = (spec.parent / "tasks.md").read_text(encoding="utf-8")
        out[spec.parent.name] = {"FR": re.findall(r"^- \*\*(FR-\d{3}) —", text, re.M), "V": re.findall(r'^<a id="(V\d{3})"></a>', text, re.M),
                                 "T": re.findall(r"^- \[[ x]\] (T\d{3}) —", tasks, re.M), "done": re.findall(r"^- \[x\] (T\d{3}) —", tasks, re.M)}
    return out


def reconcile(trace: dict, ledger: dict, specs: dict, root: Path = ROOT) -> list[str]:
    """Every problem found; empty means reconciled."""
    problems = []
    rows = trace["requirements"]
    triples = [(r["requirement"], r["task"], r["acceptance"]) for r in rows]
    if len(set(triples)) != len(triples):
        problems.append("duplicate traceability rows")
    for kind, index in (("FR", 0), ("T", 1), ("V", 2)):
        ids = [t[index] for t in triples]
        if len(set(ids)) != len(ids):
            problems.append(f"duplicate {kind} identity")
    spec_fr = [i for s in specs.values() for i in s["FR"]]; spec_t = [i for s in specs.values() for i in s["T"]]; spec_v = [i for s in specs.values() for i in s["V"]]
    for kind, found in (("FR", spec_fr), ("T", spec_t), ("V", spec_v)):
        if len(found) != len(set(found)):
            problems.append(f"an identity is defined twice in the slice specs ({kind})")
    expected = {n for n in CORE} | {n for r in AMENDMENTS.values() for n in r}
    fmt = {"FR": "FR-{:03d}", "T": "T{:03d}", "V": "V{:03d}"}
    for kind, found in (("FR", spec_fr), ("T", spec_t), ("V", spec_v)):
        if set(found) != {fmt[kind].format(n) for n in expected}:
            problems.append(f"slice-spec {kind} identities differ from 1-120 plus the recorded amendments")
    core = {fmt["V"].format(n) for n in CORE}
    if len(core & set(spec_v)) != 120:
        problems.append("the frozen core does not hold exactly 120 acceptance identities")
    if {t[2] for t in triples} != set(spec_v) or {t[0] for t in triples} != set(spec_fr) or {t[1] for t in triples} != set(spec_t):
        problems.append("traceability identities differ from the slice specs (rename/merge/split without redirect)")
    for slice_name, numbers in AMENDMENTS.items():
        if set(specs.get(slice_name, {}).get("V", [])) != {fmt["V"].format(n) for n in numbers}:
            problems.append(f"amendment slice {slice_name} does not hold its recorded identities")
    for r in rows:
        if f"V{r['requirement'][3:]}" != r["acceptance"] or f"T{r['requirement'][3:]}" != r["task"]:
            problems.append(f"{r['acceptance']}: FR/T/V numbers do not align")
        status = r["status"] if r["status"] in STATUSES else "NOT_RUN"
        done = r["task"] in specs.get(r["slice"], {}).get("done", [])
        if (status == "PASS") != done:
            problems.append(f"{r['task']}: tasks.md checkbox disagrees with status {r['status']}")
        if status == "PASS":
            if not r["evidence"]:
                problems.append(f"{r['acceptance']}: PASS without evidence")
            for e in r["evidence"]:
                if "pending" in e or not (root / e).exists():
                    problems.append(f"{r['acceptance']}: evidence {e!r} is pending or missing")
            if r.get("verified_commit") is not None and not re.fullmatch(r"[0-9a-f]{40}", r["verified_commit"]):
                problems.append(f"{r['acceptance']}: verified_commit is not a commit id")
        else:
            if any("pending" not in e and not (root / e).exists() for e in r.get("evidence", [])):
                problems.append(f"{r['acceptance']}: missing evidence path")
            if not (r.get("verification") or {}).get("reason"):
                problems.append(f"{r['acceptance']}: unresolved status without an explicit reason")
    for redirect in trace.get("redirects", []):
        if not redirect.get("adr") or not (root / redirect["adr"]).is_file():
            problems.append(f"redirect {redirect} without an ADR")
    records = {r["acceptance_id"]: r for r in ledger["records"]}
    if list(records) != [t[2] for t in triples]:
        problems.append("ledger identities differ from traceability")
    for r in rows:
        record = records.get(r["acceptance"])
        status = r["status"] if r["status"] in STATUSES else "NOT_RUN"
        if record is None or record["status"] != status:
            problems.append(f"{r['acceptance']}: ledger status differs")
        elif status not in ("PASS", "FAIL") and record.get("reason") != (r.get("verification") or {}).get("reason"):
            problems.append(f"{r['acceptance']}: ledger reason differs from the recorded reason")
    if ledger.get("verified_commit") is not None:
        problems.append("ledger claims a verified commit")
    return problems


def load():
    return json.loads(TRACE.read_text(encoding="utf-8")), json.loads(LEDGER.read_text(encoding="utf-8")), spec_identities()


def test_v116_identities_tasks_evidence_and_ledger_reconcile():
    trace, ledger, specs = load()
    problems = reconcile(trace, ledger, specs)
    assert problems == [], "\n".join(problems)
    assert len(trace["requirements"]) == 167 and sum(len(s["V"]) for s in specs.values()) == 167


def test_v116_verified_commits_exist_in_history():
    trace, _, _ = load()
    for commit in {r["verified_commit"] for r in trace["requirements"] if r.get("verified_commit")}:
        assert subprocess.run(["git", "cat-file", "-e", f"{commit}^{{commit}}"], cwd=ROOT).returncode == 0, commit


@pytest.mark.parametrize("mutate, expected", [
    (lambda t, l, s: t["requirements"][0]["evidence"].__setitem__(0, "pending-after-freeze"), "pending or missing"),
    (lambda t, l, s: t["requirements"][5].update(acceptance="V999"), "rename/merge/split"),
    (lambda t, l, s: t["requirements"].append(copy.deepcopy(t["requirements"][3])), "duplicate"),
    (lambda t, l, s: s["003-intake"]["V"].remove("V011"), "slice-spec V identities"),
    (lambda t, l, s: next(r for r in t["requirements"] if r["status"] != "PASS").update(verification={}), "explicit reason"),
    (lambda t, l, s: s["002-foundation"]["done"].remove("T001"), "checkbox"),
    (lambda t, l, s: t["redirects"].append({"from": "V010", "to": "V011"}), "without an ADR"),
    (lambda t, l, s: next(r for r in l["records"] if r["status"] == "NOT_RUN").update(reason="generic"), "ledger reason"),
])
def test_v116_negative_mutations_are_detected(mutate, expected):
    trace, ledger, specs = load()
    mutate(trace, ledger, specs)
    problems = reconcile(trace, ledger, specs)
    assert any(expected in p for p in problems), problems
