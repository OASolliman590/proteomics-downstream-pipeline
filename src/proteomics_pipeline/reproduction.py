"""Environment locks, offline reproduction and semantic comparison (packet R11, V101/V102).

Locks are qualified against what is actually installed (exact versions, no
floating pins) and against the adapter APIs the stages call.  Reproduction
re-runs a recorded configuration offline into a fresh directory and compares
scientific tables keyed by stable identifiers with the frozen tolerances;
timestamps, run IDs and cosmetic metadata are separated from scientific values.
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import socket
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

from .errors import ProteomicsError

DETERMINISTIC = (1e-10, 1e-8)
BACKEND = (1e-8, 1e-6)
VOLATILE_COLUMNS = {"run_id", "plan_hash", "started_at", "finished_at", "engine_version"}
SCIENTIFIC_TABLES = {  # relative path -> key columns
    "inputs/matrix.tsv": ["feature_id"], "inputs/observed_mask.tsv": ["feature_id"], "inputs/observations.tsv": ["observation_id"],
    "preprocessing/primary/matrix.tsv": ["feature_id"], "preprocessing/coverage/coverage.tsv": ["model_id", "contrast_id", "feature_id"],
    "designs/estimability.tsv": ["model_id", "contrast_id", "feature_id"], "dea/zero_null.tsv": ["model_id", "contrast_id", "hypothesis_type", "feature_id"],
    "dea/treat.tsv": ["model_id", "contrast_id", "hypothesis_type", "feature_id"], "dea/families.tsv": ["family_id"],
    "permanova/tests.tsv": ["feature_set_id", "analysis", "comparison", "term"], "response/descriptive.tsv": ["axis_id", "feature_id"],
    "resources/gene_mapping.tsv": ["feature_id"], "resources/gene_matrix_finite.tsv": ["feature_id"],
}
API_PROBES = {
    "limma": "all(c('inter.gene.cor') %in% names(formals(limma::camera.default))) && all(c('midp','set.statistic','nrot') %in% names(formals(limma::mroast.default))) && 'robust' %in% names(formals(limma::eBayes))",
    "DEqMS": "all(c('fit','coef_col') %in% names(formals(DEqMS::spectraCounteBayes)))",
    "proDA": "all(c('contrast','reduced_model') %in% names(formals(proDA::test_diff)))",
    "vegan": "'by' %in% names(formals(vegan::adonis2))",
    "fgsea": "all(c('pathways','stats','minSize','maxSize') %in% names(formals(fgsea::fgsea)))",
    # A-2026-10-01-17 (R14d): the classifier APIs the biomarker stage calls
    "glmnet": "all(c('alpha','lambda','family','standardize') %in% names(formals(glmnet::glmnet)))",
    "e1071": "all(c('kernel','cost','scale') %in% names(formals(e1071:::svm.default)))",
}


def parse_python_lock(path: str | Path) -> dict[str, str]:
    pins = {}
    for number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"([A-Za-z0-9_.\-]+)==([A-Za-z0-9_.+\-]+)", line)
        if not match:
            raise ProteomicsError("E_LOCK_FLOATING", f"requirements.lock line {number} is not an exact pin: {line!r}", exit_code=3)
        pins[match.group(1).lower().replace("_", "-")] = match.group(2)
    return pins


def check_python_lock(lock: str | Path, python: str = sys.executable) -> list[dict]:
    pins = parse_python_lock(lock)
    code = "import json,importlib.metadata as m; print(json.dumps({d.metadata['Name'].lower().replace('_','-'): d.version for d in m.distributions()}))"
    installed = json.loads(subprocess.run([python, "-c", code], capture_output=True, text=True, check=True, encoding="utf-8").stdout)
    return [{"package": name, "locked": version, "installed": installed.get(name)} for name, version in pins.items() if installed.get(name) != version]


def check_r_lock(lock: str | Path, rscript: str | None = None) -> dict:
    rscript = rscript or os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript")
    if not rscript:
        raise ProteomicsError("E_CAPABILITY_NOT_AVAILABLE", "Rscript is unavailable; R lock qualification is NOT_RUN", exit_code=3)
    data = json.loads(Path(lock).read_text(encoding="utf-8"))
    floating = [name for name, p in data["Packages"].items() if not p.get("Version")]
    if floating:
        raise ProteomicsError("E_LOCK_FLOATING", f"renv.lock entries without versions: {floating}", exit_code=3)
    names = sorted(data["Packages"])
    code = ("p <- commandArgs(TRUE); v <- vapply(p, function(x) if (requireNamespace(x, quietly=TRUE)) utils::packageDescription(x)$Version else NA_character_, ''); "
            "probes <- list(" + ", ".join(f"{k} = tryCatch(isTRUE({v}), error = function(e) FALSE)" for k, v in API_PROBES.items()) + "); "
            "cat(jsonlite::toJSON(list(r = paste(R.version$major, R.version$minor, sep='.'), versions = as.list(v), probes = probes), auto_unbox = TRUE, null = 'null', na = 'null'))")
    from .runtime import run_r_code
    out = run_r_code(code, names, rscript=rscript)
    if out.returncode != 0:
        raise ProteomicsError("E_LOCK_QUALIFICATION", out.stderr[-400:], exit_code=3)
    observed = json.loads(out.stdout.strip().splitlines()[-1])
    mismatches = [{"package": n, "locked": data["Packages"][n]["Version"], "installed": observed["versions"].get(n)} for n in names
                  if observed["versions"].get(n) != data["Packages"][n]["Version"]]
    r_ok = observed["r"] == data["R"]["Version"]
    return {"r_locked": data["R"]["Version"], "r_installed": observed["r"], "r_compatible": r_ok, "bioconductor": data.get("Bioconductor", {}).get("Version"),
            "mismatches": mismatches, "api_probes": observed["probes"], "qualified": r_ok and not mismatches and all(observed["probes"].values())}


def _read(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _close(a: str, b: str, tol) -> bool:
    if a == b:
        return True
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False
    if math.isnan(x) or math.isnan(y):
        return False
    return abs(x - y) <= tol[0] or abs(x - y) <= tol[1] * max(abs(x), abs(y))


def semantic_compare(left: str | Path, right: str | Path) -> dict:
    left, right = Path(left), Path(right)
    differences = []
    compared = []
    for relative, keys in SCIENTIFIC_TABLES.items():
        a, b = left / relative, right / relative
        if not a.is_file() and not b.is_file():
            continue
        if a.is_file() != b.is_file():
            differences.append({"table": relative, "kind": "presence", "left": a.is_file(), "right": b.is_file()}); continue
        ra = {tuple(r[k] for k in keys): r for r in _read(a)}; rb = {tuple(r[k] for k in keys): r for r in _read(b)}
        compared.append(relative)
        if set(ra) != set(rb):
            differences.append({"table": relative, "kind": "keys", "only_left": len(set(ra) - set(rb)), "only_right": len(set(rb) - set(ra))})
        tol = BACKEND if relative.startswith(("dea/", "permanova/", "response/")) else DETERMINISTIC
        for key in set(ra) & set(rb):
            for column, value in ra[key].items():
                if column in VOLATILE_COLUMNS:
                    continue
                if not _close(value, rb[key].get(column, ""), tol):
                    differences.append({"table": relative, "kind": "value", "key": list(key), "column": column, "left": value, "right": rb[key].get(column)})
    return {"compared_tables": compared, "n_differences": len(differences), "differences": differences[:200], "volatile_columns_ignored": sorted(VOLATILE_COLUMNS)}


@contextmanager
def offline():
    """Refuse network connections in this process and point child proxies at a dead port."""
    original = socket.create_connection
    saved = {k: os.environ.get(k) for k in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY")}
    def refuse(*_a, **_k):
        raise ProteomicsError("E_NETWORK_FORBIDDEN", "analysis attempted a network connection", exit_code=5)
    socket.create_connection = refuse
    for k in saved:
        os.environ[k] = "http://127.0.0.1:9"
    try:
        yield
    finally:
        socket.create_connection = original
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def reproduce(run: str | Path, output: str | Path) -> dict:
    from .workflow import run_command, verify_run
    run = Path(run).resolve()
    verify_run(run)
    source = json.loads((run / "provenance" / "config_source.json").read_text(encoding="utf-8"))["config_path"]
    with offline():
        payload, code = run_command(source, output)
    original = json.loads((run / "plan.json").read_text(encoding="utf-8"))["plan_hash"] if (run / "plan.json").is_file() else None
    comparison = semantic_compare(run, output)
    return {"exit_code": code, "state": payload["state"], "same_plan_hash": payload["plan_hash"] == original, "semantic": comparison,
            "reproduced": code == 0 and payload["plan_hash"] == original and comparison["n_differences"] == 0}


# --------------------------------------------------------------------------- validation ledger (V110)
ROOT_DIR = Path(__file__).resolve().parents[2]


def validate_evidence(record: dict) -> None:
    from jsonschema import Draft202012Validator
    schema = json.loads((Path(__file__).parent / "schemas" / "validation-evidence.schema.json").read_text(encoding="utf-8"))
    errors = list(Draft202012Validator(schema).iter_errors(record))
    if errors:
        raise ProteomicsError("E_EVIDENCE_SCHEMA", f"{record.get('acceptance_id')}: {errors[0].message}", exit_code=5)
    for artifact in record["artifacts"]:
        if record["status"] == "PASS" and not (ROOT_DIR / artifact["path"]).is_file():
            raise ProteomicsError("E_EVIDENCE_SCHEMA", f"{record['acceptance_id']}: evidence artifact {artifact['path']} does not exist", exit_code=5)


def build_acceptance_ledger(root: str | Path = ROOT_DIR, extra: dict | None = None) -> list[dict]:
    """Machine-readable ledger: one record per acceptance ID from traceability and packet gate evidence."""
    from .provenance import sha256_file
    root = Path(root)
    trace = json.loads((root / "specs" / "001-downstream-proteomics" / "traceability.json").read_text(encoding="utf-8"))["requirements"]
    records = []
    for row in trace:
        status = row["status"] if row["status"] in ("PASS", "FAIL", "NOT_RUN", "SKIPPED", "INAPPLICABLE") else "NOT_RUN"
        verification = row.get("verification") or {}
        gate_file = next((Path(e) for e in row.get("evidence", []) if e.endswith("gate-results.json")), None)
        commands, artifacts = [], []
        if gate_file and (root / gate_file).is_file():
            gate = json.loads((root / gate_file).read_text(encoding="utf-8"))
            commands = [{"command": c["command"], "exit_code": c["exit_code"]} for c in gate["commands"]]
            artifacts = [{"path": gate_file.as_posix(), "sha256": sha256_file(root / gate_file)}] + [{"path": c["log"], "sha256": sha256_file(root / c["log"])} for c in gate["commands"] if (root / c["log"]).is_file()]
        elif status == "PASS" and row.get("evidence"):
            artifacts = [{"path": e, "sha256": sha256_file(root / e)} for e in row["evidence"] if (root / e).is_file()]
            commands = [{"command": "historical Maintainer receipt (see artifact)", "exit_code": 0}]
        external = [e for e in row.get("evidence", []) if "/evidence/" in e and not e.endswith("gate-results.json") and (root / e).is_file()]
        artifacts += [{"path": e, "sha256": sha256_file(root / e)} for e in external if e not in {a["path"] for a in artifacts}]
        record = {"acceptance_id": row["acceptance"], "requirement": row["requirement"], "task": row["task"], "status": status,
                  "reviewed_tree": verification.get("working_source_manifest_sha256") or row.get("verified_commit"),
                  "commands": commands, "versions": {}, "artifacts": artifacts, "expected_oracle": row["criterion"],
                  "observed": (verification.get("reason") or "see packet receipt and gate logs") if status == "FAIL" else ("see packet receipt and gate logs" if commands else "no execution recorded"),
                  "reviewer_conclusion": (f"{status} by {verification.get('route', 'Maintainer receipt')}; independent audit {verification.get('independent_audit', 'n/a')}") if status == "PASS" else status,
                  "verified_commit": row.get("verified_commit")}
        if status == "FAIL":
            # Amendment A-2026-10-01-12: a FAIL record states why it failed (V110 requires a reason on every non-PASS record).
            record["reason"] = verification.get("reason") or "failed; see observed and artifacts"
        elif status != "PASS":
            record["reason"] = (extra or {}).get(row["acceptance"], "not executed in this working tree")
        validate_evidence(record)
        records.append(record)
    return records


# --------------------------------------------------------------------------- golden reference matrix (V104)
def check_reference_matrix(matrix_path: str | Path, lock_path: str | Path, root: str | Path = ROOT_DIR) -> dict:
    """Verify the frozen matrix against the lock and the test tree; differences are recorded, never hidden."""
    root = Path(root)
    matrix = json.loads(Path(matrix_path).read_text(encoding="utf-8"))
    lock = json.loads(Path(lock_path).read_text(encoding="utf-8"))["Packages"]
    findings = []
    tiers = matrix["tolerances"]
    for case in matrix["cases"]:
        if case["tolerance"] not in tiers:
            findings.append({"case_id": case["case_id"], "field": "tolerance", "frozen": case["tolerance"], "observed": "undefined tier"})
        for package, version in case["packages"].items():
            locked = (lock.get(package) or {}).get("Version")
            if locked != version:
                findings.append({"case_id": case["case_id"], "field": f"packages/{package}", "frozen": version, "observed": locked})
        path, _, token = case["oracle_test"].partition("::")
        source = root / path
        text = source.read_text(encoding="utf-8") if source.is_file() else ""
        present = (f"def {token}(" in text) if path.endswith(".py") else (f'test_that("{token}' in text)
        if not present:
            findings.append({"case_id": case["case_id"], "field": "oracle_test", "frozen": case["oracle_test"], "observed": "missing"})
    sources = "\n".join(p.read_text(encoding="utf-8") for p in list((root / "r" / "proteomicsCore" / "R").glob("*.R")) + list((root / "src" / "proteomics_pipeline").rglob("*.py")))
    for item in matrix["unsupported_combinations"]:
        if item["reason_code"] not in sources:
            findings.append({"case_id": item["combination"], "field": "reason_code", "frozen": item["reason_code"], "observed": "not raised anywhere"})
    return {"n_cases": len(matrix["cases"]), "n_unsupported": len(matrix["unsupported_combinations"]), "differences": findings, "consistent": not findings}


def write_acceptance_ledger(path: str | Path = ROOT_DIR / "docs" / "validation" / "acceptance.json", root: str | Path = ROOT_DIR) -> dict:
    root = Path(root)
    trace = json.loads((root / "specs" / "001-downstream-proteomics" / "traceability.json").read_text(encoding="utf-8"))["requirements"]
    # A-2026-10-01-21: an R12 row keeps its own recorded reason (Maintainer-only gate, operator decision or final review); the generic text is a fallback.
    reasons = {r["acceptance"]: (r.get("verification") or {}).get("reason") or "R12 private regression/release gate: Maintainer-only, NOT_RUN by an implementer"
               for r in trace if r["packet"] == "R12"}
    for r in trace:
        if r["acceptance"] not in reasons and r["status"] not in ("PASS", "FAIL"):
            reasons[r["acceptance"]] = (r.get("verification") or {}).get("reason") or "not executed or not passed in this working tree"
    records = build_acceptance_ledger(root, reasons)
    counts: dict[str, int] = {}
    for record in records:
        counts[record["status"]] = counts.get(record["status"], 0) + 1
    ledger = {"schema": "src/proteomics_pipeline/schemas/validation-evidence.schema.json", "verified_commit": None, "status_counts": counts,
              "note": "Generated from traceability.json and packet gate evidence; PASS/FAIL/NOT_RUN/SKIPPED/INAPPLICABLE stay distinct. Uncommitted Claude-route receipts record the working-source manifest hash as the reviewed tree.",
              "records": records}
    Path(path).write_text(json.dumps(ledger, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    return ledger
