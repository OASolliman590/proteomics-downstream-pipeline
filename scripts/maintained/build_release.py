#!/usr/bin/env python3
"""Release hygiene check and versioned successor export (packet R12, FR-117/FR-118, V117/V118).

`check` inventories the files a source release would contain (the tracked files of the checkout, never untracked
environments, caches, runs or private evidence) and reports every blocking finding:
  - a tracked file of a forbidden kind (private workbooks and data archives, R session files, environments,
    credentials files, private execution logs);
  - content that looks like a credential (private keys, provider tokens, hard-coded passwords/API keys);
  - an absolute private path (a user home or mounted volume) in text content;
  - a license file while no license has been selected by the operator (no license decision is recorded here:
    choosing a license is an operator decision, AGENTS.md);
  - a protected recovered-baseline file whose hash differs from the trusted registry.
The report always records ownership and license as unresolved unless an operator decision record says otherwise.

`export` writes a *new* versioned successor directory only when every gate holds: no blocking finding, an explicit
release-authorization record for that version, a destination that does not exist and lies outside the source
tree, and, for a public release, a recorded license decision and disclosure review. It never overwrites anything
and re-verifies the protected baseline after copying.

`codes` prints the typed reason-code catalogue (docs/methods/reason-codes.md) extracted from the maintained Python
and R sources, so the documented rejection messages are exactly the ones the code can raise.

Usage:
  build_release.py check [--root DIR] [--json]
  build_release.py codes [--root DIR] > docs/methods/reason-codes.md
  build_release.py export --version V --destination DIR --authorization FILE [--root DIR]
Exit codes: 0 ok, 2 refused (blocking findings or a failed gate), 5 integrity failure.
"""
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REGISTRY = "docs/validation/baseline-hashes.json"
PACKAGE_NAME = "proteomics-downstream-pipeline"
FORBIDDEN_PATTERNS = {   # tracked files of these kinds block a release (the Maintainer removes them from the source)
    "private workbook or spreadsheet": ["*.xlsx", "*.xlsm", "*.xls", "*.ods"],
    "data archive or R object": ["*.zip", "*.rds", "*.RData", "*.rda", "*.tif", "*.tiff"],
    "R session file": [".Rhistory", "*/.Rhistory"],
    "environment or cache": [".venv/*", "*/.venv/*", ".r-lib/*", "renv/library/*", "*/__pycache__/*", "*.pyc", "*.egg-info/*"],
    "credentials file": [".env", ".env.*", "*/.env", "*.pem", "*.key", "*.p12", "id_rsa", "id_ed25519", "*/id_rsa", ".netrc", "*credentials*.json"],
    "private execution output": ["runs/*", "results/*", "data/*", "private_evidence/*", "reproducibility/*", "logs/*"],
}
ALLOWED_EXCEPTIONS = [".env.example"]
SECRET_PATTERNS = {
    "private key block": re.compile(r"-----BEGIN (?:RSA |EC |DSA |OPENSSH |ENCRYPTED )?PRIVATE KEY-----"),
    "GitHub token": re.compile(r"\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{36}\b|\bgithub_pat_[A-Za-z0-9_]{40,}\b"),
    "AWS access key": re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
    "Slack token": re.compile(r"\bxox[abprs]-[A-Za-z0-9-]{10,}\b"),
    "hard-coded secret": re.compile(r"(?i)\b(?:password|passwd|api[_-]?key|secret[_-]?key|access[_-]?token)\b\s*[:=]\s*['\"][^'\"\s]{8,}['\"]"),
}
# A user home or mounted volume; the public CI runner homes (/home/runner, /Users/runner, C:\Users\runneradmin) are not private.
PRIVATE_PATH = re.compile(r"(?<![A-Za-z0-9_<])(?:/Users/(?!<|runner/)[A-Za-z0-9._-]+/|/home/(?!<|runner/)[A-Za-z0-9._-]+/|/Volumes/[A-Za-z0-9._ -]+/|[A-Za-z]:\\Users\\(?!runneradmin\\)[A-Za-z0-9._-]+\\)")
LICENSE_NAMES = re.compile(r"^(?:LICEN[CS]E|COPYING|UNLICENSE)(?:\.[A-Za-z0-9]+)?$", re.IGNORECASE)
TEXT_LIMIT = 5 * 1024 * 1024
DECISION_RECORD = "docs/validation/013-release/operator-decisions.json"   # absent until the operator decides


class ReleaseError(Exception):
    def __init__(self, code: str, message: str, exit_code: int = 2):
        super().__init__(message); self.code = code; self.message = message; self.exit_code = exit_code


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tracked_files(root: Path) -> list[str]:
    """The files a source release would contain: what git tracks in this checkout (never untracked files)."""
    try:
        out = subprocess.run(["git", "ls-files", "-z"], cwd=root, capture_output=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        raise ReleaseError("E_RELEASE_INVENTORY", f"cannot list tracked files with git: {error}") from error
    return sorted(p for p in out.decode("utf-8").split("\0") if p and (root / p).is_file())


def _forbidden_kind(path: str) -> str | None:
    if path in ALLOWED_EXCEPTIONS or path.endswith("/.env.example"):
        return None
    for kind, patterns in FORBIDDEN_PATTERNS.items():
        if any(fnmatch.fnmatchcase(path, p) or fnmatch.fnmatchcase(Path(path).name, p) for p in patterns if "/" not in p) or \
           any(fnmatch.fnmatchcase(path, p) for p in patterns if "/" in p):
            return kind
    return None


def _text(path: Path) -> str | None:
    if path.stat().st_size > TEXT_LIMIT:
        return None
    data = path.read_bytes()
    if b"\0" in data[:8192]:
        return None
    return data.decode("utf-8", errors="replace")


def operator_decisions(root: Path) -> dict:
    path = root / DECISION_RECORD
    if not path.is_file():
        return {"license": None, "ownership": None, "disclosure_review": None}
    record = json.loads(path.read_text(encoding="utf-8"))
    return {"license": record.get("license"), "ownership": record.get("ownership"), "disclosure_review": record.get("disclosure_review")}


def protected_baseline(root: Path) -> list[dict]:
    registry = root / REGISTRY
    if not registry.is_file():
        return [{"path": REGISTRY, "reason": "registry missing"}]
    problems = []
    for item in json.loads(registry.read_text(encoding="utf-8")).get("files", []):
        path = root / item["path"]
        if not path.is_file():
            problems.append({"path": item["path"], "reason": "missing"})
        elif sha256_file(path) != item["sha256"]:
            problems.append({"path": item["path"], "reason": "hash_mismatch"})
    return problems


def check(root: Path = ROOT) -> dict:
    root = Path(root).resolve()
    files = tracked_files(root)
    decisions = operator_decisions(root)
    findings = []
    for rel in files:
        kind = _forbidden_kind(rel)
        if kind:
            findings.append({"path": rel, "kind": kind, "rule": "forbidden tracked file"})
            continue
        if LICENSE_NAMES.match(Path(rel).name) and not decisions["license"]:
            findings.append({"path": rel, "kind": "license file without an operator license decision", "rule": "license"})
        text = _text(root / rel)
        if text is None:
            continue
        for name, pattern in SECRET_PATTERNS.items():
            match = pattern.search(text)
            if match:
                findings.append({"path": rel, "kind": f"credential-like content ({name})", "rule": "secret", "line": text.count("\n", 0, match.start()) + 1})
        match = PRIVATE_PATH.search(text)
        if match:
            findings.append({"path": rel, "kind": "absolute private path", "rule": "private_path", "line": text.count("\n", 0, match.start()) + 1})
    baseline = protected_baseline(root)
    for problem in baseline:
        findings.append({"path": problem["path"], "kind": f"protected baseline {problem['reason']}", "rule": "baseline"})
    blocking = bool(findings)
    return {"result_type": "ReleaseHygiene", "root_name": root.name, "n_tracked_files": len(files), "findings": findings, "blocking": blocking,
            "license": decisions["license"] or "unresolved (operator decision required, AGENTS.md)",
            "ownership": decisions["ownership"] or "unresolved (operator decision required, AGENTS.md)",
            "disclosure_review": decisions["disclosure_review"] or "not recorded",
            "public_release_allowed": (not blocking) and bool(decisions["license"]) and bool(decisions["disclosure_review"]),
            "untracked_files_included": False}


def _authorization(path: Path | None, version: str) -> dict:
    if path is None or not Path(path).is_file():
        raise ReleaseError("E_RELEASE_NOT_AUTHORIZED", "export needs an explicit release-authorization record for this version (operator decision)")
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    if record.get("release_authorized") is not True or record.get("version") != version or not record.get("authorized_by"):
        raise ReleaseError("E_RELEASE_NOT_AUTHORIZED", f"the authorization record does not authorize version {version!r} (needs release_authorized=true, version and authorized_by)")
    if record.get("scope") not in ("private_successor", "public"):
        raise ReleaseError("E_RELEASE_NOT_AUTHORIZED", "authorization scope must be 'private_successor' or 'public'")
    return record


def export(root: Path, version: str, destination: Path, authorization: Path | None) -> dict:
    root = Path(root).resolve()
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+(?:[-+][A-Za-z0-9.]+)?", version):
        raise ReleaseError("E_RELEASE_VERSION", f"version {version!r} is not a semantic version")
    record = _authorization(authorization, version)
    report = check(root)
    if report["blocking"]:
        raise ReleaseError("E_RELEASE_HYGIENE", f"{len(report['findings'])} blocking hygiene finding(s); see `check`")
    if record["scope"] == "public":
        decisions = operator_decisions(root)
        if not decisions["license"] or not decisions["disclosure_review"] or not record.get("license_decision"):
            raise ReleaseError("E_LICENSE_UNRESOLVED", "a public release needs a recorded operator license decision and disclosure review")
    destination = Path(destination).resolve()
    if destination.exists():
        raise ReleaseError("E_PATH_COLLISION", f"release destination already exists: {destination.name}")
    if destination == root or root in destination.parents:
        raise ReleaseError("E_RELEASE_DESTINATION", "the successor directory must lie outside the source tree")
    before = {p["path"] for p in protected_baseline(root)}
    target = destination / f"{PACKAGE_NAME}-{version}"
    staging = destination.parent / f".{destination.name}.staging"
    if staging.exists():
        raise ReleaseError("E_PATH_COLLISION", f"staging directory already exists: {staging.name}")
    files = tracked_files(root)
    try:
        for rel in files:
            out = staging / target.name / rel
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / rel, out)
        manifest = [f"{sha256_file(staging / target.name / rel)}  {rel}" for rel in files]
        (staging / target.name / "MANIFEST.sha256").write_text("".join(line + "\n" for line in manifest), encoding="utf-8")
        release = {"result_type": "ReleaseRecord", "package": PACKAGE_NAME, "version": version, "scope": record["scope"], "authorized_by": record["authorized_by"],
                   "license": report["license"], "ownership": report["ownership"], "n_files": len(files),
                   "manifest_sha256": sha256_file(staging / target.name / "MANIFEST.sha256"), "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")}
        (staging / target.name / "release.json").write_text(json.dumps(release, indent=2) + "\n", encoding="utf-8")
        os.replace(staging, destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    after = {p["path"] for p in protected_baseline(root)}
    if after != before:
        raise ReleaseError("E_BASELINE_CHANGED", "the protected baseline changed during export", exit_code=5)
    return release | {"path": str(target)}


CODE = re.compile(r"\b([EW]_[A-Z][A-Z0-9_]+)\b")
CODE_SOURCES = ("src/proteomics_pipeline/**/*.py", "r/proteomicsCore/R/*.R")


def reason_codes(root: Path = ROOT) -> dict[str, dict]:
    """Every typed reason code in the maintained sources, with the files that use it and one example message."""
    root = Path(root)
    codes: dict[str, dict] = {}
    for pattern in CODE_SOURCES:
        for path in sorted(root.glob(pattern)):
            rel = path.relative_to(root).as_posix()
            text = path.read_text(encoding="utf-8")
            for match in CODE.finditer(text):
                entry = codes.setdefault(match.group(1), {"files": set(), "example": None})
                entry["files"].add(rel)
                if entry["example"] is None:
                    tail = text[match.end():match.end() + 240]
                    found = (re.match(r"""["']\s*,\s*(?:reason\s*=\s*)?f?(["'])(.{8,220}?)\1""", tail) or (rel.endswith(".R") and re.match(r""":\s*(?P<q>)([^"\n]{8,220})""", tail)))
                    if found:
                        entry["example"] = found.group(2).strip()
    return {k: codes[k] for k in sorted(codes)}


def codes_markdown(root: Path = ROOT) -> str:
    lines = ["# Typed reason codes", "",
             "Generated by `python scripts/maintained/build_release.py codes` from the maintained Python and R sources (R12, V115). `E_` codes are refusals, failures or NOT_RUN reasons; `W_` codes are warnings that leave the stage COMPLETED. A test checks that this list and the sources agree in both directions. Messages are examples: `{...}` and `%s` mark values filled in at run time.",
             "", "| Code | Used in | Example message |", "|---|---|---|"]
    for code, entry in reason_codes(root).items():
        files = ", ".join(f"`{f.split('/')[-1]}`" for f in sorted(entry["files"]))
        example = (entry["example"] or "—").replace("|", "\\|").replace("`", "'")
        lines.append(f"| `{code}` | {files} | {example} |")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    c = sub.add_parser("check"); c.add_argument("--root", type=Path, default=ROOT); c.add_argument("--json", action="store_true")
    e = sub.add_parser("export"); e.add_argument("--root", type=Path, default=ROOT); e.add_argument("--version", required=True)
    e.add_argument("--destination", type=Path, required=True); e.add_argument("--authorization", type=Path)
    k = sub.add_parser("codes"); k.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    if args.command == "codes":
        sys.stdout.write(codes_markdown(args.root)); return 0
    try:
        if args.command == "check":
            report = check(args.root)
            print(json.dumps(report, indent=2))
            return 2 if report["blocking"] else 0
        print(json.dumps(export(args.root, args.version, args.destination, args.authorization), indent=2))
        return 0
    except ReleaseError as error:
        print(json.dumps({"state": "REFUSED", "error": {"code": error.code, "message": error.message}}))
        return error.exit_code


if __name__ == "__main__":
    sys.exit(main())
