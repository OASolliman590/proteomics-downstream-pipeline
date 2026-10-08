"""R12 V117 (source release hygiene) and the V118 export tooling (FR-117, FR-118; SM14, SM25).

The fixture is a scratch git tree with deliberately seeded dummy secrets, private paths, a private workbook,
untracked environments and a fake license. The oracle is an independent inventory written here (os.walk plus the
tree's own `git add` list) and an independent SHA-256 manifest of the exported directory. V118 itself stays
NOT_RUN: no real release authorization exists, so only the export tooling and its refusals are exercised.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("build_release", ROOT / "scripts" / "maintained" / "build_release.py")
BR = importlib.util.module_from_spec(spec); spec.loader.exec_module(BR)
pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="NOT_RUN: git unavailable")

# Dummy credential-like strings, assembled at run time so this test file itself never contains one.
DUMMY = {"key": "-----BEGIN " + "PRIVATE KEY-----\nMIIEvQ" + "A" * 40 + "\n-----END PRIVATE KEY-----\n",
         "github": "token = 'gh" + "p_" + "a1B2c3D4e5F6g7H8i9J0k1L2m3N4o5P6q7R8'\n",
         "aws": "aws_key = AK" + "IA" + "ABCDEFGHIJKLMNOP\n",
         "password": "password = '" + "hunter2hunter2" + "'\n"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scratch_tree(tmp_path: Path, extra: dict[str, str | bytes] | None = None, *, track_extra: bool = True) -> Path:
    root = tmp_path / "src-tree"
    files = {"README.md": "# scratch\n", "src/pkg/__init__.py": "VALUE = 1\n", "docs/usage.md": "Run `tool --help`.\n", ".env.example": "TOKEN=\n",
             "pipeline/old_script.R": "x <- 1\n"}
    for rel, content in files.items():
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text(content, encoding="utf-8")
    registry = {"schema_version": "1.0.0", "protected_roots": ["pipeline"], "files": [{"path": "pipeline/old_script.R", "sha256": sha(root / "pipeline/old_script.R")}]}
    (root / "docs/validation").mkdir(parents=True)
    (root / BR.REGISTRY).write_text(json.dumps(registry), encoding="utf-8")
    # untracked material that must never enter a release
    for rel in (".venv/lib/site.py", ".r-lib/pkg/DESCRIPTION", "runs/r1/report.html", "private_evidence/study.xlsx", "data/cohort.csv"):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text("private\n", encoding="utf-8")
    subprocess.run(["git", "init", "-q"], cwd=root, check=True)
    subprocess.run(["git", "add", "README.md", "src", "docs", ".env.example", "pipeline"], cwd=root, check=True)
    for rel, content in (extra or {}).items():
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_bytes(content if isinstance(content, bytes) else content.encode("utf-8"))
        if track_extra:
            subprocess.run(["git", "add", "-f", rel], cwd=root, check=True)
    return root


def oracle_inventory(root: Path) -> set[str]:
    """Independent: every file on disk minus the ones git does not track (git's own `ls-files --others`)."""
    on_disk = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() and ".git" not in p.relative_to(root).parts}
    others = set(subprocess.run(["git", "ls-files", "--others", "-z"], cwd=root, capture_output=True, check=True).stdout.decode().split("\0")) - {""}
    return on_disk - others


def test_v117_clean_tree_inventory_excludes_untracked_environments_and_private_files(tmp_path):
    root = scratch_tree(tmp_path)
    report = BR.check(root)
    assert report["blocking"] is False and report["findings"] == []
    assert set(BR.tracked_files(root)) == oracle_inventory(root) and report["n_tracked_files"] == 6
    assert not any(p.startswith((".venv", ".r-lib", "runs/", "private_evidence/", "data/")) for p in BR.tracked_files(root))
    assert report["license"].startswith("unresolved") and report["ownership"].startswith("unresolved") and report["public_release_allowed"] is False


@pytest.mark.parametrize("rel,content,rule", [
    ("src/keys.txt", DUMMY["key"], "secret"), ("src/config.py", DUMMY["github"], "secret"), ("docs/aws.md", DUMMY["aws"], "secret"),
    ("src/settings.py", DUMMY["password"], "secret"), ("docs/notes.md", "output in /Users/" + "alice/study/run1\n", "private_path"),
    ("docs/win.md", "C:\\Users\\" + "bob\\proteomics\\x\n", "private_path"), ("data/study_workbook.xlsx", b"PK\x03\x04dummy", "forbidden tracked file"),
    (".env", "TOKEN=x\n", "forbidden tracked file"), ("logs/private_run.log", "run log\n", "forbidden tracked file"),
    ("LICENSE", "MIT License\n\nPermission is hereby granted...\n", "license")])
def test_v117_negative_seeded_secret_private_or_license_blocks_release(tmp_path, rel, content, rule):
    root = scratch_tree(tmp_path, {rel: content})
    report = BR.check(root)
    assert report["blocking"] is True and any(f["path"] == rel and f["rule"] == rule for f in report["findings"]), report["findings"]
    with pytest.raises(BR.ReleaseError) as error:
        BR.export(root, "1.0.0", tmp_path / "out", _authorization(tmp_path, "1.0.0"))
    assert error.value.code == "E_RELEASE_HYGIENE" and not (tmp_path / "out").exists()


def test_v117_negative_changed_protected_baseline_blocks(tmp_path):
    root = scratch_tree(tmp_path)
    (root / "pipeline/old_script.R").write_text("x <- 2\n", encoding="utf-8")
    report = BR.check(root)
    assert report["blocking"] and report["findings"] == [{"path": "pipeline/old_script.R", "kind": "protected baseline hash_mismatch", "rule": "baseline"}]


def test_v117_ci_runner_paths_are_not_private(tmp_path):
    root = scratch_tree(tmp_path, {"docs/ci.log": "/home/runner/work/x/y.py and /Users/runner/work and C:\\Users\\runneradmin\\x\n"})
    assert BR.check(root)["findings"] == []


def test_v117_this_repository_has_no_blocking_hygiene_finding():
    report = BR.check(ROOT)
    assert report["findings"] == [], report["findings"]
    assert report["license"] == "PolyForm-Noncommercial-1.0.0"                                         # D-63, Maintainer decision
    assert report["public_release_allowed"] is False                                                  # no disclosure review recorded


def _authorization(tmp_path: Path, version: str, scope: str = "private_successor", **extra) -> Path:
    path = tmp_path / f"authorization-{scope}-{version}.json"
    path.write_text(json.dumps({"release_authorized": True, "version": version, "authorized_by": "synthetic test record", "scope": scope, **extra}), encoding="utf-8")
    return path


def test_v118_tooling_exports_a_new_versioned_directory_with_an_independent_manifest(tmp_path):
    root = scratch_tree(tmp_path)
    before = {p: sha(root / p) for p in BR.tracked_files(root)}
    release = BR.export(root, "1.2.3", tmp_path / "release-root", _authorization(tmp_path, "1.2.3"))
    target = tmp_path / "release-root" / "proteomics-downstream-pipeline-1.2.3"
    assert Path(release["path"]) == target.resolve()
    exported = {p.relative_to(target).as_posix(): sha(p) for p in target.rglob("*") if p.is_file()}
    assert set(exported) == set(before) | {"MANIFEST.sha256", "release.json"}
    assert all(exported[p] == before[p] for p in before)                       # byte-identical copies of tracked files only
    lines = (target / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines()
    assert {line.split("  ", 1)[1]: line.split("  ", 1)[0] for line in lines} == before
    assert {p: sha(root / p) for p in before} == before                           # the source (and its protected baseline) is unchanged
    record = json.loads((target / "release.json").read_text(encoding="utf-8"))
    assert record["license"].startswith("unresolved") and record["scope"] == "private_successor" and record["manifest_sha256"] == sha(target / "MANIFEST.sha256")


def test_v118_negative_collision_missing_authorization_and_unresolved_license_refuse(tmp_path):
    root = scratch_tree(tmp_path)
    (tmp_path / "taken").mkdir()
    cases = [(lambda: BR.export(root, "1.0.0", tmp_path / "taken", _authorization(tmp_path, "1.0.0")), "E_PATH_COLLISION"),
             (lambda: BR.export(root, "1.0.0", tmp_path / "a", None), "E_RELEASE_NOT_AUTHORIZED"),
             (lambda: BR.export(root, "1.0.0", tmp_path / "b", _authorization(tmp_path, "0.9.0")), "E_RELEASE_NOT_AUTHORIZED"),
             (lambda: BR.export(root, "1.0.0", tmp_path / "c", _authorization(tmp_path, "1.0.0", scope="public")), "E_LICENSE_UNRESOLVED"),
             (lambda: BR.export(root, "1.0.0", root / "inside", _authorization(tmp_path, "1.0.0")), "E_RELEASE_DESTINATION"),
             (lambda: BR.export(root, "v1", tmp_path / "d", _authorization(tmp_path, "v1")), "E_RELEASE_VERSION")]
    for call, code in cases:
        with pytest.raises(BR.ReleaseError) as error:
            call()
        assert error.value.code == code
    assert sorted(p.name for p in tmp_path.iterdir() if p.is_dir()) == ["src-tree", "taken"] and not any((tmp_path / "taken").iterdir())


def test_v118_cli_refusal_is_typed_json(tmp_path):
    root = scratch_tree(tmp_path)
    code = BR.main(["export", "--root", str(root), "--version", "1.0.0", "--destination", str(tmp_path / "x")])
    assert code == 2 and not (tmp_path / "x").exists()
