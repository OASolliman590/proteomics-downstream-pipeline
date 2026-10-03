"""Platform-independent plan hash (review follow-up, decisions.md open item 12; D-42).

CI run 37087074374 produced different plan hashes for the same commit and
configuration on Ubuntu (7ee163b2...) and Windows (5a774df3...).  The causes
were byte-level, not scientific: a Windows checkout gives CRLF input tables
(`* text=auto`), R and Python text-mode writers emit CRLF on Windows, and the
plan hashed the platform's R session information and solved environment.

These tests fail on the pre-fix tree (fa4d948/106fb20) and pass after it.
Oracles are hashlib and byte manipulations written here, not the adapter.
"""
from __future__ import annotations

import copy
import hashlib
import json
import shutil
from pathlib import Path

import pytest

from proteomics_pipeline import design_service, workflow
from proteomics_pipeline.errors import IntegrityError

ROOT = Path(__file__).resolve().parents[2]


def _needs_r():
    import importlib.util
    spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module.rscript() is None


pytestmark = pytest.mark.skipif(_needs_r(), reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")


def _copy_example(destination: Path, *, crlf: bool) -> Path:
    shutil.copytree(ROOT / "configs" / "examples", destination)
    if crlf:   # what a Windows checkout of `* text=auto` files looks like
        for path in (destination / "fixtures").iterdir():
            if path.suffix in (".tsv", ".json"):
                path.write_bytes(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
    return destination / "example-independent.json"


def test_same_plan_hash_from_two_roots_and_from_crlf_inputs(tmp_path):
    first = _copy_example(tmp_path / "a", crlf=False)
    second = _copy_example(tmp_path / "deeper root" / "with spaces Ω" / "b", crlf=False)
    third = _copy_example(tmp_path / "windows-checkout", crlf=True)
    assert (tmp_path / "windows-checkout" / "fixtures" / "independent-abundance.tsv").read_bytes().count(b"\r\n") > 0
    hashes = {}
    for name, config in (("a", first), ("b", second), ("crlf", third)):
        payload, code = workflow.plan_command(config, tmp_path / "plans" / name / "plan.json")
        assert code == 0
        hashes[name] = payload["plan_hash"]
    assert hashes["a"] == hashes["b"] == hashes["crlf"], hashes
    plan_a = json.loads((tmp_path / "plans" / "a" / "plan.json").read_text(encoding="utf-8"))
    plan_c = json.loads((tmp_path / "plans" / "crlf" / "plan.json").read_text(encoding="utf-8"))
    for source_a, source_c in zip(plan_a["sources"], plan_c["sources"]):
        lf = (ROOT / "configs" / "examples" / "fixtures" / source_a["file_name"]).read_bytes().replace(b"\r\n", b"\n")
        assert source_a["sha256"] == source_c["sha256"] == hashlib.sha256(lf).hexdigest()   # independent oracle
    # absolute locations never enter the hashed plan
    text = json.dumps(design_service.semantic_view(plan_a))
    assert str(tmp_path) not in text and str(tmp_path).replace("\\", "/") not in text


def test_hash_ignores_platform_bytes_of_planning_artifacts_and_environment(tmp_path):
    """Simulate a Windows R writer (CRLF JSON), a different sessionInfo() and a different solved environment."""
    config_path = _copy_example(tmp_path / "cfg", crlf=False)
    config = workflow.load(config_path)
    requested = workflow.requested_capabilities(config)
    base = tmp_path / "plan-artifacts"; base.mkdir()
    ledger = workflow.Ledger(tmp_path, "run-portability")
    plan = workflow.plan_into(config_path, config, plan_root=tmp_path, base=base, ledger=ledger, requested=requested)
    stage_dirs = {"intake": base / "inputs", "preprocessing": base / "preprocessing", "design": base / "designs"}
    # rewrite every planning JSON artifact with CRLF and change the session information, keeping stage results consistent
    for directory in stage_dirs.values():
        result_path = directory / "stage-result.json"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        for output in result["outputs"]:
            path = directory / output["relative_path"]
            if output["result_type"] == "session_info":
                path.write_bytes(b"R version 4.6.1 (2026-06-24 ucrt)\r\nPlatform: x86_64-w64-mingw32/x64\r\n")
            elif path.suffix == ".json":
                path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
            output["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        result_path.write_bytes(json.dumps(result, indent=2).encode("utf-8"))
    environment = copy.deepcopy(plan["environment"]); environment["python"] = {"version": "3.12.99"}
    rebuilt = design_service.build_plan(config, plan_root=tmp_path, stage_dirs=stage_dirs,
                                        sources=[{"artifact_id": s["artifact_id"], "path": str(config_path.parent / "fixtures" / s["file_name"]), "sha256": s["sha256"]} for s in plan["sources"]],
                                        environment=environment, code_sha256=plan["code_sha256"], capability_plan=plan["capability_plan"])
    assert any(a["sha256"] != b["sha256"] for a, b in zip(plan["artifacts"], rebuilt["artifacts"]))   # raw bytes really differ
    assert rebuilt["plan_hash"] == plan["plan_hash"]
    assert rebuilt["integrity_sha256"] != plan["integrity_sha256"]                                     # provenance still recorded


def test_verify_plan_detects_provenance_and_semantic_tampering(tmp_path):
    config = _copy_example(tmp_path / "cfg", crlf=False)
    workflow.plan_command(config, tmp_path / "o" / "plan.json")
    path = tmp_path / "o" / "plan.json"
    original = json.loads(path.read_text(encoding="utf-8"))
    assert design_service.verify_plan(path)["plan_hash"] == original["plan_hash"]
    tampered = copy.deepcopy(original); tampered["environment"]["python"] = {"version": "0.0"}
    path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(IntegrityError, match="provenance"):
        design_service.verify_plan(path)
    tampered = copy.deepcopy(original); tampered["runtime"]["seed"] = 1
    path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(IntegrityError, match="plan hash"):
        design_service.verify_plan(path)


def test_content_hash_and_code_identity_are_line_ending_independent(tmp_path):
    from proteomics_pipeline.provenance import content_sha256
    lf, crlf, binary = tmp_path / "lf.tsv", tmp_path / "crlf.tsv", tmp_path / "b.png"
    lf.write_bytes(b"a\tb\n1\t2\n"); crlf.write_bytes(b"a\tb\r\n1\t2\r\n"); binary.write_bytes(b"\x89PNG\r\n\x1a\n\x00\r\n")
    assert content_sha256(lf) == content_sha256(crlf) == hashlib.sha256(b"a\tb\n1\t2\n").hexdigest()
    assert content_sha256(binary) == hashlib.sha256(binary.read_bytes()).hexdigest()               # NUL bytes: hashed raw
    package = tmp_path / "pkg"; shutil.copytree(ROOT / "src" / "proteomics_pipeline", package, ignore=shutil.ignore_patterns("__pycache__"))
    repo = tmp_path / "repo"
    shutil.copytree(ROOT / "r" / "proteomicsCore", repo / "r" / "proteomicsCore", ignore=shutil.ignore_patterns("tests"))
    (repo / "scripts" / "maintained").mkdir(parents=True); shutil.copy(ROOT / "scripts" / "maintained" / "run_stage.R", repo / "scripts" / "maintained")
    before = workflow.code_manifest(package, repo)
    for path in [*package.rglob("*.html"), repo / "r" / "proteomicsCore" / "DESCRIPTION", repo / "r" / "proteomicsCore" / "NAMESPACE"]:
        path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
    assert workflow.code_manifest(package, repo)["sha256"] == before["sha256"]
    assert workflow.code_manifest()["sha256"] == before["sha256"]
