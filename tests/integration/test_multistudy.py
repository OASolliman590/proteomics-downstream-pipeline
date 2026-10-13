"""R15d acceptance: V185-V188 (FR-185 to FR-188; SM46; ADR 0011) and the slot adapters the slot cases need (FR-189, FR-191, FR-192).

Synthetic studies only (tests/fixtures/methods_extension/multistudy). The network layer is patched: no public dataset is
ever downloaded. Every negative case sits next to its positive case.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_pipeline import multistudy
from proteomics_pipeline.errors import ProteomicsError
from proteomics_pipeline.multistudy.manifest import file_sha256
from proteomics_pipeline.multistudy.runner import commit_status, run_study, status_record

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


F = _load("multistudy_factory", ROOT / "tests" / "fixtures" / "methods_extension" / "multistudy" / "factory.py")
FETCH = _load("fetch_public_study", ROOT / "scripts" / "maintained" / "fetch_public_study.py")
VERIFY_STATUS = _load("verify_private_status", ROOT / "scripts" / "maintained" / "verify_private_status.py")


def _manifest(tmp_path, **kwargs):
    return json.loads(F.make_study(tmp_path / "study", **kwargs).read_text(encoding="utf-8"))


# ----------------------------------------------------------------------------- V185
def test_v185_a_valid_manifest_validates_and_hashes_every_input(tmp_path):
    path = F.make_study(tmp_path / "study")
    record = multistudy.validate_manifest(json.loads(path.read_text(encoding="utf-8")), path.parent)
    assert record["state"] == "PASS" and record["files_checked"] == 4 and record["design"] == "two_group"


def test_v185_negative_a_changed_input_file_fails_the_hash_check_before_any_run(tmp_path):
    path = F.make_study(tmp_path / "study")
    (path.parent / "abundance.tsv").write_bytes((path.parent / "abundance.tsv").read_bytes() + b"\n")
    with pytest.raises(ProteomicsError) as info:
        multistudy.validate_manifest(json.loads(path.read_text(encoding="utf-8")), path.parent)
    assert info.value.code == "E_MANIFEST_HASH"


def test_v185_negative_an_unknown_design_or_format_is_refused(tmp_path):
    base = _manifest(tmp_path)
    with pytest.raises(ProteomicsError) as design_info:
        multistudy.validate_manifest({**base, "design": "nested_mystery"}, tmp_path / "study")
    assert design_info.value.code == "E_MANIFEST_DESIGN"
    with pytest.raises(ProteomicsError) as format_info:
        multistudy.validate_manifest({**base, "input": {**base["input"], "format": "excel_workbook"}}, tmp_path / "study")
    assert format_info.value.code == "E_MANIFEST_DESIGN"


def test_v185_negative_a_study_identifier_that_names_a_file_or_person_is_refused(tmp_path):
    base = _manifest(tmp_path)
    with pytest.raises(ProteomicsError) as info:
        multistudy.validate_manifest({**base, "study_id": "Dina Smith/run1"}, tmp_path / "study")
    assert info.value.code == "E_PRIVATE_IDENTIFIER"


# ----------------------------------------------------------------------------- V186
@pytest.fixture(scope="module")
def study_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("study")
    manifest = F.make_study(tmp / "study")
    result = run_study(manifest, tmp / "run")
    return {"manifest": manifest, "run": tmp / "run", "result": result}


def test_v186_validate_run_verify_and_acceptance_each_report_their_own_state(study_run):
    checks = study_run["result"]["checks"]
    assert checks == {"validate": "PASS", "run": "COMPLETED", "verify": "PASS", "acceptance": "PASS"}, checks
    assert study_run["result"]["state"] == "PASS"
    assert study_run["result"]["check_counts"] == {"passed": 4, "total": 4}


def test_v186_negative_an_artifact_changed_after_the_run_fails_verify(tmp_path, study_run, capsys):
    from proteomics_pipeline.cli import main as cli_main
    copied = tmp_path / "run"
    import shutil
    shutil.copytree(study_run["run"], copied)
    (copied / "dea" / "zero_null.tsv").write_bytes((copied / "dea" / "zero_null.tsv").read_bytes() + b"\n")
    capsys.readouterr()
    assert cli_main(["verify", "--run", str(copied), "--json"]) == 5
    out = capsys.readouterr().out
    # the verifier reports the stage-output hash mismatch as E_INTEGRITY (the existing code for changed stage outputs)
    assert '"valid":false' in out and ("E_INTEGRITY" in out or "E_ARTIFACT_HASH" in out)


# ----------------------------------------------------------------------------- V187
class _Response:
    def __init__(self, payload):
        self.payload = payload

    def read(self):
        return self.payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _approval(payload_sha: str, url: str = "https://example.invalid/public/ACC1.tsv"):
    return {"approval_ref": "approval-001", "approved_by": "Maintainer", "approved_on": "2026-10-13",
            "accessions": [{"accession": "ACC1", "url": url, "file_name": "ACC1.tsv", "sha256": payload_sha}]}


def test_v187_a_public_dataset_is_downloaded_only_against_a_matching_approval_and_recorded(tmp_path):
    payload = b"synthetic public stand-in\n"
    import hashlib
    approval = _approval(hashlib.sha256(payload).hexdigest())
    calls = []
    record = FETCH.fetch_public(approval, "ACC1", tmp_path / "outside", opener=lambda url: calls.append(url) or _Response(payload),
                                repo_root=ROOT)
    assert calls == ["https://example.invalid/public/ACC1.tsv"]
    assert record["observed_sha256"] == approval["accessions"][0]["sha256"] and record["url"] == approval["accessions"][0]["url"]
    assert (tmp_path / "outside" / "ACC1.tsv").read_bytes() == payload


def test_v187_negative_no_approval_means_no_network_call_and_no_file(tmp_path):
    def must_not_be_called(url):
        raise AssertionError("the network layer was reached without an approval")
    with pytest.raises(ProteomicsError) as info:
        FETCH.fetch_public({"approval_ref": "x", "approved_by": "M", "approved_on": "2026-10-13", "accessions": []}, "ACC1", tmp_path / "outside",
                           opener=must_not_be_called, repo_root=ROOT)
    assert info.value.code == "E_PUBLIC_DATA_UNAPPROVED"
    assert not (tmp_path / "outside").exists()


def test_v187_negative_a_changed_file_and_an_in_repository_destination_are_refused(tmp_path):
    import hashlib
    approved = _approval(hashlib.sha256(b"approved bytes").hexdigest())
    with pytest.raises(ProteomicsError) as hash_info:
        FETCH.fetch_public(approved, "ACC1", tmp_path / "outside", opener=lambda url: _Response(b"different bytes"), repo_root=ROOT)
    assert hash_info.value.code == "E_MANIFEST_HASH"
    assert not (tmp_path / "outside" / "ACC1.tsv").exists()
    with pytest.raises(ProteomicsError) as repo_info:
        FETCH.fetch_public(approved, "ACC1", ROOT / "tests" / "_scratch_public", opener=lambda url: _Response(b"approved bytes"), repo_root=ROOT)
    assert repo_info.value.code == "E_PUBLIC_DATA_UNAPPROVED"
    assert not (ROOT / "tests" / "_scratch_public").exists()


def test_v187_a_public_manifest_without_an_approval_record_fails_validation(tmp_path):
    path = F.make_study(tmp_path / "public", origin="public")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    with pytest.raises(ProteomicsError) as info:
        multistudy.validate_manifest(manifest, path.parent, approvals={})
    assert info.value.code == "E_PUBLIC_DATA_UNAPPROVED"


# ----------------------------------------------------------------------------- V188
def test_v188_a_private_study_commits_only_a_non_identifying_status_record(tmp_path, study_run):
    record = status_record(study_run["result"])
    assert set(record) == {"schema", "slot", "design", "state", "checks", "check_counts"}
    path = commit_status(study_run["result"], tmp_path / "status" / "synthetic.json")
    text = path.read_text(encoding="utf-8")
    assert "synthetic_two_group" not in text and "abundance" not in text and "study-manifest" not in text
    assert VERIFY_STATUS.check_record(path) == []


def test_v188_negative_a_status_record_that_would_carry_an_identifier_is_refused(tmp_path, study_run):
    leaking = dict(study_run["result"], _study_id="two_group")          # the design value appears in the record: a stand-in identifier
    with pytest.raises(ProteomicsError) as info:
        commit_status(leaking, tmp_path / "status" / "leak.json")
    assert info.value.code == "E_PRIVATE_IDENTIFIER"
    assert not (tmp_path / "status" / "leak.json").exists()
    bad = tmp_path / "bad"
    bad.mkdir()
    (bad / "bad.json").write_text(json.dumps({"schema": "x", "slot": "dda_label_free", "design": "two_group", "state": "PASS", "checks": {},
                                              "check_counts": {}, "input": "sample_07.raw"}), encoding="utf-8")
    assert VERIFY_STATUS.main(["--status-dir", str(bad)]) == 5


# ----------------------------------------------------------------------------- slot adapters (FR-189, FR-191, FR-192)
def test_slot_adapter_dia_columns_map_through_a_recorded_mapping_and_only_the_missing_column_is_refused():
    mapping = multistudy.manifest.map_columns("diann", ["Protein.Group", "Genes", "Precursor.Quantity"])
    assert mapping["mapping"] == {"protein_group": "Protein.Group", "genes": "Genes"} and len(mapping["mapping_sha256"]) == 64
    with pytest.raises(ProteomicsError) as info:
        multistudy.manifest.map_columns("diann", ["Protein.Group"])
    assert info.value.code == "E_MAPPING_UNMAPPED" and "genes" in info.value.message


def test_slot_adapter_maxquant_filters_remove_reverse_hits_and_contaminants_and_count_them():
    rows = [{"Protein IDs": "P1", "Reverse": "", "Potential contaminant": ""}, {"Protein IDs": "REV_P2", "Reverse": "+", "Potential contaminant": ""},
            {"Protein IDs": "CON_P3", "Reverse": "", "Potential contaminant": "+"}]
    kept, counts = multistudy.manifest.maxquant_filters(rows, {"reverse": "Reverse", "contaminant": "Potential contaminant"})
    assert [r["Protein IDs"] for r in kept] == ["P1"] and counts == {"input_rows": 3, "removed_reverse": 1, "removed_contaminant": 1, "kept": 1}


def test_slot_adapter_tmt_batch_confounding_is_an_adaptation_and_an_aliased_batch_term_is_not_fitted():
    confounded = [{"group": "A", "batch": "b1"}, {"group": "B", "batch": "b1"}, {"group": "A", "batch": "b2"}, {"group": "A", "batch": "b2"}]
    result = multistudy.manifest.tmt_batch_confounding(confounded, "group", "batch")
    assert result["confounded_batches"] == ["b2"] and result["batch_term_fitted"] is False and result["adaptation"]["reason"] == "batch_confounded_with_group"
    balanced = [{"group": g, "batch": b} for g in ("A", "B") for b in ("b1", "b2")]
    assert multistudy.manifest.tmt_batch_confounding(balanced, "group", "batch")["batch_term_fitted"] is True
