"""R12 tooling for the Maintainer-only cases V111/V112 (FR-111, FR-112; SM12).

V111 and V112 themselves stay NOT_RUN: they need the Maintainer's original archive, historical tables and a
maintained run of the private study. These tests exercise the tooling on synthetic stand-ins only: an archive and
inventory built here, and hand-made tables whose differences are computed in the test.
"""
from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

import pytest

from proteomics_pipeline import legacy_service as L

ROOT = Path(__file__).resolve().parents[2]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def synthetic_archive(tmp_path: Path, members: dict[str, bytes]) -> tuple[Path, Path]:
    archive = tmp_path / "synthetic-archive.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        for name, data in members.items():
            handle.writestr(name, data)
    inventory = {"archive": archive.name, "archive_sha256": sha(archive.read_bytes()), "file_count": len(members), "uncompressed_bytes": sum(map(len, members.values())),
                 "files": [{"path": n, "bytes": len(d), "sha256": sha(d)} for n, d in members.items()]}
    (tmp_path / "inventory.json").write_text(json.dumps(inventory), encoding="utf-8")
    return archive, tmp_path / "inventory.json"


def test_v111_tooling_verifies_archive_hashes_without_modifying_it(tmp_path):
    archive, inventory = synthetic_archive(tmp_path, {"00_INDEX/CONTENTS.csv": b"a,b\n1,2\n", "scripts/run.R": b"x <- 1\n"})
    before = sha(archive.read_bytes())
    result = L.verify_archive(archive, inventory)
    assert result["verified"] and result["file_count"] == 2 and result["archive_sha256_matches"] and result["problems"] == []
    assert sha(archive.read_bytes()) == before and sorted(p.name for p in tmp_path.iterdir()) == ["inventory.json", "synthetic-archive.zip"]   # nothing extracted


def test_v111_tooling_negative_tampered_member_and_missing_archive(tmp_path):
    archive, inventory = synthetic_archive(tmp_path, {"a.txt": b"original\n", "b.txt": b"kept\n"})
    record = json.loads(inventory.read_text(encoding="utf-8")); record["files"][0]["sha256"] = "0" * 64; inventory.write_text(json.dumps(record), encoding="utf-8")
    result = L.verify_archive(archive, inventory)
    assert not result["verified"] and {"member": "a.txt", "reason": "hash_or_size_mismatch"} in result["problems"]
    with pytest.raises(L.LegacyError) as error:
        L.verify_archive(tmp_path / "absent.zip", inventory)
    assert error.value.code == "E_LEGACY_UNAVAILABLE" and error.value.exit_code == 3   # NOT_RUN, never an invented reconstruction


def write_tsv(path: Path, rows: list[list]) -> Path:
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows), encoding="utf-8")
    return path


def test_v112_tooling_compares_with_frozen_tolerances_and_preserves_differences(tmp_path):
    legacy = write_tsv(tmp_path / "legacy.tsv", [["feature_id", "contrast_id", "logFC", "P.Value"], ["F1", "B-A", 0.5, 0.01], ["F2", "B-A", -1.25, 0.2], ["F3", "B-A", 2.0, 1e-6]])
    maintained = write_tsv(tmp_path / "maintained.tsv", [["feature_id", "contrast_id", "effect", "p_value"], ["F1", "B-A", 0.5 + 1e-12, 0.01 * (1 + 5e-7)], ["F2", "B-A", -1.2, 0.2], ["F4", "B-A", 0.1, 0.5]])
    result = L.compare_tables(legacy, maintained, keys=["feature_id", "contrast_id"], fields={"effect": {"legacy": "logFC", "kind": "effect"}, "p_value": {"legacy": "P.Value", "kind": "statistic"}})
    by = {(r["feature_id"], r["field"]): r for r in result["rows"]}
    # hand oracle: F1 effect differs by 1e-12 (<= 1e-10) and P by a relative 5e-7 (<= 1e-6): both MATCH; F2 effect differs by 0.05: DIFFERENT
    assert by[("F1", "effect")]["status"] == "MATCH" and by[("F1", "p_value")]["status"] == "MATCH"
    assert by[("F2", "effect")]["status"] == "DIFFERENT" and abs(by[("F2", "effect")]["abs_difference"] - 0.05) < 1e-12 and by[("F2", "p_value")]["status"] == "MATCH"
    assert by[("F3", None)]["status"] == "NOT_REPRODUCED" and by[("F3", None)]["side"] == "legacy only" and by[("F4", None)]["side"] == "maintained only"
    assert result["counts"] == {"MATCH": 3, "DIFFERENT": 1, "NOT_REPRODUCED": 2}


def test_v112_tooling_negative_copied_values_and_unverified_runs_are_refused(tmp_path):
    legacy = write_tsv(tmp_path / "legacy.tsv", [["feature_id", "effect"], ["F1", 0.5]])
    copied = tmp_path / "copied.tsv"; copied.write_bytes(legacy.read_bytes())
    with pytest.raises(L.LegacyError) as error:
        L.compare_tables(legacy, copied, keys=["feature_id"], fields={"effect": {"legacy": "effect", "kind": "effect"}})
    assert error.value.code == "E_LEGACY_COPIED"                                 # stored values are not evidence that R was rerun
    run = tmp_path / "run"; (run / "dea").mkdir(parents=True)
    table = write_tsv(run / "dea" / "zero_null.tsv", [["feature_id", "effect"], ["F1", 0.4]])
    (run / "dea" / "stage-result.json").write_text(json.dumps({"state": "COMPLETED", "outputs": [{"relative_path": "zero_null.tsv", "sha256": "0" * 64}]}), encoding="utf-8")
    with pytest.raises(L.LegacyError) as error:
        L.maintained_table(run, "dea/zero_null.tsv")
    assert error.value.code == "E_LEGACY_MAINTAINED_RUN"
    (run / "dea" / "stage-result.json").write_text(json.dumps({"state": "COMPLETED", "outputs": [{"relative_path": "zero_null.tsv", "sha256": sha(table.read_bytes())}]}), encoding="utf-8")
    assert L.maintained_table(run, "dea/zero_null.tsv") == table


def test_v111_tooling_private_output_never_inside_the_repository(tmp_path):
    archive, _ = synthetic_archive(tmp_path, {"a.txt": b"x\n"})
    config = tmp_path / "local-legacy.json"
    config.write_text(json.dumps({"archive": str(archive), "maintained_run": str(tmp_path), "output_dir": str(ROOT / "runs" / "legacy-private"), "tables": []}), encoding="utf-8")
    with pytest.raises(L.LegacyError) as error:
        L.run(config)
    assert error.value.code == "E_PRIVATE_OUTPUT_IN_REPO" and not (ROOT / "runs" / "legacy-private").exists()


def test_template_holds_placeholders_only():
    template = json.loads((ROOT / "configs" / "local-legacy.template.json").read_text(encoding="utf-8"))
    assert template["archive"].startswith("/absolute/path/outside/the/repository/") and template["output_dir"].startswith("/absolute/path/outside/the/repository/")
