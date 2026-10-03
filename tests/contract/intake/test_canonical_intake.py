"""V011, V012, V018, V019 — canonical intake, scale handling, immutable bundle and report."""
from __future__ import annotations

import json
import math
import shutil

import pytest

from intake_helpers import (EXAMPLES, FIXTURES, example_config, features_rows, minimal_config, observation_rows, read_tsv, sha256,
                            write_tsv)
from proteomics_pipeline.errors import CollisionError, IntegrityError, ProteomicsError
from proteomics_pipeline.intake import hierarchy, service

FIX = EXAMPLES / "fixtures"


def _wide_oracle():
    """Independent pivot of the public wide fixture by explicit keys."""
    rows = read_tsv(FIX / "independent-abundance.tsv")
    header = rows[0][1:]
    return {(r[0], o): (None if v == "NA" else float(v)) for r in rows[1:] for o, v in zip(header, r[1:])}, [r[0] for r in rows[1:]], header


# ----------------------------------------------------------------------------- V011
def test_v011_wide_and_shuffled_long_give_identical_ordered_values_and_masks(tmp_path):
    wide_cfg = example_config()
    long_cfg = example_config()
    long_cfg["input"]["format"] = "long_tsv"
    long_cfg["input"]["matrix"] = str(FIXTURES / "independent-long-shuffled.tsv")
    service.run_intake(wide_cfg, EXAMPLES, tmp_path / "wide")
    service.run_intake(long_cfg, EXAMPLES, tmp_path / "long")
    for name in ("matrix.tsv", "numeric_mask.tsv", "observed_mask.tsv", "prior_imputed_mask.tsv", "observations.tsv", "features.tsv"):
        assert (tmp_path / "wide" / name).read_bytes() == (tmp_path / "long" / name).read_bytes(), name
    oracle, features, observations = _wide_oracle()
    published = read_tsv(tmp_path / "long" / "matrix.tsv")
    assert published[0] == ["feature_id"] + observations
    assert [r[0] for r in published[1:]] == features
    for row in published[1:]:
        for observation, token in zip(observations, row[1:]):
            expected = oracle[(row[0], observation)]
            assert (token == "NA") if expected is None else float(token) == expected
    na_cells = [(r[0], o) for r in published[1:] for o, t in zip(observations, r[1:]) if t == "NA"]
    assert na_cells == [("P07", "U2"), ("P08", "T1")]
    mask = read_tsv(tmp_path / "long" / "observed_mask.tsv")
    assert sum(t == "false" for r in mask[1:] for t in r[1:]) == 2


def test_v011_negative_duplicate_long_pair_rejected_without_publication(tmp_path):
    rows = read_tsv(FIXTURES / "independent-long-shuffled.tsv")
    rows.append(list(rows[1]))
    bad = write_tsv(tmp_path / "dup-long.tsv", rows)
    config = example_config()
    config["input"].update({"format": "long_tsv", "matrix": str(bad)})
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, EXAMPLES, tmp_path / "out")
    assert error.value.code == "E_ID_DUPLICATE" and error.value.row == len(rows)
    assert not (tmp_path / "out").exists()


def test_v011_negative_duplicate_header_and_missing_long_key(tmp_path):
    rows = read_tsv(FIX / "independent-abundance.tsv")
    rows[0][2] = rows[0][1]
    config = example_config()
    config["input"]["matrix"] = str(write_tsv(tmp_path / "dup-header.tsv", rows))
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, EXAMPLES, tmp_path / "out")
    assert error.value.code == "E_ID_DUPLICATE" and error.value.field == "C1"
    long_rows = read_tsv(FIXTURES / "independent-long-shuffled.tsv")[:-1]
    config = example_config()
    config["input"].update({"format": "long_tsv", "matrix": str(write_tsv(tmp_path / "short-long.tsv", long_rows))})
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, EXAMPLES, tmp_path / "out2")
    assert error.value.code == "E_ID_MISSING"


# ----------------------------------------------------------------------------- V012
def _one_row(tmp_path, values, **kwargs):
    ids = [f"O{i}" for i in range(1, len(values) + 1)]
    matrix = write_tsv(tmp_path / "m.tsv", [["feature_id"] + ids, ["F1"] + [str(v) for v in values]])
    obs = write_tsv(tmp_path / "o.tsv", observation_rows([(o, o, "NA", "NA", "g1" if i < 3 else "g2") for i, o in enumerate(ids)]))
    feats = write_tsv(tmp_path / "f.tsv", features_rows(["F1"]))
    return minimal_config(tmp_path, matrix=matrix, observations=obs, features=feats, group_levels=["g1", "g2"], **kwargs), ids


def test_v012_log2_zero_and_negative_values_are_preserved(tmp_path):
    config, ids = _one_row(tmp_path, [1, 2, 4, 0, -1], source_scale="log2", transform="none")
    service.run_intake(config, tmp_path, tmp_path / "out")
    row = read_tsv(tmp_path / "out" / "matrix.tsv")[1]
    assert [float(x) for x in row[1:]] == [1.0, 2.0, 4.0, 0.0, -1.0]
    assert read_tsv(tmp_path / "out" / "observed_mask.tsv")[1][1:] == ["true"] * 5


def test_v012_linear_values_logged_once_and_explicit_zero_is_missing(tmp_path):
    config, ids = _one_row(tmp_path, [1, 2, 4, 0], source_scale="linear_positive", transform="log2", zero_as_missing=True)
    service.run_intake(config, tmp_path, tmp_path / "out")
    row = read_tsv(tmp_path / "out" / "matrix.tsv")[1]
    assert row[1:4] == ["0.0", "1.0", "2.0"] and row[4] == "NA"   # oracle log2([1,2,4]) = [0,1,2]
    report = json.loads((tmp_path / "out" / "intake_report.json").read_text(encoding="utf-8"))
    assert report["explicit_zero_missing_cells"] == 1 and report["scale"]["output_scale"] == "log2"
    assert read_tsv(tmp_path / "out" / "observed_mask.tsv")[1][1:] == ["true", "true", "true", "false"]


@pytest.mark.parametrize("values,zero_missing", [([1, 2, 4, 0, -1], True), ([1, 2, 4, 0], False)])
def test_v012_invalid_linear_nonpositive_is_distinct_from_encoded_zero(tmp_path, values, zero_missing):
    config, _ = _one_row(tmp_path, values, source_scale="linear_positive", transform="log2", zero_as_missing=zero_missing)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == "E_LINEAR_NONPOSITIVE"
    assert not (tmp_path / "out").exists()


def test_v012_negative_second_log_fails_before_any_transform(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(hierarchy, "log2_once", lambda v: calls.append(v))
    config, _ = _one_row(tmp_path, [1, 2, 4], source_scale="log2", transform="log2")
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == "E_SCALE_SECOND_LOG" and calls == []
    with pytest.raises(ProteomicsError) as error:
        hierarchy.scale_transition("ratio_log2", "log2")
    assert error.value.code == "E_SCALE_SECOND_LOG"
    assert hierarchy.scale_transition("linear_positive", "none")["inference_reason_code"] == "E_LOG2_REQUIRED"
    assert hierarchy.scale_transition("standardized_unknown", "none")["inference_reason_code"] == "E_SCALE_UNKNOWN_INFERENCE"


# ----------------------------------------------------------------------------- V018
def test_v018_manifest_is_self_consistent_and_roundtrips_exactly(tmp_path):
    service.run_intake(example_config(), EXAMPLES, tmp_path / "bundle")
    manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text(encoding="utf-8"))
    for artifact in manifest["artifacts"]:
        assert sha256(tmp_path / "bundle" / artifact["relative_path"]) == artifact["sha256"]   # independent hashlib oracle
    assert service.verify_bundle(tmp_path / "bundle")["feature_ids"] == [f"P0{i}" for i in range(1, 9)]
    oracle, _, _ = _wide_oracle()
    for row in read_tsv(tmp_path / "bundle" / "matrix.tsv")[1:]:
        for observation, token in zip(manifest["observation_ids"], row[1:]):
            value = oracle[(row[0], observation)]
            assert (token == "NA" and value is None) or float(token) == value
    sources = {s["artifact_id"]: s["sha256"] for s in manifest["sources"]}
    # D-42: sources are identified by their LF-normalised content hash (independent hashlib oracle; equals the raw hash for an LF file)
    import hashlib
    assert sources["source.matrix"] == hashlib.sha256((FIX / "independent-abundance.tsv").read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def test_v018_rerun_collision_and_changed_source_are_refused(tmp_path):
    config = example_config()
    service.run_intake(config, EXAMPLES, tmp_path / "bundle")
    with pytest.raises(CollisionError):
        service.run_intake(config, EXAMPLES, tmp_path / "bundle")
    copy_dir = tmp_path / "src"; shutil.copytree(FIX, copy_dir)
    config["input"]["matrix"] = str(copy_dir / "independent-abundance.tsv")
    recorded = service.source_inputs(config, EXAMPLES)
    text = (copy_dir / "independent-abundance.tsv").read_text(encoding="utf-8").replace("9.7", "9.8", 1)
    (copy_dir / "independent-abundance.tsv").write_text(text, encoding="utf-8")
    with pytest.raises(IntegrityError) as error:
        service.run_intake(config, EXAMPLES, tmp_path / "bundle2", expected_sources=recorded)
    assert error.value.code == "E_SOURCE_CHANGED" and not (tmp_path / "bundle2").exists()


@pytest.mark.parametrize("mutation", ["delete_mask", "rename_observation"])
def test_v018_negative_mutated_bundle_fails_verification(tmp_path, mutation):
    service.run_intake(example_config(), EXAMPLES, tmp_path / "bundle")
    if mutation == "delete_mask":
        (tmp_path / "bundle" / "observed_mask.tsv").unlink()
    else:
        path = tmp_path / "bundle" / "observations.tsv"
        path.write_text(path.read_text(encoding="utf-8").replace("C1\t", "C1x\t", 1), encoding="utf-8")
    with pytest.raises(IntegrityError):
        service.verify_bundle(tmp_path / "bundle")


# ----------------------------------------------------------------------------- V019
def test_v019_report_exports_declarations_counts_and_unknown_provenance(tmp_path):
    ids = ["S1a", "S1b", "S2", "S3", "S4"]
    matrix = write_tsv(tmp_path / "m.tsv", [["feature_id"] + ids, ["F1", "10", "10.5", "NA", "11", "12"], ["F2", "8", "NA", "9", "NA", "9.5"]])
    obs = write_tsv(tmp_path / "o.tsv", observation_rows([("S1a", "U1", "NA", "inj1", "a"), ("S1b", "U1", "NA", "inj2", "a"), ("S2", "U2", "NA", "NA", "a"),
                                                          ("S3", "U3", "NA", "NA", "b"), ("S4", "U4", "NA", "NA", "b")]))
    feats = write_tsv(tmp_path / "f.tsv", features_rows(["F1", "F2"]))
    config = minimal_config(tmp_path, matrix=matrix, observations=obs, features=feats, group_levels=["a", "b"],
                            technical={"mode": "aggregate", "method": "mean_log2", "biological_unit_column": "biological_unit_id"})
    config["input"]["normalization_state"] = "unknown"
    config["study"]["tissue"] = None
    service.run_intake(config, tmp_path, tmp_path / "out")
    report = json.loads((tmp_path / "out" / "intake_report.json").read_text(encoding="utf-8"))
    assert report["grain"] == "protein" and report["scale"]["source_scale"] == "log2"
    assert report["n_features"] == 2 and report["n_source_observations"] == 5 and report["n_canonical_observations"] == 4
    groups = {g["group"]: g for g in report["groups"]}
    assert groups["a"]["n_injections"] == 3 and groups["a"]["n_biological_units"] == 2 and groups["b"]["n_biological_units"] == 2
    # hand count on canonical units: F1/U2 (S2=NA) and F2/U3 (S3=NA); U1 keeps F2 via its observed injection
    assert report["missingness"]["total_cells"] == 8
    assert report["missingness"]["numeric_missing_cells"] == 2
    fields = {gap["field"] for gap in report["provenance_gaps"]}
    assert {"study.tissue", "input.normalization_state"} <= fields
    assert report["asserted_unknown_facts"] == [] and report["declared"]["tissue"] is None
    text = (tmp_path / "out" / "intake_report.json").read_text(encoding="utf-8").lower()
    assert "confirmed" not in text and "pass" not in text.split()


def test_v019_negative_unknown_mask_stays_unknown_not_all_observed(tmp_path):
    config = example_config()
    config["input"]["prior_imputation"] = "unknown"
    service.run_intake(config, EXAMPLES, tmp_path / "out")
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["original_observed_mask"]["state"] == "unknown"
    assert not (tmp_path / "out" / "observed_mask.tsv").exists()
    report = json.loads((tmp_path / "out" / "intake_report.json").read_text(encoding="utf-8"))
    assert report["missingness"]["originally_observed_cells"] is None
    assert any(g["field"] == "input.original_observed_mask" for g in report["provenance_gaps"])
