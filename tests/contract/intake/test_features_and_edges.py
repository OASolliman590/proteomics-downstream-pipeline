"""V015 feature annotations/protein groups and V020 malformed-input edge cases."""
from __future__ import annotations

import csv
import json

import pytest

from intake_helpers import features_rows, minimal_config, observation_rows, read_tsv, write_tsv
from proteomics_pipeline.errors import ProteomicsError
from proteomics_pipeline.intake import service

OBS = [("A1", "A1", "NA", "NA", "a"), ("A2", "A2", "NA", "NA", "a"), ("B1", "B1", "NA", "NA", "b"), ("B2", "B2", "NA", "NA", "b")]
IDS = [o[0] for o in OBS]


def _setup(tmp_path, feature_rows, matrix_rows=None, **kwargs):
    fids = [r[0] for r in feature_rows[1:]]
    matrix_rows = matrix_rows or [["feature_id"] + IDS] + [[f, "10", "11", "12", "13"] for f in fids]
    matrix = write_tsv(tmp_path / "m.tsv", matrix_rows)
    obs = write_tsv(tmp_path / "o.tsv", observation_rows(OBS))
    feats = write_tsv(tmp_path / "f.tsv", feature_rows)
    return minimal_config(tmp_path, matrix=matrix, observations=obs, features=feats, group_levels=["a", "b"], **kwargs)


# ----------------------------------------------------------------------------- V015
def test_v015_opaque_features_sharing_accession_and_multigene_group(tmp_path):
    header = ["feature_id", "accessions", "gene_ids", "gene_symbols", "is_decoy", "is_contaminant", "protein_group_ambiguous"]
    rows = [header,
            ["grp:0001", '["Q9X001"]', '["1001"]', '["Abc1"]', "unknown", "false", "false"],
            ["grp:0002", '["Q9X001"]', '["1001"]', '["Abc1"]', "unknown", "unknown", "false"],
            ["grp:0003", '["Q9X001", "Q9X002"]', '["1001", "1002"]', '["Abc1", "Abc2"]', "unknown", "unknown", "true"]]
    config = _setup(tmp_path, rows)
    service.run_intake(config, tmp_path, tmp_path / "out")
    published = read_tsv(tmp_path / "out" / "features.tsv")
    assert [r[0] for r in published[1:]] == ["grp:0001", "grp:0002", "grp:0003"]   # three unique rows; symbol is not a key
    by_id = {r[0]: dict(zip(published[0], r)) for r in published[1:]}
    assert json.loads(by_id["grp:0003"]["accessions"]) == ["Q9X001", "Q9X002"]
    assert json.loads(by_id["grp:0003"]["gene_ids"]) == ["1001", "1002"]
    assert by_id["grp:0003"]["protein_group_ambiguous"] == "true"
    assert all(by_id[f]["is_decoy"] == "unknown" for f in by_id)          # tri-state retained, not coerced to false
    assert by_id["grp:0001"]["is_contaminant"] == "false"


def test_v015_negative_duplicate_feature_id_with_different_annotation(tmp_path):
    header = features_rows([])[0]
    rows = [header, ["F1", '["Q1"]', "[]", '["A"]', "unknown", "unknown", "false"], ["F1", '["Q2"]', "[]", '["B"]', "unknown", "unknown", "false"]]
    matrix = [["feature_id"] + IDS, ["F1", "1", "2", "3", "4"]]
    config = _setup(tmp_path, rows, matrix_rows=matrix)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == "E_ID_DUPLICATE" and error.value.row == 3
    bad_flag = [header, ["F1", '["Q1"]', "[]", "[]", "maybe", "unknown", "false"]]
    config = _setup(tmp_path / "b", bad_flag, matrix_rows=matrix)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path / "b", tmp_path / "b" / "out")
    assert error.value.code == "E_FEATURE_FLAG" and error.value.field == "is_decoy"


# ----------------------------------------------------------------------------- V020
def _expect(tmp_path, config, code, *, row=None, field=None):
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == code, error.value.as_dict()
    if row is not None:
        assert error.value.row == row
    if field is not None:
        assert error.value.field == field
    assert not (tmp_path / "out").exists()
    return error.value


def test_v020_ragged_row(tmp_path):
    rows = [["feature_id"] + IDS, ["F1", "1", "2", "3", "4"], ["F2", "1", "2", "3"]]
    _expect(tmp_path, _setup(tmp_path, features_rows(["F1", "F2"]), rows), "E_INPUT_RAGGED", row=3)


def test_v020_duplicated_id(tmp_path):
    rows = [["feature_id"] + IDS, ["F1", "1", "2", "3", "4"], ["F1", "1", "2", "3", "5"]]
    _expect(tmp_path, _setup(tmp_path, features_rows(["F1"]), rows), "E_ID_DUPLICATE", row=3, field="feature_id")


def test_v020_mismatched_mask(tmp_path):
    config = _setup(tmp_path, features_rows(["F1"]), [["feature_id"] + IDS, ["F1", "1", "NA", "3", "4"]],
                    prior_imputation="mask_available",
                    original_mask=write_tsv(tmp_path / "mask.tsv", [["feature_id", "A1", "A2", "B1", "X9"], ["F1", "true", "false", "true", "true"]]))
    _expect(tmp_path, config, "E_ID_ALIGNMENT")
    config = _setup(tmp_path / "c", features_rows(["F1"]), [["feature_id"] + IDS, ["F1", "1", "NA", "3", "4"]],
                    prior_imputation="mask_available",
                    original_mask=write_tsv(tmp_path / "mask2.tsv", [["feature_id"] + IDS, ["F1", "true", "true", "true", "true"]]))
    _expect(tmp_path / "c", config, "E_MASK_INCONSISTENT", field="A2")


@pytest.mark.parametrize("token", ["Inf", "-inf", "NaN"])
def test_v020_nonfinite_values_rejected_not_treated_as_missing(tmp_path, token):
    rows = [["feature_id"] + IDS, ["F1", "1", token, "3", "4"]]
    _expect(tmp_path, _setup(tmp_path, features_rows(["F1"]), rows), "E_INPUT_NONFINITE", row=2, field="A2")


def test_v020_unknown_encoding(tmp_path):
    config = _setup(tmp_path, features_rows(["F1"]))
    (tmp_path / "m.tsv").write_bytes("feature_id\tA1\tA2\tB1\tB2\nF1\t1\t2\t3\t4\n# caf\xe9\n".encode("latin-1"))
    _expect(tmp_path, config, "E_INPUT_ENCODING", row=3)


def test_v020_unsupported_assay(tmp_path):
    config = _setup(tmp_path, features_rows(["F1"]), provenance_fields={"assay": "phosphoproteomics"})
    _expect(tmp_path, config, "E_UNSUPPORTED_SCOPE", field="assay")


def test_v020_quoted_tab_and_unicode_identifier_roundtrip(tmp_path):
    odd = "grp\tΩ-01 \"q\""
    header = features_rows([])[0]
    rows = [header, [odd, json.dumps(["Q1"]), "[]", json.dumps(["Ωmega"]), "unknown", "unknown", "false"]]
    config = _setup(tmp_path, rows, [["feature_id"] + IDS, [odd, "1.5", "2", "3", "4"]])
    raw = (tmp_path / "m.tsv").read_text(encoding="utf-8")
    assert '"grp\tΩ-01 ""q"""' in raw            # the fixture really is a quoted TAB
    service.run_intake(config, tmp_path, tmp_path / "out")
    with (tmp_path / "out" / "matrix.tsv").open(encoding="utf-8", newline="") as handle:
        published = list(csv.reader(handle, delimiter="\t"))
    assert published[1][0] == odd and published[1][1] == "1.5"
    assert service.verify_bundle(tmp_path / "out")["feature_ids"] == [odd]
