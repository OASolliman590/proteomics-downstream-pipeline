"""V016 — versioned vendor protein-level mappings with hand-authored expected rows."""
from __future__ import annotations

import json
import math

import pytest

from intake_helpers import FIXTURES, MAPPINGS, minimal_config, read_tsv
from proteomics_pipeline.errors import ProteomicsError
from proteomics_pipeline.intake import service

V = FIXTURES / "vendor"
L = math.log2

# Expected canonical log2 rows, authored by hand from the synthetic exports before running the adapter.
EXPECTED = {
    "diann": {"Q0001": [L(1000), L(2000), L(4000), None], "Q0002;Q0003": [L(500), None, L(250), L(125)], "Q0004": [3.0, 4.0, 5.0, 6.0]},
    "maxquant": {"Q0001": [10.0, 11.0, None, 9.0], "REV__Q0009": [2.0, None, 3.0, 4.0], "CON__Q0010;Q0011": [6.0, 7.0, 8.0, None]},
    "fragpipe": {"sp|Q0001|PRA_SYN": [1.0, 2.0, 3.0, None], "sp|Q0005|PRE_SYN": [None, 0.0, 0.0, 0.0]},
    "spectronaut": {"Q0001": [L(10), L(20), None, L(40)], "Q0007;Q0008": [None, L(3), L(6), L(12)]},
}
CASES = {
    "diann": ("diann_pg_matrix.tsv", "diann_pg_matrix.json", "diann_pg_matrix@1.8-headers-1", "diann_observations.tsv", False, ["NA", ""]),
    "maxquant": ("maxquant_proteinGroups.tsv", "maxquant_proteingroups.json", "maxquant_proteingroups@2.x-lfq-headers-1", "observations.tsv", True, ["NA", ""]),
    "fragpipe": ("fragpipe_combined_protein.tsv", "fragpipe_combined_protein.json", "fragpipe_combined_protein@20.x-maxlfq-headers-1", "observations.tsv", True, ["NA", ""]),
    "spectronaut": ("spectronaut_pg_pivot.tsv", "spectronaut_pg_pivot.json", "spectronaut_pg_pivot@18.x-pg-quantity-headers-1", "observations.tsv", False, ["NA", "", "NaN", "Filtered"]),
}


def _config(tmp_path, vendor, **override):
    matrix, mapping, profile, observations, zero, missing = CASES[vendor]
    values = dict(matrix=V / matrix, observations=V / observations, features="from_mapping", source_scale="linear_positive", transform="log2",
                  zero_as_missing=zero, fmt="mapped_protein", profile=profile, mapping=MAPPINGS / mapping, missing=missing, group_levels=["alpha", "beta"])
    values.update(override)
    config = minimal_config(tmp_path, **values)
    return config


@pytest.mark.parametrize("vendor", sorted(CASES))
def test_v016_profile_maps_only_declared_protein_columns(tmp_path, vendor):
    config = _config(tmp_path, vendor)
    service.run_intake(config, tmp_path, tmp_path / "out")
    rows = read_tsv(tmp_path / "out" / "matrix.tsv")
    assert rows[0] == ["feature_id", "A1", "A2", "B1", "B2"]
    published = {r[0]: [None if t == "NA" else float(t) for t in r[1:]] for r in rows[1:]}
    assert list(published) == list(EXPECTED[vendor])
    for feature, expected in EXPECTED[vendor].items():
        for got, want in zip(published[feature], expected):
            assert (got is None and want is None) or abs(got - want) <= 1e-10
    feature_rows = read_tsv(tmp_path / "out" / "features.tsv")
    assert all(dict(zip(feature_rows[0], r))["source_profile"] == CASES[vendor][2] for r in feature_rows[1:])


def test_v016_flags_and_members_follow_the_profile(tmp_path):
    service.run_intake(_config(tmp_path, "maxquant"), tmp_path, tmp_path / "out")
    rows = read_tsv(tmp_path / "out" / "features.tsv")
    by_id = {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}
    assert by_id["REV__Q0009"]["is_decoy"] == "true" and by_id["Q0001"]["is_decoy"] == "false"
    assert by_id["CON__Q0010;Q0011"]["is_contaminant"] == "true"
    assert json.loads(by_id["CON__Q0010;Q0011"]["gene_symbols"]) == ["Gnx", "Gny"]
    assert by_id["CON__Q0010;Q0011"]["protein_group_ambiguous"] == "true"
    service.run_intake(_config(tmp_path / "d", "diann"), tmp_path / "d", tmp_path / "d" / "out")
    drows = read_tsv(tmp_path / "d" / "out" / "features.tsv")
    assert all(dict(zip(drows[0], r))["is_decoy"] == "unknown" for r in drows[1:])   # profile has no decoy column: unknown


def test_v016_negative_peptide_grain_is_unsupported_scope(tmp_path):
    config = _config(tmp_path, "diann", matrix=V / "peptide_grain_export.tsv")
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == "E_UNSUPPORTED_SCOPE"


def test_v016_negative_unknown_header_needs_explicit_mapping(tmp_path):
    config = _config(tmp_path, "maxquant", matrix=V / "maxquant_unknown_header.tsv")
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == "E_MAPPING_HEADER" and error.value.field == "Normalized LFQ A1"


def test_v016_negative_version_is_never_inferred(tmp_path):
    config = _config(tmp_path, "maxquant", profile="maxquant_proteingroups@2.y-lfq-headers-1")
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == "E_MAPPING_VERSION"
    config = _config(tmp_path / "z", "maxquant", zero_as_missing=False)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path / "z", tmp_path / "z" / "out")
    assert error.value.code == "E_MAPPING_MISSING_ENCODING"
