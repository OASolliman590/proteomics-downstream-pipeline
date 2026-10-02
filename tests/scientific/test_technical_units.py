"""V013 biological observation hierarchy and V014 explicit technical aggregation (R02)."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "contract" / "intake"))
from intake_helpers import EXAMPLES, example_config, features_rows, minimal_config, observation_rows, read_tsv, write_tsv  # noqa: E402
from proteomics_pipeline.errors import ProteomicsError  # noqa: E402
from proteomics_pipeline.intake import service  # noqa: E402

AGG = lambda method: {"mode": "aggregate", "method": method, "biological_unit_column": "biological_unit_id"}


def _four_injections(tmp_path, technical):
    ids = ["s1_i1", "s1_i2", "s1_i3", "s1_i4", "s2", "s3", "s4"]
    obs = [("s1_i1", "spec1", "NA", "i1", "ctrl"), ("s1_i2", "spec1", "NA", "i2", "ctrl"), ("s1_i3", "spec1", "NA", "i3", "ctrl"), ("s1_i4", "spec1", "NA", "i4", "ctrl"),
           ("s2", "spec2", "NA", "NA", "ctrl"), ("s3", "spec3", "NA", "NA", "trt"), ("s4", "spec4", "NA", "NA", "trt")]
    matrix = write_tsv(tmp_path / "m.tsv", [["feature_id"] + ids, ["F1", "10", "10.2", "9.8", "10", "11", "12", "13"]])
    return minimal_config(tmp_path, matrix=matrix, observations=write_tsv(tmp_path / "o.tsv", observation_rows(obs)),
                          features=write_tsv(tmp_path / "f.tsv", features_rows(["F1"])), group_levels=["ctrl", "trt"], technical=technical)


def test_v013_four_injections_of_one_specimen_are_one_biological_unit(tmp_path):
    service.run_intake(_four_injections(tmp_path, AGG("mean_log2")), tmp_path, tmp_path / "out")
    report = json.loads((tmp_path / "out" / "intake_report.json").read_text())
    ctrl = next(g for g in report["groups"] if g["group"] == "ctrl")
    assert ctrl["n_injections"] == 5 and ctrl["n_biological_units"] == 2          # spec1 (4 injections) + spec2
    lineage = read_tsv(tmp_path / "out" / "observation_lineage.tsv")
    spec1 = [r for r in lineage[1:] if r[2] == "spec1"]
    assert len(spec1) == 4 and {r[1] for r in spec1} == {"spec1"} and all(r[4] == "4" for r in spec1)


def test_v013_paired_visits_keep_four_subjects_not_eight(tmp_path):
    config = example_config("example-paired.json")
    for key in ("matrix", "observations", "features", "source_provenance"):
        config["input"][key] = str((EXAMPLES / json.loads((EXAMPLES / "example-paired.json").read_text())["input"][key]).resolve())
    service.run_intake(config, EXAMPLES, tmp_path / "out")
    report = json.loads((tmp_path / "out" / "intake_report.json").read_text())
    assert report["n_subjects"] == 4 and report["n_canonical_observations"] == 8
    for group in report["groups"]:
        assert group["n_biological_units"] == 4 and group["n_subjects"] == 4


def test_v013_negative_unaggregated_injections_fail(tmp_path):
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(_four_injections(tmp_path, {"mode": "none", "biological_unit_column": "biological_unit_id"}), tmp_path, tmp_path / "out")
    assert error.value.code == "E_TECHNICAL_REPLICATION_UNMODELED" and not (tmp_path / "out").exists()


def _two_injections(tmp_path, method, *, second_group="before", second_unit="spec1"):
    ids = ["a_i1", "a_i2", "visit2", "other1", "other2"]
    obs = [("a_i1", "spec1", "subj1", "i1", "before"), ("a_i2", second_unit, "subj1", "i2", second_group), ("visit2", "spec1_after", "subj1", "NA", "after"),
           ("other1", "spec2", "subj2", "NA", "before"), ("other2", "spec2_after", "subj2", "NA", "after")]
    matrix = write_tsv(tmp_path / "m.tsv", [["feature_id"] + ids, ["F1", "4", "16", "64", "8", "8"]])
    return minimal_config(tmp_path, matrix=matrix, observations=write_tsv(tmp_path / "o.tsv", observation_rows(obs)),
                          features=write_tsv(tmp_path / "f.tsv", features_rows(["F1"])), group_levels=["before", "after"],
                          source_scale="linear_positive", transform="log2", technical=AGG(method))


@pytest.mark.parametrize("method,expected", [("mean_linear", math.log2(10)), ("mean_log2", 3.0)])
def test_v014_mean_linear_and_mean_log2_differ_as_specified(tmp_path, method, expected):
    service.run_intake(_two_injections(tmp_path, method), tmp_path, tmp_path / "out")
    rows = read_tsv(tmp_path / "out" / "matrix.tsv")
    values = dict(zip(rows[0][1:], rows[1][1:]))
    assert abs(float(values["spec1"]) - expected) <= 1e-10
    assert float(values["spec1_after"]) == 6.0 and "visit2" not in values        # distinct visit stays separate
    coverage = read_tsv(tmp_path / "out" / "aggregation_coverage.tsv")
    assert coverage[1][:4] == ["F1", "spec1", "2", "2"]
    report = json.loads((tmp_path / "out" / "intake_report.json").read_text())
    before = next(g for g in report["groups"] if g["group"] == "before")
    assert before["n_biological_units"] == 2 and before["n_injections"] == 3


def test_v014_negative_aggregation_across_visits_or_specimens_fails(tmp_path):
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(_two_injections(tmp_path, "mean_linear", second_group="after"), tmp_path, tmp_path / "out")
    assert error.value.code == "E_TECHNICAL_AGGREGATION_INVALID" and error.value.field == "group"
    config = _two_injections(tmp_path / "b", "mean_log2")
    obs = read_tsv(Path(config["input"]["observations"])); obs[2][2] = "subj9"
    write_tsv(Path(config["input"]["observations"]), obs)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path / "b", tmp_path / "b" / "out")
    assert error.value.code == "E_TECHNICAL_AGGREGATION_INVALID" and error.value.field == "subject_id"
