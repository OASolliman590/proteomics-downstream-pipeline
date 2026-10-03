"""V017 — historical workbook importer wraps the recovered parser read-only."""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json

import pytest

from intake_helpers import FIXTURES, MAPPINGS, ROOT, minimal_config, read_tsv, write_tsv
from proteomics_pipeline.errors import ProteomicsError
from proteomics_pipeline.intake import service

spec = importlib.util.spec_from_file_location("make_workbook", FIXTURES / "make_workbook.py")
make_workbook = importlib.util.module_from_spec(spec); spec.loader.exec_module(make_workbook)

A_NAMES = {"CTL": "CTRL", "DM": "DM", "PDM": "PDM", "DM-Treated": "DM_Treated", "PDM-Treated": "PDM_Treated"}
B_NAMES = {"CTL": "CTRL", "PDM": "PreDM", "DM": "DM", "PDM-Treated": "PreDM_Dapa", "DM-Treated": "DM_Dapa"}
SAMPLE_PREFIX = {"CTL": "CTRL", "DM": "DM", "PDM": "PDM", "DM-Treated": "DM_Treated", "PDM-Treated": "PDM_Treated"}


def _observations(tmp_path, names):
    rows = [("EXT_%s_%d" % (SAMPLE_PREFIX[label], r), "EXT_%s_%d" % (SAMPLE_PREFIX[label], r), "NA", "NA", names[label]) for label in SAMPLE_PREFIX for r in range(1, 5)]
    return write_tsv(tmp_path / "obs.tsv", [["observation_id", "biological_unit_id", "subject_id", "technical_replicate_id", "group"]] + [list(r) for r in rows])


def _config(tmp_path, crosswalk, names):
    tmp_path.mkdir(parents=True, exist_ok=True)
    workbook = tmp_path / "synthetic.xlsx"
    expected = make_workbook.write_workbook(workbook)
    cw_path = tmp_path / "crosswalk.json"; cw_path.write_text(json.dumps(crosswalk), encoding="utf-8")
    config = minimal_config(tmp_path, matrix=workbook, observations=_observations(tmp_path, names), features="from_mapping", fmt="legacy_workbook",
                            profile="historical_workbook", mapping=cw_path, group_levels=sorted(set(names.values())))
    return config, expected


def _crosswalk(name):
    return json.loads((MAPPINGS / name).read_text(encoding="utf-8"))


@pytest.mark.parametrize("name,names", [("legacy-crosswalk-method-a.json", A_NAMES), ("legacy-crosswalk-method-b.json", B_NAMES)])
def test_v017_wraps_parser_and_keeps_sheet_identity(tmp_path, name, names):
    parser_hash = hashlib.sha256((ROOT / "pipeline" / "20_external_limma_intake.py").read_bytes()).hexdigest()
    config, expected = _config(tmp_path, _crosswalk(name), names)
    service.run_intake(config, tmp_path, tmp_path / "out")
    assert hashlib.sha256((ROOT / "pipeline" / "20_external_limma_intake.py").read_bytes()).hexdigest() == parser_hash
    rows = read_tsv(tmp_path / "out" / "matrix.tsv")
    header = rows[0][1:]
    assert len(header) == 20 and [r[0] for r in rows[1:]] == ["A1", "A2", "A3"]
    for row in rows[1:]:
        for sample_id, token in zip(header, row[1:]):
            group, rep = sample_id.removeprefix("EXT_").rsplit("_", 1)
            label = {v: k for k, v in SAMPLE_PREFIX.items()}[group]
            want = expected[(row[0], f"{label} {rep}")]
            assert (token == "NA" and want is None) or float(token) == want
    source = read_tsv(tmp_path / "out" / "legacy_source_contrasts.tsv")
    sheets = {r[1] for r in source[1:]}
    assert sheets == {"DM_vs_CTL", "PDM_vs_CTL", "DM-Treated_vs_DM", "PDM-Treated_vs_PDM"}
    dm = next(dict(zip(source[0], r)) for r in source[1:] if r[1] == "DM_vs_CTL")
    assert (dm["numerator"], dm["denominator"]) == (names["DM"], "CTRL") and dm["status"] == "historical_source_value_not_refit"
    report = json.loads((tmp_path / "out" / "intake_report.json").read_text(encoding="utf-8"))
    assert report["legacy_crosswalk"]["method"] == _crosswalk(name)["method"]


def test_v017_negative_swapped_d1_d2_without_semantic_update_fails(tmp_path):
    crosswalk = _crosswalk("legacy-crosswalk-method-a.json")
    swapped = copy.deepcopy(crosswalk)
    swapped["contrasts"][0]["sheet"], swapped["contrasts"][1]["sheet"] = crosswalk["contrasts"][1]["sheet"], crosswalk["contrasts"][0]["sheet"]
    config, _ = _config(tmp_path, swapped, A_NAMES)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == "E_CONTRAST_IDENTITY"
    t_swap = copy.deepcopy(crosswalk)
    t_swap["contrasts"][2]["numerator"], t_swap["contrasts"][3]["numerator"] = crosswalk["contrasts"][3]["numerator"], crosswalk["contrasts"][2]["numerator"]
    config, _ = _config(tmp_path / "t", t_swap, A_NAMES)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path / "t", tmp_path / "t" / "out")
    assert error.value.code == "E_CONTRAST_IDENTITY"


def test_v017_negative_method_names_are_not_interchangeable(tmp_path):
    mixed = _crosswalk("legacy-crosswalk-method-a.json"); mixed["group_names"] = dict(B_NAMES)
    config, _ = _config(tmp_path, mixed, B_NAMES)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path, tmp_path / "out")
    assert error.value.code == "E_CROSSWALK_METHOD"
    config, _ = _config(tmp_path / "g", _crosswalk("legacy-crosswalk-method-b.json"), A_NAMES)
    with pytest.raises(ProteomicsError) as error:
        service.run_intake(config, tmp_path / "g", tmp_path / "g" / "out")
    assert error.value.code == "E_CROSSWALK_GROUP"
