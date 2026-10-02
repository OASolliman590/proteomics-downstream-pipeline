"""V018 (R side) — downstream R stages read only verified canonical bundles."""
from __future__ import annotations

import json
import shutil
import subprocess

import pytest

from intake_helpers import EXAMPLES, example_config, features_rows, minimal_config, observation_rows, write_tsv
from proteomics_pipeline.intake import service


def _rscript():
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("NOT_RUN: Rscript unavailable")
    probe = subprocess.run([rscript, "--vanilla", "-e", "quit(status=if(exists('read_canonical_bundle', envir=asNamespace('proteomicsCore'))) 0L else 1L)"], capture_output=True)
    if probe.returncode != 0:
        pytest.skip("NOT_RUN: installed proteomicsCore lacks the R02 bundle reader")
    return rscript


def _read(rscript, bundle):
    code = ("b <- proteomicsCore:::read_canonical_bundle(commandArgs(TRUE)[1]); "
            "cat(jsonlite::toJSON(list(ids=I(rownames(b$values)), obs=I(colnames(b$values)), na=which(is.na(b$values)), v=b$values[1,1], mask=b$mask_state, observed=sum(b$observed)), auto_unbox=TRUE, digits=NA))")
    return subprocess.run([rscript, "--vanilla", "-e", code, str(bundle)], capture_output=True, text=True)


def test_v018_r_reader_returns_verified_values(tmp_path):
    rscript = _rscript()
    service.run_intake(example_config(), EXAMPLES, tmp_path / "bundle")
    result = _read(rscript, tmp_path / "bundle")
    assert result.returncode == 0, result.stderr
    value = json.loads(result.stdout)
    assert value["ids"] == [f"P0{i}" for i in range(1, 9)] and value["obs"][0] == "C1" and value["v"] == 9.7
    assert len(value["na"]) == 2 and value["observed"] == 94 and value["mask"] == "known"


def test_v018_r_reader_handles_quoted_tab_and_unicode_ids(tmp_path):
    rscript = _rscript()
    odd = "grp\tΩ \"x\""
    ids = ["A1", "A2", "B1", "B2"]
    matrix = write_tsv(tmp_path / "m.tsv", [["feature_id"] + ids, [odd, "1.25", "NA", "3", "4"]])
    obs = write_tsv(tmp_path / "o.tsv", observation_rows([(i, i, "NA", "NA", "a" if i[0] == "A" else "b") for i in ids]))
    feats = features_rows([odd])
    config = minimal_config(tmp_path, matrix=matrix, observations=obs, features=write_tsv(tmp_path / "f.tsv", feats), group_levels=["a", "b"])
    service.run_intake(config, tmp_path, tmp_path / "bundle")
    result = _read(rscript, tmp_path / "bundle")
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["ids"] == [odd]


@pytest.mark.parametrize("mutation", ["delete_mask", "edit_value"])
def test_v018_negative_r_reader_refuses_mutated_bundle(tmp_path, mutation):
    rscript = _rscript()
    service.run_intake(example_config(), EXAMPLES, tmp_path / "bundle")
    if mutation == "delete_mask":
        (tmp_path / "bundle" / "observed_mask.tsv").unlink()
    else:
        path = tmp_path / "bundle" / "matrix.tsv"
        path.write_text(path.read_text().replace("9.7", "9.9", 1))
    result = _read(rscript, tmp_path / "bundle")
    assert result.returncode != 0 and "E_INTEGRITY" in result.stderr
