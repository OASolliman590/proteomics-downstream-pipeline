import json
from copy import deepcopy
from pathlib import Path
import pytest
from proteomics_pipeline.config import load_config
from proteomics_pipeline.errors import ConfigurationError
ROOT=Path(__file__).resolve().parents[3]
def test_loader_rejects_duplicate_ids_and_conflicting_primary(tmp_path):
    source=json.loads((ROOT/"configs/examples/example-independent.json").read_text()); source["primary_engine"]="proda"; path=tmp_path/"bad.json"; path.write_text(json.dumps(source))
    with pytest.raises(ConfigurationError) as exc:load_config(path)
    assert exc.value.code in {"E_CONFIG_SCHEMA", "E_PRIMARY_MODEL_CONFLICT"}
def test_loader_rejects_nan(tmp_path):
    path=tmp_path/"nan.json"; path.write_text('{"schema_version":"1.2.0","value":NaN}')
    with pytest.raises(ConfigurationError) as exc:load_config(path)
    assert exc.value.code=="E_CONFIG_SCHEMA" and exc.value.pointer=="/value"

@pytest.mark.parametrize("payload,pointer",[("{\"schema_version\":\"1.2.0\",\"nested\":{\"value\":Infinity}}","/nested/value"),("{\"schema_version\":\"1.2.0\",\"nested\":[0,-Infinity]}","/nested/1")])
def test_loader_rejects_nested_nonfinite_json_with_pointer(tmp_path,payload,pointer):
    path=tmp_path/"nonfinite.json"; path.write_text(payload)
    with pytest.raises(ConfigurationError) as exc:load_config(path)
    assert exc.value.code=="E_CONFIG_SCHEMA" and exc.value.pointer==pointer

def test_loader_rejects_nonfinite_yaml_with_pointer(tmp_path):
    path=tmp_path/"nonfinite.yaml"; path.write_text("schema_version: '1.2.0'\nnested:\n  value: .nan\n")
    with pytest.raises(ConfigurationError) as exc:load_config(path)
    assert exc.value.code=="E_CONFIG_SCHEMA" and exc.value.pointer=="/nested/value"
def test_only_four_defaults_are_inserted(tmp_path):
    source=json.loads((ROOT/"configs/examples/example-independent.json").read_text())
    for key in ("primary_engine","primary_hypothesis","response_mode","score_test"):source.pop(key)
    path=tmp_path/"defaults.json"; path.write_text(json.dumps(source)); result=load_config(path)
    assert result["primary_engine"]=="limma"; assert result["primary_hypothesis"]=="zero_null"; assert result["response_mode"]=="descriptive_only"; assert result["score_test"]=="off"

@pytest.mark.parametrize("mutator", [
    lambda value: value.update(schema_version="1.0.0"),
    lambda value: value["contrasts"][0].update(design_id="missing-design"),
    lambda value: value["multiplicity_families"][0]["model_ids"].append("missing-model"),
    lambda value: value["contrasts"].append(value["contrasts"][0].copy()),
])
def test_frozen_config_negatives_are_typed(tmp_path, mutator):
    source=json.loads((ROOT/"configs/examples/example-independent.json").read_text())
    mutator(source)
    path=tmp_path/"bad.json"; path.write_text(json.dumps(source))
    with pytest.raises(ConfigurationError) as exc: load_config(path)
    assert exc.value.code in {"E_CONFIG_SCHEMA","E_REFERENCE_UNKNOWN","E_ID_DUPLICATE"}
    assert exc.value.pointer is not None

def test_invalid_explicit_score_test_is_rejected_not_defaulted(tmp_path):
    source=json.loads((ROOT/"configs/examples/example-independent.json").read_text()); source["score_test"]="bogus"
    path=tmp_path/"bad.json"; path.write_text(json.dumps(source))
    with pytest.raises(ConfigurationError) as exc: load_config(path)
    assert exc.value.code=="E_CONFIG_SCHEMA"

def test_yaml_executable_tag_is_rejected(tmp_path):
    path=tmp_path/"bad.yaml"; path.write_text("!!python/object/apply:os.system ['echo unsafe']")
    with pytest.raises(ConfigurationError) as exc: load_config(path)
    assert exc.value.code=="E_CONFIG_PARSE"

@pytest.mark.parametrize("mutator",[
    lambda value: value["preprocessing"]["sensitivities"].append({"id":"sensitivity-1","method":"complete_case","model_id":"missing-model"}),
    lambda value: value.setdefault("omnibus_tests",[]).append({"id":"omnibus-1","model_id":"missing-model","coefficient_names":["x"]}),
    lambda value: value["multiplicity_families"][0].update(axis_ids=["missing-axis"]),
    lambda value: value["multiplicity_families"][0].update(collection_ids=["missing-collection"]),
    lambda value: value["pathways"].update(mapping_resource_id="missing-resource"),
    lambda value: value["models"].append({**deepcopy(value["models"][0]),"id":"sensitivity-model","role":"sensitivity","execution_requirement":"optional","design_id":"missing-design"}),
])
def test_schema_valid_relational_references_are_typed(tmp_path,mutator):
    source=json.loads((ROOT/"configs/examples/example-independent.json").read_text()); mutator(source)
    path=tmp_path/"bad-relation.json"; path.write_text(json.dumps(source))
    with pytest.raises(ConfigurationError) as exc: load_config(path)
    assert exc.value.code=="E_REFERENCE_UNKNOWN" and exc.value.pointer is not None
