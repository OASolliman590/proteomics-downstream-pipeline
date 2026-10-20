from __future__ import annotations
import json, math
from pathlib import Path
from typing import Any
from .errors import ConfigurationError
_DEFAULTS = {"primary_engine":"limma","primary_hypothesis":"zero_null","response_mode":"descriptive_only","score_test":"off"}
def _pointer(parts) -> str:
    escaped = [str(p).replace("~","~0").replace("/","~1") for p in parts]
    return "/" + "/".join(escaped) if escaped else ""
def _read(path: Path) -> dict[str, Any]:
    try: text = path.read_text(encoding="utf-8")
    except OSError as exc: raise ConfigurationError("E_CONFIG_READ", str(exc)) from exc
    try:
        if path.suffix.lower() in {".yaml",".yml"}:
            try: import yaml
            except ImportError as exc: raise ConfigurationError("E_CONFIG_PARSE","PyYAML is required for YAML configuration") from exc
            value = yaml.safe_load(text)
        else: value = json.loads(text)
    except ConfigurationError:
        raise
    except Exception as exc: raise ConfigurationError("E_CONFIG_PARSE", str(exc)) from exc
    if not isinstance(value, dict): raise ConfigurationError("E_CONFIG_SCHEMA","configuration root must be an object","")
    _reject_nonfinite(value)
    return value
def _reject_nonfinite(value: Any, parts=()):
    if isinstance(value,float) and not math.isfinite(value):
        raise ConfigurationError("E_CONFIG_SCHEMA","non-finite numeric values are not allowed",_pointer(parts))
    if isinstance(value,dict):
        for key,item in value.items(): _reject_nonfinite(item,parts+(key,))
    elif isinstance(value,list):
        for index,item in enumerate(value): _reject_nonfinite(item,parts+(index,))
def _schema(): return json.loads((Path(__file__).parent/"schemas"/"analysis.schema.json").read_text(encoding="utf-8"))
def _validate_schema(value):
    try:
        from jsonschema import Draft202012Validator, FormatChecker
    except ImportError as exc: raise ConfigurationError("E_CONFIG_VALIDATOR","jsonschema is required for configuration validation",exit_code=3) from exc
    errors = sorted(Draft202012Validator(_schema(), format_checker=FormatChecker()).iter_errors(value), key=lambda e:list(e.absolute_path))
    if errors:
        error=errors[0]; raise ConfigurationError("E_CONFIG_SCHEMA", error.message, _pointer(error.absolute_path))
def _unique(items, path, key="id"):
    seen={}
    for index,item in enumerate(items):
        ident=item.get(key)
        if ident in seen: raise ConfigurationError("E_ID_DUPLICATE",f"duplicate id: {ident}",f"{path}/{index}/{key}")
        seen[ident]=index
def _semantic(value):
    _unique([value["design"], *value.get("additional_designs",[])],"/designs"); _unique(value["models"],"/models"); _unique(value["contrasts"],"/contrasts"); _unique(value.get("multiplicity_families",[]),"/multiplicity_families"); _unique(value.get("resources",[]),"/resources"); _unique(value.get("omnibus_tests",[]),"/omnibus_tests")
    design_ids={value["design"]["id"],*(d["id"] for d in value.get("additional_designs",[]))}; model_ids={m["id"] for m in value["models"]}; contrast_ids={c["id"] for c in value["contrasts"]}; resource_ids={r["id"] for r in value.get("resources",[])}
    primary_models=[m for m in value["models"] if m.get("role")=="primary"]
    if len(primary_models)!=1: raise ConfigurationError("E_PRIMARY_MODEL_CONFLICT","exactly one model must have role=primary","/models")
    primary=primary_models[0]; pindex=value["models"].index(primary)
    if primary["design_id"]!=value["design"]["id"]: raise ConfigurationError("E_PRIMARY_MODEL_CONFLICT","primary model design does not match design.id",f"/models/{pindex}/design_id")
    if primary["engine"]!=value["primary_engine"]: raise ConfigurationError("E_PRIMARY_MODEL_CONFLICT","primary model engine does not match primary_engine",f"/models/{pindex}/engine")
    if primary["hypothesis"]!=value["primary_hypothesis"]: raise ConfigurationError("E_PRIMARY_MODEL_CONFLICT","primary model hypothesis does not match primary_hypothesis",f"/models/{pindex}/hypothesis")
    if primary.get("execution_requirement")!="required": raise ConfigurationError("E_PRIMARY_MODEL_CONFLICT","primary model must be required",f"/models/{pindex}/execution_requirement")
    # Amendment A-2026-10-01-44 (ADR 0011, R15a): welch exists only when declared under methods_extension.welch (FR-169).
    welch_block=(value.get("methods_extension") or {}).get("welch") or {}
    welch_models=[(index,m) for index,m in enumerate(value["models"]) if m["engine"]=="welch"]
    if (welch_models or value["primary_engine"]=="welch") and not welch_block.get("enabled"):
        raise ConfigurationError("E_WELCH_NOT_DECLARED","welch is requested but methods_extension.welch.enabled is not true; no welch result is produced and no limma result is relabelled",
                                 "/methods_extension/welch/enabled")
    pathway_block=(value.get("methods_extension") or {}).get("pathways") or {}
    if pathway_block.get("enabled"): _unique(pathway_block.get("libraries",[]),"/methods_extension/pathways/libraries")   # A-2026-10-01-46 (R15b)
    if value["primary_engine"]=="welch" and any(m["engine"]=="limma" for m in value["models"]):
        # Amendment A-2026-10-01-46: a welch primary writes dea/, which a limma model would also write; the two cannot share a run.
        raise ConfigurationError("E_PRIMARY_MODEL_CONFLICT","a welch primary writes the dea/ contract and cannot share a run with limma models","/models")
    if welch_block.get("enabled") and not welch_models:
        raise ConfigurationError("E_WELCH_NOT_DECLARED","methods_extension.welch is enabled but no model declares engine welch","/methods_extension/welch")
    for index,model in welch_models:
        # A primary Welch model maps to role primary; a secondary Welch model is a non-primary model (role sensitivity) with its own family.
        expected="primary" if model["role"]=="primary" else "secondary"
        if welch_block.get("role","secondary")!=expected: raise ConfigurationError("E_PRIMARY_MODEL_CONFLICT","methods_extension.welch.role must be primary for the primary welch model and secondary otherwise","/methods_extension/welch/role")
    for index,contrast in enumerate(value["contrasts"]):
        if contrast["design_id"] not in design_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","contrast design does not exist",f"/contrasts/{index}/design_id")
    for index,model in enumerate(value["models"]):
        if model["design_id"] not in design_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","model design does not exist",f"/models/{index}/design_id")
    for index,sensitivity in enumerate(value["preprocessing"].get("sensitivities",[])):
        if sensitivity["model_id"] not in model_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","sensitivity model does not exist",f"/preprocessing/sensitivities/{index}/model_id")
        # Audit 2026-10-02 (A-2026-10-01-10): an imputed sensitivity matrix may never feed the primary model.
        if sensitivity["model_id"]==primary["id"]: raise ConfigurationError("E_SENSITIVITY_PRIMARY_MODEL","a sensitivity matrix cannot be assigned to the primary model; declare a separate role=sensitivity model",f"/preprocessing/sensitivities/{index}/model_id")
        # Amendment A-2026-10-01-46 (ADR 0011, F4): the welch engine reads the observed primary matrix only; a sensitivity matrix on a welch model would be ignored silently.
        if any(m["id"]==sensitivity["model_id"] and m["engine"]=="welch" for m in value["models"]): raise ConfigurationError("E_WELCH_MATRIX_UNSUPPORTED","the welch engine reads the observed primary matrix only; a sensitivity matrix cannot be assigned to a welch model",f"/preprocessing/sensitivities/{index}/model_id")
    for index,omnibus in enumerate(value.get("omnibus_tests",[])):
        if omnibus["model_id"] not in model_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","omnibus model does not exist",f"/omnibus_tests/{index}/model_id")
        if omnibus.get("reduced_design_id") is not None and omnibus["reduced_design_id"] not in design_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","omnibus reduced design does not exist",f"/omnibus_tests/{index}/reduced_design_id")
    pathways=value.get("pathways",{})
    for key in ("mapping_resource_id",):
        if pathways.get(key) is not None and pathways[key] not in resource_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","pathways resource does not exist",f"/pathways/{key}")
    if pathways.get("gene_set_resource_ids"):
        for resource_index,resource_id in enumerate(pathways["gene_set_resource_ids"]):
            if resource_id not in resource_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","pathways gene-set resource does not exist",f"/pathways/gene_set_resource_ids/{resource_index}")
    if pathways.get("linear_sensitivity_model_id") is not None and pathways["linear_sensitivity_model_id"] not in model_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","pathways sensitivity model does not exist","/pathways/linear_sensitivity_model_id")
    response=value.get("response",{}); axes=response.get("axes",[]); scores=response.get("scores",[]); axis_ids={a["id"] for a in axes}; score_ids={s["id"] for s in scores}; _unique(axes,"/response/axes"); _unique(scores,"/response/scores")
    for index,axis in enumerate(axes):
        if axis["model_id"] not in model_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","response axis model does not exist",f"/response/axes/{index}/model_id")
        for key in ("disease_contrast","treatment_contrast","residual_contrast"):
            if axis[key] not in contrast_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","response axis contrast does not exist",f"/response/axes/{index}/{key}")
    if response.get("direction_resource_id") is not None and response["direction_resource_id"] not in resource_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","response direction resource does not exist","/response/direction_resource_id")
    for index,score in enumerate(scores):
        if score["model_id"] not in model_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","response score model does not exist",f"/response/scores/{index}/model_id")
        if score["contrast_id"] not in contrast_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","response score contrast does not exist",f"/response/scores/{index}/contrast_id")
        if score["resource_id"] not in resource_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN","response score resource does not exist",f"/response/scores/{index}/resource_id")
    omnibus={x["id"] for x in value.get("omnibus_tests",[])}
    for index,family in enumerate(value["multiplicity_families"]):
        for model in family["model_ids"]:
            if model not in model_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN",f"unknown model id: {model}",f"/multiplicity_families/{index}/model_ids")
        for contrast in family["contrast_ids"]:
            if contrast not in contrast_ids and contrast not in omnibus: raise ConfigurationError("E_REFERENCE_UNKNOWN",f"unknown contrast id: {contrast}",f"/multiplicity_families/{index}/contrast_ids")
        for axis in family.get("axis_ids",[]):
            if axis not in axis_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN",f"unknown axis id: {axis}",f"/multiplicity_families/{index}/axis_ids")
        for score in family.get("score_ids",[]):
            if score not in score_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN",f"unknown score id: {score}",f"/multiplicity_families/{index}/score_ids")
        for collection in family.get("collection_ids",[]):
            if collection not in resource_ids: raise ConfigurationError("E_REFERENCE_UNKNOWN",f"unknown collection/resource id: {collection}",f"/multiplicity_families/{index}/collection_ids")
def load_config(path: str | Path) -> dict[str, Any]:
    value=_read(Path(path))
    for key,default in _DEFAULTS.items(): value.setdefault(key,default)
    _validate_schema(value); _semantic(value); return value
def defaults(): return dict(_DEFAULTS)
