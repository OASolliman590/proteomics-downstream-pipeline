"""Canonical protein input service (R02: FR-011 to FR-020).

``run_intake`` reads one declared protein-level input, aligns it with
observation and feature metadata by key, applies the one-time scale
transition and any explicit technical aggregation, and publishes an
immutable canonical bundle (values, separate masks, metadata, lineage,
intake report and a SHA-256 manifest).  It performs no normalization, QC,
design or model fitting.
"""
from __future__ import annotations

import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from ..errors import CollisionError, IntegrityError, ProteomicsError
from ..provenance import canonical_json_bytes, content_sha256, sha256_bytes, sha256_file
from . import hierarchy as hier
from . import legacy_workbook, vendor
from .wide_long import (InputError, RawMatrix, format_number, parse_bool_mask, parse_long, parse_wide,
                        read_text_utf8, read_tsv, write_tsv)

BUNDLE_SCHEMA = "1.2.0"
FEATURE_COLUMNS = ("feature_id", "accessions", "gene_ids", "gene_symbols", "is_decoy", "is_contaminant", "protein_group_ambiguous")
FLAG_COLUMNS = ("is_decoy", "is_contaminant", "protein_group_ambiguous")
SUPPORTED_ASSAYS = ("lfq_dda", "lfq_dia", "tmt", "processed_protein")
FROM_MAPPING = "from_mapping"


def capabilities():
    return [{"id": "intake", "dependencies": [], "required_r_packages": []}]


# --------------------------------------------------------------------------- features
def _json_array(token: str, *, file: str, row: int, field: str) -> list[str]:
    if token in ("", "NA"):
        return []
    try:
        value = json.loads(token)
    except json.JSONDecodeError:
        raise InputError("E_FEATURE_ANNOTATION", f"{field} must be a JSON array of strings", file=file, row=row, field=field) from None
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise InputError("E_FEATURE_ANNOTATION", f"{field} must be a JSON array of strings", file=file, row=row, field=field)
    return value


def read_features(path: str) -> tuple[list[str], list[dict]]:
    header, body = read_tsv(path)
    for column in FEATURE_COLUMNS:
        if column not in header:
            raise InputError("E_METADATA_COLUMN", f"feature metadata is missing column {column!r}", file=path, row=1, field=column)
    rows: list[dict] = []
    seen: dict[str, int] = {}
    for number, raw in enumerate(body, start=2):
        row = dict(zip(header, raw))
        feature_id = row["feature_id"]
        if not feature_id:
            raise InputError("E_ID_MISSING", "feature_id is empty", file=path, row=number, field="feature_id")
        if feature_id in seen:
            raise InputError("E_ID_DUPLICATE", f"duplicate feature_id {feature_id!r} (first at row {seen[feature_id]}); differing annotations do not make it unique", file=path, row=number, field="feature_id")
        seen[feature_id] = number
        for column in ("accessions", "gene_ids", "gene_symbols"):
            row[column] = json.dumps(_json_array(row[column], file=path, row=number, field=column), ensure_ascii=False)
        for column in FLAG_COLUMNS:
            value = row[column]
            if value in ("", "NA"):
                value = "unknown"
            if value not in ("true", "false", "unknown"):
                raise InputError("E_FEATURE_FLAG", f"{column} must be true, false or unknown (got {value!r})", file=path, row=number, field=column)
            row[column] = value
        rows.append(row)
    return list(header), rows


# --------------------------------------------------------------------------- helpers
def _resolve(base: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    return path if path.is_absolute() else (base / path)


def _provenance(path: Path, assay: str) -> dict:
    try:
        value = json.loads(read_text_utf8(path))
    except json.JSONDecodeError as exc:
        raise InputError("E_PROVENANCE", f"source provenance is not valid JSON: {exc}", file=str(path)) from exc
    if not isinstance(value, dict):
        raise InputError("E_PROVENANCE", "source provenance must be a JSON object", file=str(path))
    grain = value.get("grain")
    if grain is not None and grain != "protein":
        raise InputError("E_UNSUPPORTED_SCOPE", f"source provenance declares grain {grain!r}; only protein-level input is supported", file=str(path), field="grain")
    declared = value.get("assay")
    if declared is not None and declared not in SUPPORTED_ASSAYS:
        raise InputError("E_UNSUPPORTED_SCOPE", f"source provenance declares unsupported assay {declared!r}", file=str(path), field="assay")
    if declared is not None and declared != assay:
        raise InputError("E_ASSAY_MISMATCH", f"source provenance assay {declared!r} differs from configured assay {assay!r}", file=str(path), field="assay")
    return value


def _align(kind: str, matrix_ids: list[str], metadata_ids: list[str], file: str) -> None:
    extra = [i for i in matrix_ids if i not in set(metadata_ids)]
    absent = [i for i in metadata_ids if i not in set(matrix_ids)]
    if extra or absent:
        raise InputError("E_ID_ALIGNMENT", f"{kind} keys do not align: in matrix only {extra[:5]}, in metadata only {absent[:5]}", file=file)


def source_inputs(config: dict, config_dir: Path) -> list[dict]:
    """List every source file the intake reads, with its current content SHA-256.

    D-42: the content hash normalises CRLF to LF for text files, so the same table checked out on Windows
    (CRLF) and on POSIX (LF) has one identity; any other byte change is still detected (E_SOURCE_CHANGED)."""
    spec = config["input"]
    roles = [("matrix", spec["matrix"]), ("observations", spec["observations"]), ("source_provenance", spec["source_provenance"])]
    if spec["features"] != FROM_MAPPING:
        roles.append(("features", spec["features"]))
    for role in ("original_observed_mask", "prior_imputed_mask", "mapping"):
        if spec.get(role):
            roles.append((role, spec[role]))
    items = []
    for role, value in roles:
        path = _resolve(config_dir, value)
        if not path.is_file():
            raise InputError("E_INPUT_READ", f"{role} file does not exist", file=str(path))
        items.append({"artifact_id": f"source.{role}", "path": str(path.resolve()), "sha256": content_sha256(path)})
    return items


# --------------------------------------------------------------------------- main intake
def run_intake(config: dict, config_dir: str | Path, output_dir: str | Path, *, expected_sources: list[dict] | None = None) -> dict:
    """Validate and publish the canonical bundle into a new ``output_dir``.

    ``expected_sources`` (from a plan or stage request) makes intake refuse
    source bytes that changed after they were recorded.
    """
    config_dir = Path(config_dir).resolve()
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise CollisionError(f"canonical bundle destination already exists: {output_dir}")
    spec = config["input"]
    preprocessing = config["preprocessing"]
    group_column = config["design"]["group_column"]
    technical = preprocessing["technical_replicates"]
    unit_column = technical["biological_unit_column"]
    missing_tokens = list(spec["missing_values"])
    scale = hier.scale_transition(config["source_scale"], preprocessing["transform"])
    sources = source_inputs(config, config_dir)
    if expected_sources is not None:
        expected = {item["artifact_id"]: item["sha256"] for item in expected_sources}
        for item in sources:
            if item["artifact_id"] in expected and expected[item["artifact_id"]] != item["sha256"]:
                raise IntegrityError(f"source bytes changed since they were recorded: {item['artifact_id']}", code="E_SOURCE_CHANGED")
    paths = {item["artifact_id"].removeprefix("source."): Path(item["path"]) for item in sources}
    provenance = _provenance(paths["source_provenance"], config["assay"])
    observations = hier.read_observations(str(paths["observations"]), group_column=group_column, unit_column=unit_column)
    source_columns = None
    if "source_column" in observations.header:
        source_columns = {row["source_column"]: row["observation_id"] for row in observations.rows if row["source_column"]}

    legacy = None
    fmt = spec["format"]
    if fmt == "mapped_protein":
        if not spec.get("mapping"):
            raise InputError("E_MAPPING_PROFILE", "format mapped_protein requires input.mapping (a versioned mapping profile)")
        profile = vendor.load_profile(paths["mapping"], spec["profile"])
        vendor.check_scale_declaration(profile, source_scale=config["source_scale"], zero_as_missing=spec["linear_zero_as_missing"], missing_tokens=missing_tokens)
        matrix, feature_rows = vendor.map_vendor_table(paths["matrix"], profile, missing_tokens=missing_tokens, source_columns=source_columns)
        feature_header = list(FEATURE_COLUMNS) + ["source_profile"]
        if spec["features"] != FROM_MAPPING:
            raise InputError("E_METADATA_COLUMN", f"mapped_protein derives features from the mapping; set input.features to {FROM_MAPPING!r}")
    elif fmt == "legacy_workbook":
        if not spec.get("mapping"):
            raise InputError("E_CROSSWALK", "format legacy_workbook requires input.mapping (the Method A/B semantic crosswalk)")
        crosswalk = legacy_workbook.load_crosswalk(paths["mapping"])
        matrix, feature_rows, resolved_contrasts, source_contrast_rows = legacy_workbook.import_workbook(paths["matrix"], crosswalk, missing_tokens=missing_tokens)
        legacy_workbook.check_observation_groups(crosswalk, observations, group_column)
        legacy = {"method": crosswalk["method"], "crosswalk_version": crosswalk["crosswalk_version"], "contrasts": resolved_contrasts, "rows": source_contrast_rows}
        feature_header = list(FEATURE_COLUMNS) + ["source_profile"]
        if spec["features"] != FROM_MAPPING:
            raise InputError("E_METADATA_COLUMN", f"legacy_workbook derives features from the workbook; set input.features to {FROM_MAPPING!r}")
    else:
        if spec["features"] == FROM_MAPPING:
            raise InputError("E_METADATA_COLUMN", "wide/long input requires a feature metadata table")
        matrix = parse_wide(paths["matrix"], missing_tokens) if fmt == "wide_tsv" else parse_long(paths["matrix"], missing_tokens)
        feature_header, feature_rows = read_features(str(paths["features"]))

    feature_by_id = {row["feature_id"]: row for row in feature_rows}
    _align("feature", matrix.feature_ids, list(feature_by_id), str(paths.get("features", paths["matrix"])))
    _align("observation", matrix.observation_ids, [row["observation_id"] for row in observations.rows], str(paths["observations"]))
    # canonical order: metadata order, so equivalent wide and long inputs give identical output
    feature_order = [row["feature_id"] for row in feature_rows]
    source_order = [row["observation_id"] for row in observations.rows]

    # one-time missing encoding on the source scale
    source_values: dict[tuple[str, str], float | None] = {}
    encoded_zero = 0
    for feature_id in feature_order:
        for observation_id in source_order:
            value, zero_missing = hier.convert_source_value(matrix.values[(feature_id, observation_id)], config["source_scale"], spec["linear_zero_as_missing"],
                                                            file=str(paths["matrix"]), row=matrix.source_rows.get(feature_id, 0), field=observation_id)
            encoded_zero += zero_missing
            source_values[(feature_id, observation_id)] = value

    observed, prior_imputed, mask_state = _masks(spec, paths, feature_order, source_order, source_values)
    hierarchy = hier.build_hierarchy(observations, source_order, mode=technical["mode"], method=technical.get("method"),
                                     unit_column=unit_column, group_column=group_column, file=str(paths["observations"]))
    canonical = hier.canonicalize(feature_order, hierarchy, source_values, observed, prior_imputed, scale)

    staging = output_dir.parent / f".{output_dir.name}.intake-staging-{os.getpid()}"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    try:
        artifacts = _write_bundle(staging, config=config, scale=scale, feature_order=feature_order, feature_header=feature_header,
                                  feature_by_id=feature_by_id, hierarchy=hierarchy, observations=observations, canonical=canonical,
                                  mask_state=mask_state, sources=sources, provenance=provenance, encoded_zero=encoded_zero,
                                  legacy=legacy, source_order=source_order)
        if output_dir.exists():
            output_dir.rmdir()  # only an empty pre-created stage directory reaches here
        os.replace(staging, output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return artifacts


def _masks(spec, paths, feature_order, source_order, source_values):
    numeric = {key: value is not None for key, value in source_values.items()}
    state = spec["prior_imputation"]
    if state == "none_documented":
        if spec.get("original_observed_mask"):
            supplied = parse_bool_mask(paths["original_observed_mask"], feature_order, source_order, "original_observed_mask")
            if supplied != numeric:
                raise InputError("E_MASK_INCONSISTENT", "with prior_imputation=none_documented the original-observed mask must equal numeric availability", file=str(paths["original_observed_mask"]))
        return dict(numeric), {key: False for key in numeric}, {"state": "known", "source": "derived_from_numeric_availability_none_documented"}
    if state == "mask_available":
        if not spec.get("original_observed_mask"):
            raise InputError("E_ORIGINAL_MASK_REQUIRED", "prior_imputation=mask_available requires input.original_observed_mask")
        observed = parse_bool_mask(paths["original_observed_mask"], feature_order, source_order, "original_observed_mask")
        for key, flag in observed.items():
            if flag and not numeric[key]:
                raise InputError("E_MASK_INCONSISTENT", f"cell {key!r} is marked originally observed but has no numeric value", file=str(paths["original_observed_mask"]), field=key[1])
        if spec.get("prior_imputed_mask"):
            imputed = parse_bool_mask(paths["prior_imputed_mask"], feature_order, source_order, "prior_imputed_mask")
            for key in imputed:
                if imputed[key] and (observed[key] or not numeric[key]):
                    raise InputError("E_MASK_INCONSISTENT", f"cell {key!r} cannot be both observed and prior-imputed, or imputed without a value", file=str(paths["prior_imputed_mask"]), field=key[1])
                if numeric[key] and not observed[key] and not imputed[key]:
                    raise InputError("E_MASK_INCONSISTENT", f"cell {key!r} has a value that is neither observed nor prior-imputed", file=str(paths["prior_imputed_mask"]), field=key[1])
        else:
            imputed = {key: numeric[key] and not observed[key] for key in numeric}
        return observed, imputed, {"state": "known", "source": "supplied_original_observed_mask"}
    return None, None, {"state": "unknown", "reason": "prior_imputation=unknown: genuine observation status cannot be distinguished from numeric availability",
                        "consequence": "E_ORIGINAL_MASK_REQUIRED for primary observed-coverage inference"}


def _mask_rows(feature_order, observation_ids, mask):
    for feature_id in feature_order:
        yield [feature_id] + ["true" if mask[(feature_id, o)] else "false" for o in observation_ids]


def _write_bundle(root: Path, *, config, scale, feature_order, feature_header, feature_by_id, hierarchy, observations, canonical,
                  mask_state, sources, provenance, encoded_zero, legacy, source_order) -> dict:
    units = hierarchy.units
    written: list[tuple[str, str, str]] = []

    def record(artifact_id: str, relative: str, result_type: str):
        written.append((artifact_id, relative, result_type))

    write_tsv(root / "matrix.tsv", ["feature_id"] + units, ([f] + [format_number(canonical.values[(f, o)]) for o in units] for f in feature_order))
    record("matrix", "matrix.tsv", "MatrixArtifact")
    write_tsv(root / "numeric_mask.tsv", ["feature_id"] + units, _mask_rows(feature_order, units, canonical.numeric))
    record("numeric_mask", "numeric_mask.tsv", "NumericAvailabilityMask")
    if canonical.observed is not None:
        write_tsv(root / "observed_mask.tsv", ["feature_id"] + units, _mask_rows(feature_order, units, canonical.observed))
        record("observed_mask", "observed_mask.tsv", "OriginalObservedMask")
        write_tsv(root / "prior_imputed_mask.tsv", ["feature_id"] + units, _mask_rows(feature_order, units, canonical.prior_imputed))
        record("prior_imputed_mask", "prior_imputed_mask.tsv", "PriorImputedMask")
    obs_header = [c for c in observations.header if c != "source_column"] + ["n_injections"]
    write_tsv(root / "observations.tsv", obs_header, ([hierarchy.canonical_rows[u].get(c) for c in obs_header] for u in units))
    record("observations", "observations.tsv", "ObservationMetadata")
    lineage_header = ["source_observation_id", "observation_id", "biological_unit_id", "aggregation", "n_injections"]
    write_tsv(root / "observation_lineage.tsv", lineage_header, ([row[c] for c in lineage_header] for row in hierarchy.lineage))
    record("observation_lineage", "observation_lineage.tsv", "ObservationHierarchy")
    if canonical.coverage:
        coverage_header = ["feature_id", "observation_id", "n_injections", "n_numeric", "n_observed", "n_contributing"]
        write_tsv(root / "aggregation_coverage.tsv", coverage_header, ([row[c] for c in coverage_header] for row in canonical.coverage))
        record("aggregation_coverage", "aggregation_coverage.tsv", "AggregationCoverage")
    write_tsv(root / "features.tsv", feature_header, ([feature_by_id[f].get(c, "") for c in feature_header] for f in feature_order))
    record("features", "features.tsv", "FeatureMetadata")
    if legacy is not None:
        legacy_header = ["contrast_id", "sheet", "historical_label", "numerator", "denominator", "feature_id", "source_log2fc", "source_p_value", "source_adjusted_p", "status"]
        write_tsv(root / "legacy_source_contrasts.tsv", legacy_header, ([row[c] for c in legacy_header] for row in legacy["rows"]))
        record("legacy_source_contrasts", "legacy_source_contrasts.tsv", "HistoricalSourceContrasts")

    # missingness tallies (numeric availability and, when known, genuine observation)
    per_observation = []
    for o in units:
        per_observation.append({"observation_id": o, "group": hierarchy.canonical_rows[o][hierarchy.group_column],
                                "n_features": len(feature_order),
                                "n_numeric_missing": sum(not canonical.numeric[(f, o)] for f in feature_order),
                                "n_observed": None if canonical.observed is None else sum(canonical.observed[(f, o)] for f in feature_order),
                                "n_prior_imputed": None if canonical.prior_imputed is None else sum(canonical.prior_imputed[(f, o)] for f in feature_order)})
    miss_header = list(per_observation[0])
    write_tsv(root / "missingness_by_observation.tsv", miss_header, ([row[c] for c in miss_header] for row in per_observation))
    record("missingness_by_observation", "missingness_by_observation.tsv", "MissingnessSummary")

    report = intake_report(config=config, scale=scale, feature_order=feature_order, hierarchy=hierarchy, canonical=canonical,
                           mask_state=mask_state, sources=sources, provenance=provenance, encoded_zero=encoded_zero,
                           n_source_observations=len(source_order), legacy=legacy)
    (root / "scale_decision.json").write_bytes(canonical_json_bytes(scale))
    record("scale_decision", "scale_decision.json", "ScaleDecision")
    (root / "intake_report.json").write_bytes((json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True) + "\n").encode("utf-8"))   # D-42: LF on every platform
    record("intake_report", "intake_report.json", "InputValidation")
    group_header = ["group", "n_injections", "n_biological_units", "n_subjects", "small_n"]
    write_tsv(root / "sample_n.tsv", group_header, ([row[c] if c != "small_n" else str(row[c]).lower() for c in group_header] for row in hierarchy.counts))
    record("sample_n", "sample_n.tsv", "ObservationHierarchy")

    artifacts = [{"artifact_id": a, "relative_path": r, "sha256": sha256_file(root / r), "result_type": t} for a, r, t in written]
    manifest = {
        "schema_version": BUNDLE_SCHEMA,
        "bundle": "canonical_protein_input",
        "artifacts": artifacts,
        "feature_ids": feature_order,
        "observation_ids": units,
        "feature_order_sha256": sha256_bytes(canonical_json_bytes(feature_order)),
        "observation_order_sha256": sha256_bytes(canonical_json_bytes(units)),
        "sources": [{"artifact_id": s["artifact_id"], "file_name": Path(s["path"]).name, "sha256": s["sha256"]} for s in sources],
        "original_observed_mask": mask_state,
        "scale": scale,
        "assay": config["assay"],
        "group_column": hierarchy.group_column,
        "unit_column": hierarchy.unit_column,
    }
    (root / "manifest.json").write_bytes(canonical_json_bytes(manifest) + b"\n")
    return {"manifest": manifest, "report": report}


def intake_report(*, config, scale, feature_order, hierarchy, canonical, mask_state, sources, provenance, encoded_zero, n_source_observations, legacy) -> dict:
    """InputValidation report: declarations and counts, never asserted unknown facts (FR-019)."""
    total = len(feature_order) * len(hierarchy.units)
    study = config["study"]
    gaps = []
    if study.get("tissue") in (None, ""):
        gaps.append({"field": "study.tissue", "state": "unknown", "consequence": "tissue identity is not asserted; tissue-dependent interpretation is unavailable"})
    if config["input"]["normalization_state"] == "unknown":
        gaps.append({"field": "input.normalization_state", "state": "unknown", "consequence": "upstream normalization is not asserted"})
    if mask_state["state"] == "unknown":
        gaps.append({"field": "input.original_observed_mask", "state": "unknown", "consequence": mask_state["consequence"]})
    if study.get("design_evidence") == "unknown":
        gaps.append({"field": "study.design_evidence", "state": "unknown", "consequence": "biological design facts are declarations, not confirmed evidence"})
    for key in ("upstream_software", "search_engine", "quantification_method", "acquisition"):
        if provenance.get(key) in (None, "", "unknown"):
            gaps.append({"field": f"source_provenance.{key}", "state": "unknown", "consequence": "upstream processing is not asserted"})
    warnings = [{"code": "W_SMALL_N", "message": f"group {row['group']!r} has {row['n_biological_units']} biological unit(s)"} for row in hierarchy.counts if row["small_n"]]
    return {
        "result_type": "InputValidation",
        "status": "published",
        "grain": "protein",
        "assay": config["assay"],
        "input_format": config["input"]["format"],
        "input_profile": config["input"]["profile"],
        "scale": scale,
        "explicit_zero_missing_cells": encoded_zero,
        "n_features": len(feature_order),
        "n_source_observations": n_source_observations,
        "n_canonical_observations": len(hierarchy.units),
        "n_subjects": hier.total_subjects(hierarchy),
        "technical_replicates": {"mode": hierarchy.mode, "method": hierarchy.method, "biological_unit_column": hierarchy.unit_column},
        "groups": hierarchy.counts,
        "missingness": {
            "total_cells": total,
            "numeric_missing_cells": sum(not v for v in canonical.numeric.values()),
            "originally_observed_cells": None if canonical.observed is None else sum(canonical.observed.values()),
            "prior_imputed_cells": None if canonical.prior_imputed is None else sum(canonical.prior_imputed.values()),
            "original_observed_mask": mask_state,
        },
        "declared": {"study_id": study["id"], "experimental_unit": study["experimental_unit"], "design_evidence": study["design_evidence"],
                     "tissue": study.get("tissue"), "normalization_state": config["input"]["normalization_state"],
                     "prior_imputation": config["input"]["prior_imputation"], "organism_taxonomy_id": config["organism"]["taxonomy_id"]},
        "provenance_gaps": gaps,
        "asserted_unknown_facts": [],
        "legacy_crosswalk": None if legacy is None else {"method": legacy["method"], "crosswalk_version": legacy["crosswalk_version"], "contrasts": legacy["contrasts"],
                                                          "source_statistics_status": "historical_source_value_not_refit"},
        "sources": [{"role": s["artifact_id"], "file_name": Path(s["path"]).name, "sha256": s["sha256"]} for s in sources],
        "warnings": warnings,
        "errors": [],
    }


# --------------------------------------------------------------------------- verification
def _keys_of(path: Path, kind: str) -> tuple[list[str], list[str]]:
    header, body = read_tsv(path)
    return header[1:], [row[0] for row in body]


def verify_bundle(bundle_dir: str | Path) -> dict:
    """Recompute every hash and cross-check matrix/mask/metadata keys.

    Raises IntegrityError (exit 5) on any mismatch so no downstream fit can
    start from a modified bundle.
    """
    root = Path(bundle_dir)
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise IntegrityError("canonical bundle manifest is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for artifact in manifest["artifacts"]:
        path = root / artifact["relative_path"]
        if not path.is_file():
            raise IntegrityError(f"canonical artifact is missing: {artifact['relative_path']}")
        if sha256_file(path) != artifact["sha256"]:
            raise IntegrityError(f"canonical artifact hash mismatch: {artifact['relative_path']}")
    ids = {a["artifact_id"] for a in manifest["artifacts"]}
    required = {"matrix", "numeric_mask", "observations", "features", "scale_decision", "intake_report"}
    if manifest["original_observed_mask"]["state"] == "known":
        required |= {"observed_mask", "prior_imputed_mask"}
    if not required <= ids:
        raise IntegrityError(f"canonical bundle lacks required artifacts {sorted(required - ids)}")
    features, observations = manifest["feature_ids"], manifest["observation_ids"]
    for artifact in ("matrix", "numeric_mask", "observed_mask", "prior_imputed_mask"):
        if artifact not in ids:
            continue
        columns, rows = _keys_of(root / f"{artifact}.tsv", artifact)
        if columns != observations or rows != features:
            raise IntegrityError(f"{artifact} keys do not match the manifest order")
    _, observation_rows = _keys_of(root / "observations.tsv", "observations")
    _, feature_rows = _keys_of(root / "features.tsv", "features")
    if observation_rows != observations or feature_rows != features:
        raise IntegrityError("metadata keys do not match the manifest order")
    return manifest


# --------------------------------------------------------------------------- stage protocol
def execute(request: dict) -> dict:
    """Stage-protocol handler: publish the bundle into ``output_temp_dir``."""
    from ..config import load_config
    from ..runtime import stage_result, validate_stage_request

    validate_stage_request(request)
    config_path = request["config_path"]
    output_dir = Path(request["output_temp_dir"])
    started = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    try:
        config = load_config(config_path)
        built = run_intake(config, Path(config_path).parent, output_dir, expected_sources=request["inputs"])
    except ProteomicsError as error:
        result = stage_result(request["run_id"], request["stage_id"], "intake", "FAILED", plan_hash=request["plan_hash"],
                              exit_code=error.exit_code, reason_code=error.code, message=error.message)
        result["started_at"] = started
        return result
    result = stage_result(request["run_id"], request["stage_id"], "intake", "COMPLETED", plan_hash=request["plan_hash"], exit_code=0,
                          message=f"canonical bundle: {len(built['manifest']['feature_ids'])} features x {len(built['manifest']['observation_ids'])} observations")
    result["started_at"] = started
    result["outputs"] = [dict(a) for a in built["manifest"]["artifacts"]] + [{"artifact_id": "manifest", "relative_path": "manifest.json",
                                                                               "sha256": sha256_file(output_dir / "manifest.json"), "result_type": "CanonicalManifest"}]
    result["warnings"] = [{"code": w["code"], "stage_id": request["stage_id"], "severity": "warning", "message": w["message"], "artifact_ref": "intake_report.json"}
                          for w in built["report"]["warnings"]]
    return result
