"""Observation hierarchy, scale transitions and explicit technical aggregation.

Implements SM02 (one-time scale transitions) and SM03 (biological versus
technical units) for FR-012, FR-013 and FR-014.  Technical injections never
increase biological n; paired visits share a subject but stay distinct
specimens; aggregation only merges injections of the same specimen.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from statistics import fmean

from .wide_long import InputError, read_tsv

NULL_TOKENS = ("NA", "")
REQUIRED_OBSERVATION_COLUMNS = ("observation_id", "biological_unit_id", "subject_id", "technical_replicate_id")
# Columns allowed to differ between technical injections of one specimen.
TECHNICAL_COLUMNS = frozenset({"observation_id", "technical_replicate_id"})

# SM02 scale-state table: (source_scale, transform) -> (output_scale, inference_reason or None)
_TRANSITIONS = {
    ("linear_positive", "log2"): ("log2", None),
    ("linear_positive", "none"): ("linear_positive", "E_LOG2_REQUIRED"),
    ("log2", "none"): ("log2", None),
    ("ratio_log2", "none"): ("ratio_log2", None),
    ("standardized_unknown", "none"): ("standardized_unknown", "E_SCALE_UNKNOWN_INFERENCE"),
}
_REJECTED = {
    ("log2", "log2"): ("E_SCALE_SECOND_LOG", "source is already log2; a second logarithm is forbidden"),
    ("ratio_log2", "log2"): ("E_SCALE_SECOND_LOG", "source is already a log2 ratio; a second logarithm is forbidden"),
    ("standardized_unknown", "log2"): ("E_SCALE_TRANSFORM_UNJUSTIFIED", "a log2 transform of an unknown standardized scale is not justified"),
}


def scale_transition(source_scale: str, transform: str) -> dict:
    """Return the frozen ScaleDecision or raise before any value is transformed."""
    key = (source_scale, transform)
    if key in _REJECTED:
        code, message = _REJECTED[key]
        raise InputError(code, message)
    if key not in _TRANSITIONS:
        raise InputError("E_SCALE_UNKNOWN", f"unsupported scale transition {source_scale!r} -> {transform!r}")
    output_scale, reason = _TRANSITIONS[key]
    return {
        "source_scale": source_scale,
        "requested_transform": transform,
        "applied_transform": "log2" if transform == "log2" else "none",
        "output_scale": output_scale,
        "inference_eligible": reason is None,
        "inference_reason_code": reason,
        "estimand_note": "log2-ratio difference relative to the declared reference" if output_scale == "ratio_log2" else None,
    }


def convert_source_value(value: float | None, source_scale: str, zero_as_missing: bool, *, file: str, row: int, field: str) -> tuple[float | None, bool]:
    """Apply missing-zero encoding on the source scale.

    Returns (value, encoded_zero_missing).  Valid log2 zero/negative values are
    preserved; linear zero is missing only when explicitly declared so.
    """
    if value is None or source_scale != "linear_positive":
        return value, False
    if value == 0:
        if zero_as_missing:
            return None, True
        raise InputError("E_LINEAR_NONPOSITIVE", "linear zero is not declared as a missing encoding (linear_zero_as_missing=false)", file=file, row=row, field=field)
    if value < 0:
        raise InputError("E_LINEAR_NONPOSITIVE", f"linear abundance {value!r} is negative", file=file, row=row, field=field)
    return value, False


def log2_once(value: float | None) -> float | None:
    if value is None:
        return None
    if value <= 0:
        raise InputError("E_LINEAR_NONPOSITIVE", f"cannot take log2 of nonpositive value {value!r}")
    return math.log2(value)


def _null(token: str) -> str | None:
    return None if token in NULL_TOKENS else token


@dataclass
class Observations:
    header: list[str]
    rows: list[dict[str, str | None]]
    by_id: dict[str, dict[str, str | None]] = field(default_factory=dict)


def read_observations(path: str, *, group_column: str, unit_column: str) -> Observations:
    header, body = read_tsv(path)
    for column in REQUIRED_OBSERVATION_COLUMNS + (group_column, unit_column):
        if column not in header:
            raise InputError("E_METADATA_COLUMN", f"observation metadata is missing required column {column!r}", file=path, row=1, field=column)
    observations = Observations(header=list(header), rows=[])
    for number, raw in enumerate(body, start=2):
        row = {name: _null(token) for name, token in zip(header, raw)}
        observation_id = row["observation_id"]
        if observation_id is None:
            raise InputError("E_ID_MISSING", "observation_id is empty", file=path, row=number, field="observation_id")
        if observation_id in observations.by_id:
            raise InputError("E_ID_DUPLICATE", f"duplicate observation_id {observation_id!r}", file=path, row=number, field="observation_id")
        for column in (unit_column, group_column):
            if row[column] is None:
                raise InputError("E_METADATA_VALUE", f"{column} is required for every observation", file=path, row=number, field=column)
        row["_row"] = str(number)
        observations.rows.append(row)
        observations.by_id[observation_id] = row
    if not observations.rows:
        raise InputError("E_INPUT_EMPTY", "observation metadata has no rows", file=path)
    return observations


@dataclass
class Hierarchy:
    mode: str
    method: str | None
    unit_column: str
    group_column: str
    units: list[str]                      # canonical observation order
    members: dict[str, list[str]]         # canonical id -> source observation ids
    canonical_rows: dict[str, dict[str, str | None]]
    counts: list[dict]
    lineage: list[dict]


def build_hierarchy(observations: Observations, order: list[str], *, mode: str, method: str | None,
                    unit_column: str, group_column: str, file: str) -> Hierarchy:
    """Group source observations into biological units and validate replication semantics."""
    if mode == "none" and method is not None:
        raise InputError("E_TECHNICAL_CONFIG", "technical_replicates.method is only allowed with mode=aggregate")
    if mode == "aggregate" and method not in ("mean_linear", "mean_log2"):
        raise InputError("E_TECHNICAL_CONFIG", "mode=aggregate requires method mean_linear or mean_log2")
    members: dict[str, list[str]] = {}
    for observation_id in order:
        unit = observations.by_id[observation_id][unit_column]
        members.setdefault(unit, []).append(observation_id)
    repeated = {unit: ids for unit, ids in members.items() if len(ids) > 1}
    if repeated and mode == "none":
        unit, ids = next(iter(repeated.items()))
        raise InputError("E_TECHNICAL_REPLICATION_UNMODELED", f"biological unit {unit!r} has {len(ids)} observations ({', '.join(ids)}) but technical_replicates.mode=none; declare an explicit aggregation", file=file, field=unit_column)
    canonical_rows: dict[str, dict[str, str | None]] = {}
    lineage: list[dict] = []
    for unit, ids in members.items():
        rows = [observations.by_id[i] for i in ids]
        if len(ids) > 1:
            technical = [row["technical_replicate_id"] for row in rows]
            if any(t is None for t in technical) or len(set(technical)) != len(technical):
                raise InputError("E_TECHNICAL_ID", f"injections of {unit!r} need unique technical_replicate_id values", file=file, row=int(rows[0]["_row"]), field="technical_replicate_id")
            for column in observations.header:
                if column in TECHNICAL_COLUMNS:
                    continue
                distinct = {row[column] for row in rows}
                if len(distinct) > 1:
                    raise InputError("E_TECHNICAL_AGGREGATION_INVALID", f"aggregation unit {unit!r} spans different {column!r} values {sorted(str(v) for v in distinct)}; only injections of the same specimen/condition/visit may be averaged", file=file, field=column)
        canonical_id = unit if mode == "aggregate" else ids[0]
        base = {key: value for key, value in rows[0].items() if key != "_row"}
        base["observation_id"] = canonical_id
        if mode == "aggregate":
            base["technical_replicate_id"] = None
        base["n_injections"] = str(len(ids))
        canonical_rows[canonical_id] = base
        for source_id in ids:
            lineage.append({"source_observation_id": source_id, "observation_id": canonical_id, "biological_unit_id": unit,
                            "aggregation": method if (mode == "aggregate") else "none", "n_injections": len(ids)})
    units = list(canonical_rows)
    counts = []
    groups: list[str] = []
    for unit in units:
        group = canonical_rows[unit][group_column]
        if group not in groups:
            groups.append(group)
    for group in groups:
        group_units = [u for u in units if canonical_rows[u][group_column] == group]
        subjects = {canonical_rows[u]["subject_id"] or f"unit:{canonical_rows[u][unit_column]}" for u in group_units}
        counts.append({"group": group,
                       "n_injections": sum(len(members[canonical_rows[u][unit_column]]) for u in group_units),
                       "n_biological_units": len(group_units),
                       "n_subjects": len(subjects),
                       "small_n": len(group_units) < 2})
    return Hierarchy(mode=mode, method=method, unit_column=unit_column, group_column=group_column, units=units,
                     members={canonical: members[canonical_rows[canonical][unit_column]] for canonical in units},
                     canonical_rows=canonical_rows, counts=counts, lineage=lineage)


def total_subjects(hierarchy: Hierarchy) -> int:
    return len({row["subject_id"] or f"unit:{row[hierarchy.unit_column]}" for row in hierarchy.canonical_rows.values()})


@dataclass
class CanonicalValues:
    values: dict[tuple[str, str], float | None]
    numeric: dict[tuple[str, str], bool]
    observed: dict[tuple[str, str], bool] | None       # None = original mask unknown
    prior_imputed: dict[tuple[str, str], bool] | None
    coverage: list[dict]


def canonicalize(feature_ids: list[str], hierarchy: Hierarchy, source_values: dict[tuple[str, str], float | None],
                 observed: dict[tuple[str, str], bool] | None, prior_imputed: dict[tuple[str, str], bool] | None,
                 scale: dict) -> CanonicalValues:
    """Aggregate declared technical injections and apply the one-time transform.

    mean_linear averages positive linear values then takes log2 once;
    mean_log2 takes log2 once then averages.  When the original-observed mask
    is known only genuinely observed injections contribute to an observed
    specimen value; a specimen with only prior-imputed injections stays
    marked prior-imputed.
    """
    transform = scale["applied_transform"]
    method = hierarchy.method if hierarchy.mode == "aggregate" else None
    if method == "mean_linear" and scale["source_scale"] != "linear_positive":
        raise InputError("E_SCALE_REVERSE_UNSUPPORTED", "mean_linear aggregation requires linear_positive source values; reverse-scale aggregation of log values is not exposed")
    values: dict[tuple[str, str], float | None] = {}
    numeric: dict[tuple[str, str], bool] = {}
    obs_out: dict[tuple[str, str], bool] | None = {} if observed is not None else None
    imp_out: dict[tuple[str, str], bool] | None = {} if prior_imputed is not None else None
    coverage: list[dict] = []
    for feature_id in feature_ids:
        for canonical_id in hierarchy.units:
            sources = hierarchy.members[canonical_id]
            key = (feature_id, canonical_id)
            available = [s for s in sources if source_values[(feature_id, s)] is not None]
            if observed is not None:
                genuine = [s for s in available if observed[(feature_id, s)]]
                contributing = genuine if genuine else available
                is_observed = bool(genuine)
            else:
                contributing, is_observed = available, None
            if method == "mean_linear":
                raw = [source_values[(feature_id, s)] for s in contributing]
                value = fmean(raw) if raw else None
                value = log2_once(value) if transform == "log2" else value
            else:
                raw = [source_values[(feature_id, s)] for s in contributing]
                if transform == "log2":
                    raw = [log2_once(v) for v in raw]
                value = (fmean(raw) if len(raw) > 1 else raw[0]) if raw else None
            values[key] = value
            numeric[key] = value is not None
            if obs_out is not None:
                obs_out[key] = bool(is_observed)
            if imp_out is not None:
                imp_out[key] = value is not None and not is_observed if observed is not None else any(prior_imputed[(feature_id, s)] for s in contributing)
            if len(sources) > 1:
                coverage.append({"feature_id": feature_id, "observation_id": canonical_id, "n_injections": len(sources),
                                 "n_numeric": len(available), "n_observed": None if observed is None else sum(observed[(feature_id, s)] for s in available),
                                 "n_contributing": len(contributing)})
    return CanonicalValues(values=values, numeric=numeric, observed=obs_out, prior_imputed=imp_out, coverage=coverage)
