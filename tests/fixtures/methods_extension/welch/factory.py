"""Synthetic fixtures for the Welch engine (packet R15a, V168-V173).

Every value is drawn from seeded standard-library Gaussian generators with planted structure; nothing
here is study data. The builders reuse the design builders of tests/fixtures/design/builders.py, so
every packet writes the same dataset layout and the same analysis configuration shape.
"""
from __future__ import annotations

import copy
import csv
import importlib.util
import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
_spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def two_group(n_num: int = 6, n_den: int = 9, features: int = 200, *, seed: int = 168, planted: dict | None = None, sd_num: float = 0.4,
              sd_den: float = 1.6, base: float = 10.0) -> tuple[dict, list[dict]]:
    """Two groups with unequal variances; ``planted`` maps feature index (0-based) to a shift of the denominator group."""
    rng = random.Random(seed)
    observations = [{"observation_id": f"A{i + 1}", "group": "A"} for i in range(n_num)] + [{"observation_id": f"B{i + 1}", "group": "B"} for i in range(n_den)]
    values = {}
    for f in range(features):
        shift = (planted or {}).get(f, 0.0)
        row = []
        for o in observations:
            sd = sd_num if o["group"] == "A" else sd_den
            mean = base + 0.02 * f + (shift if o["group"] == "B" else 0.0)
            row.append(mean + rng.gauss(0.0, sd))
        values[f"F{f + 1:03d}"] = row
    return values, observations


def three_group(per_group: int = 4, features: int = 60, *, seed: int = 15, shifts: dict | None = None, sd: float = 0.35) -> tuple[dict, list[dict]]:
    """Three groups (control, disease, treated) with planted up and down features (shifts maps feature index to (group, shift))."""
    rng = random.Random(seed)
    groups = ["control", "disease", "treated"]
    observations = [{"observation_id": f"{g[:2].upper()}{i + 1}", "group": g} for g in groups for i in range(per_group)]
    values = {}
    for f in range(features):
        row = []
        for o in observations:
            shift = (shifts or {}).get(f, {}).get(o["group"], 0.0)
            row.append(10.0 + 0.05 * f + shift + rng.gauss(0.0, sd))
        values[f"G{f + 1:03d}"] = row
    return values, observations


def paired(subjects: int = 8, features: int = 200, *, seed: int = 171, planted: dict | None = None, sd: float = 0.3) -> tuple[dict, list[dict]]:
    """Subjects measured under conditions A and B, paired within subject; planted within-subject differences (B minus A)."""
    rng = random.Random(seed)
    observations = []
    for s in range(subjects):
        observations.append({"observation_id": f"S{s + 1}A", "group": "A", "subject_id": f"subj{s + 1}", "biological_unit_id": f"S{s + 1}A"})
        observations.append({"observation_id": f"S{s + 1}B", "group": "B", "subject_id": f"subj{s + 1}", "biological_unit_id": f"S{s + 1}B"})
    values = {}
    subject_effect = {s: rng.gauss(0.0, 1.0) for s in range(subjects)}
    for f in range(features):
        shift = (planted or {}).get(f, 0.0)
        row = []
        for o in observations:
            s = int(o["subject_id"][4:]) - 1
            row.append(10.0 + 0.01 * f + subject_effect[s] + (shift if o["group"] == "B" else 0.0) + rng.gauss(0.0, sd))
        values[f"P{f + 1:03d}"] = row
    return values, observations


def contrast(cid: str, numerator: str, denominator: str, *, design_id: str = "joint", role: str = "primary") -> dict:
    return B.contrast(cid, numerator, denominator, role=role, design_id=design_id)


def study_config(directory: Path, values: dict, observations: list[dict], *, groups, contrasts, primary: str = "welch", methods: dict | None = None,
                 design_overrides: dict | None = None, extra_models: list | None = None, families: list | None = None, seed: int = 20261008,
                 extra_columns=(), phase: int = 5, genes: bool = False) -> Path:
    """An analysis configuration whose primary model uses ``primary`` (welch or limma) on the observed primary matrix.

    ``methods`` is the methods_extension block (None omits it). ``extra_models`` appends secondary models.
    ``phase`` is 5 for every Welch configuration (ADR 0011); Phase 1-4 configurations cannot enable the extension.
    """
    files = B.dataset(directory, values, observations, extra_columns=tuple(extra_columns))
    if genes:   # synthetic gene identifiers (G0001 for feature G001), so the offline libraries can map features to genes
        rows = [["feature_id", "accessions", "gene_ids", "gene_symbols", "is_decoy", "is_contaminant", "protein_group_ambiguous"]]
        for feature in values:
            gene = "G%04d" % int(feature[1:])
            rows.append([feature, json.dumps([feature]), json.dumps([gene]), json.dumps([gene]), "unknown", "unknown", "false"])
        B.write_tsv(directory / "features.tsv", rows)

    def mutate(c: dict) -> None:
        c["runtime"]["seed"] = seed
        c["runtime"]["phase"] = phase
        c["primary_engine"] = primary
        c["models"][0]["engine"] = primary
        c["models"][0]["id"] = "welch-main" if primary == "welch" else "limma-main"
        if primary == "welch":
            c["models"][0].pop("trend", None)
            c["models"][0].pop("robust", None)
        if methods is not None:
            c["methods_extension"] = copy.deepcopy(methods)
        if extra_models:
            c["models"].extend(copy.deepcopy(extra_models))
        if families is not None:
            c["multiplicity_families"] = copy.deepcopy(families)
        else:
            primary_id = c["models"][0]["id"]
            family_id = ((methods or {}).get("welch") or {}).get("family", "welch_zero_null") if primary == "welch" else "primary-family"
            c["multiplicity_families"] = [{"id": family_id, "hypothesis_type": "protein_zero_null", "model_ids": [primary_id],
                                           "contrast_ids": [x["id"] for x in contrasts], "adjustment": "BH", "denominator": "finite_eligible",
                                           "q_cutoff": 0.05, "role": "primary"}]

    return B.config(directory, files, groups=groups, contrasts=contrasts, design_overrides=design_overrides,
                    model_overrides={"id": "welch-main" if primary == "welch" else "limma-main", "engine": primary, "role": "primary"},
                    families=[], mutate=mutate)


def welch_methods(*, paired_value="auto", paper=None, role="primary", family="welch_zero_null", contrasts=None, ci_level=None) -> dict:
    block = {"welch": {"enabled": True, "execution_requirement": "required", "role": role, "family": family, "paired": paired_value}}
    if paper is not None:
        block["welch"]["paper_compatible"] = paper
    if contrasts is not None:
        block["welch"]["contrasts"] = list(contrasts)
    if ci_level is not None:
        block["welch"]["ci_level"] = ci_level
    return block


def write_json(path: Path, value) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")
    return path
