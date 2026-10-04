"""One synthetic design factory for every post-DE packet (R14a-R14f, V131-V167).

Values come from seeded standard-library Gaussian generators with planted structure; nothing is
study data.  The same factory produces every V167 generality cell (two-group, three-group, paired,
repeated with subject blocks, continuous exposure, 3-vs-3, 15-vs-80, 200 units, and human, mouse and
rat identifier conventions), so every packet tests the same design cells.
"""
from __future__ import annotations

import copy
import csv
import importlib.util
import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(B)

# Identifier conventions only (synthetic accessions/symbols in each species' style; no real mapping is implied).
SPECIES = {"human": {"taxonomy_id": 9606, "name": "Homo sapiens", "gene": lambda i: f"HSG{i:04d}", "accession": lambda i: f"P9{i:04d}"},
           "mouse": {"taxonomy_id": 10090, "name": "Mus musculus", "gene": lambda i: f"Msg{i:04d}", "accession": lambda i: f"Q8{i:04d}"},
           "rat": {"taxonomy_id": 10116, "name": "Rattus norvegicus", "gene": lambda i: f"Rng{i:04d}", "accession": lambda i: f"O7{i:04d}"}}


def write_dataset(directory: Path, values: dict, observations: list[dict], *, species: str = "rat", shared_genes: dict | None = None, extra_columns=()) -> dict:
    """Like builders.dataset, but features carry species-style accessions and gene identifiers.
    shared_genes maps feature id -> gene index, so several protein groups can share one gene."""
    directory.mkdir(parents=True, exist_ok=True)
    files = B.dataset(directory, values, observations, extra_columns=tuple(extra_columns))
    style = SPECIES[species]
    rows = [["feature_id", "accessions", "gene_ids", "gene_symbols", "is_decoy", "is_contaminant", "protein_group_ambiguous"]]
    for i, feature in enumerate(values, start=1):
        gene = (shared_genes or {}).get(feature, i)
        rows.append([feature, json.dumps([style["accession"](i)]), json.dumps([style["gene"](gene)]), json.dumps([style["gene"](gene)]), "unknown", "unknown", "false"])
    B.write_tsv(directory / "features.tsv", rows)
    return files


def config(directory: Path, files: dict, *, groups, contrasts, post_de: dict, families=None, design_overrides=None, model_overrides=None, mutate=None,
           species: str = "rat", seed: int = 20261003) -> Path:
    def change(c):
        c["runtime"]["phase"] = 4; c["runtime"]["execution_profile"] = "smoke_test"; c["runtime"]["seed"] = seed
        c["organism"] = {"taxonomy_id": SPECIES[species]["taxonomy_id"], "scientific_name": SPECIES[species]["name"]}
        c["post_de"] = copy.deepcopy(post_de)
        if mutate:
            mutate(c)
    return B.config(directory, files, groups=groups, contrasts=contrasts, families=families, design_overrides=design_overrides,
                    model_overrides=model_overrides, mutate=change)


def gaussian_groups(groups: dict[str, int], n_features: int, *, effects: dict | None = None, seed: int = 1, sd: float = 0.3,
                    base=lambda f: 10 + 0.1 * f, prefix: str = "F") -> tuple[dict, list[dict]]:
    """Independent units; effects maps (group, feature index 0-based) -> shift in log2 units."""
    rng = random.Random(seed)
    observations = [{"observation_id": f"{g}{i + 1}", "group": g} for g, n in groups.items() for i in range(n)]
    values = {}
    for f in range(n_features):
        values[f"{prefix}{f + 1:02d}"] = [base(f) + (effects or {}).get((o["group"], f), 0.0) + rng.gauss(0, sd) for o in observations]
    return values, observations


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def family(fid, contrasts, *, hypothesis="protein_zero_null", models=("limma-main",), role="secondary"):
    return {"id": fid, "hypothesis_type": hypothesis, "model_ids": list(models), "contrast_ids": list(contrasts), "adjustment": "BH",
            "denominator": "finite_eligible", "q_cutoff": 0.05, "role": role}


def leaf(family_id, contrast_id, criterion="family_q", threshold=0.05, direction="any", **extra):
    return {"family_id": family_id, "contrast_id": contrast_id, "criterion": criterion, "threshold": threshold, "direction": direction, **extra}


# ----------------------------------------------------------------------------- R14a fixture
SETS_PLANTED = {"concordant_up": ["F01", "F02"], "concordant_down": ["F05", "F06"], "discordant": ["F09", "F10"],
                "only_b": ["F03", "F04", "F07", "F08"], "only_c": ["F11", "F12"], "d_vs_c_up": ["F03", "F04", "F13", "F14", "F15"]}


def four_group_sets(seed: int = 8, shift: float = 1.5):
    """Groups A, B, C, D (5 units each); contrasts B-A and C-A share A, while B-A and D-C share no units.
    Planted (hand-labelled in SETS_PLANTED): B-A up F01-F04, F09-F10 and down F05-F08; C-A up F01-F02, F11-F12 and
    down F05-F06, F09-F10; D = C plus +shift on F03-F04 and F13-F15.  F01 and F02 share one gene."""
    effects = {}
    for i in (0, 1, 2, 3, 8, 9): effects[("B", i)] = shift
    for i in (4, 5, 6, 7): effects[("B", i)] = -shift
    for i in (0, 1, 10, 11): effects[("C", i)] = shift
    for i in (4, 5, 8, 9): effects[("C", i)] = -shift
    for (g, i), v in list(effects.items()):
        if g == "C": effects[("D", i)] = v
    for i in (2, 3, 12, 13, 14): effects[("D", i)] = effects.get(("D", i), 0.0) + shift
    values, obs = gaussian_groups({"A": 5, "B": 5, "C": 5, "D": 5}, 30, effects=effects, seed=seed)
    return values, obs, {"F01": 1, "F02": 1}


# ----------------------------------------------------------------------------- R14b fixtures
def sex_three_groups(*, variant: str = "partial", n: int = 8, seed: int = 21, n_features: int = 30):
    """Groups A, B, C x sex.  variant 'partial': A 2M/6F, B 6M/2F, C 4M/4F, so part of the B-A effect is a sex effect;
    'balanced': 4M/4F everywhere; 'confounded': C entirely male; 'one_male': as partial but A has a single male.
    Sex shifts F01-F06 by +1.2 (M); group B shifts F01-F03 and F07-F09 by +1.0; group C shifts F10-F12 by -1.2."""
    males = {"partial": {"A": 2, "B": 6, "C": 4}, "balanced": {"A": 4, "B": 4, "C": 4}, "confounded": {"A": 4, "B": 4, "C": 8},
             "one_male": {"A": 1, "B": 6, "C": 4}}[variant]
    rng = random.Random(seed)
    observations = [{"observation_id": f"{g}{i + 1}", "group": g, "sex": "M" if i < males[g] else "F"} for g in ("A", "B", "C") for i in range(n)]
    values = {}
    for f in range(n_features):
        column = []
        for o in observations:
            shift = (1.2 if o["sex"] == "M" and f < 6 else 0.0) + (1.0 if o["group"] == "B" and (f < 3 or 6 <= f < 9) else 0.0) + (-1.2 if o["group"] == "C" and 9 <= f < 12 else 0.0)
            column.append(10 + 0.1 * f + shift + rng.gauss(0, 0.3))
        values[f"F{f + 1:02d}"] = column
    return values, observations


def unbalanced_with_injections(*, n_a: int = 31, n_b: int = 11, n_features: int = 25, seed: int = 31, leverage: float = 3.0):
    """31-vs-11 independent units with known B-A effects on F01-F08 (+1.0); unit B1 is measured by two technical
    injections (aggregated to one biological unit) and is a planted high-leverage unit: +leverage on F01-F04."""
    rng = random.Random(seed)
    observations = []
    for g, n in (("A", n_a), ("B", n_b)):
        for i in range(n):
            unit = f"{g}{i + 1}"
            injections = 2 if unit == "B1" else 1
            for j in range(injections):
                observations.append({"observation_id": f"{unit}_i{j + 1}" if injections > 1 else unit, "biological_unit_id": unit,
                                     "technical_replicate_id": f"i{j + 1}" if injections > 1 else "NA", "group": g})
    unit_noise = {}
    values = {}
    for f in range(n_features):
        column = []
        for o in observations:
            key = (o["biological_unit_id"], f)
            if key not in unit_noise:
                unit_noise[key] = rng.gauss(0, 0.4)
            shift = (1.0 if o["group"] == "B" and f < 8 else 0.0) + (leverage if o["biological_unit_id"] == "B1" and f < 4 else 0.0)
            column.append(10 + 0.1 * f + shift + unit_noise[key] + rng.gauss(0, 0.05))
        values[f"F{f + 1:02d}"] = column
    return values, observations


# ----------------------------------------------------------------------------- R14c fixtures
def phenotype_design(*, n: int = 10, seed: int = 41, n_features: int = 30, missing: int = 0, kind: str = "linear"):
    """Groups A and B (n units each), continuous covariate age and phenotype score.
    kind 'linear': F01-F10 rise 0.4 log2 per phenotype unit; group B shifts F11-F15 by +1; age adds 0.02 per year everywhere.
    kind 'one_group': the phenotype is recorded only in group B.  kind 'simpson': phenotype ranges overlap (A 1-7, B 4-10);
    F01-F05 fall 0.4 per unit within each group while group B is 6 log2 higher, so the pooled slope is positive.
    missing: the phenotype is missing ("NA") for the first `missing` units."""
    rng = random.Random(seed)
    observations = []
    for g in ("A", "B"):
        for i in range(n):
            if kind == "simpson":
                score = rng.uniform(1, 7) if g == "A" else rng.uniform(4, 10)
            else:
                score = rng.gauss(5, 2)
            observations.append({"observation_id": f"{g}{i + 1}", "group": g, "age": f"{rng.uniform(20, 60):.3f}", "score": f"{score:.6f}"})
    for o in observations[:missing]:
        o["score"] = "NA"
    if kind == "one_group":
        for o in observations:
            if o["group"] == "A":
                o["score"] = "NA"
    means = {g: sum(float(o["score"]) for o in observations if o["group"] == g and o["score"] != "NA") / max(1, sum(1 for o in observations if o["group"] == g and o["score"] != "NA"))
             for g in ("A", "B")}
    values = {}
    for f in range(n_features):
        column = []
        for o in observations:
            s = float(o["score"]) if o["score"] != "NA" else means[o["group"]]
            if kind == "simpson":
                shift = (6.0 * (o["group"] == "B") - 0.4 * (s - means[o["group"]])) if f < 5 else 0.0
            else:
                shift = (0.4 * s if f < 10 else 0.0) + (1.0 if o["group"] == "B" and 10 <= f < 15 else 0.0)
            column.append(10 + 0.1 * f + shift + 0.02 * float(o["age"]) + rng.gauss(0, 0.3))
        values[f"F{f + 1:02d}"] = column
    return values, observations
