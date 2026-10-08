"""Deterministic synthetic PERMANOVA fixtures with known structure (R13, V121-V130).

Values come from a seeded standard-library Gaussian generator; nothing is
study data.  Structures:
  null        - every group drawn from the same distribution (no separation);
  location    - group "B" shifted by +2.5 SD on the first half of features;
  dispersion  - identical centroids; group "B" deviations are the mirrored,
                four-fold-scaled deviations of group "A" (so group means are
                exactly equal and only the spread differs).
"""
from __future__ import annotations

import random


def matrix(kind: str, *, n_per_group: int = 6, n_features: int = 12, groups=("A", "B"), seed: int = 20261001):
    rng = random.Random(seed)
    observations = [{"observation_id": f"{g}{i + 1}", "group": g} for g in groups for i in range(n_per_group)]
    values = {}
    for f in range(n_features):
        base = 10 + f * 0.5
        if kind == "dispersion":
            deviations = [rng.gauss(0, 0.3) for _ in range(n_per_group)]
            mean_dev = sum(deviations) / n_per_group
            deviations = [d - mean_dev for d in deviations]
            column = [base + d for d in deviations] + [base - 4 * d for d in deviations]
            for g in groups[2:]:
                column += [base + rng.gauss(0, 0.3) for _ in range(n_per_group)]
        else:
            column = []
            for g in groups:
                shift = 2.5 * 0.3 if (kind == "location" and g == "B" and f < n_features // 2) else 0.0
                column += [base + shift + rng.gauss(0, 0.3) for _ in range(n_per_group)]
        values[f"F{f + 1:02d}"] = column
    return values, observations


def pure_location(*, n_per_group: int = 6, n_features: int = 12, seed: int = 20261001, shift: float = 1.0):
    """Audit 2026-10-02: group B = group A deviations + a constant shift, so spreads are identical
    (PERMDISP F = 0) and the only difference is the centroid: the exact label is location_shift."""
    rng = random.Random(seed)
    observations = [{"observation_id": f"{g}{i + 1}", "group": g} for g in ("A", "B") for i in range(n_per_group)]
    values = {}
    for f in range(n_features):
        deviations = [rng.gauss(0, 0.3) for _ in range(n_per_group)]
        mean_dev = sum(deviations) / n_per_group
        deviations = [d - mean_dev for d in deviations]
        base = 10 + f * 0.5
        values[f"F{f + 1:02d}"] = [base + d for d in deviations] + [base + shift + d for d in deviations]
    return values, observations


def subject_matrix(kind: str, *, n_subjects: int = 8, reps: int = 2, n_features: int = 10, shift: float = 1.5, seed: int = 20261002):
    """Audit 2026-10-02 subject-blocked fixtures.  between/between_null: group constant within each subject
    (reps observations per subject); within/within_null: each subject observed once in A and once in B;
    mixed: half the subjects in both groups, half in one group only.  A subject random effect is shared by a
    subject's observations; *_null kinds have no group shift."""
    rng = random.Random(seed)
    subjects = [f"S{i + 1}" for i in range(n_subjects)]
    if kind == "between_unbalanced":   # D-59: group constant within subject, 1-3 observations per subject
        layout = [(s, "A" if i < n_subjects // 2 else "B") for i, s in enumerate(subjects) for _ in range(1 + i % 3)]
    elif kind.startswith("between"):
        layout = [(s, "A" if i < n_subjects // 2 else "B") for i, s in enumerate(subjects) for _ in range(reps)]
    elif kind.startswith("within"):
        layout = [(s, g) for s in subjects for g in ("A", "B")]
    elif kind == "mixed":
        half = n_subjects // 2
        layout = [(s, g) for s in subjects[:half] for g in ("A", "B")] + [(s, "A" if i % 2 == 0 else "B") for i, s in enumerate(subjects[half:]) for _ in range(2)]
    else:
        raise ValueError(kind)
    counter: dict[str, int] = {}
    observations = []
    for subject, group in layout:
        counter[subject] = counter.get(subject, 0) + 1
        observations.append({"observation_id": f"{subject}_{counter[subject]}", "biological_unit_id": f"{subject}_{counter[subject]}", "subject_id": subject, "group": group})
    effect = shift if kind in ("between", "within", "mixed", "between_unbalanced") else 0.0
    subject_effect = {s: [rng.gauss(0, 0.3) for _ in range(n_features)] for s in subjects}
    values = {f"F{f + 1:02d}": [10 + 0.2 * f + subject_effect[o["subject_id"]][f] + rng.gauss(0, 0.3) + (effect if (o["group"] == "B" and f < n_features // 2) else 0.0)
                                for o in observations] for f in range(n_features)}
    return values, observations
