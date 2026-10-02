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
