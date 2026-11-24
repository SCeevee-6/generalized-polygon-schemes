"""Check five-line neighbourhoods using the 21 opposite-pair coordinates."""
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path

import numpy as np

from theta_profiles import expand, profiles, seed_distances, verify


def line_types(points):
    lookup = {tuple(row): i for i, row in enumerate(points.tolist())}
    lines = []
    for i, a in enumerate(points):
        candidates = np.where((np.abs(points[i:] - a) <= 1).all(axis=1))[0] + i
        third = np.where(points[candidates] == a, a - 1, np.maximum(points[candidates], a))
        for j, c in zip(candidates, third.tolist()):
            k = lookup.get(tuple(c))
            if k is not None and j <= k:
                lines.append((i, j, k))
    return np.asarray(lines, dtype=int)


def unsupported(points):
    lines = line_types(points)
    minimum = points[lines].min(axis=1)
    bits = 1 << np.arange(points.shape[1])
    rejected = []
    for i, p in enumerate(points):
        incident = np.where((lines == i).any(axis=1))[0]
        active = (p > 0) & (p < 4)
        masks = set(((minimum[incident] < p) & active) @ bits)
        # On every line a coordinate at four has one neighbour at three.
        assert np.all(minimum[incident][:, p == 4] == 3)

        @lru_cache(None)
        def cover(remaining, slots):
            if not remaining:
                return slots == 0 or 0 in masks
            if not slots:
                return False
            bit = remaining & -remaining
            return any(
                cover(remaining ^ m, slots - 1) for m in masks
                if m & bit and m & remaining == m
            )

        if not cover(int(active @ bits), 5):
            rejected.append(i)
    return rejected


def main():
    dist = seed_distances()
    reps = profiles(dist)
    points, labels = expand(reps)
    data = json.loads(Path(__file__).with_name("theta_no_extension.json").read_text())
    counts = Counter({tuple(p): n for p, n in data})
    weights = np.array([counts[tuple(reps[i])] for i in labels])
    verify(dist, points, weights)
    excluded = unsupported(points)
    assert len(excluded) == 3 and len(set(labels[excluded])) == 1
    print("unconditionally unsupported under no extension", len(excluded))
    occupied = np.flatnonzero(weights)
    rejected = occupied[unsupported(points[occupied])]
    assert (len(set(labels[rejected])), int(weights[rejected].sum())) == (14, 249)
    print("saved support: rejected profiles", len(rejected), "orbits", len(set(labels[rejected])),
          "point multiplicity", weights[rejected].sum())


if __name__ == "__main__":
    main()
