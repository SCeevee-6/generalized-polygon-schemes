"""Check the opposite-line criterion on every opposite pair in the ATLAS model."""

from collections import Counter
from itertools import permutations
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "flags"))
from o2_build import octagon
from theta_profiles import seed_distances


def only(mask):
    ids = np.flatnonzero(mask)
    assert len(ids) == 1
    return int(ids[0])


def check_suboctagon(lines, chosen):
    points, positions = np.unique(lines[sorted(chosen)], return_inverse=True)
    assert len(chosen) == 30 and len(points) == 45
    incidence = csr_matrix(
        (np.ones(90, dtype=int), (np.repeat(np.arange(30), 3), positions.ravel())), shape=(30, 45)
    )
    assert np.all(incidence.sum(axis=0) == 2)
    # Suppress the degree-two points. The cubic graph has the GQ(2,2)
    # intersection array, which also certifies girth eight and diameter four.
    graph = (incidence @ incidence.T).toarray()
    np.fill_diagonal(graph, 0)
    assert graph.max() == 1 and np.all(graph.sum(1) == 3)
    d = shortest_path(csr_matrix(graph), unweighted=True)
    assert np.isfinite(d).all() and d.max() == 4
    for h, (back, forward) in enumerate([(0, 3), (1, 2), (1, 2), (1, 2), (3, 0)]):
        for x, y in zip(*np.where(d == h)):
            shell = d[x, graph[y].astype(bool)]
            assert np.count_nonzero(shell == h - 1) == back
            assert np.count_nonzero(shell == h + 1) == forward


def main():
    gens = np.load(ROOT / "data/generators.npy")
    pd, lines = octagon(gens)
    assert lines.shape == (2925, 3)
    lookup = {tuple(triple): j for j, triple in enumerate(lines)}
    orbit = [0]
    seen = {0}
    for j in orbit:
        for gen in gens:
            k = lookup[tuple(sorted(gen[lines[j]]))]
            if k not in seen:
                orbit.append(k)
                seen.add(k)
    assert len(orbit) == len(lines)
    inc = csr_matrix(
        (np.ones(lines.size), (np.repeat(np.arange(len(lines)), 3), lines.ravel())),
        shape=(len(lines), len(pd)),
    )
    ld = shortest_path(inc @ inc.T, unweighted=True).astype(int)
    near = pd[:, lines].min(axis=2) == 1
    opposite = np.flatnonzero(ld[0] == 4)
    assert len(opposite) == 2048
    through_pair = {
        tuple(sorted((x, y))): j
        for j, (a, b, c) in enumerate(lines)
        for x, y in ((a, b), (a, c), (b, c))
    }
    expected_seed = seed_distances()
    by_copy = Counter()
    for partner in opposite:
        cols = np.flatnonzero((ld[0] == 2) & (ld[partner] == 2))
        assert len(cols) == 3
        rows = np.flatnonzero((ld[cols] == 2).all(axis=0))
        assert len(rows) == 3 and 0 in rows and partner in rows
        third_row = next(r for r in rows if r not in (0, partner))
        tips = np.empty((3, 3), dtype=int)
        chosen = set(rows) | set(cols)
        core = np.empty((2, 3, 2), dtype=int)
        for i, r in enumerate((0, partner, third_row)):
            for j, c in enumerate(cols):
                link = only((ld[r] == 1) & (ld[c] == 1))
                chosen.add(link)
                if i < 2:
                    core[i, j] = (
                        next(p for p in lines[link] if p in lines[r]),
                        next(p for p in lines[link] if p in lines[c]),
                    )
                tips[i, j] = next(
                    p for p in lines[link] if p not in lines[r] and p not in lines[c]
                )
        seed = []
        for j, c in enumerate(cols):
            a, l = core[0, j]
            b, r = core[1, j]
            middle = next(p for p in lines[c] if p not in (l, r))
            seed.extend((a, l, r, b, tips[0, j], middle, tips[1, j]))
        assert np.array_equal(pd[np.ix_(seed, seed)], expected_seed)
        extensions = np.flatnonzero((pd[tips[:2].ravel()] == 3).all(axis=0))
        assert set(extensions) == set(tips[2])
        for perm in permutations(range(3)):
            triple = tips[np.arange(3), perm]
            centre = only(near[triple].all(axis=0))
            chosen.add(centre)
            for tip in triple:
                p = lines[centre][only(pd[tip, lines[centre]] == 1)]
                chosen.add(through_pair[tuple(sorted((tip, p)))])
        key = frozenset(chosen)
        if key not in by_copy:
            check_suboctagon(lines, key)
        by_copy[key] += 1
    assert set(by_copy.values()) == {8}
    print("line orbit", len(orbit), "opposite partners checked", len(opposite))
    print("three extensions and six centered transversals", sum(by_copy.values()))
    print("verified thin copies through line 0", len(by_copy), "opposite partners per copy", 8)
    print(
        "unordered opposite pairs",
        len(orbit) * len(opposite) // 2,
        "thin copies",
        len(orbit) * len(by_copy) // 30,
    )


if __name__ == "__main__":
    main()
