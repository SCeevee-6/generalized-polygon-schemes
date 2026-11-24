"""Verify the explicit nine-root incidence structure from its line triples."""

from itertools import combinations, permutations
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, eye
from scipy.sparse.csgraph import shortest_path

from grid_directions import skeleton
from theta_profiles import COS, distance_polynomials


def walk_counts(inc):
    """Count nonbacktracking paths from each point in batches, using integer arithmetic."""
    n, m = inc.shape
    identity = eye(n, dtype=np.int64, format="csr")
    deficits = []
    opposite_paths = 0
    for start in range(0, n, 64):
        unit = identity[start : start + 64]
        one = unit @ inc
        two = one @ inc.T - 5 * unit
        three = two @ inc - 2 * one
        four = three @ inc.T - 4 * two
        # Unique paths through radius four at every point exclude incidence cycles of length <=8.
        assert (unit + two + four).max() == 1
        assert (one + three).max() == 1
        five = four @ inc - 2 * three
        assert five.min() >= 0
        opposite_paths += int((five.data * (five.data - 1)).sum() // 2)
        six = five @ inc.T - 4 * four
        seven = six @ inc - 2 * five
        walks = (one + three + five + seven).toarray()
        assert walks.min() >= 0 and np.all(walks.sum(axis=1) == m)
        deficits.extend(np.count_nonzero(walks == 0, axis=1))
    # A 10-cycle gives five pairs of opposite vertices, one of each type.
    assert opposite_paths > 0 and opposite_paths % 5 == 0
    return np.array(deficits), opposite_paths // 5


def verify(lines):
    n, m = 1755, 2925
    assert lines.shape == (m, 3) and np.all(np.diff(lines, axis=1) > 0)
    assert lines.min() == 0 and lines.max() == n - 1
    ids = {tuple(line): j for j, line in enumerate(lines)}
    assert len(ids) == m
    inc = csr_matrix(
        (np.ones(3 * m, dtype=np.int64), (lines.ravel(), np.repeat(np.arange(m), 3))),
        shape=(n, m),
    )
    assert np.all(inc.sum(axis=1) == 5)
    graph = inc @ inc.T
    graph.setdiag(0)
    graph.eliminate_zeros()
    assert graph.max() == 1
    distances = shortest_path(graph, unweighted=True)
    assert np.isfinite(distances).all()
    dist = distances.astype(np.int8)
    tips = dist[:9]
    profiles, seed_lines = skeleton((0,) * 9)
    assert np.array_equal(tips[:, :81].T, profiles)
    assert all(tuple(sorted(line)) in ids for line in seed_lines)
    distances_on_lines = np.sort(tips[:, lines], axis=2)
    assert np.all(distances_on_lines[:, :, 1] == distances_on_lines[:, :, 0] + 1)
    assert np.all(distances_on_lines[:, :, 2] == distances_on_lines[:, :, 1])
    shells = [1, 10, 80, 640, 1024]
    assert all(np.bincount(row).tolist() == shells for row in tips)
    for i, j in combinations(range(9), 2):
        shell_pairs = np.bincount(5 * tips[i] + tips[j], minlength=25).reshape(5, 5)
        assert np.array_equal(shell_pairs, [p[3] for p in distance_polynomials()])
    cos = COS[tips]
    for a in cos:
        assert not ((cos * a) @ cos.T).any()
    theory = {tuple(range(3 * i, 3 * i + 3)) for i in range(3)}
    theory.update(tuple(range(j, 9, 3)) for j in range(3))
    near = distances_on_lines[:, :, 0] == 1
    centers = {
        t: np.flatnonzero(near[list(t)].all(axis=0)).tolist() for t in combinations(range(9), 3)
    }
    assert {t for t, l in centers.items() if l} == theory
    assert all(len(centers[t]) == 1 for t in theory)
    for perm in permutations(range(3)):
        assert not centers[tuple(3 * i + j for i, j in enumerate(perm))]
    selected = [ids[tuple(sorted(l))] for l in seed_lines[:6]]
    dual = shortest_path(inc.T @ inc, unweighted=True, indices=selected)[:, selected]
    expected = np.full((6, 6), 2)
    expected[:3, :3] = 4
    expected[3:, 3:] = 4
    np.fill_diagonal(expected, 0)
    assert np.array_equal(dual, expected)
    deficit, cycles = walk_counts(inc)
    assert np.array_equal(np.flatnonzero(deficit == 0), np.arange(9))
    # Independently count the unreachable lines by shortest distances.
    via_distances = np.zeros(n, dtype=int)
    diameter = 2 * int(dist.max())
    for batch in np.array_split(lines, 40):
        mixed = 2 * dist[:, batch].min(axis=2) + 1
        via_distances += np.count_nonzero(mixed > 7, axis=1)
        # For distinct lines, cross incidence distance is one plus the minimum
        # from a point of the first to the second. Diagonal placeholders equal two.
        diameter = max(diameter, int(mixed.max()), 1 + int(mixed[lines].min(axis=1).max()))
    assert np.array_equal(deficit, via_distances)

    print("points", n, "lines", m, "degrees", 5, 3)
    print("marked shells", shells, "required seed lines", len(seed_lines))
    print("centered rows and columns", len(theory), "centered transversals", 0)
    print(
        "point roots with a unique path of length at most seven to every line",
        int((deficit == 0).sum()),
    )
    print(
        "incidence girth", 10, "cycles of length ten", cycles, "point diameter", int(dist.max())
    )
    print(
        "incidence diameter",
        diameter,
        "point-line pairs at distance above seven",
        int(deficit.sum()),
    )


if __name__ == "__main__":
    verify(np.loadtxt(Path(__file__).with_name("rooted_grid_lines.txt"), dtype=int))
