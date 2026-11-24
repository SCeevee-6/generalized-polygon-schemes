"""Enumerate opposite-point selector parity solutions in the Ree–Tits model."""
from collections import Counter
from itertools import product
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "octagon"))
from geom import octagon


def reduce_binary(matrix):
    matrix = matrix.copy()
    pivots = []
    for column in range(matrix.shape[1]):
        row = len(pivots)
        candidates = np.flatnonzero(matrix[row:, column])
        if not len(candidates):
            continue
        other = row + candidates[0]
        matrix[[row, other]] = matrix[[other, row]]
        clear = matrix[:, column].copy()
        clear[row] = False
        matrix[clear] ^= matrix[row]
        pivots.append(column)
    return matrix[:len(pivots)], pivots


def main():
    _, adjacency, distance, lines = octagon()
    root = 0
    opposite = np.flatnonzero(distance[root] == 4)
    outer_lines = lines[distance[root, lines].min(axis=1) == 1]
    root_lines = lines[(lines == root).any(axis=1)]
    parents = [line[distance[root, line] == 1].item() for line in outer_lines]
    branch = np.array([np.flatnonzero((root_lines == p).any(axis=1)).item()
                       for p in parents])
    matrix = distance[opposite[:, None, None], outer_lines].min(axis=2) == 2
    assert matrix.shape == (1024, 40)
    assert np.all(matrix.sum(axis=1) == 5)
    assert np.all(np.bincount(branch) == 8)
    reduced, pivots = reduce_binary(matrix)
    outer = csr_matrix(adjacency[np.ix_(opposite, opposite)])
    components, labels = connected_components(outer)
    ranks = [len(reduce_binary(matrix[labels == j])[1]) for j in range(components)]
    assert len(pivots) == 29 and ranks == [29, 29]
    print("Full rank", len(pivots), "component ranks", ranks)

    free = np.setdiff1d(np.arange(40), pivots)
    one_branch = branch == 0
    assert np.all(matrix[:, one_branch].sum(axis=1) == 1)
    seen = Counter()
    for bits in product((0, 1), repeat=len(free)):
        vector = np.zeros(40, dtype=np.int64)
        vector[free] = bits
        vector[pivots] = reduced[:, free] @ vector[free] % 2
        vector ^= one_branch
        assert np.all(matrix @ vector % 2 == 1)
        counts = tuple(sorted(np.bincount(branch[vector == 1], minlength=5)))
        if sum(counts) in (10, 20, 30, 40):
            seen[counts] += 1
    assert seen == {
        (2, 2, 2, 2, 2): 16,
        (0, 4, 4, 4, 8): 240,
        (4, 4, 4, 4, 4): 288,
        (6, 6, 6, 6, 6): 48,
        (8, 8, 8, 8, 8): 1,
    }
    for counts, number in sorted(seen.items(), key=lambda item: (sum(item[0]), item[0])):
        print(tuple(map(int, counts)), number)


if __name__ == "__main__":
    main()
