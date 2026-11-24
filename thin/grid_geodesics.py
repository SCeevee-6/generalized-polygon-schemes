"""Necessary identifications of geodesics between named points and lines."""

from itertools import combinations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path


def check_geodesics(point_distance, lines):
    """Verify the tree identifications for every named triangle of perimeter < 16."""
    lines = np.asarray(lines)
    mixed = 2 * point_distance[:, lines].min(axis=2) + 1
    dual = mixed[lines, :].min(axis=1) + 1
    np.fill_diagonal(dual, 0)
    distance = np.block([[2 * point_distance, mixed], [mixed.T, dual]])
    n = len(distance)
    assert distance.max() == 8
    parent = list(range(n))
    paths = {}

    def root(i):
        while i != parent[i]:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def identify(i, j):
        a, b = sorted((root(i), root(j)))
        if a != b:
            assert b >= n, (a, b, "distinct named vertices identified")
            parent[b] = a

    for i, j in combinations(range(n), 2):
        d = int(distance[i, j])
        if d == 8:
            continue
        start = len(parent)
        middle = tuple(range(start, start + d - 1))
        parent.extend(middle)
        paths[i, j] = (i, *middle, j)

    for (i, j), ij in paths.items():
        for k in np.flatnonzero(
            (np.arange(n) > j) & (distance[i] + distance[j] + distance[i, j] < 16)
        ):
            ik, jk = paths[i, k], paths[j, k]
            shared = (len(ij) + len(ik) - len(jk) - 1) // 2
            for s in range(1, shared + 1):
                identify(ij[s], ik[s])
            for s in range(1, len(ij) - shared):
                identify(ij[-1 - s], jk[s])
            for s in range(1, len(ik) - shared):
                identify(ik[-1 - s], jk[-1 - s])

    _, labels = np.unique([root(i) for i in range(len(parent))], return_inverse=True)
    edges = np.array(
        [(labels[a], labels[b]) for path in paths.values() for a, b in zip(path, path[1:])]
    )
    graph = csr_matrix((np.ones(len(edges)), edges.T), shape=(int(labels.max()) + 1,) * 2)
    actual = shortest_path(graph, directed=False, unweighted=True, indices=labels[:n])
    assert np.all(actual[:, labels[:n]] >= distance)
