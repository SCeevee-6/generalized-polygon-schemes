"""Exact boundary reconstruction and mixed-moment checks in the Ree--Tits model."""
from itertools import combinations
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "flags"))
from o2_build import octagon


def five_cliques(ds):
    def extend(q, pool):
        if len(q) == 5:
            yield tuple(q)
            return
        while len(pool) >= 5 - len(q):
            v = int(pool[0])
            pool = pool[1:]
            yield from extend(q + [v], pool[ds[v, pool] == 2])
    return sorted(extend([], np.arange(len(ds))))


def main():
    d, _ = octagon(np.load(ROOT / "data/generators.npy"))
    s = np.flatnonzero(d[0] == 3)
    o = np.flatnonzero(d[0] == 4)
    ds = d[np.ix_(s, s)]
    qs = five_cliques(ds)
    actual = {tuple(np.flatnonzero(d[s, v] == 1)): v for v in o}
    assert len(qs) == len(actual) == 1024 and set(qs) == set(actual)
    v = np.array([actual[q] for q in qs])
    h = np.zeros((640, 1024), dtype=np.int64)
    for i, q in enumerate(qs):
        h[list(q), i] = 1
    assert np.all(h.sum(axis=0) == 5) and np.all(h.sum(axis=1) == 8)
    a = np.zeros((1024, 1024), dtype=np.int64)
    for u in range(640):
        for i, j in combinations(np.flatnonzero(h[u]), 2):
            assert set(qs[i]) & set(qs[j]) == {u}
            distances = {int(ds[b, c]) for b in qs[i] if b != u
                         for c in qs[j] if c != u}
            assert distances in ({3}, {4})
            a[i, j] = a[j, i] = distances == {3}
    assert np.array_equal(a, d[np.ix_(v, v)] == 1)
    assert np.all(a.sum(axis=1) == 5)
    components, _ = connected_components(csr_matrix(a))
    g = np.array([64, 32, 8, 0, -1])
    b = g[ds]
    f = g[d[np.ix_(v, s)]]
    go = g[d[np.ix_(v, v)]]
    eye = np.eye(1024, dtype=np.int64)
    lap = 5 * eye - a
    assert np.array_equal(lap @ f, h.T @ b)
    assert np.array_equal(lap @ go, h.T @ f.T)
    assert np.array_equal(lap @ go @ lap, h.T @ b @ h)
    assert np.all(f.sum(axis=0) == 0)
    k = h.T @ h
    r = k - 5 * eye - a
    assert np.all((r == 0) | (r == 1))
    assert np.all(r.sum(axis=1) == 30)
    dist = shortest_path(csr_matrix(a), unweighted=True)
    assert np.all(dist[r == 1] >= 6)
    # Sparse multiplication keeps the walk-moment check inexpensive.
    asp = csr_matrix(a)
    power = eye.copy()
    traces, mixed = [], []
    for j in range(8):
        traces.append(int(np.trace(power)))
        if j <= 5:
            mixed.append(int(np.sum(k * power)))
        power = asp @ power
    assert traces == [1024, 0, 5120, 0, 46080, 0, 496640, 0]
    assert mixed == [5120, 5120, 25600, 46080, 230400, 496640]
    assert components == 2
    print("1024 cliques, 2560 outer edges, 2 components")
    print("All reconstruction, integer Laplacian and mixed-moment checks passed")


if __name__ == "__main__":
    main()
