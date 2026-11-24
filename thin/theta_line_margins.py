"""Exact rational witness for point and line margins on a grid-free seed."""

from itertools import combinations_with_replacement, product
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from theta_profiles import COS, GROUP, LINES, distance_polynomials, expand, profiles, seed_distances
from theta_incident_lines import line_types, unsupported


def quotient(source):
    """Radial adjacency quotient of the incidence graph about a given type (0=P)."""
    degree = np.array([5 if (source + i) % 2 == 0 else 3 for i in range(9)])
    q = np.zeros((9, 9), dtype=np.int64)
    for i in range(9):
        inward = 0 if i == 0 else degree[i] if i == 8 else 1
        if i:
            q[i, i - 1] = inward
        if i < 8:
            q[i, i + 1] = degree[i] - inward
    w = [np.eye(9, dtype=np.int64), q]
    w.append(q @ q - np.diag(degree))
    for i in range(3, 9):
        w.append(q @ w[-1] - (degree[:, None] - 1) * w[-2])
    assert not (w[8] % degree[:, None]).any()
    w[8] //= degree[:, None]
    return np.asarray(w)


def geometry():
    dist = seed_distances()
    lines = np.asarray(LINES)
    mixed = 2 * dist[:, lines].min(axis=2) + 1
    dual = mixed[lines, :].min(axis=1) + 1
    np.fill_diagonal(dual, 0)
    D = np.block([[2 * dist, mixed], [mixed.T, dual]])
    P, labels = expand(profiles(dist))
    # Eliminate profiles that cannot have five incident lines, iteratively.
    while True:
        bad = unsupported(P)
        if not bad:
            break
        keep = np.ones(len(P), dtype=bool)
        keep[bad] = False
        P, labels = P[keep], labels[keep]
    _, labels = np.unique(labels, return_inverse=True)
    L = line_types(P)
    # Line orbits under the same arm and endpoint symmetries.
    plook = {tuple(p): i for i, p in enumerate(P)}
    group = np.asarray([[plook[tuple(p[list(g)])] for p in P] for g in GROUP])
    keys = [min(tuple(sorted(g[line])) for g in group) for line in L]
    _, ll = np.unique(keys, axis=0, return_inverse=True)
    mixed = 2 * P[L].min(axis=1) + 1
    dual = mixed[:, lines].min(axis=2) + 1
    for i, l in enumerate(lines):
        named = sorted(plook[tuple(dist[j])] for j in l)
        dual[np.all(np.sort(L, axis=1) == named, axis=1), i] = 0
    pp = np.hstack([2 * P, 2 * P[:, lines].min(axis=2) + 1])
    lp = np.hstack([mixed, dual])
    assert np.all(pp.max(axis=1) <= 8) and np.all(lp.max(axis=1) <= 8)
    return D, P, L, labels, ll, pp, lp


def system():
    D, P, L, pl, ll, pp, lp = geometry()
    n = int(pl.max()) + 1
    m = int(ll.max()) + 1
    assert (len(P), n, len(L), m) == (753, 78, 8982, 824)
    print("profiles", len(P), n, "line types", len(L), m, flush=True)
    rows = {}

    def add(a, b, rhs):
        key = tuple((int(i), int(v)) for i, v in enumerate(np.r_[a, b]) if v)
        assert rows.get(key, rhs) == rhs, (key, rows.get(key), rhs)
        rows[key] = int(rhs)

    def pcount(weights):
        c = np.zeros(n, dtype=np.int64)
        np.add.at(c, pl, weights.astype(np.int64))
        return c

    def lcount(weights):
        c = np.zeros(m, dtype=np.int64)
        np.add.at(c, ll, weights.astype(np.int64))
        return c

    pzero, lzero = np.zeros(n, dtype=int), np.zeros(m, dtype=int)
    Q = [quotient(0), quotient(1)]
    # Cross-check the point-point-point margins against the independent
    # point-graph recurrence used for the earlier 21-point certificate.
    for k, matrix in enumerate(distance_polynomials()):
        assert np.array_equal(Q[0][2 * k, ::2, ::2], matrix)
    assert np.array_equal(Q[0][:, 0, :].sum(axis=1), [1, 5, 10, 40, 80, 320, 640, 2560, 1024])
    assert np.array_equal(Q[1][:, 0, :].sum(axis=1), [1, 3, 12, 24, 96, 192, 768, 1536, 2048])
    for i, j in combinations_with_replacement(range(len(D)), 2):
        si, sj = int(i >= 21), int(j >= 21)
        for typ, prof, count in ((0, pp, pcount), (1, lp, lcount)):
            for h in range(si ^ typ, 9, 2):
                for k in range(sj ^ typ, 9, 2):
                    rhs = Q[si][k, D[i, j], h]
                    c = count((prof[:, i] == h) & (prof[:, j] == k))
                    add(c if typ == 0 else pzero, c if typ == 1 else lzero, rhs)
    # Cubic point identities, with all repetitions.
    cos = COS[P]
    for i, j, k in combinations_with_replacement(range(21), 3):
        add(pcount(cos[:, i] * cos[:, j] * cos[:, k]), lzero, 0)
    # Incidence degree and inward line counts at every profile.
    for i, p in enumerate(P):
        incident = (L == i).sum(axis=1)
        a = np.zeros(n, dtype=int)
        a[pl[i]] = -5
        add(a, lcount(incident), 0)
        for j in range(21):
            if 0 < p[j] < 4:
                a[pl[i]] = -1
                add(a, lcount(incident * (P[L, j].min(axis=1) < p[j])), 0)
    I, J, V = [], [], []
    for i, key in enumerate(rows):
        for j, v in key:
            I.append(i)
            J.append(j)
            V.append(v)
    A = coo_matrix((V, (I, J)), shape=(len(rows), n+m), dtype=np.int64).tocsr()
    b = np.array(list(rows.values()), dtype=np.int64)
    print("matrix", A.shape, A.nnz, flush=True)
    sizes = np.r_[np.bincount(pl), np.bincount(ll)]
    assert not (12 % sizes).any()
    return A, b, sizes


def main():
    A, b, sizes = system()
    data = json.loads(Path(__file__).with_name("theta_line_margins.json").read_text())
    denominator = data["denominator"]
    assert type(denominator) is int and denominator > 0
    weights = [0] * A.shape[1]
    for key, value in data["weights"].items():
        j = int(key)
        assert 0 <= j < len(weights) and type(value) is int and value > 0
        assert weights[j] == 0
        weights[j] = value
    # Use Python integers, since these products do not fit in int64.
    for i in range(A.shape[0]):
        start, end = A.indptr[i:i + 2]
        actual = sum(int(v) * weights[j] for j, v in zip(A.indices[start:end], A.data[start:end]))
        assert actual == denominator * int(b[i]), i
    assert sum(n * int(s) for n, s in zip(weights[:78], sizes[:78])) == 1755 * denominator
    assert sum(n * int(s) for n, s in zip(weights[78:], sizes[78:])) == 2925 * denominator
    print("exact rational witness", sum(w > 0 for w in weights), "occupied orbits")
    print("all equations satisfied, total points 1755, total lines 2925")


if __name__ == '__main__':
    main()
