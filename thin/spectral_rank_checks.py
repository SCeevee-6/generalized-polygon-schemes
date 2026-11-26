"""Spectral and rank checks recorded in spectral_rank_checks.md."""

import os
import random
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "flags"))
import go24_common as gc  # noqa: E402
import indexed_twelve_metrics as itm  # noqa: E402
from grid_path_closure import LineConstraints  # noqa: E402
from opposite_link_surfaces import failed_metrics  # noqa: E402

SHELLS = (1, 10, 80, 640, 1024)
B, C = (10, 8, 8, 8, 0), (0, 1, 1, 1, 5)


def cosines(theta):
    w = [Fraction(1), Fraction(theta, 10)]
    for i in range(1, 4):
        a = 10 - B[i] - C[i]
        w.append((theta*w[i] - C[i]*w[i-1] - a*w[i]) / B[i])
    assert C[4]*w[3] + 5*w[4] == theta*w[4]
    return w


def gegenbauer(k, lam, x):
    p0, p1 = Fraction(1), 2*lam*x
    if k == 0:
        return p0
    for n in range(1, k):
        p0, p1 = p1, (2*(n+lam)*x*p1 - (n+2*lam-1)*p0) / (n+1)
    return p1


def design_strength():
    w = cosines(-5)
    mult = Fraction(1755) / sum(k*x*x for k, x in zip(SHELLS, w))
    assert mult == 78
    lam = Fraction(78-2, 2)
    sums = [sum(k*gegenbauer(d, lam, x) for k, x in zip(SHELLS, w)) for d in range(1, 6)]
    assert sums[:3] == [0, 0, 0] and sums[3] > 0
    return sums


def psd(matrix):
    a, previous, rank = matrix.astype(object), 1, 0
    while a.size:
        d = np.diag(a)
        if (d < 0).any() or a[d == 0].any():
            return False, rank
        a = a[np.ix_(d > 0, d > 0)]
        if not a.size:
            break
        p = a[0, 0]
        a, previous, rank = (p*a[1:, 1:] - np.outer(a[1:, 0], a[0, 1:])) // previous, p, rank+1
    return True, rank


def rank2(matrix):
    basis, rank = {}, 0
    for row in matrix:
        v = int("".join(map(str, (row % 2)[::-1])), 2)
        while v:
            h = v.bit_length() - 1
            if h not in basis:
                basis[h], rank = v, rank + 1
                break
            v ^= basis[h]
    return rank


def witness_gram():
    scaled = {th: [int(x*320) for x in cosines(th)] for th in (5, 1, -3, -5)}
    mult = {5: 351, 1: 650, -3: 675, -5: 78}
    constraints = LineConstraints()
    metrics = (1 << failed_metrics(constraints)).astype(np.uint8)
    rows = []
    for word, choices in sorted(itm.WITNESSES.items()):
        masks, _ = itm.check_case(word, choices, metrics, constraints)
        keep = [i for i in range(len(masks)) if not any(masks[i, j] == 1 for j in range(i))]
        masks = masks[np.ix_(keep, keep)].astype(int)
        single = (masks & (masks - 1)) == 0
        rng, clique = random.Random(0), []
        for _ in range(30):  # greedy cliques of determined pairs, seeded orders
            order, current = list(range(len(masks))), []
            rng.shuffle(order)
            for v in order:
                if all(single[v, u] for u in current):
                    current.append(v)
            clique = max(clique, current, key=len)
        dist = np.log2(masks[np.ix_(clique, clique)]).astype(int)
        ranks = {}
        for th, seq in scaled.items():
            ok, r = psd(np.array(seq, dtype=object)[dist])
            assert ok and r <= mult[th]
            ranks[th] = r
        rows.append((word, len(keep), len(clique), ranks[-5]))
    return rows


def ree_tits():
    gens = [np.asarray(g) for g in np.load("data/generators.npy", allow_pickle=True)]
    d = gc.point_distances(gens)
    lines = np.array(sorted(map(sorted, gc.lines_from_distances(d))))
    to_point = np.min(np.stack([d[lines[:, k]] for k in range(3)]), axis=0)
    dl = np.min(np.stack([to_point[:, lines[:, k]] for k in range(3)]), axis=0) + 1
    np.fill_diagonal(dl, 0)
    assert np.bincount(dl[0]).tolist() == [1, 12, 96, 768, 2048]
    cycles = Counter()
    for c in np.nonzero(dl[0] == 2)[0]:
        bs = np.nonzero((dl[0] == 4) & (dl[c] == 2))[0]
        ds = np.nonzero((dl[c] == 4) & (dl[0] == 2))[0]
        m = (dl[np.ix_(bs, ds)] == 2).astype(int)
        assert m.shape == (64, 64) and (m.sum(0) == 2).all() and (m.sum(1) == 2).all()
        # A 2-regular bipartite graph is a union of 4-cycles iff M M^T = 2(I + P).
        mm = m @ m.T
        cycles[int((mm == 2).sum() - 64) // 64, rank2(m)] += 1
    eye = np.eye(len(dl), dtype=np.int64)
    a1, a2 = (dl == 1).astype(np.int64), (dl == 2).astype(np.int64)
    line_ranks = {"A2": rank2(a2), "A2+I": rank2(a2 + eye),
                  "A2+A1": rank2(a2 + a1), "A2+A1+I": rank2(a2 + a1 + eye)}
    point_ranks = {"A+I": rank2((d == 1).astype(np.int64) + np.eye(1755, dtype=np.int64)),
                   "A4": rank2((d == 4).astype(np.int64))}
    return cycles, line_ranks, point_ranks


def main():
    print("Gegenbauer sums k=1..5:", [str(x) for x in design_strength()])
    for row in witness_gram():
        print("witness %s: %d points, determined clique %d, rank at -5: %d" % row)
    cycles, line_ranks, point_ranks = ree_tits()
    assert dict(cycles) == {(1, 32): 96}
    assert min(line_ranks.values()) >= 48
    print("Ree-Tits: all 96 indexed systems at a line are 32 four-cycles (2-rank 32)")
    print("line 2-ranks", line_ranks, "point 2-ranks", point_ranks)


if __name__ == "__main__":
    main()
