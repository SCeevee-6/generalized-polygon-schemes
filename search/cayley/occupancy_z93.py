"""Occupancy model for {56,42,20;1,6,28} on Z9^3 (a test of a short proof; it does not give one).

On G = Z9^3 the order-9 atoms of S correspond to points of PG(2,3) = subgroups of order 3 of 3G
(the atom of s goes to <3s>), and the index-3 subgroups M of G to lines. Constraints used:
  - an atom [s] is invariant under translation by <3s>, so for 0 != h in L = <3s> the atoms over L
    give 6 n_L common neighbours of 0 and h, and -h is one more if h is in S: n_L <= 2 if L\\0 is in
    S, n_L <= 1 if it is at distance 2, n_L = 0 at distance 3;
  - |S cap 3G| = 56 - 6a, where a is the number of order-9 atoms in S;
  - a character with kernel M takes the value 56 - 9 off(M), off(M) = atoms outside M, so
    off(M) is 4, 6 or 7 (eigenvalues 20, 2, -7);
  - MacWilliams on the subgroup 3G: the number of characters of order 3 in E_j is
    (1/27) sum_i x_i Q_ij, with x_i = |S_i cap 3G|.
Prints the number of configurations that pass.
"""
import itertools
import sys

import numpy as np

P = np.array([[1, 56, 392, 280], [1, 20, 14, -35], [1, 2, -13, 10], [1, -7, 14, -8]], dtype=float)
Q = 729 * np.linalg.inv(P)


def points():
    pts = []
    for v in itertools.product(range(3), repeat=3):
        if any(v):
            f = next(x for x in v if x)
            v = tuple((x * (1 if f == 1 else 2)) % 3 for x in v)
            if v not in pts:
                pts.append(v)
    return pts


def main(out):
    pts = points()
    inc = np.array([[sum(p[i] * l[i] for i in range(3)) % 3 == 0 for p in pts] for l in pts], dtype=np.int64)
    alln = np.array(list(itertools.product(range(3), repeat=13)), dtype=np.int64)
    tot = alln.sum(1)
    off_to_j = {4: 1, 6: 2, 7: 3}
    count = 0
    by_a = {}
    for a in range(5, 10):
        n1 = 28 - 3 * a                       # number of order-3 subgroups of 3G inside S
        if not 0 <= n1 <= 13:
            continue
        X = alln[tot == a]
        off = a - X @ inc.T
        ok = np.isin(off, [4, 6, 7]).all(1)
        for xv, ov in zip(X[ok], off[ok]):
            if (xv == 2).sum() > n1:
                continue
            lines = {j: int((ov == o).sum()) for o, j in off_to_j.items()}
            for x2 in range(0, 27 - 2 * n1, 2):
                x3 = 26 - 2 * n1 - x2
                xs = [1, 2 * n1, x2, x3]
                N = [sum(xs[i] * Q[i, j] for i in range(4)) / 27 for j in range(4)]
                if not all(abs(N[j] - 2 * lines[j]) < 1e-6 for j in (1, 2, 3)):
                    continue
                zeros = int((xv == 0).sum())
                if x3 // 2 <= zeros and 13 - zeros <= n1 + x2 // 2:
                    count += 1
                    by_a[a] = by_a.get(a, 0) + 1
    msg = f"configurations passing: {count}; by a: {by_a}"
    open(out, "w", newline="\n").write(msg + "\n")
    print(msg)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "occupancy_z93_output.txt")
