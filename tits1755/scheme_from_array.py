"""Orbital algebra of the Tits group T = 2F4(2)' acting on the 1755 points of GO(2,4).

Reproduces in exact arithmetic:
  1. valencies k_i and all 125 intersection numbers p^k_ij from the intersection array
  2. eigenvalues, multiplicities, eigenmatrices P, Q and Krein parameters q^k_ij
  3. Q-polynomial test over all 24 orderings of E_1..E_4
  4. Bose-Mesner closure test for all 15 set partitions of {1,2,3,4}
  5. (--atlas) direct verification against the ATLAS v3 permutation representation
     TF42G1-p1755B0 (standard generators, GAP format)

Independent of compute_tits1755.py: derives everything from the intersection array alone
(stdlib, exact Fractions); --atlas then checks the derived table against the permutation action.

Usage:  python scheme_from_array.py [--atlas] [--out DIR]      (default DIR: out_array)
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import urllib.request
from collections import deque
from fractions import Fraction as Fr
from itertools import permutations

# Collinearity graph of GO(s,t): b = (s(t+1), st, st, st), c = (1, 1, 1, t+1)
S, T = 2, 4
D = 4
B = [Fr(S * (T + 1)), Fr(S * T), Fr(S * T), Fr(S * T), Fr(0)]
C = [Fr(0), Fr(1), Fr(1), Fr(1), Fr(T + 1)]
N = (S + 1) * (S * T + 1) * (S * S * T * T + 1)  # 1755
ATLAS = "https://brauer.maths.qmul.ac.uk/Atlas/exc/TF42/gap/TF42G1-p1755B0.g{}"
R = range(D + 1)


def parameters():
    a = [B[0] - B[i] - C[i] for i in R]
    k = [Fr(1)]
    for i in range(D):
        k.append(k[i] * B[i] / C[i + 1])
    assert sum(k) == N
    return a, k


def intersection_numbers(a):
    """p[i][j][h] = p^h_ij from A_1 A_i = b_{i-1} A_{i-1} + a_i A_i + c_{i+1} A_{i+1}."""
    g = lambda v, x: v[x] if 0 <= x <= D else Fr(0)
    p = [[[Fr(0)] * (D + 1) for _ in R] for _ in R]
    for j in R:
        rows = [[Fr(0)] * (D + 1) for _ in R]
        rows[0][j] = Fr(1)
        if j < D:
            rows[1][j + 1] += C[j + 1]
        rows[1][j] += a[j]
        if j > 0:
            rows[1][j - 1] += B[j - 1]
        for i in range(1, D):
            rows[i + 1] = [(C[h] * g(rows[i], h - 1) + (a[h] - a[i]) * g(rows[i], h)
                            + B[h] * g(rows[i], h + 1) - B[i - 1] * g(rows[i - 1], h)) / C[i + 1]
                           for h in R]
        for i in R:
            p[i][j] = rows[i][:]
    for i in R:
        for j in R:
            assert all(x.denominator == 1 and x >= 0 for x in p[i][j])
    return p


def check_p(p, k):
    for i in R:
        for j in R:
            assert p[i][j] == p[j][i]
            assert sum(p[i][j][h] * k[h] for h in R) == k[i] * k[j]


def spectrum(a, k):
    """Integer eigenvalues of L1 (tridiagonal), first eigenmatrix P[l][i], multiplicities, Q[i][l]."""
    found = []
    for th in range(int(B[0]), -int(B[0]) - 1, -1):
        v = [Fr(1), Fr(th)]
        for i in range(1, D):
            v.append(((th - a[i]) * v[i] - B[i - 1] * v[i - 1]) / C[i + 1])
        if B[D - 1] * v[D - 1] + a[D] * v[D] == th * v[D]:  # last row of A_1 A_D
            found.append((th, v))
    assert len(found) == D + 1, "non-integral eigenvalues"
    thetas = [th for th, _ in found]
    P = [v for _, v in found]
    m = [N / sum(v[i] ** 2 / k[i] for i in R) for v in P]
    assert sum(m) == N and all(x.denominator == 1 for x in m)
    Q = [[m[l] * P[l][i] / k[i] for l in R] for i in R]
    return thetas, P, m, Q


def krein(P, m, k):
    return [[[m[i] * m[j] / N * sum(P[i][t] * P[j][t] * P[h][t] / k[t] ** 2 for t in R)
              for h in R] for j in R] for i in R]


def q_polynomial_orderings(q):
    out = []
    for perm in permutations(range(1, D + 1)):
        o = (0,) + perm
        if all((q[o[1]][o[x]][o[y]] == 0) == (abs(x - y) > 1) for x in R for y in R):
            out.append(o)
    return out


def set_partitions(xs):
    if len(xs) == 1:
        yield [xs]
        return
    for sm in set_partitions(xs[1:]):
        for n, sub in enumerate(sm):
            yield sm[:n] + [[xs[0]] + sub] + sm[n + 1:]
        yield [[xs[0]]] + sm


def fusion_violations(p, part):
    """Number of (block,block,block) triples where sum_{i in Bm, j in Bn} p^h_ij depends on h in Bk."""
    blocks = [[0]] + part
    return sum(len({sum(p[i][j][h] for i in bm for j in bn) for h in bk}) > 1
               for bm in blocks for bn in blocks for bk in blocks)


# ---------------- ATLAS verification ----------------
def _inv(f):
    g = [0] * len(f)
    for x, y in enumerate(f):
        g[y] = x
    return g


def _parse_gap_perm(txt):
    body = re.sub(r"\s+", "", txt.split(":=", 1)[1]).rstrip(";")
    perm = list(range(N + 1))
    for cyc in re.findall(r"\(([^()]+)\)", body):
        pts = [int(x) for x in cyc.split(",")]
        for u, v in zip(pts, pts[1:] + pts[:1]):
            perm[u] = v
    return perm


def atlas_check(p):
    gens = [_parse_gap_perm(urllib.request.urlopen(ATLAS.format(i), timeout=60).read().decode())
            for i in (1, 2)]
    tr = {1: list(range(N + 1))}           # transversal: tr[y](1) = y
    dq = deque([1])
    while dq:
        y = dq.popleft()
        for g in gens:
            if g[y] not in tr:
                tr[g[y]] = [g[tr[y][x]] for x in range(N + 1)]
                dq.append(g[y])
    assert len(tr) == N
    rinv = {y: _inv(r) for y, r in tr.items()}
    seen, orbits = set(), []                # orbits of Stab(1) via Schreier generators
    for st in range(1, N + 1):
        if st in seen:
            continue
        orb, dq = {st}, deque([st])
        while dq:
            x = dq.popleft()
            for y in range(1, N + 1):
                for g in gens:
                    z = rinv[g[y]][g[tr[y][x]]]
                    if z not in orb:
                        orb.add(z)
                        dq.append(z)
        seen |= orb
        orbits.append(orb)
    sizes = sorted(len(o) for o in orbits)
    assert len(set(sizes)) == len(sizes)    # distinct subdegrees -> label by size
    label = [None] * (N + 1)
    for o in orbits:
        for x in o:
            label[x] = sizes.index(len(o))
    typ = lambda u, v: label[rinv[u][v]]
    reps = {}
    for x in range(1, N + 1):
        reps.setdefault(label[x], x)
    pc = [[[0] * (D + 1) for _ in R] for _ in R]
    for h, y in reps.items():
        for z in range(1, N + 1):
            pc[label[z]][typ(z, y)][h] += 1
    match = all(pc[i][j][h] == p[i][j][h] for i in R for j in R for h in R)
    dist, dq = {1: 0}, deque([1])
    while dq:
        u = dq.popleft()
        for v in range(1, N + 1):
            if v not in dist and typ(u, v) == 1:
                dist[v] = dist[u] + 1
                dq.append(v)
    return sizes, match, all(dist[v] == label[v] for v in range(1, N + 1))


def _write(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", action="store_true", help="verify against the ATLAS permutation representation")
    ap.add_argument("--out", default="out_array")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    a, k = parameters()
    p = intersection_numbers(a)
    check_p(p, k)
    thetas, P, m, Q = spectrum(a, k)
    q = krein(P, m, k)
    assert all(x >= 0 for qi in q for qij in qi for x in qij)
    qzeros = [(i, j, h) for i in range(1, D + 1) for j in range(i, D + 1) for h in range(1, D + 1) if q[i][j][h] == 0]
    qpoly = q_polynomial_orderings(q)
    fus = [(pt, len(pt) + 1, fusion_violations(p, pt)) for pt in set_partitions([1, 2, 3, 4])]

    o = args.out
    _write(f"{o}/orbital_valencies.csv", ["i", "k_i", "a_i", "b_i", "c_i"], [[i, k[i], a[i], B[i], C[i]] for i in R])
    _write(f"{o}/intersection_numbers.csv", ["i", "j", "k", "p^k_ij"],
           [[i, j, h, p[i][j][h]] for i in R for j in R for h in R])
    _write(f"{o}/P_matrix.csv", ["E_l", "theta", "m_l"] + [f"A_{i}" for i in R], [[l, thetas[l], m[l]] + P[l] for l in R])
    _write(f"{o}/Q_matrix.csv", ["A_i"] + [f"E_{l}" for l in R], [[i] + [str(x) for x in Q[i]] for i in R])
    _write(f"{o}/krein_parameters.csv", ["i", "j", "k", "q^k_ij"],
           [[i, j, h, str(q[i][j][h])] for i in R for j in R for h in R])
    _write(f"{o}/fusion_results.csv", ["partition", "classes", "violations", "scheme"],
           [[str(pt), c, v, v == 0] for pt, c, v in fus])

    print(f"n={N} k={[int(x) for x in k]} theta={thetas} m={[int(x) for x in m]}")
    print("p^k_ij: 125 nonneg integers, symmetric, sum_k p^k_ij k_k = k_i k_j: OK")
    print(f"Krein >= 0: OK; nontrivial zeros (i<=j): {qzeros}; Q-polynomial orderings: {qpoly}")
    print(f"closing fusions: {[pt for pt, _, v in fus if v == 0]} (of {len(fus)})")
    if args.atlas:
        sizes, match, dist_ok = atlas_check(p)
        print(f"ATLAS: subdegrees={sizes} p^k_ij match={match} R1 distance == orbital label: {dist_ok}")


if __name__ == "__main__":
    main()
