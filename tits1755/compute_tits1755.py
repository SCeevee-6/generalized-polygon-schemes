#!/usr/bin/env python3
"""Orbital algebra of T = 2F4(2)' acting on 1755 points (ATLAS v3, TF42G1-p1755B0).

Computed from the raw permutation generators only:
  - suborbits of the point stabiliser (Schreier generators + union-find)
  - 125 intersection numbers p^k_ij (transporter-based count)
  - R1 graph check: BFS distance from 0 equals the suborbit label of every point
  - eigenmatrices P, Q and multiplicities (exact, sympy)
  - 25 Krein parameters q^k_ij (exact) and the vanishing non-trivial ones
  - Q-polynomial test over the 24 orderings (necessary triangle-zero condition)
  - fusion closure test over all 15 partitions of {1,2,3,4}

Outputs in ./out: intersection_numbers.csv, krein_parameters.csv, fusion_test.csv,
summary.json.  Generators are cached in ./cache.
Needs: sympy.  Network: brauer.maths.qmul.ac.uk (ATLAS generator files).
"""
import csv
import json
import math
import os
import re
import urllib.request
from collections import defaultdict, deque
from itertools import permutations

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
CACHE = os.path.join(HERE, "cache")
BASE = "https://brauer.maths.qmul.ac.uk/Atlas/exc/TF42/gap/TF42G1-p1755B0.g"
N = 1755


def fetch(i):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "TF42G1-p1755B0.g%d" % i)
    if not os.path.exists(path):
        with urllib.request.urlopen(BASE + str(i), timeout=60) as r:
            data = r.read()
        with open(path, "wb") as f:
            f.write(data)
    with open(path, encoding="utf-8") as f:
        return f.read()


def parse_perm(text, n=N):
    body = re.sub(r"\s+", "", text.split(":=", 1)[1].strip().rstrip(";"))
    perm = list(range(n))
    for cyc in re.findall(r"\(([^)]*)\)", body):
        idx = [int(t) - 1 for t in cyc.split(",")]
        for k, x in enumerate(idx):
            perm[x] = idx[(k + 1) % len(idx)]
    if sorted(perm) != list(range(n)):
        raise ValueError("parsed object is not a permutation")
    return perm


def perm_order(p):
    seen = [False] * len(p)
    o = 1
    for x in range(len(p)):
        if not seen[x]:
            length, y = 0, x
            while not seen[y]:
                seen[y] = True
                y = p[y]
                length += 1
            o = math.lcm(o, length)
    return o


def compose(f, g):
    """(f o g)(x) = f(g(x))."""
    return [f[y] for y in g]


def inverse(p):
    inv = [0] * len(p)
    for x, y in enumerate(p):
        inv[y] = x
    return inv


def transversal(gens, n=N):
    """trans[y] is a group element mapping 0 to y (BFS tree)."""
    trans = [None] * n
    trans[0] = list(range(n))
    q = deque([0])
    while q:
        p = q.popleft()
        for g in gens:
            y = g[p]
            if trans[y] is None:
                trans[y] = compose(g, trans[p])
                q.append(y)
    if any(t is None for t in trans):
        raise ValueError("action not transitive")
    return trans


def suborbits(gens, trans, n=N):
    """Orbits of the stabiliser of 0, from Schreier generators (union-find)."""
    inv_trans = [inverse(t) for t in trans]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in range(n):
        for g in gens:
            s = compose(inv_trans[g[p]], compose(g, trans[p]))  # s fixes 0
            for x in range(n):
                a, b = find(x), find(s[x])
                if a != b:
                    parent[a] = b
    orbs = defaultdict(list)
    for x in range(n):
        orbs[find(x)].append(x)
    return sorted(orbs.values(), key=len)


def intersection_numbers(R, trans):
    """p[(i,j,k)] = |R_i(0) cap R_j(y)| for (0,y) in R_k."""
    d = len(R)
    sets = [set(r) for r in R]
    p = {}
    for k in range(d):
        g = trans[min(R[k])]  # g(0) = y, so g(R_j(0)) = R_j(y)
        for j in range(d):
            gRj = {g[x] for x in R[j]}
            for i in range(d):
                p[(i, j, k)] = len(sets[i] & gRj)
    return p


def bfs_R1(R, trans, n=N):
    dist = [-1] * n
    dist[0] = 0
    q = deque([0])
    while q:
        u = q.popleft()
        for r in R[1]:
            v = trans[u][r]
            if dist[v] < 0:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def spectral(p, ks):
    d = len(ks)
    n = sum(ks)
    B = [sp.Matrix(d, d, lambda j, h, i=i: p[(i, j, h)]) for i in range(d)]
    if not all((B[i] * B[j] - B[j] * B[i]).is_zero_matrix
               for i in range(d) for j in range(d)):
        raise ValueError("intersection matrices do not commute")
    thetas = sorted(B[1].eigenvals().keys(), key=lambda x: -x)
    if len(thetas) != d:
        raise ValueError("B_1 does not have d distinct eigenvalues")
    P = sp.zeros(d, d)
    for l, th in enumerate(thetas):
        v = (B[1] - th * sp.eye(d)).nullspace()[0]
        v = v / v[0]
        for i in range(d):
            P[l, i] = (B[i] * v)[0]
    mult = [sp.Integer(n) / sum(P[l, i] ** 2 / ks[i] for i in range(d)) for l in range(d)]
    Q = sp.Matrix(d, d, lambda i, l: mult[l] * P[l, i] / ks[i])
    if not (P * Q - n * sp.eye(d)).is_zero_matrix:
        raise ValueError("P Q != n I")
    return thetas, mult, P, Q


def krein_all(Q, d):
    """Solve Q x = (Q_ai Q_aj)_a for x = (q^k_ij)_k."""
    Qi = Q.inv()
    kr = {}
    for i in range(d):
        for j in range(d):
            w = sp.Matrix([Q[a, i] * Q[a, j] for a in range(d)])
            x = Qi * w
            for k in range(d):
                kr[(i, j, k)] = x[k]
    return kr


def qpoly_orderings(kr, d):
    """Orderings (eigen-index sequences starting at 0) passing the triangle-zero test."""
    good = []
    for perm in permutations(range(1, d)):
        order = (0,) + perm
        bad = 0
        for a in range(d):
            for b in range(d):
                for c in range(d):
                    if abs(a - b) <= c <= a + b:
                        continue
                    if kr[(order[a], order[b], order[c])] != 0:
                        bad += 1
        if bad == 0:
            good.append(order)
    return good


def partitions(xs):
    if not xs:
        yield []
        return
    first, rest = xs[0], xs[1:]
    for smaller in partitions(rest):
        for i in range(len(smaller)):
            yield smaller[:i] + [[first] + smaller[i]] + smaller[i + 1:]
        yield [[first]] + smaller


def fusion_test(p, d):
    out = []
    for part in partitions(list(range(1, d))):
        blocks = [[0]] + part
        viol = 0
        for bi in blocks:
            for bj in blocks:
                for bk in blocks:
                    if len(bk) < 2:
                        continue
                    sums = {sum(p[(i, j, k)] for i in bi for j in bj) for k in bk}
                    if len(sums) > 1:
                        viol += 1
        out.append((part, viol))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    gens = [parse_perm(fetch(1)), parse_perm(fetch(2))]
    orders = [perm_order(g) for g in gens]
    assert orders == [2, 3], orders
    trans = transversal(gens)
    R = suborbits(gens, trans)
    ks = [len(r) for r in R]
    assert ks == [1, 10, 80, 640, 1024], ks
    d = len(R)
    cls = [0] * N
    for k, r in enumerate(R):
        for x in r:
            cls[x] = k
    p = intersection_numbers(R, trans)
    dist = bfs_R1(R, trans)
    graph_ok = all(dist[x] == cls[x] for x in range(N))
    r1_symmetric = all(0 in {trans[u][r] for r in R[1]} for u in R[1])
    p_sym = all(p[(i, j, k)] == p[(j, i, k)] for i in range(d) for j in range(d) for k in range(d))
    rowsum_ok = all(sum(p[(i, j, k)] * ks[k] for k in range(d)) == ks[i] * ks[j]
                    for i in range(d) for j in range(d))
    thetas, mult, P, Q = spectral(p, ks)
    kr = krein_all(Q, d)
    kr_nonneg = all(v >= 0 for v in kr.values())
    kr_zero_nontrivial = sorted((i, j, k) for (i, j, k), v in kr.items()
                                if v == 0 and min(i, j, k) >= 1)
    kr_diag_ok = all(kr[(i, i, 0)] == mult[i] for i in range(d))
    qord = qpoly_orderings(kr, d)
    fus = fusion_test(p, d)

    with open(os.path.join(OUT, "intersection_numbers.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["i", "j", "k", "p_k_ij"])
        for (i, j, k), v in sorted(p.items()):
            w.writerow([i, j, k, v])
    with open(os.path.join(OUT, "krein_parameters.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["i", "j", "k", "q_k_ij"])
        for (i, j, k), v in sorted(kr.items()):
            w.writerow([i, j, k, str(v)])
    with open(os.path.join(OUT, "fusion_test.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["partition_of_1_to_4", "violated_closure_conditions"])
        for part, viol in fus:
            w.writerow([str(part), viol])

    summary = {
        "suborbit_sizes": ks,
        "eigenvalues": [int(t) for t in thetas],
        "multiplicities": [str(m) for m in mult],
        "R1_BFS_matches_suborbits_at_all_points": graph_ok,
        "R1_symmetric": r1_symmetric,
        "p_symmetric_in_ij": p_sym,
        "sum_k p_ij k_k = k_i k_j": rowsum_ok,
        "krein_all_nonnegative": kr_nonneg,
        "krein_q0_ii_equals_m_i": kr_diag_ok,
        "krein_zero_nontrivial_(i,j,k)_eigen_indices": kr_zero_nontrivial,
        "Q_polynomial_orderings_found": len(qord),
        "fusions_zero_violations": [str(part) for part, viol in fus if viol == 0],
        "fusion_violation_counts": [viol for _, viol in fus],
    }
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
