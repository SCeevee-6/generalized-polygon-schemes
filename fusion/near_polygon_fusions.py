"""Fusions of the distance schemes of regular near hexagons and near octagons with symbolic parameters.

A regular near 2d-gon with parameters (s, t_2, ..., t_{d-1}, t) has intersection numbers
c_i = t_i + 1, b_i = s(t - t_i), a_i = (s - 1)(t_i + 1), with t_1 = 0 and t_d = t. The p^k_ij are computed as
rational functions of the parameters from the array (A_i as polynomials in A_1, reduced modulo the minimal
polynomial). For every partition of {1,...,d} the closure conditions are cleared of denominators and a Groebner
basis over QQ (grevlex) is computed; the factored basis is reported. As an independent check every partition is
also tested at all integer points of a box (s <= S, 0 <= t_2 <= ... <= t <= T, integral valencies).

Usage: python near_polygon_fusions.py [S T]     (writes near_polygon_fusions_output.txt next to the script)
"""
import os, sys, itertools, time
from fractions import Fraction as Fr
import sympy as sp

s, t, t2, t3, x = sp.symbols('s t t2 t3 x')
HERE = os.path.dirname(os.path.abspath(__file__))


def array(d):
    ts = {3: [0, t2, t], 4: [0, t2, t3, t]}[d]  # t_1 .. t_d
    b = [s * (t + 1)] + [s * (t - ts[i - 1]) for i in range(1, d)]  # b_0 .. b_{d-1}, b_i = s(t - t_i)
    c = [None] + [ti + 1 for ti in ts]  # c_1 .. c_d
    return b, c


def scheme(d):
    """p[k, i, j] = p^k_ij as rational functions: B_1 is the tridiagonal intersection matrix and
    B_{i+1} = (B_1 B_i - a_i B_i - b_{i-1} B_{i-1}) / c_{i+1}; then p^k_ij = (B_i)[k, j]."""
    b, c = array(d)
    k = b[0]
    bb = b + [0]
    cc = [0] + c[1:]
    a = [sp.expand(k - bb[i] - cc[i]) for i in range(d + 1)]
    n = d + 1
    L = sp.zeros(n, n)
    for j in range(n):
        if j > 0:
            L[j - 1, j] = bb[j - 1]
        L[j, j] = a[j]
        if j < d:
            L[j + 1, j] = cc[j + 1]
    B = [sp.eye(n), L]
    for i in range(1, d):
        M = (L * B[i] - a[i] * B[i] - bb[i - 1] * B[i - 1]) / cc[i + 1]
        B.append(M.applyfunc(sp.cancel))
    p = {}
    for i in range(n):
        for j in range(n):
            for kk in range(n):
                p[kk, i, j] = sp.factor(B[i][kk, j])
    return p


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def conditions(p, part):
    bl = [[0]] + part
    eqs = set()
    for A in bl:
        for B in bl:
            for C in bl:
                vals = [sp.together(sum(p[h, i, j] for i in A for j in B)) for h in C]
                for v_ in vals[1:]:
                    num = sp.expand(sp.numer(sp.together(v_ - vals[0])))
                    if num != 0:
                        eqs.add(num)
    return sorted(eqs, key=sp.default_sort_key)


def int_array(d, sv, ts, tv):
    tl = [0] + list(ts) + [tv]
    b = [sv * (tv + 1)] + [sv * (tv - tl[i - 1]) for i in range(1, d)]
    c = [ti + 1 for ti in tl]
    return b, c


def p_numeric(b, c):
    """exact Fractions p[k][i][j] from an integer array (tridiagonal recurrence)."""
    d = len(b)
    bb = b + [0]; cc = [0] + c
    k = b[0]
    a = [k - bb[i] - cc[i] for i in range(d + 1)]
    n = d + 1
    L = [[Fr(0)] * n for _ in range(n)]
    for j in range(n):
        if j > 0: L[j - 1][j] = Fr(bb[j - 1])
        L[j][j] = Fr(a[j])
        if j < d: L[j + 1][j] = Fr(cc[j + 1])
    B = [[[Fr(int(r == q)) for q in range(n)] for r in range(n)], L]
    def mm(X, Y):
        return [[sum(X[r][u] * Y[u][q] for u in range(n)) for q in range(n)] for r in range(n)]
    for i in range(1, d):
        Ai = mm(L, B[i])
        B.append([[(Ai[r][q] - a[i] * B[i][r][q] - bb[i - 1] * B[i - 1][r][q]) / cc[i + 1] for q in range(n)]
                  for r in range(n)])
    return lambda kk, i, j: B[i][kk][j]


def closed_num(pf, part):
    bl = [[0]] + part
    for A in bl:
        for Bk in bl:
            for C in bl:
                if len({sum(pf(h, i, j) for i in A for j in Bk) for h in C}) > 1:
                    return False
    return True


def scan(d, S, T):
    hits = {}
    for sv in range(1, S + 1):
        for tv in range(0, T + 1):
            for ts in itertools.combinations_with_replacement(range(0, tv + 1), d - 2):
                b, c = int_array(d, sv, ts, tv)
                if any(x <= 0 for x in b):
                    continue
                kv = [Fr(1)]
                for i in range(d):
                    kv.append(kv[-1] * b[i] / c[i])
                if any(x.denominator != 1 for x in kv):
                    continue
                pf = p_numeric(b, c)
                for part in set_partitions(list(range(1, d + 1))):
                    if len(part) in (1, d):
                        continue
                    part = sorted(sorted(q) for q in part)
                    if closed_num(pf, part):
                        hits.setdefault(str(part), []).append((sv,) + tuple(ts) + (tv,))
    return hits


def main(S=12, T=40):
    out = []
    for d in (3, 4):
        t0 = time.time()
        p = scheme(d)
        vars_ = [s, t2, t] if d == 3 else [s, t2, t3, t]
        out.append(f'## d = {d}, parameters {vars_}; valencies ' + ', '.join(str(sp.factor(p[0, i, i])) for i in range(d + 1)))
        for part in set_partitions(list(range(1, d + 1))):
            if len(part) in (1, d):
                continue
            part = sorted(sorted(q) for q in part)
            eqs = conditions(p, part)
            G = sp.groebner(eqs, *vars_, order='grevlex', domain=sp.QQ)
            fac = [sp.factor(g) for g in G.exprs]
            out.append(f'{part}: ' + ('no solution (basis {1})' if G.exprs == [1] else '; '.join(str(f) for f in fac)))
        out.append(f'(symbolic part {time.time() - t0:.0f}s)')
        hits = scan(d, S, T)
        out.append(f'integer scan s <= {S}, t <= {T} (t_i nondecreasing, valencies integral, b_i > 0):')
        for k_, v_ in sorted(hits.items()):
            out.append(f'  {k_}: {len(v_)} points, e.g. {v_[:12]}')
        out.append('')
    txt = '\n'.join(out)
    open(os.path.join(HERE, 'near_polygon_fusions_output.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main(*(int(a) for a in sys.argv[1:3]))
