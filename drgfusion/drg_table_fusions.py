"""Closed fusions of the distance schemes of all intersection arrays in Brouwer's DRG tables.

For each array {b_0..b_{d-1}; c_1..c_d} the intersection numbers p^k_ij are computed exactly (Fractions)
from the tridiagonal recurrence A_1 A_i = b_{i-1} A_{i-1} + a_i A_i + c_{i+1} A_{i+1}. Every set partition
of {1..d} with at least two blocks is tested for closure (block sums of p^h_ij constant on each block).
Closed fusions with exactly two blocks give strongly regular graphs; their parameters are recorded.

Usage: python drg_table_fusions.py <table_dir> <out.csv>

Arrays of diameter above 7 are skipped (the number of set partitions grows as the Bell numbers).
"""
import os, re, sys, csv
from fractions import Fraction as F
from itertools import product


def parse(line):
    m = re.fullmatch(r'\s*i\(([\d,\s]+);([\d,\s]+)\)\s*', line)
    b = [int(x) for x in m.group(1).split(',')]
    c = [int(x) for x in m.group(2).split(',')]
    return b, c


def intersection_numbers(b, c):
    d = len(b)
    bb = b + [0]
    cc = [0] + c
    k = b[0]
    a = [k - bb[i] - cc[i] for i in range(d + 1)]
    # A_i as polynomials in A_1: v_0 = 1, v_1 = x, c_{i+1} v_{i+1} = (x - a_i) v_i - b_{i-1} v_{i-1}
    v = [[F(1)], [F(0), F(1)]]
    for i in range(1, d):
        nxt = [F(0)] * (i + 2)
        for e, co in enumerate(v[i]):
            nxt[e + 1] += co
            nxt[e] -= a[i] * co
        for e, co in enumerate(v[i - 1]):
            nxt[e] -= bb[i - 1] * co
        v.append([co / cc[i + 1] for co in nxt])
    # tridiagonal matrix L (A_1 on basis A_0..A_d): A_1 A_j = b_{j-1} A_{j-1} + a_j A_j + c_{j+1} A_{j+1}
    n = d + 1
    L = [[F(0)] * n for _ in range(n)]
    for j in range(n):
        if j > 0: L[j - 1][j] = F(bb[j - 1])
        L[j][j] = F(a[j])
        if j < d: L[j + 1][j] = F(cc[j + 1])
    # B_i = v_i(L); then p^k_ij = (B_i)[k][j]  (A_i A_j = sum_k p^k_ij A_k)
    def matmul(X, Y):
        return [[sum(X[r][t] * Y[t][s] for t in range(n)) for s in range(n)] for r in range(n)]
    I = [[F(int(r == s)) for s in range(n)] for r in range(n)]
    powers = [I]
    for _ in range(d):
        powers.append(matmul(powers[-1], L))
    B = []
    for i in range(n):
        M = [[sum(v[i][e] * powers[e][r][s] for e in range(len(v[i]))) for s in range(n)] for r in range(n)]
        B.append(M)
    p = [[[B[i][kk][j] for kk in range(n)] for j in range(n)] for i in range(n)]  # p[i][j][k]
    kv = [1]
    for i in range(d):
        kv.append(kv[-1] * F(b[i]) / F(c[i]))
    return p, kv


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def closed(p, blocks):
    blocks = [[0]] + blocks
    for A in blocks:
        for Bk in blocks:
            for C in blocks:
                vals = {sum(p[i][j][h] for i in A for j in Bk) for h in C}
                if len(vals) > 1:
                    return False
    return True


def main(tdir, out, maxd=7):
    rows = []
    for fn in sorted(os.listdir(tdir)):
        for line in open(os.path.join(tdir, fn)).read().splitlines():
            if not line.strip():
                continue
            b, c = parse(line)
            d = len(b)
            if d > maxd:
                continue
            p, kv = intersection_numbers(b, c)
            n = sum(kv)
            assert all(x.denominator == 1 for x in kv), line
            for part in set_partitions(list(range(1, d + 1))):
                if len(part) in (1, d):
                    continue
                part = sorted(sorted(bl) for bl in part)
                if not closed(p, part):
                    continue
                info = ''
                if len(part) == 2:
                    S = part[0] if 1 in part[0] else part[1]
                    S = part[0]
                    kS = sum(kv[i] for i in S)
                    lam = sum(p[i][j][S[0]] for i in S for j in S)
                    other = part[1]
                    mu = sum(p[i][j][other[0]] for i in S for j in S)
                    info = f'srg({n},{kS},{lam},{mu})'
                    S2 = part[1]; k2 = sum(kv[i] for i in S2)
                    lam2 = sum(p[i][j][S2[0]] for i in S2 for j in S2); mu2 = sum(p[i][j][S[0]] for i in S2 for j in S2)
                    info += f' / srg({n},{k2},{lam2},{mu2})'
                rows.append((fn, line.strip(), d, int(n), '|'.join(''.join(map(str, bl)) for bl in part), info))
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['table', 'array', 'd', 'v', 'fusion', 'srg_parameters'])
        w.writerows(rows)
    print(f'{len(rows)} closed proper fusions written to {out}')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
