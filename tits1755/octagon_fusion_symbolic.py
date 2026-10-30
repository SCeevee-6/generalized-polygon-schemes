#!/usr/bin/env python3
"""Symbolic fusion test for the point-distance scheme of a generalized octagon GO(s,t).

Intersection array of a thick GO(s,t) point graph:
  b = (s(t+1), st, st, st), c = (1, 1, 1, t+1), a_i = k - b_i - c_i, valency k = s(t+1).
B_i[j,h] = p^h_{ij}; built by the three-term recursion A_1 A_i = b_{i-1}A_{i-1} + a_i A_i + c_{i+1}A_{i+1}.
At (s,t) = (2,4) the matrices are cross-checked against ./out/intersection_numbers.csv
(produced by compute_tits1755.py) when that file exists.
For each partition of {1,2,3,4} the closure differences that must vanish are listed.
Usage: python octagon_fusion_symbolic.py [--solve]   (--solve adds sympy.solve per partition; slow)
Output: ./out/octagon_fusion_symbolic.txt
"""
import csv
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
D = 5
s, t = sp.symbols("s t")
bv = [s * (t + 1), s * t, s * t, s * t, sp.Integer(0)]
cv = [sp.Integer(0), sp.Integer(1), sp.Integer(1), sp.Integer(1), t + 1]
av = [sp.Integer(0), s - 1, s - 1, s - 1, (s - 1) * (t + 1)]


def build():
    B1 = sp.zeros(D, D)
    for i in range(D):
        if i > 0:
            B1[i, i - 1] = bv[i - 1]
        B1[i, i] = av[i]
        if i < D - 1:
            B1[i, i + 1] = cv[i + 1]
    Bs = [sp.eye(D), B1]
    for i in range(1, D - 1):
        nxt = (B1 * Bs[i] - av[i] * Bs[i] - bv[i - 1] * Bs[i - 1]) / cv[i + 1]
        Bs.append(nxt.applyfunc(sp.cancel))
    return Bs


def cross_check(Bs):
    path = os.path.join(OUT, "intersection_numbers.csv")
    if not os.path.exists(path):
        return "cross-check skipped (no intersection_numbers.csv)"
    n_ok = 0
    with open(path) as f:
        for row in csv.DictReader(f):
            i, j, k, v = int(row["i"]), int(row["j"]), int(row["k"]), int(row["p_k_ij"])
            got = Bs[i].subs({s: 2, t: 4})[j, k]
            if got != v:
                raise AssertionError("mismatch at %s" % ((i, j, k),))
            n_ok += 1
    return "cross-check at (s,t)=(2,4): %d of 125 values agree" % n_ok


def partitions(xs):
    if not xs:
        yield []
        return
    first, rest = xs[0], xs[1:]
    for smaller in partitions(rest):
        for i in range(len(smaller)):
            yield smaller[:i] + [[first] + smaller[i]] + smaller[i + 1:]
        yield [[first]] + smaller


def closure_diffs(Bs, blocks):
    out = []
    for bi in blocks:
        for bj in blocks:
            for bk in blocks:
                if len(bk) < 2:
                    continue
                exprs = [sp.cancel(sum(Bs[i][j, k] for i in bi for j in bj)) for k in bk]
                for e in exprs[1:]:
                    diff = sp.factor(sp.cancel(exprs[0] - e))
                    if diff != 0:
                        out.append(diff)
    return out


def main():
    do_solve = "--solve" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    Bs = build()
    lines = [cross_check(Bs), ""]
    for part in partitions([1, 2, 3, 4]):
        blocks = [[0]] + part
        diffs = closure_diffs(Bs, blocks)
        uniq = sorted({str(d) for d in diffs})
        lines.append("partition %s : %d violated closure conditions" % (part, len(diffs)))
        for u in uniq:
            lines.append("    " + u)
        if do_solve and diffs:
            lines.append("    solutions (s,t): %s" % sp.solve(diffs, [s, t], dict=True))
        lines.append("")
    with open(os.path.join(OUT, "octagon_fusion_symbolic.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(lines[0])
    print("partitions processed: 15; output written to out/octagon_fusion_symbolic.txt")


if __name__ == "__main__":
    main()
