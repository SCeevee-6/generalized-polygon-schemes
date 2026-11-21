"""Atom bound for distance-regular Cayley graphs on abelian groups.

If Cay(G,S) is distance-regular with integral eigenvalues, S is a union of atoms
[s] = {us : u a unit mod n}, n = ord(s). If p^2 | n and h = (n/p)s, then
[s] + h = [s] (u + n/p is a unit mod n), so every element of [s] is a common
neighbour of 0 and h, and phi(n) <= |S cap (S - h)| <= max(a1, c2).
If the Sylow p-subgroup of G has exponent p^e, e >= 2, the elements whose p-part
has order < p^e form a proper subgroup, so a connected S contains such an s and
phi(p^e) <= max(a1, c2).

This script checks integrality of the eigenvalues and applies the bound to every
abelian group of order v for the five arrays that pass the abelian Krein test.
"""
import sys
from fractions import Fraction
import sympy as sp

ARRAYS = ["{52,42,16;1,6,28}", "{56,42,20;1,6,28}", "{49,40,22;1,5,28}",
          "{54,42,11;1,3,44}", "{63,48,24;1,6,28}"]


def parse(a):
    b, c = a.strip("{}").split(";")
    return [int(x) for x in b.split(",")], [int(x) for x in c.split(",")]


def params(b, c):
    d = len(b)
    k = b[0]
    bb = b + [0]
    cc = [0] + c
    a = [k - bb[i] - cc[i] for i in range(d + 1)]
    ki = [1]
    for i in range(d):
        ki.append(Fraction(ki[-1] * b[i], c[i]))
    L = sp.zeros(d + 1, d + 1)          # tridiagonal intersection matrix
    for i in range(d + 1):
        L[i, i] = a[i]
        if i > 0:
            L[i - 1, i] = cc[i]
        if i < d:
            L[i + 1, i] = bb[i]
    ev = L.charpoly().all_roots()
    return k, a, ki, ev


def partitions(n, m=None):
    if m is None:
        m = n
    if n == 0:
        yield ()
        return
    for i in range(min(n, m), 0, -1):
        for rest in partitions(n - i, i):
            yield (i,) + rest


def abelian_groups(v):
    f = sp.factorint(v)
    pieces = [[(p, lam) for lam in partitions(e)] for p, e in f.items()]
    from itertools import product
    for combo in product(*pieces):
        yield combo


def name(combo):
    out = []
    for p, lam in combo:
        for x in lam:
            out.append(f"Z{p**x}")
    return "x".join(out)


def main(out):
    lines = []
    for A in ARRAYS:
        b, c = parse(A)
        k, a, ki, ev = params(b, c)
        v = int(sum(ki))
        integral = all(e.is_integer for e in ev)
        bound = max(a[1], c[1])
        lines.append(f"{A}  v={v}  eigenvalues={[int(e) for e in ev]}  integral={integral}  max(a1,c2)={bound}")
        for combo in abelian_groups(v):
            hits = []
            for p, lam in combo:
                e = lam[0]
                if e >= 2:
                    ph = p ** (e - 1) * (p - 1)
                    hits.append((p, e, ph, ph > bound))
            killed = [h for h in hits if h[3]]
            if killed:
                p, e, ph, _ = killed[0]
                verdict = f"excluded: phi({p}^{e}) = {ph} > {bound}"
            elif hits:
                verdict = "not excluded: " + ", ".join(f"phi({p}^{e}) = {ph}" for p, e, ph, _ in hits)
            else:
                verdict = "not excluded: exponent squarefree"
            lines.append(f"  {name(combo):28s} {verdict}")
    open(out, "w", newline="\n").write("\n".join(lines) + "\n")
    return lines


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "atom_bound_output.txt")
