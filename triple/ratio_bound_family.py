"""Ratio bound for sets of pairwise non-opposite points in a generalized octagon of order (s, t).

Opposite = distance 4 in the collinearity graph. The eigenvalues of A_4 follow from the standard sequences
of A_1 (eigenvalues s(t+1), s-1+sqrt(2st), s-1, s-1-sqrt(2st), -(t+1)). For t = s^2 (Ree-Tits orders) the
script compares the ratio bound n(-l_min)/(k_4 - l_min) with the size (s+1)(st+1) of a line neighbourhood
(the points at distance <= 1 from a line). Usage: python triple/ratio_bound_family.py
"""
import sympy as sp

s = sp.symbols('s', positive=True)
for label, t in (('t = s^2', s ** 2),):
    b = [s * (t + 1), s * t, s * t, s * t]; c = [1, 1, 1, t + 1]
    bb = b + [0]; cc = [0] + c; a = [b[0] - bb[i] - cc[i] for i in range(5)]
    kv = [sp.Integer(1)]
    for i in range(4):
        kv.append(sp.simplify(kv[-1] * b[i] / c[i]))
    n = sp.factor(sum(kv))
    r = sp.sqrt(2 * s * t)
    ths = [s * (t + 1), s - 1 + r, s - 1, s - 1 - r, -(t + 1)]
    ev = []
    for th in ths:
        u = [sp.Integer(1), th / b[0]]
        for i in range(1, 4):
            u.append(((th - a[i]) * u[i] - cc[i] * u[i - 1]) / bb[i])
        ev.append(sp.simplify(sp.expand(kv[4] * u[4])))
    print(label, ': eigenvalues of A_4:', ev)
    for q in (2, 8, 32, 128):
        vals = [sp.nsimplify(sp.simplify(e.subs(s, q))) for e in ev]
        lmin = min(vals[1:])
        bound = sp.simplify(n.subs(s, q) * (-lmin) / (vals[0] - lmin))
        print(f'  s = {q}: A_4 spectrum {vals}; ratio bound {bound}; (s+1)(st+1) = {(q + 1) * (q ** 3 + 1)}')
    # symbolic: assume the minimum is the eigenvalue at theta = s-1+sqrt(2st) or s-1-sqrt(2st)
    for idx in (1, 3, 4):
        bnd = sp.factor(sp.simplify(n * (-ev[idx]) / (ev[0] - ev[idx])))
        print(f'  ratio bound if l_min = eigenvalue #{idx}: {bnd}')
