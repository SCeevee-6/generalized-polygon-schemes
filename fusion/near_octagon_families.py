"""Spectra and feasibility of the parameter families from near_octagon_note.md.

Family A: (s, m, s(m+1), (m+1)(s^2+1) - 1), fusion {1,2}|{3,4}.
Family C: (2, i, 2i, 2i+1), fusions {1,2,4}|{3} and {1,4}|{2}|{3}.
Prints the factored cubic factor of the characteristic polynomial, the multiplicities for family C as rational
functions of i, and for family A every (s, m) with s <= S, m <= M whose valencies and multiplicities are all
integers (exact arithmetic; multiplicities of irrational eigenvalues are computed in Q[X]/(f)).

Usage: python near_octagon_families.py [S M]      (defaults 12 200; writes near_octagon_families_output.txt)
"""
import os, sys
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
s, t, t2, t3, x, m, i, X = sp.symbols('s t t2 t3 x m i X')


def symbolic():
    ts = [0, t2, t3, t]
    b = [s * (t + 1)] + [s * (t - ts[j - 1]) for j in range(1, 4)]
    c = [u + 1 for u in ts]
    a = [sp.expand(b[0] - b[j] - (c[j - 1] if j > 0 else 0)) for j in range(4)] + [b[0] - c[3]]
    v = [sp.Integer(1), x]
    for j in range(1, 4):
        v.append(sp.cancel(((x - a[j]) * v[j] - b[j - 1] * v[j - 1]) / c[j]))
    mp = sp.numer(sp.together((x - a[4]) * v[4] - b[3] * v[3]))
    g = sp.cancel(mp / ((x - s * (t + 1)) * (x + t + 1)))
    g = sp.Poly(g, x)
    return v, sp.expand(g.as_expr() / g.LC())


def exact_mults(b, c):
    d = len(b)
    k = [sp.Integer(1)]
    for j in range(d):
        k.append(k[-1] * sp.Rational(b[j], c[j]))
    a = [b[0] - b[j] - (c[j - 1] if j > 0 else 0) for j in range(d)] + [b[0] - c[-1]]
    u = [sp.Integer(1), X / b[0]]
    for j in range(1, d):
        u.append(sp.expand(((X - a[j]) * u[j] - c[j - 1] * u[j - 1]) / b[j]))
    mp = sp.expand((X - a[d]) * u[d] - c[d - 1] * u[d - 1])
    S = sp.Poly(sp.expand(sum(k[j] * u[j] ** 2 for j in range(d + 1))), X)
    out = []
    for f, _ in sp.factor_list(mp, X)[1]:
        fp = sp.Poly(f, X)
        if fp.degree() == 1:
            th = -fp.all_coeffs()[1] / fp.all_coeffs()[0]
            out.append((th, sum(k) / S.eval(th)))
        else:
            inv = sp.invert(S.rem(fp), fp)
            out.append((fp.as_expr(), sp.Poly(sum(k) * inv.as_expr(), X).rem(fp).as_expr()))
    return k, out


def main(S=12, M=200):
    out = []
    v, g = symbolic()
    subA = {t2: m, t3: s * (m + 1), t: (m + 1) * (s ** 2 + 1) - 1}
    subC = {s: 2, t2: i, t3: 2 * i, t: 2 * i + 1}
    out.append('family A: g = ' + str(sp.factor(g.subs(subA, simultaneous=True))))
    out.append('family C: g = ' + str(sp.factor(g.subs(subC, simultaneous=True))))
    ks = [sp.cancel(w.subs(x, s * (t + 1)).subs(subC, simultaneous=True)) for w in v]
    out.append('family C: valencies ' + str([sp.factor(k_) for k_ in ks]))
    for th in (2 * i + 3, i + 1, -1, -2 * i - 2):
        den = sum(sp.cancel(w.subs(subC, simultaneous=True).subs(x, th)) ** 2 / ks[j] for j, w in enumerate(v))
        out.append(f'family C: multiplicity of {th}: {sp.factor(sp.cancel(sum(ks) / den))}')
    good = []
    n = 0
    for sv in range(2, S + 1):
        for mv in range(0, M + 1):
            tt = (mv + 1) * (sv * sv + 1) - 1
            ts = [0, mv, sv * (mv + 1), tt]
            b = [sv * (tt + 1)] + [sv * (tt - ts[j - 1]) for j in range(1, 4)]
            c = [u + 1 for u in ts]
            k = [sp.Integer(1)]
            for j in range(4):
                k.append(k[-1] * sp.Rational(b[j], c[j]))
            if any(not kk.is_Integer for kk in k):
                continue
            n += 1
            _, mm = exact_mults(b, c)
            if all(sp.sympify(q).is_Integer for _, q in mm):
                good.append((sv, mv, int(sum(k)), [str(th) for th, _ in mm], [int(q) for _, q in mm]))
    out.append(f'family A, 2 <= s <= {S}, 0 <= m <= {M}: {n} parameter sets with integral valencies, '
               f'{len(good)} with integral multiplicities:')
    for gd in good:
        out.append('  (s, m) = (%d, %d): v = %d, eigenvalues %s, multiplicities %s' % (gd[0], gd[1], gd[2], gd[3], gd[4]))
    txt = '\n'.join(out)
    open(os.path.join(HERE, 'near_octagon_families_output.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main(*(int(q) for q in sys.argv[1:3]))
