"""Krein test for normal Cayley graphs (connection set closed under conjugation) on any group.

If S_1 is a union of conjugacy classes, so is every distance class, and the eigenspaces are E_k = sum of the
central idempotents e_phi over a set of irreducible characters phi. Comparing coefficients of phi in E_i o E_j gives
    q^k_ij * phi(1) = sum_{chi in E_i, psi in E_j} chi(1) psi(1) <chi psi, phi>,   an integer, for every phi in E_k.
So for every k and every phi in E_k, phi(1) is a multiple of D_k := lcm of the denominators of q^k_ij (i, j),
and every q^k_ij is rational. Also phi(1) divides v, and sum_{phi in E_k} phi(1)^2 = m_k.

An array is excluded on all groups of order v if for some k >= 1 no admissible degree exists, i.e. no multiple f of
D_k with f | v and f^2 <= m_k (or some q^k_ij is irrational).
Output: normal_krein_filter.csv
"""
import csv, os
import sys
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from krein_filter import ARRAYS  # noqa: E402


def krein_all(s):
    b, c = s.split(';'); b = [int(x) for x in b.split(',')]; c = [int(x) for x in c.split(',')]
    d = len(b); k = b[0]; B = b + [0]; C = [0] + c; a = [k - B[i] - C[i] for i in range(d + 1)]
    L = sp.zeros(d + 1)
    for i in range(d + 1):
        L[i, i] = a[i]
        if i < d: L[i + 1, i] = B[i]
        if i > 0: L[i - 1, i] = C[i]
    ks = [1]
    for i in range(d): ks.append(ks[-1] * B[i] // C[i + 1])
    v = sum(ks); th = sorted(L.eigenvals().keys(), key=lambda z: -float(z)); U = []
    for t in th:
        u = [sp.Integer(1), t / k]
        for i in range(1, d): u.append(sp.simplify(((t - a[i]) * u[i] - C[i] * u[i - 1]) / B[i]))
        U.append(u)
    m = [sp.nsimplify(sp.simplify(v / sum(ks[i] * U[j][i] ** 2 for i in range(d + 1)))) for j in range(d + 1)]
    q = {(i, j, h): sp.nsimplify(sp.simplify(m[i] * m[j] / v * sum(ks[l] * U[i][l] * U[j][l] * U[h][l]
                                                                   for l in range(d + 1))))
         for i in range(d + 1) for j in range(d + 1) for h in range(d + 1)}
    return v, th, [int(x) for x in m], q, d


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    rows = []
    for s in ARRAYS:
        v, th, m, q, d = krein_all(s)
        verdict = 'passes'; info = []
        for h in range(1, d + 1):
            vals = [q[i, j, h] for i in range(d + 1) for j in range(d + 1)]
            if not all(x.is_rational for x in vals):
                info.append(f'E{h}: irrational'); verdict = f'excluded (E{h} has irrational Krein parameters)'
                continue
            D = sp.ilcm(*[sp.fraction(x)[1] for x in vals])
            degs = [f for f in sp.divisors(v) if f % D == 0 and f * f <= m[h]]
            info.append(f'E{h}: D={D}, degrees {degs}')
            if not degs and verdict == 'passes':
                verdict = f'excluded (no degree for E{h})'
        lin = [h for h in range(1, d + 1) if 'D=1,' in info[h - 1]]
        if verdict == 'passes' and not lin:
            verdict = 'excluded unless G is perfect (no E_k can hold a linear character)'
        rows.append(['{' + s + '}', v, ' '.join(map(str, m)), '; '.join(info), verdict])
    with open(os.path.join(here, 'normal_krein_filter.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['array', 'v', 'multiplicities', 'per_eigenspace', 'verdict']); w.writerows(rows)
    print(sum(r[-1] == 'passes' for r in rows), 'pass;', sum(r[-1].startswith('excluded') for r in rows), 'excluded')


if __name__ == '__main__':
    main()
