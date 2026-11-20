"""Krein integrality test for the 18 arrays with no known status (search/drg/drg_note.md).

If a distance-regular graph is a Cayley graph on an abelian group, its distance scheme is a translation scheme,
the dual scheme lives on the character group, and the Krein parameters q^k_ij are its intersection numbers.
So every q^k_ij must be a nonnegative integer. This is a known lemma (see e.g. the paper on Fourier analysis
on distance-regular Cayley graphs over abelian groups, Electron. J. Combin. 33(1) (2026) #P1.14).

Exact arithmetic with sympy. Output: krein_filter.csv (array, v, eigenvalues, multiplicities, number of
non-integral Krein parameters, abelian Cayley possible yes/no).
"""
import csv, os
import sympy as sp

ARRAYS = ["18,14,5;1,2,14", "44,36,5;1,9,40", "27,16,9;1,3,12", "48,35,9;1,7,40", "42,32,9;1,6,24",
          "78,50,9;1,15,60", "39,32,16;1,4,24", "32,28,9;1,2,28", "52,42,16;1,6,28", "56,42,20;1,6,28",
          "80,63,11;1,9,70", "49,40,22;1,5,28", "63,50,8;1,5,56", "21,20,10;1,1,12", "40,33,19;1,3,22",
          "54,42,11;1,3,44", "63,48,24;1,6,28", "99,80,12;1,10,88"]


def scheme(s):
    b, c = s.split(';'); b = [int(x) for x in b.split(',')]; c = [int(x) for x in c.split(',')]
    d = len(b); k = b[0]; B = b + [0]; C = [0] + c; a = [k - B[i] - C[i] for i in range(d + 1)]
    L = sp.zeros(d + 1)
    for i in range(d + 1):
        L[i, i] = a[i]
        if i < d: L[i + 1, i] = B[i]
        if i > 0: L[i - 1, i] = C[i]
    ks = [1]
    for i in range(d): ks.append(ks[-1] * B[i] // C[i + 1])
    v = sum(ks)
    th = sorted(L.eigenvals().keys(), key=lambda z: -float(z))
    U = []
    for t in th:
        u = [sp.Integer(1), t / k]
        for i in range(1, d): u.append(sp.simplify(((t - a[i]) * u[i] - C[i] * u[i - 1]) / B[i]))
        U.append(u)
    m = [sp.nsimplify(sp.simplify(v / sum(ks[i] * U[j][i] ** 2 for i in range(d + 1)))) for j in range(d + 1)]
    bad = 0
    for i in range(d + 1):
        for j in range(d + 1):
            for h in range(d + 1):
                # q^h_ij = (m_i m_j / v) sum_l k_l u_l(theta_i) u_l(theta_j) u_l(theta_h)
                q = sp.nsimplify(sp.simplify(m[i] * m[j] / v *
                                             sum(ks[l] * U[i][l] * U[j][l] * U[h][l] for l in range(d + 1))))
                if not q.is_integer: bad += 1
    return v, th, m, bad


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, 'krein_filter.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['array', 'v', 'eigenvalues', 'multiplicities', 'nonintegral_krein', 'abelian_cayley_possible'])
        for s in ARRAYS:
            v, th, m, bad = scheme(s)
            w.writerow(['{' + s + '}', v, ' '.join(str(t) for t in th), ' '.join(str(x) for x in m), bad,
                        'yes' if bad == 0 else 'no'])
    print('wrote krein_filter.csv')


if __name__ == '__main__':
    main()
