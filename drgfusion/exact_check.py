"""Exact (sympy) Krein parameters and absolute bound for given intersection arrays.

Usage: python exact_check.py "14,10,3;1,5,12" "21,16,10,1;1,2,16,21"
q^h_ij = (m_i m_j / v) sum_l k_l u_i(l) u_j(l) u_h(l), u = cosine sequences.
"""
import sys
import sympy as sp


def scheme(b, c):
    d = len(b); k = b[0]
    bb = list(b) + [0]; cc = [0] + list(c)
    a = [k - bb[i] - cc[i] for i in range(d + 1)]
    L = sp.zeros(d + 1)
    for i in range(d + 1):
        L[i, i] = a[i]
        if i < d: L[i, i + 1] = bb[i]
        if i > 0: L[i, i - 1] = cc[i]
    kk = [sp.Integer(1)]
    for i in range(d):
        kk.append(kk[-1] * bb[i] / cc[i + 1])
    v = sum(kk)
    th = sorted(sp.roots(L.charpoly().as_expr()).keys(), key=lambda x: -float(x))
    U = []
    for t in th:
        u = [sp.Integer(1), t / k]
        for i in range(1, d):
            u.append(sp.simplify(((t - a[i]) * u[i] - cc[i] * u[i - 1]) / bb[i]))
        U.append(u)
    m = [sp.simplify(v / sum(kk[l] * U[j][l] ** 2 for l in range(d + 1))) for j in range(d + 1)]
    return v, kk, th, U, m


def main(arrays):
    for s in arrays:
        b, c = [tuple(int(x) for x in p.split(',')) for p in s.split(';')]
        v, kk, th, U, m = scheme(b, c)
        d = len(b)
        print('array {%s}  v=%s  eigenvalues %s  multiplicities %s' % (s, v, th, m))
        for i in range(1, d + 1):
            for j in range(i, d + 1):
                q = [sp.nsimplify(sp.simplify(m[i] * m[j] / v * sum(kk[l] * U[i][l] * U[j][l] * U[h][l]
                                                                    for l in range(d + 1)))) for h in range(d + 1)]
                supp = [h for h in range(d + 1) if q[h] != 0]
                tot = sum(m[h] for h in supp)
                bound = m[i] * (m[i] + 1) / 2 if i == j else m[i] * m[j]
                flag = []
                if any(x < 0 for x in q): flag.append('KREIN FAILS')
                if tot > bound: flag.append('ABSOLUTE BOUND FAILS (%s > %s)' % (tot, bound))
                if flag:
                    print('  q^h_{%d%d} = %s  support sum %s  %s' % (i, j, q, tot, '; '.join(flag)))


if __name__ == '__main__':
    main(sys.argv[1:])
