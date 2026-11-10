"""Triple intersection numbers of the GO(2,4) array {10,8,8,8;1,1,1,5} with the vanishing Krein parameter q^4_44.

For vertices u, v, w with d(u,v), d(u,w), d(v,w) = (U, V, W), [i j k] counts x with d(u,x) = i, d(v,x) = j,
d(w,x) = k. The marginals are intersection numbers, and q^4_44 = 0 gives (Coolsaet-Jurisic)
sum_{i,j,k} Q_{i4} Q_{j4} Q_{k4} [i j k] = 0. For each of the 17 triple types the script computes, by LP,
the range of every [i j k] with and without the Krein equation, lists values that the Krein equation forces,
and checks them on the actual octagon (ATLAS permutation representation TF42G1-p1755B0).
Usage: python go24_triple_lp.py   (writes go24_triple_lp_output.txt)
"""
import os, sys, re, itertools
import numpy as np, sympy as sp, requests
from scipy.optimize import linprog
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'drgfusion'))
sys.path.insert(0, os.path.join(HERE, '..', 'tits1755'))
import drg_table_fusions as dtf

b, c = [10, 8, 8, 8], [1, 1, 1, 5]
p, kv = dtf.intersection_numbers(b, c)
n = sum(kv)
th = [10, 5, 1, -3, -5]
bb = b + [0]; cc = [0] + c
a = [b[0] - bb[i] - cc[i] for i in range(5)]
U = []
for t_ in th:
    u = [sp.Integer(1), sp.Rational(t_, 10)]
    for i in range(1, 4):
        u.append(((t_ - a[i]) * u[i] - cc[i] * u[i - 1]) / bb[i])
    U.append(u)
m = [sp.nsimplify(n / sum(kv[i] * U[j][i] ** 2 for i in range(5))) for j in range(5)]
Q = sp.Matrix(5, 5, lambda i, j: m[j] * U[j][i])
krein = lambda i, j, k: sp.nsimplify(sp.Rational(m[i] * m[j], n) * sum(kv[l] * U[i][l] * U[j][l] * U[k][l] for l in range(5)))
P = lambda k, i, j: int(p[i][j][k])


def system(Uu, Vv, Ww):
    idx = {}
    for i, j, k in itertools.product(range(5), repeat=3):
        if P(Uu, i, j) > 0 and P(Vv, i, k) > 0 and P(Ww, j, k) > 0:
            idx[(i, j, k)] = len(idx)
    rows, rhs = [], []
    for (A, R) in ((lambda i, j, k: (i, j), Uu), (lambda i, j, k: (i, k), Vv), (lambda i, j, k: (j, k), Ww)):
        for x, y in itertools.product(range(5), repeat=2):
            r = [0] * len(idx)
            for key in idx:
                if A(*key) == (x, y):
                    r[idx[key]] = 1
            if any(r) or P(R, x, y):
                rows.append(r); rhs.append(P(R, x, y))
    rows.append([Q[i, 4] * Q[j, 4] * Q[k, 4] for (i, j, k) in idx]); rhs.append(0)
    return list(idx), np.array(sp.Matrix(rows).tolist(), dtype=float), np.array(rhs, dtype=float)


def ranges(A, bv, nv):
    lo, hi = [], []
    for v in range(nv):
        cvec = np.zeros(nv); cvec[v] = 1
        r1 = linprog(cvec, A_eq=A, b_eq=bv, bounds=(0, None), method='highs')
        r2 = linprog(-cvec, A_eq=A, b_eq=bv, bounds=(0, None), method='highs')
        lo.append(r1.fun); hi.append(-r2.fun)
    return lo, hi


def octagon_distances():
    gens = []
    for g in ('g1', 'g2'):
        t = requests.get(f'https://brauer.maths.qmul.ac.uk/Atlas/exc/TF42/gap/TF42G1-p1755B0.{g}', timeout=60).text
        perm = list(range(1755))
        for cyc in re.findall(r'\(([\d,\s]+)\)', t[t.find('('):]):
            pts = [int(x) - 1 for x in cyc.replace('\n', '').split(',')]
            for x, y in zip(pts, pts[1:] + pts[:1]):
                perm[x] = y
        gens.append(np.array(perm))
    from coset_schemes import orbitals
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path
    L, val, reps, _ = orbitals(gens, 1755)
    k1 = list(val).index(10)
    return shortest_path(csr_matrix((L == k1).astype(np.int8)), unweighted=True).astype(np.int8)


def main():
    out = [f'multiplicities {m}; vanishing Krein parameters (i<=j): ' +
           str([(i, j, k) for i in range(1, 5) for j in range(i, 5) for k in range(1, 5) if krein(i, j, k) == 0])]
    D = octagon_distances()
    types = [(x, y, z) for x in range(1, 5) for y in range(x, 5) for z in range(y, 5) if P(x, y, z) > 0]
    for ty in types:
        keys, A, bv = system(*ty)
        lo, hi = ranges(A, bv, len(keys))
        lo0, hi0 = ranges(A[:-1], bv[:-1], len(keys))
        forced = [(keys[i], round(lo[i], 6)) for i in range(len(keys)) if hi[i] - lo[i] < 1e-7 < hi0[i] - lo0[i]]
        ntight = sum((hi0[i] - lo0[i]) - (hi[i] - lo[i]) > 1e-6 for i in range(len(keys)))
        # one triple of this type in the octagon
        Uu, Vv, Ww = ty
        u = 0
        vs = np.nonzero(D[u] == Uu)[0]
        trip = None
        for v in vs[:50]:
            ws = np.nonzero((D[u] == Vv) & (D[v] == Ww))[0]
            if len(ws):
                trip = (u, v, ws[0]); break
        cnt = np.zeros((5, 5, 5), int)
        for x in range(1755):
            cnt[D[trip[0], x], D[trip[1], x], D[trip[2], x]] += 1
        ok = all(abs(cnt[k] - val) < 1e-6 for k, val in forced)
        inrange = all(lo[i] - 1e-6 <= cnt[keys[i]] <= hi[i] + 1e-6 for i in range(len(keys)))
        out.append(f'type {ty}: {len(keys)} unknowns, free dimension {len(keys) - np.linalg.matrix_rank(A)}; '
                   f'Krein tightens {ntight} ranges; forced only with Krein: {forced}; '
                   f'octagon triple agrees with forced values: {ok}; all values within LP ranges: {inrange}')
    txt = '\n'.join(out)
    open(os.path.join(HERE, 'go24_triple_lp_output.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
