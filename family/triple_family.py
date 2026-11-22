"""Triple intersection numbers of generalized octagons of order (s, s^2): what the LP, the Krein equation and the
general geometric constraints give at (s, t) = (2, 4) and (8, 64).

For each triple type (U, V, W) (U <= V <= W, p^W_UV > 0) of the point graph, intersection array
{s(t+1), st, st, st; 1, 1, 1, t+1}:
  * LP: marginals equal intersection numbers, plus the Krein equation q^4_44 = 0 (E_4 = eigenvalue -(t+1));
  * G: the constraints proved for every GO(s, t) in triple/triple_relations_note.md: zeros from the geodesic-tree
    test (triple_zero_proofs.feasible, run with the given s, t), [313] = t in types (2,3,4) and (2,4,4), and
    [133] = t in type (2,4,4);
  * free dimension of LP and of LP + G (exact rank over Q);
  * implicit zeros: unknowns whose maximum under the LP with nonnegativity is 0 (scaled HiGHS LP, threshold 1e-9,
    with the smallest positive scaled maximum printed as a margin); adding them, the exact free dimension and the
    number of unknowns determined (e_j in the row space, exact).
At (2, 4) the results are compared with the Ree-Tits octagon (triple/out_orbits/triple_orbit_vectors.csv): for the
types where the octagon has one vector, does LP + G already determine it?
Usage: python family/triple_family.py   (writes family/out/triple_family_output.txt, triple_family.csv)
"""
import os, sys, itertools, csv, re
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'triple'))
import numpy as np
from scipy.optimize import linprog
import triple_zero_proofs as tz

OUT = os.path.join(HERE, 'out'); os.makedirs(OUT, exist_ok=True)
log = []
def say(*x):
    s_ = ' '.join(str(y) for y in x); log.append(s_); print(s_, flush=True)

def isqrt_exact(n):
    r = int(round(n ** 0.5))
    assert r * r == n, n
    return r

def scheme(s, t):
    b = [s * (t + 1), s * t, s * t, s * t]; c = [1, 1, 1, t + 1]
    d = 4; bb = b + [0]; cc = [0] + c; k = b[0]
    a = [k - bb[i] - cc[i] for i in range(5)]
    kv = [1]
    for i in range(4): kv.append(kv[-1] * b[i] // c[i])
    # p^h_ij from the tridiagonal recurrence: A_1 A_j = b_{j-1} A_{j-1} + a_j A_j + c_{j+1} A_{j+1}
    n = 5
    L = [[F(0)] * n for _ in range(n)]
    for j in range(n):
        if j > 0: L[j - 1][j] = F(bb[j - 1])
        L[j][j] = F(a[j])
        if j < d: L[j + 1][j] = F(cc[j + 1])
    mm = lambda X, Y: [[sum(X[r_][q] * Y[q][s_] for q in range(n)) for s_ in range(n)] for r_ in range(n)]
    # B_i = matrix of multiplication by A_i in basis A_0..A_4: B_0 = I, B_1 = L, c_{i+1} B_{i+1} = (L - a_i) B_i - b_{i-1} B_{i-1}
    I = [[F(int(r_ == s_)) for s_ in range(n)] for r_ in range(n)]
    B = [I, L]
    for i in range(1, d):
        X = mm(L, B[i])
        B.append([[(X[r_][s_] - a[i] * B[i][r_][s_] - bb[i - 1] * B[i - 1][r_][s_]) / cc[i + 1] for s_ in range(n)] for r_ in range(n)])
    p = [[[B[i][h][j] for j in range(n)] for i in range(n)] for h in range(n)]   # p[h][i][j] = p^h_ij
    for h in range(n):
        for i in range(n):
            for j in range(n):
                assert p[h][i][j].denominator == 1 and p[h][i][j] >= 0
    rt = isqrt_exact(2 * s * t)
    th = [k, s - 1 + rt, s - 1, s - 1 - rt, -t - 1]
    U = []
    for x in th:
        u = [F(1), F(x, k)]
        for i in range(1, 4):
            u.append(((x - a[i]) * u[i] - cc[i] * u[i - 1]) / bb[i])
        U.append(u)
    v = sum(kv)
    m = [F(v) / sum(kv[i] * U[j][i] ** 2 for i in range(5)) for j in range(5)]
    Q = [[m[j] * U[j][i] for j in range(5)] for i in range(5)]
    q444 = m[4] * m[4] / v * sum(kv[l] * U[4][l] ** 3 for l in range(5))
    return p, Q, q444, v

def system(p, Q, Uu, Vv, Ww):
    P = lambda h, i, j: p[h][i][j]
    idx = {}
    for i, j, k in itertools.product(range(5), repeat=3):
        if P(Uu, i, j) > 0 and P(Vv, i, k) > 0 and P(Ww, j, k) > 0:
            idx[(i, j, k)] = len(idx)
    rows, rhs = [], []
    for (Af, R) in ((lambda i, j, k: (i, j), Uu), (lambda i, j, k: (i, k), Vv), (lambda i, j, k: (j, k), Ww)):
        for x, y in itertools.product(range(5), repeat=2):
            r_ = [F(0)] * len(idx)
            for key in idx:
                if Af(*key) == (x, y):
                    r_[idx[key]] = F(1)
            if any(r_) or P(R, x, y):
                rows.append(r_); rhs.append(F(P(R, x, y)))
    rows.append([Q[i][4] * Q[j][4] * Q[k][4] for (i, j, k) in idx]); rhs.append(F(0))
    return list(idx), rows, rhs

def rank(rows):
    M = [list(r_) for r_ in rows]
    rk, ncol = 0, len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        pv = M[rk][c]
        M[rk] = [x / pv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk

def implicit_zeros(keys, rows, rhs, P_ub):
    """unknowns forced to 0 by the LP with nonnegativity: max x_j = 0. Variables are scaled by an upper bound
    (a marginal), rows by their largest entry, so the LP is well conditioned; a coordinate counts as 0 when its
    scaled maximum is below 1e-9 (the others are far above this)."""
    ub = np.array([float(u) for u in P_ub])
    A = np.array([[float(x) for x in r_] for r_ in rows]) * ub
    bv = np.array([float(x) for x in rhs])
    nr = np.abs(A).max(axis=1); nr[nr == 0] = 1
    A = A / nr[:, None]; bv = bv / nr
    zeros, worst = [], 0.0
    for j in range(len(keys)):
        cvec = np.zeros(len(keys)); cvec[j] = -1
        res = linprog(cvec, A_eq=A, b_eq=bv, bounds=(0, 1), method='highs')
        if res.status != 0:
            return None, None
        if -res.fun < 1e-9:
            zeros.append(keys[j])
        else:
            worst = max(worst, 0) if worst else -res.fun
            worst = min(worst, -res.fun)
    return zeros, worst

def determined(keys, rows):
    rk = rank(rows)
    return [kk for c, kk in enumerate(keys) if rank(rows + [[F(int(c == c2)) for c2 in range(len(keys))]]) == rk]

HAND = {(2, 3, 4): ['313'], (2, 4, 4): ['133', '313']}
obs = {}
with open(os.path.join(HERE, '..', 'triple', 'out_orbits', 'triple_orbit_vectors.csv')) as fh:
    for row in csv.DictReader(fh):
        ty = tuple(int(x) for x in re.findall(r'\d', row['type']))
        vv = {m_.group(1): int(m_.group(2)) for m_ in re.finditer(r'\[(\d\d\d)\]=(\d+)', row['vector'])}
        obs.setdefault(ty, set()).add(tuple(sorted(vv.items())))

table = []
for (s, t) in [(2, 4), (8, 64)]:
    p, Q, q444, v = scheme(s, t)
    say(f'(s, t) = ({s}, {t}): v = {v}, q^4_44 = {q444}')
    assert q444 == 0
    tz.s_, tz.t_ = s, t
    types = [(x, y, z) for x in range(1, 5) for y in range(x, 5) for z in range(y, 5) if p[z][x][y] > 0]
    for ty in types:
        keys, rows, rhs = system(p, Q, *ty)
        n = len(keys)
        rk = rank(rows)
        G = {}
        for kk in keys:
            i, j, k = kk
            if 1 <= min(kk) and max(kk) <= 3 and not tz.feasible(*ty, i, j, k):
                G[kk] = F(0)
        for h in HAND.get(ty, []):
            G[tuple(map(int, h))] = F(t)
        rowsG = rows + [[F(int(keys[c] == kk)) for c in range(n)] for kk in G]
        rhsG = rhs + list(G.values())
        rkG = rank(rowsG)
        P = lambda h, i, j: p[h][i][j]
        ub = [max(1, min(P(ty[0], i, j), P(ty[1], i, k), P(ty[2], j, k))) for (i, j, k) in keys]
        out = {}
        for tag, R_, b_ in (('LP', rows, rhs), ('LP_G', rowsG, rhsG)):
            z, smallest = implicit_zeros(keys, R_, b_, ub)
            if z is None:
                out[tag] = ('infeasible', None, None); continue
            Rz = R_ + [[F(int(keys[c] == kk)) for c in range(n)] for kk in z]
            det = determined(keys, Rz)
            out[tag] = (len(z), n - rank(Rz), len(det), smallest)
        line = {'s': s, 't': t, 'type': ty, 'unknowns': n, 'free_LP': n - rk, 'free_LP_G': n - rkG, 'n_G': len(G),
                'implicit_zeros_LP': out['LP'][0], 'free_after_zeros_LP': out['LP'][1], 'determined_LP': out['LP'][2],
                'implicit_zeros_LP_G': out['LP_G'][0], 'free_after_zeros_LP_G': out['LP_G'][1],
                'determined_LP_G': out['LP_G'][2], 'min_positive_scaled_max': out['LP_G'][3]}
        if s == 2 and ty in obs:
            line['RT_vectors'] = len(obs[ty])
        table.append(line)
        say('  ' + ', '.join(f'{k_}={v_}' for k_, v_ in line.items() if k_ not in ('s', 't')))
with open(os.path.join(OUT, 'triple_family.csv'), 'w', newline='\n') as fh:
    keys_ = ['s', 't', 'type', 'unknowns', 'free_LP', 'free_LP_G', 'n_G', 'implicit_zeros_LP', 'free_after_zeros_LP',
             'determined_LP', 'implicit_zeros_LP_G', 'free_after_zeros_LP_G', 'determined_LP_G', 'RT_vectors']
    fh.write(','.join(keys_) + '\n')
    for line in table:
        fh.write(','.join(f'"{line.get(k_, "")}"' for k_ in keys_) + '\n')
open(os.path.join(OUT, 'triple_family_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
