"""Character-trace test for Cayley realisations of a distance-regular graph with arbitrary inverse-closed
connection set. Independent reimplementation. Input: lines 'G n i [ratclass sizes] [galois-orbit sizes] [[T_r(K)]]'.
Variables: w[K][i] = |K cap S_i| (i=1..D), y[r][j] = multiplicity of nonprincipal eigenvalue theta_j in one character of
Galois orbit r. Equations: sum_K T_r(K) w[K][i] = n_r sum_j P[j][i] y[r][j];  sum_r n_r d_r y[r][j] = m_j;
sum_i w[K][i] = |K|; sum_K w[K][i] = k_i;  sum_j y[r][j] = d_r."""
import sys, json, ast, time
from fractions import Fraction
import numpy as np, sympy as sp
def parse_tables(path):
    for raw in open(path):
        if not raw.startswith('G '): continue
        _, n, i, rest = raw.split(' ', 3)
        items, lev, st = [], 0, 0
        for pos, ch in enumerate(rest):
            if ch == '[':
                if lev == 0: st = pos
                lev += 1
            elif ch == ']':
                lev -= 1
                if lev == 0: items.append(ast.literal_eval(rest[st:pos+1]))
        yield int(n), int(i), items[0], items[1], items[2]
def eigenmatrix(b, c):
    """P[j][i] = p_i(theta_j) for intersection array {b0..b_{D-1}; c1..cD}, theta_0 = k first, exact integers required."""
    D = len(b); k = b[0]
    a = [0] + [k - b[i] - c[i-1] for i in range(1, D)] + [k - c[D-1]]
    L = sp.zeros(D+1, D+1)
    for i in range(D+1):
        L[i, i] = a[i]
        if i > 0: L[i, i-1] = b[i-1]
        if i < D: L[i, i+1] = c[i]
    x = sp.symbols('x'); ev = sorted(sp.Poly(L.charpoly(x), x).all_roots(), reverse=True)
    assert all(e.is_integer for e in ev), ev
    P = []
    for th in ev:
        p = [sp.Integer(1), th]
        for i in range(1, D):
            p.append((th*p[i] - a[i]*p[i] - b[i-1]*p[i-1]) / c[i])
        P.append([int(v) for v in p])
    kk = [row for row in P][0]
    mult = []
    # multiplicities m_j = v / sum_i p_i(theta_j)^2 / k_i
    v = sum(kk)
    for row in P:
        s = sum(Fraction(row[i]**2, kk[i]) for i in range(D+1)); mult.append(Fraction(v)/s)
    assert all(m.denominator == 1 for m in mult)
    return P, [int(m) for m in mult]
def model_status_cpsat(sizes, orb, T, P, mult, tlimit=60, workers=1):
    from ortools.sat.python import cp_model
    D = len(P) - 1; nK = len(sizes); v = sum(sizes)
    deg = [T[r][0] // orb[r] for r in range(len(T))]
    assert sum(orb[r]*deg[r]**2 for r in range(len(T))) == v
    for r in range(len(T)):          # rational orthogonality
        for s in range(len(T)):
            val = sum(sizes[K]*T[r][K]*T[s][K] for K in range(nK))
            assert val == (v*orb[r] if r == s else 0)
    m = cp_model.CpModel()
    w = [[m.new_int_var(0, sizes[K], f'w{K}_{i}') for i in range(1, D+1)] for K in range(nK)]
    triv = [r for r in range(len(T)) if deg[r] == 1 and orb[r] == 1 and all(t == 1 for t in T[r])]
    assert len(triv) == 1
    y = {r: [m.new_int_var(0, deg[r], f'y{r}_{j}') for j in range(1, D+1)] for r in range(len(T)) if r not in triv}
    m.add(sum(w[0]) == 0)      # identity class (first) lies in S_0
    assert sizes[0] == 1
    for K in range(1, nK): m.add(sum(w[K]) == sizes[K])
    for i in range(D): m.add(sum(w[K][i] for K in range(nK)) == P[0][i+1])
    for r, yr in y.items(): m.add(sum(yr) == deg[r])
    for j in range(D): m.add(sum(orb[r]*deg[r]*y[r][j] for r in y) == mult[j+1])
    for r, yr in y.items():
        for i in range(D):
            m.add(sum(T[r][K]*w[K][i] for K in range(1, nK)) == orb[r]*sum(P[j+1][i+1]*yr[j] for j in range(D)))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tlimit; s.parameters.num_workers = workers
    st = s.solve(m)
    return s.status_name(st)
def model_status(sizes, orb, T, P, mult, tlimit=60, backend='highs', relax=False):
    if backend == 'cpsat': return model_status_cpsat(sizes, orb, T, P, mult, tlimit)
    from scipy.optimize import milp, LinearConstraint, Bounds
    D = len(P) - 1; nK = len(sizes); v = sum(sizes)
    deg = [T[r][0] // orb[r] for r in range(len(T))]
    assert sum(orb[r]*deg[r]**2 for r in range(len(T))) == v and sizes[0] == 1
    for r in range(len(T)):
        for s in range(len(T)):
            assert sum(sizes[K]*T[r][K]*T[s][K] for K in range(nK)) == (v*orb[r] if r == s else 0)
    triv = [r for r in range(len(T)) if all(t == 1 for t in T[r])]; assert len(triv) == 1
    R = [r for r in range(len(T)) if r not in triv]
    idx = {}
    for K in range(1, nK):
        for i in range(D): idx[('w', K, i)] = len(idx)
    for r in R:
        for j in range(D): idx[('y', r, j)] = len(idx)
    nv = len(idx); rows = []; lo = []; hi = []
    def eq(coef, rhs):
        a = np.zeros(nv)
        for k_, c_ in coef: a[idx[k_]] += c_
        rows.append(a); lo.append(rhs); hi.append(rhs)
    for K in range(1, nK): eq([(('w', K, i), 1) for i in range(D)], sizes[K])
    for i in range(D): eq([(('w', K, i), 1) for K in range(1, nK)], P[0][i+1])
    for r in R: eq([(('y', r, j), 1) for j in range(D)], deg[r])
    for j in range(D): eq([(('y', r, j), orb[r]*deg[r]) for r in R], mult[j+1])
    for r in R:
        for i in range(D):
            eq([(('w', K, i), T[r][K]) for K in range(1, nK)] + [(('y', r, j), -orb[r]*P[j+1][i+1]) for j in range(D)], 0)
    ub = np.zeros(nv)
    for (t_, a_, b_), q in idx.items(): ub[q] = sizes[a_] if t_ == 'w' else deg[a_]
    res = milp(np.zeros(nv), constraints=LinearConstraint(np.array(rows), lo, hi), integrality=(np.zeros(nv) if relax else np.ones(nv)),
               bounds=Bounds(np.zeros(nv), ub), options=dict(time_limit=tlimit))
    return {0: 'FEASIBLE', 2: 'INFEASIBLE', 1: 'UNKNOWN'}.get(res.status, 'STATUS%d' % res.status)
if __name__ == '__main__':
    tables, b, c, out = sys.argv[1], json.loads(sys.argv[2]), json.loads(sys.argv[3]), sys.argv[4]
    tl = float(sys.argv[5]) if len(sys.argv) > 5 else 60
    P, mult = eigenmatrix(b, c)
    v = sum(P[0])
    res = []
    for n, gid, sizes, orb, T in parse_tables(tables):
        if n != v: continue
        t = time.time(); st = model_status(sizes, orb, T, P, mult, tl)
        res.append(dict(n=n, id=gid, status=st, sec=round(time.time()-t, 2), classes=len(sizes)))
        with open(out, 'w') as f: json.dump(dict(array=[b, c], P=P, mult=mult, results=res), f)
        print(n, gid, st, round(time.time()-t, 1), flush=True)
