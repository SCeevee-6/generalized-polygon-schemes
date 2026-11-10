"""Exhaustive search for strongly regular and distance-regular graphs that are unions of orbitals.

Input: intersection numbers p^k_ij of a homogeneous coherent configuration (commutative or not), as a CSV
with columns k,i,j,p_k_ij (zero rows may be omitted); class 0 is the diagonal.

Any undirected graph on which a transitive group G acts is a union S of G-orbitals with S = S^T. The
search runs through all such S with a Gray code over the transposition-closed units, keeps the S for
which sum_{i,j in S} p^k_ij is constant on S (constant a_1), and for those builds the distance partition
D_0 = {0}, D_1 = S, D_2, ... from the structure constants. For (u,v) in class k the coefficient of A_k in
A_{D_i} A_S is |D_i(u) cap Gamma(v)|; the graph is distance-regular iff this is constant on D_{i-1}
(b_{i-1}), on D_i (a_i) and on D_{i+1} (c_{i+1}) for every i. Integer arithmetic throughout.

Usage: python srg_drg_search.py <intersection_numbers.csv> <out.json> [gray|dfs|cpsat|cuts] [time_limit_s] [nomp]

With 'nomp' the cuts model also excludes the complete multipartite graphs (complements of closed subsets).

The dfs method (search_dfs) prunes with the monotonicity of the walk counts and is kept for comparison; cpsat (search_cpsat) adds the condition that walk counts off S lie in {0, c_2} and is used for the large ranks.
"""
import sys, json, csv, time
import numpy as np


def load_p(path):
    rows = []
    with open(path) as f:
        rd = csv.reader(f)
        next(rd)
        for r in rd:
            if len(r) >= 4:
                rows.append(tuple(int(x) for x in r[:4]))
    r = 1 + max(max(a[:3]) for a in rows)
    p = np.zeros((r, r, r), dtype=np.int64)  # p[i, j, k] = p^k_ij
    for k, i, j, v in rows:
        p[i, j, k] = v
    return p


def structure(p):
    r = p.shape[0]
    val = p[:, :, 0].sum(axis=1)
    pair = [int(np.nonzero(p[i, :, 0])[0][0]) for i in range(r)]
    units, seen = [], set()
    for i in range(1, r):
        if i not in seen:
            u = sorted({i, pair[i]})
            seen.update(u)
            units.append(u)
    return val, pair, units


def walk(p, X, S):
    r = p.shape[0]
    x = np.zeros(r, dtype=np.int64); x[list(X)] = 1
    y = np.zeros(r, dtype=np.int64); y[list(S)] = 1
    return np.einsum('a,abk,b->k', x, p, y)


def distance_regular(p, S):
    """Intersection array if the union S is a connected DRG; 'disconnected'; or None."""
    r = p.shape[0]
    layers = [[0], sorted(S)]
    seen = {0} | set(S)
    b, c = [], [1]
    b.append(int(walk(p, S, S)[0]))  # b_0 = k
    while True:
        Di = layers[-1]
        co = walk(p, Di, S)
        nxt = sorted(k for k in range(r) if co[k] > 0 and k not in seen)
        for L in (layers[-2], Di, nxt):
            if L and len({int(co[k]) for k in L}) != 1:
                return None
        if len(layers) > 2:
            b.append(int(co[layers[-2][0]]))  # b_{i-1}
        if not nxt:
            break
        c.append(int(co[nxt[0]]))  # c_{i+1}
        layers.append(nxt)
        seen.update(nxt)
    if len(seen) < r:
        return 'disconnected'
    return {'diameter': len(layers) - 1, 'b': b, 'c': c, 'layers': layers}


def search(p):
    r = p.shape[0]
    val, pair, units = structure(p)
    U = len(units)
    Q = np.zeros((U, U, r), dtype=np.int64)
    for a, ua in enumerate(units):
        for bi, ub in enumerate(units):
            Q[a, bi] = p[np.ix_(ua, ub)].sum(axis=(0, 1))
    Rsym = Q + Q.transpose(1, 0, 2)
    Qd = np.array([Q[a, a] for a in range(U)])
    cls_mask = np.zeros((U, r), dtype=bool)
    for a, u in enumerate(units):
        cls_mask[a, u] = True
    inS = np.zeros(U, dtype=bool)
    row = np.zeros((U, r), dtype=np.int64)
    cvec = np.zeros(r, dtype=np.int64)
    cand, total = [], 0
    for g in range(1, 1 << U):
        x = (g & -g).bit_length() - 1
        if not inS[x]:
            cvec += row[x] + Qd[x]; inS[x] = True; row += Rsym[:, x]
        else:
            inS[x] = False; row -= Rsym[:, x]; cvec -= row[x] + Qd[x]
        if inS.all():
            continue
        total += 1
        v = cvec[cls_mask[inS].any(axis=0)]
        if v.min() == v.max():
            cand.append(np.nonzero(inS)[0].tolist())
    return val, pair, units, cand, total


def search_dfs(p, order=None):
    """All unions S (transposition-closed, nonempty, proper) with sum_{i,j in S} p^k_ij constant on S.

    Depth-first over units. With S' the units chosen so far and W = S' plus the undecided units,
    every final S satisfies S' <= S <= W, and c_k(X) = sum_{i,j in X} p^k_ij is monotone in X.
    So a branch is cut when max_{k in S'} c_k(S') > min_{k in S'} c_k(W)."""
    r = p.shape[0]
    val, pair, units = structure(p)
    U = len(units)
    Q = np.zeros((U, U, r), dtype=np.int64)
    for a, ua in enumerate(units):
        for bi, ub in enumerate(units):
            Q[a, bi] = p[np.ix_(ua, ub)].sum(axis=(0, 1))
    Rsym = Q + Q.transpose(1, 0, 2)
    Qd = np.array([Q[a, a] for a in range(U)])
    cls = [np.array(u) for u in units]
    if order is None:
        order = list(range(U))
    cand = []
    stats = {'nodes': 0}
    full = np.ones(U, dtype=bool)
    cW0 = Q.sum(axis=(0, 1))
    rowW0 = Rsym.sum(axis=1)

    def rec(t, inS, cS, rowS, cW, rowW, smask):
        stats['nodes'] += 1
        if smask.any():
            if cS[smask].max() > cW[smask].min():
                return
        if t == U:
            if smask.any() and not inS.all():
                if cS[smask].min() == cS[smask].max():
                    cand.append(np.nonzero(inS)[0].tolist())
            return
        x = order[t]
        # include x
        inS2 = inS.copy(); inS2[x] = True
        sm2 = smask.copy(); sm2[cls[x]] = True
        rec(t + 1, inS2, cS + rowS[x] + Qd[x], rowS + Rsym[:, x], cW, rowW, sm2)
        # exclude x
        rec(t + 1, inS, cS, rowS, cW - (rowW[x] - Qd[x]), rowW - Rsym[:, x], smask)

    import sys as _s
    _s.setrecursionlimit(10000)
    rec(0, np.zeros(U, dtype=bool), np.zeros(r, dtype=np.int64), np.zeros((U, r), dtype=np.int64),
        cW0.copy(), rowW0.copy(), np.zeros(r, dtype=bool))
    return val, pair, units, cand, stats['nodes']


def search_cpsat(p, workers=8, time_limit=None):
    """All transposition-closed S with constant a_1 on S and, off S, walk counts in {0, c_2} with at
    least one class of count c_2. Every connected distance-regular union S (diameter >= 2) satisfies
    this. Enumerated exhaustively with OR-Tools CP-SAT (y_ab = x_a and x_b linearises the quadratic walk
    counts c_k(S) = sum_{i,j in S} p^k_ij)."""
    from ortools.sat.python import cp_model
    r = p.shape[0]
    val, pair, units = structure(p)
    U = len(units)
    Q = np.zeros((U, U, r), dtype=np.int64)
    for a, ua in enumerate(units):
        for bi, ub in enumerate(units):
            Q[a, bi] = p[np.ix_(ua, ub)].sum(axis=(0, 1))
    n = int(val.sum())
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f'x{a}') for a in range(U)]
    z = [m.NewBoolVar(f'z{a}') for a in range(U)]
    y = {}
    for a in range(U):
        for b in range(a + 1, U):
            v = m.NewBoolVar(f'y{a}_{b}')
            m.AddBoolAnd([x[a], x[b]]).OnlyEnforceIf(v)
            m.AddBoolOr([x[a].Not(), x[b].Not()]).OnlyEnforceIf(v.Not())
            y[a, b] = v
    a1 = m.NewIntVar(0, n, 'a1')
    c2 = m.NewIntVar(1, n, 'c2')
    for u in range(U):
        k = units[u][0]
        expr = sum(int(Q[a, a, k]) * x[a] for a in range(U) if Q[a, a, k]) + \
            sum(int(Q[a, b, k] + Q[b, a, k]) * y[a, b] for (a, b) in y if Q[a, b, k] + Q[b, a, k])
        m.Add(expr == a1).OnlyEnforceIf(x[u])
        m.Add(expr == c2).OnlyEnforceIf(z[u])
        m.Add(expr == 0).OnlyEnforceIf([x[u].Not(), z[u].Not()])
        m.AddImplication(z[u], x[u].Not())
    m.Add(sum(x) >= 1)
    m.Add(sum(z) >= 1)
    sols = []

    class CB(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
        def on_solution_callback(self):
            sols.append([a for a in range(U) if self.Value(x[a])])

    s = cp_model.CpSolver()
    s.parameters.enumerate_all_solutions = True
    s.parameters.num_workers = 1
    if time_limit:
        s.parameters.max_time_in_seconds = time_limit
    st = s.Solve(m, CB())
    status = s.StatusName(st)
    uniq = sorted({tuple(t) for t in sols})
    return val, pair, units, [list(t) for t in uniq], status


def closed_subsets(p):
    """All closed subsets (unions of classes containing 0 and closed under the product) other than the whole
    set, as sorted tuples of class indices. A union S of orbitals gives a connected graph iff S lies in
    no proper closed subset."""
    r = p.shape[0]
    supp = p.sum(axis=(0, 1)) * 0
    prod = [[set(np.nonzero(p[i, j])[0].tolist()) for j in range(r)] for i in range(r)]
    def closure(X):
        X = set(X) | {0}
        while True:
            new = set()
            for i in X:
                for j in X:
                    new |= prod[i][j]
            if new <= X:
                return frozenset(X)
            X |= new
    seen = set()
    frontier = [closure([])]
    full = frozenset(range(r))
    while frontier:
        C = frontier.pop()
        if C in seen or C == full:
            continue
        seen.add(C)
        for k in range(r):
            if k not in C:
                D = closure(C | {k})
                if D not in seen and D != full:
                    frontier.append(D)
    return sorted(tuple(sorted(C)) for C in seen)


def maximal_closed_subsets(p):
    cs = [set(c) for c in closed_subsets(p)]
    return [tuple(sorted(c)) for c in cs if not any(c < d for d in cs)]


def search_cpsat_cuts(p, workers=6, time_limit=3600.0, log=None, exclude_multipartite=False):
    """Same model as search_cpsat plus connectivity (S in no maximal closed subset), solved repeatedly in feasibility mode with all workers; after each
    solution a no-good cut on the x variables is added. Ends with status INFEASIBLE (complete) or UNKNOWN."""
    from ortools.sat.python import cp_model
    r = p.shape[0]
    val, pair, units = structure(p)
    U = len(units)
    Q = np.zeros((U, U, r), dtype=np.int64)
    for a, ua in enumerate(units):
        for bi, ub in enumerate(units):
            Q[a, bi] = p[np.ix_(ua, ub)].sum(axis=(0, 1))
    n = int(val.sum())
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f'x{a}') for a in range(U)]
    z = [m.NewBoolVar(f'z{a}') for a in range(U)]
    y = {}
    for a in range(U):
        for b in range(a + 1, U):
            v = m.NewBoolVar(f'y{a}_{b}')
            m.AddImplication(v, x[a]); m.AddImplication(v, x[b])
            m.AddBoolOr([x[a].Not(), x[b].Not(), v])
            y[a, b] = v
    a1 = m.NewIntVar(0, n, 'a1')
    c2 = m.NewIntVar(1, n, 'c2')
    for u in range(U):
        k = units[u][0]
        expr = sum(int(Q[a, a, k]) * x[a] for a in range(U) if Q[a, a, k]) + \
            sum(int(Q[a, b, k] + Q[b, a, k]) * y[a, b] for (a, b) in y if Q[a, b, k] + Q[b, a, k])
        m.Add(expr == a1).OnlyEnforceIf(x[u])
        m.Add(expr == c2).OnlyEnforceIf(z[u])
        m.Add(expr == 0).OnlyEnforceIf([x[u].Not(), z[u].Not()])
        m.AddImplication(z[u], x[u].Not())
    m.Add(sum(x) >= 1)
    m.Add(sum(z) >= 1)
    # connectivity: S must not lie in a proper closed subset (a block of imprimitivity relation)
    for C in maximal_closed_subsets(p):
        out = [a for a in range(U) if units[a][0] not in C]
        m.AddBoolOr([x[a] for a in out])
    if exclude_multipartite:
        # S = complement of a proper closed subset C != {0} gives a complete multipartite graph; cut these
        for C in closed_subsets(p):
            if len(C) > 1:
                pat = [units[a][0] not in C for a in range(U)]
                m.AddBoolOr([x[a].Not() if pat[a] else x[a] for a in range(U)])
    sols, status, t_used = [], None, 0.0
    while True:
        s = cp_model.CpSolver()
        s.parameters.num_workers = workers
        s.parameters.max_time_in_seconds = max(1.0, time_limit - t_used)
        st = s.Solve(m)
        t_used += s.WallTime()
        status = s.StatusName(st)
        if log:
            log(f'{status} {round(t_used,1)}s sols={len(sols)}')
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            break
        sol = [a for a in range(U) if s.Value(x[a])]
        sols.append(sol)
        if log:
            S = sorted(i for a in sol for i in units[a])
            log(f'solution units={sol} classes={S} drg={distance_regular(p, S)}')
        m.AddBoolOr([x[a].Not() if a in sol else x[a] for a in range(U)])
    return val, pair, units, sols, status


def run(inp, method='gray', **kw):
    t0 = time.time()
    p = load_p(inp)
    val, pair, units, cand, total = {'gray': search, 'dfs': search_dfs, 'cpsat': search_cpsat,
                                     'cuts': search_cpsat_cuts}[method](p, **kw)
    n = int(val.sum())
    drg, ndisc = [], 0
    for us in cand:
        S = sorted(i for a in us for i in units[a])
        res = distance_regular(p, S)
        if res == 'disconnected':
            ndisc += 1
        elif res is not None:
            e = {'S': S, 'k': int(val[S].sum()), **res}
            if res['diameter'] == 2:
                k, b1, c2 = res['b'][0], res['b'][1], res['c'][1]
                e['srg'] = [n, k, k - b1 - 1, c2]
            drg.append(e)
    return {'input': inp, 'n': n, 'rank': int(p.shape[0]), 'units': units, 'valencies': val.tolist(),
            'method': method, 'subsets_or_nodes': total, 'a1_constant': len(cand), 'a1_constant_disconnected': ndisc,
            'connected_drg': drg, 'seconds': round(time.time() - t0, 1)}


def summary(info):
    srg = [e['srg'] for e in info['connected_drg'] if e['diameter'] == 2]
    big = [(e['b'], e['c']) for e in info['connected_drg'] if e['diameter'] >= 3]
    return (f"n={info['n']} rank={info['rank']} units={len(info['units'])} {info['method']}={info['subsets_or_nodes']} "
            f"a1const={info['a1_constant']} disc={info['a1_constant_disconnected']} SRG={srg} DRG(d>=3)={big} "
            f"t={info['seconds']}s")


if __name__ == '__main__':
    kw = {}
    if len(sys.argv) > 4:
        kw = {'time_limit': float(sys.argv[4]), 'log': lambda m: print(m, flush=True),
              'exclude_multipartite': len(sys.argv) > 5 and sys.argv[5] == 'nomp'}
    info = run(sys.argv[1], sys.argv[3] if len(sys.argv) > 3 else 'gray', **kw)
    with open(sys.argv[2], 'w') as f:
        json.dump(info, f, indent=1)
    print(summary(info))
