"""Quotient-and-lift search for distance-regular Cayley graphs on abelian groups.

For N = d_1 Z_{n_1} + ... + d_r Z_{n_r} <= G = Z_{n_1} + ... + Z_{n_r}, put w_i(g+N) = |S_i cap (g+N)|.
Characters of G/N are characters of G trivial on N, so
  (1) sum_{g in G/N} chi(g) w_i(g) = P_{j i}  for chi in E_j cap N^perp,
  (2) |N| sum_{chi in E_j cap N^perp} chi(g) = sum_i Q_{i j} w_i(g) + Q_{0 j}[g = 0],
and w_i is constant on the atoms of G/N. With N = 0 this is the model of abelian_cayley.py.
--mode joint: the full model on G plus the quotient models as extra constraints, solved once.
--mode split: enumerate all solutions of the first quotient model, then solve the joint model for each.

usage: python quotient_lift.py "63,48,24;1,6,28" 5,5,5,4,2 --quot 5,5,5,1,1 --mode split --out f.json
"""
import sys, os, json, time, argparse
try:
    from ortools.sat.python import cp_model  # before numpy on Windows
except OSError:
    cp_model = None
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction
import abelian_cayley as ac


def tables(ns):
    E = ac.group(ns); lab, reps, e, idx = ac.atoms(E, ns)
    na = len(reps); size = np.bincount(lab, minlength=na)
    M = np.zeros((len(E), na)); M[np.arange(len(E)), lab] = 1
    P = ((E[reps] * (e // np.array(ns))) @ E.T) % e
    Kf = np.cos(2 * np.pi * P / e) @ M; K = np.rint(Kf).astype(np.int64)
    assert np.abs(Kf - K).max() < 1e-6
    Kd = (K.T * size[None, :]) / size[:, None]; Kd2 = np.rint(Kd).astype(np.int64)
    assert np.abs(Kd - Kd2).max() < 1e-6
    return dict(E=E, lab=lab, reps=reps, idx=idx, na=na, size=size, K=K, Kd=Kd2, ns=ns)


def PQ(b, c):
    v, th, ms = ac.spectrum(b, c); k = b[0]; d = len(b); B_ = b + [0]; C_ = [0] + c
    a_ = [k - B_[i] - C_[i] for i in range(d + 1)]
    ks = [1]
    for i in range(d): ks.append(ks[-1] * B_[i] // C_[i + 1])
    U = []
    for t in th:
        u = [Fraction(1), Fraction(t, k)]
        for i in range(1, d): u.append(((t - a_[i]) * u[i] - C_[i] * u[i - 1]) / B_[i])
        U.append(u)
    Pm = [[ks[i] * U[j][i] for i in range(d + 1)] for j in range(d + 1)]
    Qm = [[ms[j] * U[j][i] for j in range(d + 1)] for i in range(d + 1)]
    assert all(z.denominator == 1 for r in Pm + Qm for z in r)
    return v, th, ms, ks, [[int(z) for z in r] for r in Pm], [[int(z) for z in r] for r in Qm]


def add_level(m, T, f, ks, ms, Pm, Qm, tag):
    """Model on group T with fibre size f. Returns w[(atom, i)] (i = 1..d) and y[(atom, j)]."""
    d = len(ks) - 1; R = range(1, d + 1); na = T['na']; size = T['size']; K = T['K']; Kd = T['Kd']
    w = {}
    for a in range(na):
        for i in R:
            hi = min(f - (1 if a == 0 else 0), ks[i])
            w[a, i] = m.NewBoolVar(f'{tag}w{a}_{i}') if f == 1 else m.NewIntVar(0, hi, f'{tag}w{a}_{i}')
    for a in range(na):
        m.Add(sum(w[a, i] for i in R) == f - (1 if a == 0 else 0))
    for i in R: m.Add(sum(int(size[a]) * w[a, i] for a in range(na)) == ks[i])
    y = {(q, j): m.NewBoolVar(f'{tag}y{q}_{j}') for q in range(1, na) for j in R}
    for q in range(1, na): m.AddExactlyOne([y[q, j] for j in R])
    for j in R: m.Add(sum(int(size[q]) * y[q, j] for q in range(1, na)) <= ms[j])
    for q in range(1, na):
        for i in R:
            m.Add(sum(int(K[q, a]) * w[a, i] for a in range(na) if K[q, a]) == sum(Pm[j][i] * y[q, j] for j in R))
    for a in range(na):
        for j in R:
            m.Add(f * sum(int(Kd[a, q]) * y[q, j] for q in range(1, na) if Kd[a, q])
                  == sum(Qm[i][j] * w[a, i] for i in R) + (Qm[0][j] if a == 0 else 0) - (f * 1 if False else 0))
    return w, y


def link(m, G, Tq, ds, wG, wq, d):
    """sum over atoms a of G mapping into quotient atom aq of |a| x[a,i] = |aq| w_q[aq,i]."""
    R = range(1, d + 1); ds = np.array(ds)
    img = {}
    for a in range(1, G['na']):
        g = G['E'][G['reps'][a]] % ds
        aq = int(Tq['lab'][Tq['idx'][tuple(g)]])
        img.setdefault(aq, []).append(a)
    for aq in range(Tq['na']):
        for i in R:
            m.Add(sum(int(G['size'][a]) * wG[a, i] for a in img.get(aq, [])) == int(Tq['size'][aq]) * wq[aq, i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('array'); ap.add_argument('ns'); ap.add_argument('--quot', action='append', default=[])
    ap.add_argument('--mode', default='joint'); ap.add_argument('--tl', type=float, default=3600)
    ap.add_argument('--qtl', type=float, default=3600); ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--maxq', type=int, default=100000); ap.add_argument('--out', required=True)
    ap.add_argument('--fullmodel', default='yes'); ap.add_argument('--cases', default=None)
    A = ap.parse_args()
    b, c = ac.parse(A.array); ns = [int(z) for z in A.ns.split(',')]; d = len(b); R = range(1, d + 1)
    v, th, ms, ks, Pm, Qm = PQ(b, c)
    G = tables(ns); quots = [[int(z) for z in q.split(',')] for q in A.quot]
    Ts = [tables(q) for q in quots]; fs = [v // int(np.prod(q)) for q in quots]
    t0 = time.time(); log = dict(array=A.array, ns=ns, quot=quots, mode=A.mode)

    def joint(fix=None):
        m = cp_model.CpModel()
        wG, yG = add_level(m, G, 1, ks, ms, Pm, Qm, 'G')
        for i in R: m.Add(wG[0, i] == 0)
        for j in R: m.Add(sum(int(G['size'][q]) * yG[q, j] for q in range(1, G['na'])) == ms[j])
        for t, (T, f, ds) in enumerate(zip(Ts, fs, quots)):
            wq, yq = add_level(m, T, f, ks, ms, Pm, Qm, f'q{t}')
            link(m, G, T, ds, wG, wq, d)
            if t == 0 and fix is not None:
                for key, val in fix.items(): m.Add(wq[key] == val)
        return m, wG

    def check(s, wG):
        sel = [a for a in range(1, G['na']) if s.Value(wG[a, 1])]
        Sidx = np.where(np.isin(G['lab'], sel))[0]
        return ac.is_drg(G['E'], ns, Sidx, G['idx'], b, c)[0], sel

    if A.mode == 'joint':
        m, wG = joint(); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = A.tl
        s.parameters.num_search_workers = A.workers; st = s.Solve(m); log['status'] = s.StatusName(st)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE): log['drg'], log['atoms'] = check(s, wG)
    else:
        T, f = Ts[0], fs[0]
        mq = cp_model.CpModel(); wq, yq = add_level(mq, T, f, ks, ms, Pm, Qm, 'q')
        for t in range(1, len(Ts)):
            pass
        sols = []; keys = sorted(wq.keys()); fsol = open(A.out + '.sols', 'w') if A.mode == 'count' else None
        log['sol_keys'] = keys

        class CB(cp_model.CpSolverSolutionCallback):
            def __init__(s2): super().__init__(); s2.n = 0
            def on_solution_callback(s2):
                s2.n += 1
                tup = tuple(s2.Value(wq[k_]) for k_ in keys)
                if fsol is not None: fsol.write(','.join(map(str, tup)) + '\n')
                else: sols.append(dict(zip(keys, tup)))
                if s2.n >= A.maxq: s2.StopSearch()
        s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = A.qtl
        if A.maxq == 1: s.parameters.num_search_workers = A.workers   # feasibility only
        else: s.parameters.enumerate_all_solutions = True; s.parameters.num_search_workers = 1
        cb = CB(); st = s.Solve(mq, cb); log['quotient_status'] = s.StatusName(st); log['n_quotient'] = cb.n
        if fsol is not None: fsol.close()
        log['quotient_complete'] = (s.StatusName(st) in ('OPTIMAL', 'INFEASIBLE')) and cb.n < A.maxq
        uniq = {tuple(sorted(z.items())): z for z in sols}; log['n_quotient_distinct_w'] = len(uniq)
        res = []
        rng = None if A.cases is None else range(*[int(z_) for z_ in A.cases.split(':')])
        for n_, z in enumerate(uniq.values() if A.mode == 'split' else []):
            if rng is not None and n_ not in rng: continue
            m, wG = joint(fix=z); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = A.tl
            s.parameters.num_search_workers = A.workers; st = s.Solve(m); nm = s.StatusName(st)
            r = {'case': n_, 'status': nm, 'sec': round(s.WallTime(), 1)}
            if st in (cp_model.OPTIMAL, cp_model.FEASIBLE): r['drg'], r['atoms'] = check(s, wG)
            res.append(r)
            with open(A.out + '.progress', 'a') as fp: fp.write(json.dumps(r) + '\n')
        log['lifts'] = res
        log['status'] = ('INFEASIBLE' if rng is None and log['quotient_complete'] and all(r['status'] == 'INFEASIBLE' for r in res)
                         else 'FOUND' if any(r.get('drg') for r in res) else 'OPEN')
    log['seconds'] = round(time.time() - t0, 1)
    open(A.out, 'w').write(json.dumps(log) + '\n')


if __name__ == '__main__':
    main()
