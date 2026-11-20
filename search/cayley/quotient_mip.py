"""SCIP version of the quotient model in quotient_lift.py (second solver for quotient infeasibility).

usage: python quotient_mip.py "63,48,24;1,6,28" 5,5,5 8 --out f.json
(quotient group Z_{d_1} + ... , fibre size f = v / |G/N|)
"""
import sys, os, json, time, argparse
try:
    from ortools.linear_solver import pywraplp
except OSError:
    pywraplp = None
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import abelian_cayley as ac
from quotient_lift import tables, PQ


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('array'); ap.add_argument('qs'); ap.add_argument('f', type=int)
    ap.add_argument('--out', required=True); ap.add_argument('--tl', type=float, default=3600)
    A = ap.parse_args(); b, c = ac.parse(A.array); d = len(b); R = range(1, d + 1)
    v, th, ms, ks, Pm, Qm = PQ(b, c); qs = [int(z) for z in A.qs.split(',')]; f = A.f
    assert f * int(np.prod(qs)) == v
    T = tables(qs); na = T['na']; size = T['size']; K = T['K']; Kd = T['Kd']
    s = pywraplp.Solver.CreateSolver('SCIP'); s.SetTimeLimit(int(A.tl * 1000))
    w = {(a, i): s.IntVar(0, min(f - (a == 0), ks[i]), f'w{a}_{i}') for a in range(na) for i in R}
    y = {(q, j): s.BoolVar(f'y{q}_{j}') for q in range(1, na) for j in R}
    for a in range(na): s.Add(sum(w[a, i] for i in R) == f - (1 if a == 0 else 0))
    for i in R: s.Add(sum(int(size[a]) * w[a, i] for a in range(na)) == ks[i])
    for q in range(1, na): s.Add(sum(y[q, j] for j in R) == 1)
    for j in R: s.Add(sum(int(size[q]) * y[q, j] for q in range(1, na)) <= ms[j])
    for q in range(1, na):
        for i in R:
            s.Add(sum(int(K[q, a]) * w[a, i] for a in range(na)) == sum(Pm[j][i] * y[q, j] for j in R))
    for a in range(na):
        for j in R:
            s.Add(f * sum(int(Kd[a, q]) * y[q, j] for q in range(1, na))
                  == sum(Qm[i][j] * w[a, i] for i in R) + (Qm[0][j] if a == 0 else 0))
    t0 = time.time(); st = s.Solve()
    nm = {s.OPTIMAL: 'FEASIBLE', s.FEASIBLE: 'FEASIBLE', s.INFEASIBLE: 'INFEASIBLE'}.get(st, 'UNKNOWN')
    open(A.out, 'w').write(json.dumps(dict(array=A.array, quotient=qs, f=f, solver='SCIP', status=nm,
                                           seconds=round(time.time() - t0, 1))) + '\n')


if __name__ == '__main__':
    main()
