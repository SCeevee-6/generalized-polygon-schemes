"""Second solver for abelian_cayley.py: the same primal-dual atom model written as a MIP and solved with
SCIP (through OR-Tools' pywraplp), so that CP-SAT infeasibility answers get an independent check.

usage: python abelian_cayley_mip.py "52,42,16;1,6,28" 25,25 --out file.json [--tl 1800]
"""
import sys, json, argparse, time
try:
    from ortools.linear_solver import pywraplp  # must load before numpy on Windows
except OSError:
    pywraplp = None
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction
import abelian_cayley as ac


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('array'); ap.add_argument('ns'); ap.add_argument('--out', required=True)
    ap.add_argument('--tl', type=float, default=1800)
    A = ap.parse_args()
    b, c = ac.parse(A.array); ns = [int(z) for z in A.ns.split(',')]
    v, th, ms = ac.spectrum(b, c); k = b[0]; d = len(b)
    E = ac.group(ns); lab, reps, e, idx = ac.atoms(E, ns)
    na = len(reps); size = np.bincount(lab, minlength=na)
    M = np.zeros((v, na)); M[np.arange(v), lab] = 1
    P = ((E[reps] * (e // np.array(ns))) @ E.T) % e
    Kf = np.cos(2 * np.pi * P / e) @ M; K = np.rint(Kf).astype(np.int64)
    assert np.abs(Kf - K).max() < 1e-6
    B_ = b + [0]; C_ = [0] + c; a_ = [k - B_[i] - C_[i] for i in range(d + 1)]
    ks = [1]
    for i in range(d): ks.append(ks[-1] * B_[i] // C_[i + 1])
    U = []
    for t in th:
        u = [Fraction(1), Fraction(t, k)]
        for i in range(1, d): u.append(((t - a_[i]) * u[i] - C_[i] * u[i - 1]) / B_[i])
        U.append(u)
    Pm = [[int(ks[i] * U[j][i]) for i in range(d + 1)] for j in range(d + 1)]
    Qm = [[int(ms[j] * U[j][i]) for j in range(d + 1)] for i in range(d + 1)]
    Kd = np.rint((K.T * size[None, :]) / size[:, None]).astype(np.int64)
    s = pywraplp.Solver.CreateSolver('SCIP')
    s.SetTimeLimit(int(A.tl * 1000))
    R = range(1, d + 1); NA = range(1, na)
    x = {(a, i): s.BoolVar(f'x{a}_{i}') for a in NA for i in R}
    y = {(q, j): s.BoolVar(f'y{q}_{j}') for q in NA for j in R}
    for a in NA: s.Add(sum(x[a, i] for i in R) == 1)
    for q in NA: s.Add(sum(y[q, j] for j in R) == 1)
    for i in R: s.Add(sum(int(size[a]) * x[a, i] for a in NA) == ks[i])
    for j in R: s.Add(sum(int(size[q]) * y[q, j] for q in NA) == ms[j])
    for q in NA:
        for i in R:
            s.Add(sum(int(K[q, a]) * x[a, i] for a in NA if K[q, a]) == sum(Pm[j][i] * y[q, j] for j in R))
    for a in NA:
        for j in R:
            s.Add(sum(int(Kd[a, q]) * y[q, j] for q in NA if Kd[a, q]) == sum(Qm[i][j] * x[a, i] for i in R))
    t0 = time.time(); st = s.Solve()
    names = {s.OPTIMAL: 'FEASIBLE', s.FEASIBLE: 'FEASIBLE', s.INFEASIBLE: 'INFEASIBLE', s.NOT_SOLVED: 'UNKNOWN'}
    res = dict(array=A.array, ns=ns, natoms=na, solver='SCIP', status=names.get(st, str(st)),
               seconds=round(time.time() - t0, 1))
    if res['status'] == 'FEASIBLE':
        sel = [a for a in NA if x[a, 1].solution_value() > 0.5]
        Sidx = np.where(np.isin(lab, sel))[0]
        res['drg'] = ac.is_drg(E, ns, Sidx, idx, b, c)[0]; res['atoms'] = sel
    open(A.out, 'w').write(json.dumps(res) + '\n')


if __name__ == '__main__':
    main()
