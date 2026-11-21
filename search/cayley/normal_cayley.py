"""Normal Cayley graphs (S a union of conjugacy classes) with a given intersection array, from rational character
data exported by export_normal.g.

A normal Cayley graph with integral eigenvalues has S power-closed (Godsil-Spiga), i.e. a union of rational
classes R, and the eigenspaces are unions of Galois orbits O of irreducibles. With T[O][R] the rationalized character
value, |O| its norm and chi_O(1) = T[O][1]/|O|:
  (P)  sum_R w_O(R) x[R,i] = sum_j P_ji y[O,j],    w_O(R) = |R| T[O][R] / T[O][1]   (an integer),
  (Q)  sum_O chi_O(1) T[O][R] y[O,j] = sum_i Q_ij x[R,i],
  sum_R |R| x[R,i] = k_i,   sum_O |O| chi_O(1)^2 y[O,j] = m_j,
and chi_O(1) must be a multiple of the Krein denominator D_j (normal_krein_filter.py). Solved with CP-SAT.
A solution is only a candidate (the group itself is not in the data) and has to be checked in GAP.

usage: python normal_cayley.py tables.txt "56,42,20;1,6,28" --out results.jsonl [--tl 600]
"""
import sys, os, json, time, argparse, re, ast
try:
    from ortools.sat.python import cp_model  # before numpy on Windows
except OSError:
    cp_model = None
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from math import lcm
import abelian_cayley as ac
from quotient_lift import PQ


def top_level_lists(text):
    """Split '[..] [..] [[..],..]' into its three top-level lists (GAP may print '[ [' with a space)."""
    out, depth, start = [], 0, None
    for k, ch in enumerate(text):
        if ch == '[':
            if depth == 0: start = k
            depth += 1
        elif ch == ']':
            depth -= 1
            if depth == 0: out.append(text[start:k + 1])
    assert len(out) == 3, len(out)
    return out


def parse_tables(path):
    for line in open(path):
        line = line.strip()
        if not line.startswith('G '): continue
        _, n, i, rest = line.split(' ', 3)
        # rest = "[sizes] [norms] [[values]]"
        parts = top_level_lists(rest)
        sizes, norms, vals = (ast.literal_eval(q) for q in parts)
        yield int(n), int(i), sizes, norms, vals


def krein_den(b, c):
    v, th, ms, ks, Pm, Qm = PQ(b, c); d = len(b)
    # q^h_ij = (1/v) sum_l Q_li Q_lj Q_lh / k_l ... via standard formula with P, Q
    D = {}
    for h in range(1, d + 1):
        dens = []
        for i in range(d + 1):
            for j in range(d + 1):
                q = Fraction(ms[i] * ms[j], v) * sum(Fraction(ks[l]) * Fraction(Qm[l][i], ms[i]) * Fraction(Qm[l][j], ms[j])
                                                    * Fraction(Qm[l][h], ms[h]) for l in range(d + 1))
                dens.append(q.denominator)
        D[h] = lcm(*dens)
    return D


def solve(n, gid, sizes, norms, vals, arr, tl, workers, limit=20):
    b, c = ac.parse(arr); d = len(b); R = range(1, d + 1)
    v, th, ms, ks, Pm, Qm = PQ(b, c); assert v == n == sum(sizes)
    D = krein_den(b, c)
    nr = len(sizes); no = len(norms)
    deg = [Fraction(vals[o][0], norms[o]) for o in range(no)]
    assert all(x.denominator == 1 for x in deg) and deg[0] == 1 and norms[0] == 1
    deg = [int(x) for x in deg]
    assert sum(norms[o] * deg[o] ** 2 for o in range(no)) == n
    W = [[Fraction(sizes[r] * vals[o][r], vals[o][0]) for r in range(nr)] for o in range(no)]
    assert all(x.denominator == 1 for row in W for x in row)
    W = [[int(x) for x in row] for row in W]
    m = cp_model.CpModel()
    x = {(r, i): m.NewBoolVar(f'x{r}_{i}') for r in range(1, nr) for i in R}
    y = {(o, j): m.NewBoolVar(f'y{o}_{j}') for o in range(1, no) for j in R}
    for r in range(1, nr): m.AddExactlyOne([x[r, i] for i in R])
    for o in range(1, no): m.AddExactlyOne([y[o, j] for j in R])
    for i in R: m.Add(sum(sizes[r] * x[r, i] for r in range(1, nr)) == ks[i])
    for j in R: m.Add(sum(norms[o] * deg[o] ** 2 * y[o, j] for o in range(1, no)) == ms[j])
    for o in range(1, no):
        for j in R:
            if deg[o] % D[j]: m.Add(y[o, j] == 0)
        for i in R:
            m.Add(sum(W[o][r] * x[r, i] for r in range(1, nr) if W[o][r]) == sum(Pm[j][i] * y[o, j] for j in R))
    for r in range(1, nr):
        for j in R:
            m.Add(sum(deg[o] * vals[o][r] * y[o, j] for o in range(1, no) if vals[o][r])
                  == sum(Qm[i][j] * x[r, i] for i in R))
    sols = []; status = None; t0 = time.time()
    while len(sols) < limit:
        s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_search_workers = workers
        st = s.Solve(m); status = s.StatusName(st)
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): break
        sel = [r for r in range(1, nr) if s.Value(x[r, 1])]
        sols.append(sel)
        m.AddBoolOr([x[r, 1].Not() for r in sel] + [x[r, 1] for r in range(1, nr) if r not in sel])
    return dict(n=n, id=gid, array=arr, nratcl=nr, norbits=no, status=status, n_candidates=len(sols),
                candidates_S1_ratclasses=sols, seconds=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('tables'); ap.add_argument('array')
    ap.add_argument('--out', required=True); ap.add_argument('--tl', type=float, default=600)
    ap.add_argument('--workers', type=int, default=8)
    A = ap.parse_args()
    done = set()
    if os.path.exists(A.out):
        for line in open(A.out): z = json.loads(line); done.add((z['n'], z['id'], z['array']))
    b, c = ac.parse(A.array); v = PQ(b, c)[0]
    for n, gid, sizes, norms, vals in parse_tables(A.tables):
        if n != v or (n, gid, A.array) in done: continue
        res = solve(n, gid, sizes, norms, vals, A.array, A.tl, A.workers)
        with open(A.out, 'a') as f: f.write(json.dumps(res) + '\n')


if __name__ == '__main__':
    main()
