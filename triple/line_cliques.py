"""27-cliques of the 'not opposite' graph of GO(2,4) (points adjacent when at distance 1, 2 or 3).

For a line L, the 27 points at distance <= 1 from L (L itself plus the 24 points collinear with a point of
L) are pairwise at distance <= 3. This script checks that
  1. the ratio (Hoffman) bound for cocliques of the opposition graph Gamma_4 is 27, so omega <= 27
     (exact, from the eigenvalues of A_4);
  2. every 27-clique through point 0 is one of the 45 line neighbourhoods containing 0
     (CP-SAT: maximise the clique size with a no-good cut for each of the 45; the optimum must be < 27).
By transitivity of T, 2 then gives exactly 1755 * 45 / 27 = 2925 cliques of size 27, one per line.
Usage: python triple/line_cliques.py        (writes triple/out_orbits/line_cliques_output.txt)
"""
import ortools  # noqa: F401  (import before numpy on this host)
from ortools.sat.python import cp_model
import os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'flags'))
import numpy as np
import go24_common as gc

OUT = os.path.join(HERE, 'out_orbits'); os.makedirs(OUT, exist_ok=True)
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)

# 1. ratio bound, exact
b, c = [10, 8, 8, 8], [1, 1, 1, 5]
kv = [1]
for i in range(4): kv.append(kv[-1] * b[i] // c[i])
n = sum(kv)
bb = b + [0]; cc = [0] + c; a = [b[0] - bb[i] - cc[i] for i in range(5)]
ev4 = []
for th in [10, 5, 1, -3, -5]:
    u = [Fraction(1), Fraction(th, 10)]
    for i in range(1, 4):
        u.append(((th - a[i]) * u[i] - cc[i] * u[i - 1]) / bb[i])
    ev4.append(kv[4] * u[4])
lmin = min(ev4)
bound = n * (-lmin) / (kv[4] - lmin)
say(f'eigenvalues of A_4 (for A_1 = 10, 5, 1, -3, -5): {[str(x) for x in ev4]}; ratio bound n(-l)/(k-l) = {bound}')

# 2. CP-SAT
Tg, _ = gc.generators()
D = gc.point_distances(Tg)
lines = gc.lines_from_distances(D)
nbhd = []
for L in lines:
    S = set(L)
    for x in L:
        S |= set(np.nonzero(D[x] == 1)[0].tolist())
    assert len(S) == 27
    assert all(D[x, y] <= 3 for x in S for y in S)
    if 0 in S:
        nbhd.append(sorted(S))
say(f'line neighbourhoods: 2925 of size 27, all cliques; {len(nbhd)} contain point 0')
cand = [int(v) for v in np.nonzero((D[0] >= 1) & (D[0] <= 3))[0]]
pos = {v: i for i, v in enumerate(cand)}
opp = [(i, j) for i in range(len(cand)) for j in range(i + 1, len(cand)) if D[cand[i], cand[j]] == 4]
say(f'candidates {len(cand)}, opposite pairs among them {len(opp)}')

def solve(cuts, tl=1800):
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f'x{i}') for i in range(len(cand))]
    for i, j in opp:
        m.AddBoolOr([x[i].Not(), x[j].Not()])
    for S in cuts:
        m.Add(sum(x[pos[v]] for v in S if v != 0) <= 25)
    m.Maximize(sum(x))
    s = cp_model.CpSolver(); s.parameters.num_workers = 8; s.parameters.max_time_in_seconds = tl
    st = s.Solve(m)
    return s.StatusName(st), int(s.ObjectiveValue()) + 1, int(s.BestObjectiveBound()) + 1, s.WallTime()

st, val, bd, wt = solve([])
say(f'max clique through 0, no cuts: status {st}, size {val}, bound {bd}, {wt:.0f} s')
st, val, bd, wt = solve(nbhd)
say(f'max clique through 0 avoiding the 45 line neighbourhoods: status {st}, size {val}, bound {bd}, {wt:.0f} s')
open(os.path.join(OUT, 'line_cliques_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
