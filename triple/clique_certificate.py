"""Second, solver-independent check of the 27-clique result for GO(2,4) (see line_cliques.py for the CP-SAT run).

Graph: points of GO(2,4), adjacent when not opposite (distance 1, 2 or 3 in the collinearity graph).
  1. T = 2F4(2)' is transitive on the 1755 points (orbit computation from the ATLAS generators), so it is
     enough to find all maximum cliques through point 0.
  2. Method A: igraph (Cliquer back end) - clique number and all largest cliques of the graph on the 730
     points not opposite 0.
  3. Method B: our own branch and bound (bitsets, greedy-colouring bound, Tomita-style) enumerating every
     clique of size 26 in the same graph. No external solver.
  4. Both lists are compared with the 45 line neighbourhoods through 0; 1755 * 45 / 27 = 2925.
  5. Hoffman equality A_4 1_C = 16 (1 - 1_C) is checked for every one of the 2925 line neighbourhoods.
Usage: python triple/clique_certificate.py     (writes triple/out_orbits/clique_certificate_output.txt)
"""
import os, sys, time, platform
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'flags'))
import numpy as np
import igraph as ig
import go24_common as gc

OUT = os.path.join(HERE, 'out_orbits'); os.makedirs(OUT, exist_ok=True)
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)

say(f'python {platform.python_version()}, numpy {np.__version__}, python-igraph {ig.__version__}')
Tg, _ = gc.generators()
lab = gc.orbits(Tg, 1755)
say(f'T-orbits on points: {len(set(lab.tolist()))} (sizes {sorted(np.bincount(lab).tolist())})')
assert len(set(lab.tolist())) == 1
D = gc.point_distances(Tg)
say(f'distance distribution from point 0: {np.bincount(D[0]).tolist()}')
lines = gc.lines_from_distances(D)
n = D.shape[0]

# line neighbourhoods
C = np.zeros((n, len(lines)), dtype=np.int64)
for j, L in enumerate(lines):
    S = set(L)
    for x in L:
        S |= set(np.nonzero(D[x] == 1)[0].tolist())
    C[list(S), j] = 1
assert (C.sum(0) == 27).all()
nb = {frozenset(np.nonzero(C[:, j])[0].tolist()) for j in range(len(lines))}
say(f'lines {len(lines)}, distinct line neighbourhoods {len(nb)}')

# 5. Hoffman equality on every line
A4 = (D == 4).astype(np.int64)
lhs = A4 @ C
ok = np.array_equal(lhs, 16 * (1 - C))
say(f'Hoffman equality A_4 1_C = 16(1 - 1_C) on all {len(lines)} line neighbourhoods: {ok}')
assert ok

# graph on the points not opposite 0
cand = np.nonzero((D[0] >= 1) & (D[0] <= 3))[0]
sub = D[np.ix_(cand, cand)] <= 3
np.fill_diagonal(sub, False)
m = len(cand)
say(f'candidates {m}, edges {int(sub.sum()) // 2}')
through0 = {S for S in nb if 0 in S}
say(f'line neighbourhoods containing 0: {len(through0)}')

# Method A: igraph / Cliquer
iu = np.triu_indices(m, 1); msk = sub[iu]
G = ig.Graph(n=m, edges=list(zip(iu[0][msk].tolist(), iu[1][msk].tolist())))
t = time.time()
LA = G.largest_cliques()
tA = time.time() - t
setA = {frozenset([0] + cand[list(c)].tolist()) for c in LA}
say(f'A (igraph largest_cliques): size {len(LA[0]) + 1} through 0, {len(LA)} cliques, '
    f'{len(setA)} distinct, equal to the line neighbourhoods through 0: {setA == through0}, {tA:.0f} s')

# Method B: own branch and bound, enumerate all cliques of size TARGET in the candidate graph
TARGET = 26
nbr = [0] * m
for i in range(m):
    b = 0
    for j in np.nonzero(sub[i])[0].tolist():
        b |= 1 << j
    nbr[i] = b

def colour_bound(P):
    """Greedy colouring of P; returns vertices in colour order with their colour numbers."""
    order, cols, k, U = [], [], 0, P
    while U:
        k += 1
        Q = U
        while Q:
            v = (Q & -Q).bit_length() - 1
            Q &= ~(1 << v) & ~nbr[v]
            U &= ~(1 << v)
            order.append(v); cols.append(k)
    return order, cols

found = []
nodes = [0]
def expand(R, P):
    nodes[0] += 1
    order, cols = colour_bound(P)
    for idx in range(len(order) - 1, -1, -1):
        if len(R) + cols[idx] < TARGET:
            return
        v = order[idx]
        R2 = R + [v]
        if len(R2) == TARGET:
            found.append(tuple(R2))
        else:
            P2 = P & nbr[v]
            if P2:
                expand(R2, P2)
        P &= ~(1 << v)

t = time.time()
expand([], (1 << m) - 1)
tB = time.time() - t
setB = {frozenset([0] + cand[list(c)].tolist()) for c in found}
say(f'B (own branch and bound): {len(found)} cliques of size {TARGET} in the candidate graph, {len(setB)} distinct, '
    f'equal to the line neighbourhoods through 0: {setB == through0}, {nodes[0]} nodes, {tB:.0f} s')
# a clique of size 27 in the candidate graph would contain 27 cliques of size 26, so B also shows omega = 26 there
say(f'conclusion: maximum non-opposite sets have size 27; those through 0 are exactly the {len(through0)} line '
    f'neighbourhoods; by transitivity there are {n * len(through0) // 27} in total, one per line')
open(os.path.join(OUT, 'clique_certificate_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
