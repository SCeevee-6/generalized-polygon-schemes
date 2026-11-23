"""Independent check of the one-component recognition of O(2).
Uses the stored generators by default. ATLAS_DIR can select GAP permutations instead."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geom, recog
from scipy.sparse.csgraph import connected_components
G, A, D, lines = geom.octagon()
assert np.bincount(D[0]).tolist() == [1, 10, 80, 640, 1024]
O = np.flatnonzero(D[0] == 4); _, lab = connected_components(A[np.ix_(O, O)])
assert np.bincount(lab).tolist() == [512, 512]
def ancestors(b):
    p2 = [p for p in np.flatnonzero(A[b]) if D[0, p] == 2][0]; sib = np.flatnonzero(A[b] & A[p2])[0]
    p1 = [p for p in np.flatnonzero(A[p2]) if D[0, p] == 1][0]; p2s = np.flatnonzero(A[p2] & A[p1])[0]
    p1s = np.flatnonzero(A[p1] & A[0])[0]
    return [frozenset((b, sib)), p2, frozenset((p2, p2s)), p1, frozenset((p1, p1s)), 0]
def closure(seed, allowed, cap=60):
    cur = set(seed)
    while True:
        cl = np.array(sorted(cur)); sub = D[np.ix_(cl, cl)]; new = set(cur)
        for i, j in zip(*np.nonzero(np.triu((sub >= 1) & (sub <= 3)))):
            a, b = cl[i], cl[j]; d = D[a, b]
            new.update(map(int, np.flatnonzero(D[a] + D[b] == d)))
            if d == 1: new.update(map(int, np.flatnonzero(A[a] & A[b])))
        if not new <= allowed or len(new) > cap: return None
        if new == cur: return cur
        cur = new
def thin(P):
    import networkx as nx
    Ls = [l for l in lines if set(l) <= P]
    if len(P) != 45 or len(Ls) != 30: return False
    g = nx.Graph(); g.add_edges_from((int(p), ('L', k)) for k, l in enumerate(Ls) for p in l)
    return all(d == 2 for n, d in g.degree() if not isinstance(n, tuple)) and nx.girth(g) == 16 and nx.diameter(g) == 8
rng = np.random.default_rng(1)
for c in (0, 1):
    Y = O[lab == c]; E = [tuple(map(int, e)) for e in np.argwhere(np.triu(A[np.ix_(Y, Y)]))]
    r = recog.recognise(E)                       # uses only the uncoloured edge list
    assert r['ok']
    third = np.array([np.flatnonzero(A[Y[u]] & A[Y[v]])[0] for u, v in E])
    pairing_ok = all(third[i] == third[r['mate'][i]] for i in range(len(E)))
    bt = {r['bid'][i]: third[i] for i in range(len(E))}
    levels_ok = all(len({ancestors(bt[j])[lev] for j in g}) == 1 for lev, gs in enumerate(r['levels']) for g in gs)
    assert pairing_ok and levels_ok
    assert [len(groups) for groups in r['levels']] == [320, 80, 40, 10, 5, 1]
    # Two neighbours per line and four neighbours from this component at each boundary point.
    boundary = np.flatnonzero(D[0] == 3)
    assert np.all(A[np.ix_(boundary, Y)].sum(axis=1) == 4)
    allowed = set(np.flatnonzero(D[0] <= 3).tolist()) | set(Y.tolist()); found = None
    for _ in range(5000):
        s = [0] + list(map(int, rng.choice(Y, 3, replace=False))); cl = closure(s, allowed)
        if cl is not None and len(cl) == 45 and thin(cl): found = s; break
    assert found is not None
    # Every selected line belongs to the forced half. A line reaching shell
    # four has exactly its two outer endpoints there, joined by a K-edge.
    for line in lines:
        if set(line) <= cl and (D[0, line] == 4).any():
            outer = line[D[0, line] == 4]
            assert len(outer) == 2 and set(outer) <= set(Y)
    print('component', c, 'mates forced', True, 'pairing correct', pairing_ok, 'six merges correct', levels_ok, 'GO(2,1) seed', found)
