"""Recover centres and collinearity of O(2) from one component of the 13-overlap graph of its thin suboctagons.

Builds O(2) from generators of 2F4(2) on 1755 points (generators.npy, ATLAS TF42d2 c,d), takes the orbit of one
GO(2,1), joins two copies when they share a 13-point ball, and keeps one component (12480 vertices). The centre labels
are then used only to check the intrinsic rules:
  rule 1  opposite edges of a K4 have the same centre;
  rule 2  in a chordless 4-cycle B-B1-C1-C with |N(B)&N(C1)| = |N(B1)&N(C)| = 4, opposite edges have the same centre
          iff there are 8 edges among the vertices outside the cycle adjacent to exactly two of its vertices;
  rule 3  two classes from rules 1-2 have the same centre iff no edge joins their vertex sets;
  rule 4  two centres are collinear iff they label two edges of a triangle.
Run: python overlap_recognition.py generators.npy   (about 40 minutes, mostly rule 2)
"""
import sys, random
from collections import defaultdict, Counter
import numpy as np, networkx as nx
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
def octagon(gens, n=1755):
    idx = np.arange(n*n); a, b = idx//n, idx % n
    r = np.concatenate([idx]*len(gens)); c = np.concatenate([g[a]*n+g[b] for g in gens])
    _, lab = connected_components(csr_matrix((np.ones(len(r), np.int8), (r, c)), shape=(n*n,)*2), directed=True, connection='weak')
    L = lab.reshape(n, n); col = [l for l in np.unique(L[0]) if (L[0] == l).sum() == 10][0]
    A = L == col; D = shortest_path(csr_matrix(A), unweighted=True).astype(np.int8)
    return A, D
def closure(seed, A, D, cap=50):
    cur = set(seed)
    while True:
        cl = np.array(sorted(cur)); sub = D[np.ix_(cl, cl)]; new = set(cur)
        for i, j in zip(*np.nonzero(np.triu((sub >= 1) & (sub <= 3)))):
            x, y = cl[i], cl[j]; d = D[x, y]
            new.update(map(int, np.flatnonzero(D[x] + D[y] == d)))
            if d == 1: new.update(map(int, np.flatnonzero(A[x] & A[y])))
        if len(new) > cap: return None
        if new == cur: return cur
        cur = new
def main(path):
    gens = [g.astype(np.int64) for g in np.load(path)]
    A, D = octagon(gens)
    rng = np.random.default_rng(3)
    while True:
        s = [0, int(rng.choice(np.flatnonzero(D[0] == 4)))] + list(map(int, rng.choice(1755, 2, replace=False)))
        c = closure(s, A, D)
        if c is not None and len(c) == 45: break
    B0 = frozenset(c); seen = {B0}; copies = [B0]; q = [B0]
    while q:
        B = q.pop()
        for g in gens:
            Bn = frozenset(int(g[p]) for p in B)
            if Bn not in seen: seen.add(Bn); copies.append(Bn); q.append(Bn)
    assert len(copies) == 24960
    balls = defaultdict(list)
    for i, B in enumerate(copies):
        Bs = np.array(sorted(B)); sub = D[np.ix_(Bs, Bs)]
        for j, p in enumerate(Bs): balls[frozenset(Bs[sub[j] <= 2].tolist())].append((i, int(p)))
    G = nx.Graph(); cen = {}
    for v in balls.values():
        assert len(v) == 2 and v[0][1] == v[1][1]
        (x, p), (y, _) = v; G.add_edge(x, y); cen[(x, y)] = cen[(y, x)] = p
    comps = sorted(nx.connected_components(G), key=min); assert [len(c) for c in comps] == [12480, 12480]
    H = G.subgraph(comps[0]); nbr = {v: set(H[v]) for v in H}
    E = list(H.edges()); eid = {}
    for i, (x, y) in enumerate(E): eid[(x, y)] = eid[(y, x)] = i
    par = list(range(len(E)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    def union(x, y): par[find(x)] = find(y)
    bad = 0
    for (x, y) in E:                               # rule 1
        cm = nbr[x] & nbr[y]
        for c in cm:
            for d in cm:
                if c < d and d in nbr[c]:
                    union(eid[(x, y)], eid[(c, d)]); union(eid[(x, c)], eid[(y, d)]); union(eid[(x, d)], eid[(y, c)])
    for (B, B1) in E:                              # rule 2
        for C in nbr[B] - {B1}:
            for C1 in (nbr[B1] & nbr[C]) - {B}:
                if C1 in nbr[B] or C in nbr[B1]: continue
                if len(nbr[B] & nbr[C1]) != 4 or len(nbr[B1] & nbr[C]) != 4: continue
                cyc = (B, B1, C1, C); cnt = Counter(z for v in cyc for z in nbr[v])
                two = [z for z, k in cnt.items() if k == 2 and z not in cyc]
                if sum(1 for i in range(len(two)) for j in range(i+1, len(two)) if two[j] in nbr[two[i]]) == 8:
                    union(eid[(B, B1)], eid[(C, C1)]); union(eid[(B, C)], eid[(B1, C1)])
    cls = defaultdict(list)
    for i in range(len(E)): cls[find(i)].append(i)
    assert all(len({cen[E[i]] for i in c}) == 1 for c in cls.values()), 'rules 1-2 mixed two centres'
    print('rules 1-2: classes', Counter(len(c) for c in cls.values()))
    keys = list(cls); verts = {k: {v for i in cls[k] for v in E[i]} for k in keys}
    owner = defaultdict(set)
    for k in keys:
        for v in verts[k]: owner[v].add(k)
    partner = {}
    for k in keys:                                 # rule 3
        touching = Counter(k2 for v in verts[k] for w in nbr[v] for k2 in owner[w])
        cand = [k2 for k2 in keys if k2 != k and touching[k2] == 0 and not (verts[k] & verts[k2])]
        assert len(cand) == 1; partner[k] = cand[0]
    label = {}
    for k in keys:
        root = min(k, partner[k]); label[k] = root
    ok = all((cen[E[cls[k][0]]] == cen[E[cls[partner[k]][0]]]) for k in keys)
    print('rule 3: centres', len(set(label.values())), 'all correct', ok)
    col = set()                                    # rule 4 on recovered labels
    lab_of = {i: label[find(i)] for i in range(len(E))}
    for (x, y) in E:
        for c in nbr[x] & nbr[y]: col.add(tuple(sorted((lab_of[eid[(x, y)]], lab_of[eid[(x, c)]]))))
    truth = {k: cen[E[cls[k][0]]] for k in keys}
    pairs = {tuple(sorted((truth[a], truth[b]))) for a, b in col}
    print('rule 4: collinear pairs', len(pairs), 'all collinear', all(D[a, b] == 1 for a, b in pairs), 'expected', int((D == 1).sum()//2))
if __name__ == '__main__': main(sys.argv[1])
