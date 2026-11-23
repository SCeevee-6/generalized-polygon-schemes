"""Independent check of one-component recognition (own BFS code, ground truth from geom.py)."""
import numpy as np, json, sys
from collections import deque
from itertools import combinations
def bfs_all(adj, n, sources=None):
    sources = range(n) if sources is None else sources
    out = np.full((len(list(sources)), n), 10**6, dtype=np.int32)
    for k, s in enumerate(sources):
        d = out[k]; d[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if d[v] > d[u]+1: d[v] = d[u]+1; q.append(v)
    return out
def recognise(edges, n=512):
    """edges: list of (u,v) of an uncoloured 5-regular graph. Returns reconstruction or reason for failure."""
    adj = [[] for _ in range(n)]
    for u,v in edges: adj[u].append(v); adj[v].append(u)
    dist = bfs_all(adj, n)
    E = np.array(edges); m = len(E)
    dd = np.minimum.reduce([dist[E[:,a]][:,E[:,b]] for a in (0,1) for b in (0,1)])  # delta(e,f)
    np.fill_diagonal(dd, -1)
    mates = [np.flatnonzero(dd[i] >= 6) for i in range(m)]
    cnt = np.array([len(x) for x in mates])
    if not (cnt == 1).all(): return dict(ok=False, why='mate counts', hist=np.bincount(cnt).tolist())
    mate = np.array([x[0] for x in mates]); assert (mate[mate] == np.arange(m)).all()
    # boundary ids
    bid = -np.ones(m, int); nb = 0
    for i in range(m):
        if bid[i] < 0: bid[i] = bid[mate[i]] = nb; nb += 1
    assert nb == m//2
    # H: points 0..n-1, lines n..n+m-1, boundary n+m..n+m+nb-1
    N = n+m+nb; H = [[] for _ in range(N)]
    for i,(u,v) in enumerate(E):
        L = n+i; b = n+m+bid[i]
        for p in (u,v,b): H[L].append(p); H[p].append(L)
    hb = bfs_all(H, N, [n+m+j for j in range(nb)])[:, n+m:]   # boundary-boundary incidence distances in H
    groups = [[j] for j in range(nb)]; levels = []
    for L, width in zip((2,4,6,8,10,12), (2,4,2,4,2,5)):
        k = len(groups); thr = 16-L
        g = np.array(groups)
        dmin = hb[g.ravel()][:, g.ravel()].reshape(k, g.shape[1], k, g.shape[1]).min(axis=(1,3))
        ok = dmin >= thr; np.fill_diagonal(ok, True)
        # must be disjoint union of cliques of size width
        seen = np.zeros(k, bool); new = []
        for i in range(k):
            if seen[i]: continue
            cl = np.flatnonzero(ok[i])
            if len(cl) != width or not ok[np.ix_(cl, cl)].all(): return dict(ok=False, why='level %d'%L, row=ok.sum(1).tolist())
            seen[cl] = True; new.append(sorted(x for c in cl for x in groups[c]))
        levels.append(new); groups = new
    return dict(ok=True, mate=mate, bid=bid, levels=levels, dist=dist, adj=adj)
