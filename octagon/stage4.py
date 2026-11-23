"""Stage 4 for the covers left after stage 3. By the boundary count every boundary point carries exactly two K-edges,
which are mates; so the mate graph of each category has a perfect matching that is the boundary pairing.
1. Unit propagation: a lifted edge with exactly one remaining mate is paired with it (both leave the graph); a lifted
   edge with no remaining mate rejects the cover.
2. With the forced boundary points build the partial incidence graph H (points, 1280 lines, forced boundary points).
   A forced boundary point needs a sibling at H-distance >= 14: another forced point of its category, or a remaining
   mate pair (distance 1 + the smaller H-distance to its two lines). None rejects the cover.
Usage: python stage4.py stage3_fp0.json"""
import sys, os, json, numpy as np
from collections import deque
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *
import stage3
def test(vol):
    E = np.array(stage3.edges_of(vol)); adj = cover(vol); dist = np.array([bfs(adj, s) for s in range(512)])
    catv = np.repeat([stage3.S.CATS.index(CAT[k]) for k in range(40)], 32)
    forced = []; options = {}
    for c in range(10):
        ids = np.flatnonzero(catv == c); X = E[ids]; flat = X.ravel()
        dd = dist[np.ix_(flat, flat)].reshape(len(X), 2, len(X), 2).min((1, 3)); np.fill_diagonal(dd, 0)
        nb = {i: set(np.flatnonzero(dd[i] >= 6)) for i in range(len(X))}
        alive = set(range(len(X))); q = deque(i for i in alive if len(nb[i]) <= 1)
        while q:
            i = q.popleft()
            if i not in alive: continue
            if len(nb[i]) == 0: return 'reject: lifted edge without mate', c
            if len(nb[i]) == 1:
                j = next(iter(nb[i])); forced.append((c, int(ids[i]), int(ids[j])))
                for x in (i, j):
                    alive.discard(x)
                    for y in nb[x]:
                        nb[y].discard(x)
                        if y in alive and len(nb[y]) <= 1: q.append(y)
                    nb[x] = set()
        options[c] = [(int(ids[i]), int(ids[j])) for i in alive for j in nb[i] if i < j]
    nl = len(E); rows = []; cols = []
    for i, (u, v) in enumerate(E): rows += [u, v]; cols += [512+i, 512+i]
    for j, (c, e1, e2) in enumerate(forced): rows += [512+e1, 512+e2]; cols += [512+nl+j]*2
    N = 512 + nl + len(forced)
    H = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N)); H = H + H.T
    if not forced: return 'no forced pair', None
    dH = shortest_path(H, indices=np.arange(512+nl, N), unweighted=True)
    for j, (c, e1, e2) in enumerate(forced):
        cand = sum(1 for j2, (c2, _, _) in enumerate(forced) if c2 == c and j2 != j and dH[j, 512+nl+j2] >= 14)
        cand += sum(1 for (a, b) in options[c] if 1 + min(dH[j, 512+a], dH[j, 512+b]) >= 14)
        if cand == 0: return 'reject: forced boundary point without sibling', c
    return 'survives', len(forced)
if __name__ == '__main__':
    S3 = json.load(open(sys.argv[1])); rem = [o for o in S3 if not o['rejected']]; out = []
    for o in rem:
        r = test(o['vol']); out.append(dict(vol=o['vol'], result=list(map(str, r)))); print(r, flush=True)
    json.dump(out, open(sys.argv[1].replace('stage3', 'stage4'), 'w'))
