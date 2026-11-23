"""Base graph (Clebsch = even 5-words, Hamming distance 4), cycles, voltage covers, mate/sibling tests."""
import numpy as np
from collections import deque
VERTS = [x for x in range(32) if bin(x).count('1') % 2 == 0]
VI = {v: i for i, v in enumerate(VERTS)}
EDGES = sorted({tuple(sorted((u, u ^ (31 ^ (1 << i))))) for u in VERTS for i in range(5)})
EI = {e: k for k, e in enumerate(EDGES)}
COLOR = [31 ^ a ^ b for a, b in EDGES]          # single bit = coordinate fixed along the edge
BIT = [a & c for (a, b), c in zip(EDGES, COLOR)]  # value of that coordinate
CAT = [(c, b) for c, b in zip(COLOR, BIT)]
NB = {u: [u ^ (31 ^ (1 << i)) for i in range(5)] for u in VERTS}
def tree_and_chords():
    par = {0: None}; q = [0]
    for u in q:
        for v in NB[u]:
            if v not in par: par[v] = u; q.append(v)
    tree = {EI[tuple(sorted((u, p)))] for u, p in par.items() if p is not None}
    chords = [k for k in range(40) if k not in tree]
    return par, q, sorted(tree), chords
PAR, BFSORD, TREE, CHORDS = tree_and_chords()
def simple_cycles(maxlen=7):
    """All simple cycles of length 4..maxlen as 40-bit edge masks (each once)."""
    out = {}
    def dfs(path, mask):
        u = path[-1]
        for v in NB[u]:
            e = 1 << EI[tuple(sorted((u, v)))]
            if v == path[0] and len(path) >= 3:
                if path[1] < u:  # orientation dedupe
                    out[mask | e] = len(path)
            elif v > path[0] and v not in path and len(path) < maxlen:
                dfs(path + [v], mask | e)
    for s in VERTS: dfs([s], 0)
    return out
def chord_coords(mask):
    return sum(1 << j for j, k in enumerate(CHORDS) if (mask >> k) & 1)
def cover(vol):
    """vol: length-40 ints in [0,32). Returns adjacency list of the 512-vertex derived graph, vertex 32*i+t."""
    adj = [[] for _ in range(512)]
    for k, (a, b) in enumerate(EDGES):
        ia, ib, w = VI[a], VI[b], vol[k]
        for t in range(32):
            adj[32*ia+t].append(32*ib+(t ^ w)); adj[32*ib+(t ^ w)].append(32*ia+t)
    return adj
def bfs(adj, s, n=512):
    d = np.full(n, 99, np.int16); d[s] = 0; q = deque([s])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if d[v] == 99: d[v] = d[u]+1; q.append(v)
    return d
def fibre_dist(adj):
    """F[i, j, t] = distance from (i,0) to (j,t). Valid for abelian deck XOR action."""
    return np.array([bfs(adj, 32*i).reshape(16, 32) for i in range(16)])
def dist(F, i, s, j, t): return F[i, j, s ^ t]
def mate_counts(vol, F):
    """Number of eligible mates (same category, delta>=6) of a lift of each base edge."""
    out = []
    for k, (a, b) in enumerate(EDGES):
        ia, ib = VI[a], VI[b]; cnt = 0
        for k2, (c, d) in enumerate(EDGES):
            if CAT[k2] != CAT[k]: continue
            ic, id_ = VI[c], VI[d]
            for t in range(32):
                if k2 == k and t == 0: continue
                m = min(F[ia, ic, t], F[ia, id_, t ^ vol[k2]], F[ib, ic, t ^ vol[k]], F[ib, id_, t ^ vol[k] ^ vol[k2]])
                cnt += m >= 6
        out.append(cnt)
    return out
