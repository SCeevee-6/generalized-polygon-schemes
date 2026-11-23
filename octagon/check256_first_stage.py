"""Independent rerun of the first stage of the regular 256-component exclusion.
All 2^25 switching classes of signings of the Clebsch graph (chord bits on a BFS tree), kept if every eigenvalue of
the signed adjacency lies in [-4,-2] or [-1, 69/20]; orbits under the 1920 base automorphisms; then the edge-tree
capacity test (each vertex of the 32-vertex double cover receives at most 8 of the 170 tree visits)."""
import sys, os, time, numpy as np, numba as nb
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *
import elem_setup as S
ci = {k: i for i, k in enumerate(CHORDS)}
ea = np.array([VI[a] for a, b in EDGES]); eb = np.array([VI[b] for a, b in EDGES])
def eigs(xs):
    H = np.zeros((len(xs), 16, 16))
    for k in range(40):
        s = 1 - 2*((xs >> ci[k]) & 1) if k in ci else np.ones(len(xs))
        H[:, ea[k], eb[k]] = s; H[:, eb[k], ea[k]] = s
    return np.linalg.eigvalsh(H)
mats = []
for g in S.AUT:
    h = S.ginv(g); cols = []
    for j in range(25):
        m = S.FUND[j]; cols.append(chord_coords(sum(1 << S.eimg(h, kk) for kk in range(40) if m >> kk & 1)))
    mats.append(cols)
M = np.array(mats, np.int64)
@nb.njit
def canon(xs, M):
    out = xs.copy()
    for g in range(M.shape[0]):
        for t in range(xs.shape[0]):
            x = xs[t]; y = 0
            for j in range(25):
                v = x & M[g, j]
                v ^= v >> 16; v ^= v >> 8; v ^= v >> 4; v ^= v >> 2; v ^= v >> 1
                y |= (v & 1) << j
            if y < out[t]: out[t] = y
    return out
def rank2(cols):
    piv = {}; r = 0
    for v in cols:
        while v:
            h = v.bit_length()-1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; r += 1; break
    return r
burn = sum(2**(25 - rank2([sum(((int(Mg[j]) >> i) & 1) << j for j in range(25)) ^ (1 << i) for i in range(25)])) for Mg in M) // len(M)
print('orbits of all signings (Burnside)', burn, flush=True)
eps = 1e-9; keep = []
for st in range(0, 1 << 25, 1 << 19):
    xs = np.arange(st, st + (1 << 19), dtype=np.int64); ev = eigs(xs)
    ok = (((ev >= -4-eps) & (ev <= -2+eps)) | ((ev >= -1-eps) & (ev <= 69/20+eps))).all(1)
    keep.append(xs[ok])
keep = np.concatenate(keep); can = canon(keep, M); reps = np.unique(can); size = Counter(can.tolist())
print('spectral: signings', len(keep), 'orbits', len(reps), flush=True)
def cap_ok(x):
    s = np.array([1 - 2*((int(x) >> ci[k]) & 1) if k in ci else 1 for k in range(40)])
    for k, (a, b) in enumerate(EDGES):
        load = Counter()
        def walk(v, prev, sh, d):
            load[(v, sh)] += 1
            if d == 3: return
            for w in NB[v]:
                if w != prev: walk(w, v, sh*s[EI[tuple(sorted((v, w)))]], d+1)
        walk(a, b, 1, 0); walk(b, a, s[k], 0)
        if max(load.values()) > 8: return False
    return True
good = [int(r) for r in reps if cap_ok(r)]
print('capacity: orbits', len(good), 'signings', sum(size[r] for r in good), flush=True)
