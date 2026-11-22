"""Enumerate all inclusion-maximal pairwise non-opposite sets of size >= MIN that contain vertex 0.

Bron-Kerbosch with Tomita pivoting and a greedy-colouring bound, on bitsets (numba). The search runs in the
graph on the vertices non-opposite to 0, with 'adjacent' = 'not opposite'; a maximal clique there plus 0 is
maximal in the whole graph. Nothing about lines is given to the search. Afterwards each set is tested against
all line neighbourhoods, and the sets that lie in none are classified and split into orbits of the stabiliser
of their central line (when they have one).

Usage: python triple/escape_bk.py --geom points --min 15 [--gens atlas|file.npy] [--out DIR]
"""
import argparse, json, os, time
from collections import Counter
import numpy as np
from numba import njit
from escape_common import build

M1 = np.uint64(0x5555555555555555); M2 = np.uint64(0x3333333333333333)
M4 = np.uint64(0x0F0F0F0F0F0F0F0F); H01 = np.uint64(0x0101010101010101)


@njit(cache=True)
def popc(w):
    w = w - ((w >> np.uint64(1)) & M1)
    w = (w & M2) + ((w >> np.uint64(2)) & M2)
    w = (w + (w >> np.uint64(4))) & M4
    return int((w * H01) >> np.uint64(56))


@njit(cache=True)
def count(S):
    c = 0
    for i in range(S.shape[0]):
        c += popc(S[i])
    return c


@njit(cache=True)
def first(S):
    for i in range(S.shape[0]):
        w = S[i]
        if w != 0:
            b = 0
            while (w >> np.uint64(b)) & np.uint64(1) == 0:
                b += 1
            return 64 * i + b
    return -1


@njit(cache=True)
def colour_bound(P, adj, need):
    U = P.copy(); Q = np.empty_like(P); col = 0
    while count(U) > 0 and col < need:
        col += 1
        Q[:] = U
        while True:
            v = first(Q)
            if v < 0:
                break
            wi, bi = v // 64, np.uint64(v % 64)
            U[wi] &= ~(np.uint64(1) << bi)
            Q &= ~adj[v]
            Q[wi] &= ~(np.uint64(1) << bi)
    return col


@njit(cache=True)
def enumerate_bk(adj, n, minc, cap, maxd):
    W = adj.shape[1]
    P = np.zeros((maxd + 1, W), np.uint64); X = np.zeros((maxd + 1, W), np.uint64)
    B = np.zeros((maxd + 1, W), np.uint64); R = np.zeros(maxd + 1, np.int64)
    sizes = np.zeros(maxd + 2, np.int64)
    out = -np.ones((cap, maxd + 1), np.int16); nout = 0; nodes = 0
    for v in range(n):
        P[0, v // 64] |= np.uint64(1) << np.uint64(v % 64)
    d = 0
    init = True
    while d >= 0:
        if init:
            init = False; nodes += 1
            B[d, :] = 0
            cp = count(P[d])
            if cp == 0:
                if count(X[d]) == 0 and d >= minc:
                    sizes[d] += 1
                    if nout < cap:
                        out[nout, :d] = R[:d]
                    nout += 1
            elif cp + d >= minc and colour_bound(P[d], adj, minc - d) >= minc - d:
                best = -1; piv = -1
                for wi in range(W):
                    w = P[d, wi] | X[d, wi]
                    while w != 0:
                        b = 0
                        while (w >> np.uint64(b)) & np.uint64(1) == 0:
                            b += 1
                        u = 64 * wi + b
                        w &= ~(np.uint64(1) << np.uint64(b))
                        sc = 0
                        for j in range(W):
                            sc += popc(adj[u, j] & P[d, j])
                        if sc > best:
                            best = sc; piv = u
                for j in range(W):
                    B[d, j] = P[d, j] & ~adj[piv, j]
        v = first(B[d])
        if v < 0 or count(P[d]) + d < minc:
            d -= 1
            continue
        wi, bi = v // 64, np.uint64(v % 64)
        B[d, wi] &= ~(np.uint64(1) << bi)
        R[d] = v
        for j in range(W):
            P[d + 1, j] = P[d, j] & adj[v, j]
            X[d + 1, j] = X[d, j] & adj[v, j]
        P[d, wi] &= ~(np.uint64(1) << bi)
        X[d, wi] |= np.uint64(1) << bi
        d += 1; init = True
    return sizes, out[:min(nout, cap)], nodes


def schreier_orbits(sets, seed, gens, act):
    """Orbits of Stab(seed) on `sets` (all with the same seed object), via Schreier generators on the orbit of seed.
    act(g, obj) maps an object (tuple) under permutation g."""
    tr = {seed: np.arange(len(gens[0]))}; seq = [seed]
    for x in seq:
        for g in gens:
            y = act(g, x)
            if y not in tr:
                tr[y] = g[tr[x]]; seq.append(y)
    inv = {x: np.argsort(p) for x, p in tr.items()}
    index = {x: i for i, x in enumerate(sets)}
    parent = list(range(len(sets)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    ident = np.arange(len(gens[0]))
    for x in seq:
        for g in gens:
            s = inv[act(g, x)][g[tr[x]]]
            if np.array_equal(s, ident):
                continue
            for i, S in enumerate(sets):
                a, b = find(i), find(index[tuple(sorted(s[list(S)]))])
                if a != b:
                    parent[a] = b
    return len(seq), sorted(Counter(find(i) for i in range(len(sets))).values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geom', default='points', choices=['points', 'lines'])
    ap.add_argument('--min', type=int, default=15)
    ap.add_argument('--gens', default='atlas', help="'atlas', 'gap:DIR' or a .npy file")
    ap.add_argument('--tag', default='')
    ap.add_argument('--cap', type=int, default=2000000)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_escape'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    D, nbhd, G = build(a.geom, a.gens)
    n_all = D.shape[0]
    cand = np.nonzero((D[0] > 0) & (D[0] < 4))[0]
    n = len(cand); W = (n + 63) // 64
    adj = np.zeros((n, W), np.uint64)
    sub = (D[np.ix_(cand, cand)] > 0) & (D[np.ix_(cand, cand)] < 4)
    for i in range(n):
        for j in np.nonzero(sub[i])[0]:
            adj[i, j // 64] |= np.uint64(1) << np.uint64(j % 64)
    K = max(len(N) for N in nbhd)
    print(f'{a.geom}: {n_all} vertices, {n} candidates through 0, {len(nbhd)} line neighbourhoods of size {K}', flush=True)
    sizes, cl, nodes = enumerate_bk(adj, n, a.min - 1, a.cap, K + 1)
    t1 = time.time()
    counts = {int(i) + 1: int(c) for i, c in enumerate(sizes) if c}
    print(f'maximal sets through 0 by size: {counts}; {nodes} search nodes; {t1 - t0:.0f} s', flush=True)
    sets = [tuple(sorted([0] + cand[r[r >= 0]].tolist())) for r in cl]
    in_nb = [any(set(S) <= N for N in nbhd) for S in sets]
    res = {'geom': a.geom, 'gens': a.gens if not a.gens.startswith('gap:') else 'gap:TF42G1-p1755B0', 'min': a.min, 'counts': counts, 'nodes': int(nodes), 'seconds': round(t1 - t0),
           'in_line_nbhd_by_size': dict(Counter(len(S) for S, f in zip(sets, in_nb) if f)),
           'escaping_by_size': dict(Counter(len(S) for S, f in zip(sets, in_nb) if not f))}
    esc = [S for S, f in zip(sets, in_nb) if not f]
    # structure of escaping sets: collinear pairs, full lines, central line, projection profile
    from collections import defaultdict
    prof = Counter(); by_centre = defaultdict(list); overlap = defaultdict(Counter)
    inc = np.zeros((len(nbhd), n_all), bool)
    for i, N in enumerate(nbhd):
        inc[i, list(N)] = True
    for S in esc:
        d = D[np.ix_(S, S)]
        deg = (d == 1).sum(axis=0)
        cen = tuple(np.array(S)[deg == deg.max()]) if deg.max() > 0 else ()
        cen = tuple(int(c) for c in cen)
        full = bool(len(cen) > 1 and (D[np.ix_(cen, cen)] <= 1).all())
        key = (tuple(sorted(Counter(deg.tolist()).items())), full)
        if full:
            rest = [x for x in S if x not in cen]
            dl = D[np.ix_(cen, rest)]
            key = key + (int(dl.min(axis=0).max()), tuple(sorted((dl == dl.min(axis=0)).sum(axis=1).tolist())))
            if 0 in cen:
                by_centre[cen].append(S)
        prof[key] += 1
        overlap[key][int(inc[:, list(S)].sum(axis=1).max())] += 1
    res['escaping_structure'] = {str(k): v for k, v in prof.items()}
    res['max_overlap_with_a_line_nbhd'] = {str(k): dict(v) for k, v in overlap.items()}
    print('escaping structure (degree multiset, central line full, max dist to it, projection counts):', flush=True)
    for k, v in prof.items():
        print('  ', k, v, 'largest overlap with a line neighbourhood:', dict(overlap[k]), flush=True)
    if by_centre:
        cen = sorted(by_centre)[0]
        fam = by_centre[cen]
        act = lambda g, x: tuple(sorted(g[list(x)].tolist()))
        nlines, orb = schreier_orbits(fam, cen, G, act)
        res['central_line'] = list(map(int, cen)); res['sets_on_central_line'] = len(fam)
        res['stabiliser_orbits_on_them'] = orb; res['lines_in_orbit'] = nlines
        res['global_orbit_sizes'] = [nlines * o for o in orb]
        print(f'central line {cen}: {len(fam)} escaping sets, Stab orbits {orb}, '
              f'{nlines} lines, global {[nlines * o for o in orb]}', flush=True)
    with open(os.path.join(a.out, f'escape_bk_{a.geom}_min{a.min}{a.tag}.json'), 'w', newline='\n') as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(a.out, f'escape_bk_{a.geom}_min{a.min}{a.tag}_sets.txt'), 'w', newline='\n') as f:
        for S, fl in zip(sets, in_nb):
            f.write(('N ' if fl else 'E ') + ' '.join(map(str, S)) + '\n')


if __name__ == '__main__':
    main()
