"""Fast numba kernel (v2): per-node base precomputation, XOR-translated kill masks, lex-leader symmetry pruning."""
import numpy as np, numba as nb
MJ = np.array([0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF, 0x0000FFFF], np.int64)
@nb.njit(cache=True, inline='always')
def popc(x):
    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    return ((x * 0x01010101) & 0xFFFFFFFF) >> 24
@nb.njit(cache=True, inline='always')
def xortrans(m, c, MJ):
    for j in range(5):
        if (c >> j) & 1:
            s = 1 << j; M = MJ[j]
            m = ((m & M) << s) | ((m >> s) & M)
    return m
@nb.njit(cache=True)
def sval_below(y, s, lev):
    # XOR of y over bits of s strictly below lev
    acc = 0; i = 0
    s = s & ((1 << lev) - 1)
    while s:
        if s & 1: acc ^= y[i]
        s >>= 1; i += 1
    return acc
@nb.njit(cache=True)
def cap_ok(al, edgepairs):
    for c in range(10):
        tot = 0
        for e in range(4):
            mc = popc(al[edgepairs[c, e, 0]]) + popc(al[edgepairs[c, e, 1]]) + popc(al[edgepairs[c, e, 2]]) + popc(al[edgepairs[c, e, 3]])
            if mc == 0: tot += 12
            elif mc == 1: tot += 6
            elif mc == 2: tot += 4
            else: tot += 3
        if tot > 24: return False
    return True
@nb.njit(cache=True)
def sym_step(y, lev, gsupp, gtop, st_done, st_rank, st_piv, st_pivc, st_dead, ng, ncol):
    """Advance every live symmetry g as far as its columns are determined (support top <= lev).
    Returns -1 if some g gives a lexicographically smaller canonical sequence (prune), else 0. Updates states in place."""
    for g in range(ng):
        if st_dead[g]: continue
        p = st_done[g]
        while p <= lev and gtop[g, p] <= lev:
            # value of (g.q) at position p
            s = gsupp[g, p]; v = 0; i = 0
            while s:
                if s & 1: v ^= y[i]
                s >>= 1; i += 1
            # canonicalise against g's echelon basis
            vv = v; cc = 0
            for h in range(4, -1, -1):
                if (vv >> h) & 1 and st_piv[g, h] != 0:
                    vv ^= st_piv[g, h]; cc ^= st_pivc[g, h]
            if vv != 0:
                r = st_rank[g]; cv = 1 << r
                h = 0
                for hh in range(4, -1, -1):
                    if (vv >> hh) & 1:
                        h = hh; break
                st_piv[g, h] = vv; st_pivc[g, h] = cc ^ (1 << r); st_rank[g] = r + 1
            else:
                cv = cc
            if cv < y[p]: return -1
            if cv > y[p]:
                st_dead[g] = True; break
            p += 1
        st_done[g] = p
    return 0
@nb.njit(cache=True)
def search(y0, r0, lev0, al0, gstart, gmask, pstart, ppair, pmask, edgepairs, ncol, leaves, counts,
           gsupp, gtop, ng, MJ, maxnodes, stop, s2, ord_chord, eA, eB, pair_k, pair_k2, pair_cat, edge_slot, cel, stat):
    npairs = al0.shape[0]
    y = y0.copy(); r = np.zeros(ncol+1, np.int64); r[lev0] = r0
    al = np.zeros((ncol+1, npairs), np.int64); al[lev0] = al0
    idx = np.zeros(ncol+1, np.int64)
    forb = np.zeros(ncol+1, np.int64)
    bm = np.zeros((ncol+1, npairs), np.int64)
    touched = np.zeros((ncol+1, npairs), np.int64); ntouch = np.zeros(ncol+1, np.int64)
    # symmetry states per level
    S_done = np.zeros((ncol+1, ng), np.int64); S_rank = np.zeros((ncol+1, ng), np.int64)
    S_piv = np.zeros((ncol+1, ng, 5), np.int64); S_pivc = np.zeros((ncol+1, ng, 5), np.int64)
    S_dead = np.zeros((ncol+1, ng), np.bool_)
    nleaf = 0; nodes = 0
    T = np.zeros((3, 512), np.int64)
    assert ncol <= 27
    lev = lev0
    # initialise symmetry state on the prefix (levels < lev0) - caller passes prefix; recompute from scratch
    if ng > 0 and lev0 > 0:
        if sym_step(y, lev0-1, gsupp, gtop, S_done[lev0], S_rank[lev0], S_piv[lev0], S_pivc[lev0], S_dead[lev0], ng, ncol) < 0:
            return -1, 0
    enter = True
    while lev >= lev0:
        if enter:
            # per-node precomputation for level lev
            # chunk XOR tables over positions < lev
            for k in range(3):
                T[k, 0] = 0
                for x in range(1, 512):
                    lb = x & (-x); b = 0
                    while (lb >> b) != 1: b += 1
                    pos_ = 9*k + b
                    T[k, x] = T[k, x & (x-1)] ^ (y[pos_] if pos_ < lev else 0)
            f = 0
            for q in range(gstart[lev], gstart[lev+1]):
                s_ = gmask[q]
                f |= 1 << (T[0, s_ & 511] ^ T[1, (s_ >> 9) & 511] ^ T[2, (s_ >> 18) & 511])
            forb[lev] = f
            nt = 0; last = -1
            for q in range(pstart[lev], pstart[lev+1]):
                p = ppair[q]
                if p != last:
                    touched[lev, nt] = p; nt += 1; bm[lev, p] = 0; last = p
                s_ = pmask[q]; bm[lev, p] |= 1 << (T[0, s_ & 511] ^ T[1, (s_ >> 9) & 511] ^ T[2, (s_ >> 18) & 511])
            ntouch[lev] = nt
            idx[lev] = 0; enter = False
        rr = r[lev]; left = ncol - lev
        forced = (5 - rr) == left
        nc = 1 if forced else ((1 << rr) - 1 + (1 if rr < 5 else 0))
        i = idx[lev]
        if i >= nc:
            lev -= 1; continue
        idx[lev] = i + 1
        if forced or i == (1 << rr) - 1:
            c = 1 << rr; new = True
        else:
            c = i + 1; new = False
        if (forb[lev] >> c) & 1: continue
        y[lev] = c
        for p in range(npairs): al[lev+1, p] = al[lev, p]
        for t in range(ntouch[lev]):
            p = touched[lev, t]
            al[lev+1, p] &= ~xortrans(bm[lev, p], c, MJ)
        if not cap_ok(al[lev+1], edgepairs): continue
        if ng > 0:
            for g in range(ng):
                S_done[lev+1, g] = S_done[lev, g]; S_rank[lev+1, g] = S_rank[lev, g]; S_dead[lev+1, g] = S_dead[lev, g]
                for h in range(5):
                    S_piv[lev+1, g, h] = S_piv[lev, g, h]; S_pivc[lev+1, g, h] = S_pivc[lev, g, h]
            if sym_step(y, lev, gsupp, gtop, S_done[lev+1], S_rank[lev+1], S_piv[lev+1], S_pivc[lev+1], S_dead[lev+1], ng, ncol) < 0:
                continue
        counts[lev+1] += 1; nodes += 1
        if nodes > maxnodes: return -2, nleaf
        r[lev+1] = rr + (1 if new else 0)
        if lev + 1 == stop:
            if s2:
                st2 = stage2(y, al[lev+1], ord_chord, eA, eB, pair_k, pair_k2, pair_cat, edge_slot, cel)
                stat[st2] += 1
                if st2 == 0 or st2 == 3: continue
            if nleaf < leaves.shape[0]:
                for j in range(ncol): leaves[nleaf, j] = y[j]
                for p in range(npairs): leaves[nleaf, ncol+p] = al[lev+1, p]
            nleaf += 1
            continue
        lev += 1; enter = True
    return 0, nleaf

@nb.njit(cache=True)
def stage2(y, al, ord_chord, eA, eB, pair_k, pair_k2, pair_cat, edge_slot, cel):
    """Leaf filter. Returns 3 = disconnected (chord voltages rank<5), 0 = rejected by H-sibling test,
    1 = survivor (every category tested), 2 = survivor with some untested category."""
    ncol = y.shape[0]
    vol = np.zeros(40, np.int64)
    piv = np.zeros(5, np.int64)
    for p in range(ncol):
        k = ord_chord[p]
        if k >= 0:
            vol[k] = y[p]
            v = y[p]
            for h in range(4, -1, -1):
                if (v >> h) & 1 and piv[h] != 0: v ^= piv[h]
            if v != 0:
                for h in range(4, -1, -1):
                    if (v >> h) & 1:
                        piv[h] = v; break
    rk = 0
    for h in range(5):
        if piv[h] != 0: rk += 1
    if rk < 5: return 3
    # mate graphs per category on 128 lifts: index slot*32+s ; adjacency bitsets (2 words of 64)
    madj = np.zeros((10, 128, 2), np.uint64)
    for pi in range(al.shape[0]):
        m = al[pi]
        if m == 0: continue
        c = pair_cat[pi]; j1 = edge_slot[pair_k[pi]]; j2 = edge_slot[pair_k2[pi]]
        for t in range(32):
            if (m >> t) & 1:
                for s in range(32):
                    u = j1*32 + s; w = j2*32 + (s ^ t)
                    if u == w: continue
                    madj[c, u, w >> 6] |= np.uint64(1) << np.uint64(w & 63)
                    madj[c, w, u >> 6] |= np.uint64(1) << np.uint64(u & 63)
    # per category: triangle-free? forced pairs = vertices of degree 1 whose neighbour has degree 1
    forced_of = -np.ones((10, 128), np.int64)   # partner if forced
    tested = np.zeros(10, np.bool_)
    for c in range(10):
        deg = np.zeros(128, np.int64)
        for u in range(128):
            deg[u] = popc64(madj[c, u, 0]) + popc64(madj[c, u, 1])
        tri = False
        for u in range(128):
            for w in range(u+1, 128):
                if (madj[c, u, w >> 6] >> np.uint64(w & 63)) & np.uint64(1):
                    if (madj[c, u, 0] & madj[c, w, 0]) != 0 or (madj[c, u, 1] & madj[c, w, 1]) != 0:
                        tri = True; break
            if tri: break
        if tri: continue
        nf = 0
        for u in range(128):
            if deg[u] == 1:
                w = -1
                for x in range(128):
                    if (madj[c, u, x >> 6] >> np.uint64(x & 63)) & np.uint64(1):
                        w = x; break
                if deg[w] == 1:
                    forced_of[c, u] = w; nf += 1
        if nf > 0: tested[c] = True
    # H partial: points 0..511, lines 512 + 32k + s, boundary nodes 1792 + id
    bid = -np.ones(1280, np.int64); nb_ = 0
    cat_of_b = np.zeros(640, np.int64); bl1 = np.zeros(640, np.int64); bl2 = np.zeros(640, np.int64)
    for c in range(10):
        if not tested[c]: continue
        for u in range(128):
            w = forced_of[c, u]
            if w > u:
                # global line ids
                lu = cat_lines(c, u, cel); lw = cat_lines(c, w, cel)
                bid[lu] = nb_; bid[lw] = nb_; cat_of_b[nb_] = c; bl1[nb_] = lu; bl2[nb_] = lw; nb_ += 1
    N = 1792 + nb_
    # point -> lines
    pl = np.zeros((512, 5), np.int64); pc = np.zeros(512, np.int64)
    lp = np.zeros((1280, 2), np.int64)
    for k in range(40):
        for s in range(32):
            a = eA[k]*32 + s; b = eB[k]*32 + (s ^ vol[k]); L = 32*k + s
            lp[L, 0] = a; lp[L, 1] = b
            pl[a, pc[a]] = L; pc[a] += 1; pl[b, pc[b]] = L; pc[b] += 1
    dist = np.empty(N, np.int64); q = np.empty(N, np.int64)
    for b0 in range(nb_):
        c = cat_of_b[b0]
        for i in range(N): dist[i] = 999
        dist[1792 + b0] = 0; q[0] = 1792 + b0; qh = 0; qt = 1
        while qh < qt:
            x = q[qh]; qh += 1; dx = dist[x]
            if dx >= 14: break
            if x >= 1792:
                bb = x - 1792
                for L in (bl1[bb], bl2[bb]):
                    if dist[512+L] > dx+1: dist[512+L] = dx+1; q[qt] = 512+L; qt += 1
            elif x >= 512:
                L = x - 512
                for j in range(2):
                    pnt = lp[L, j]
                    if dist[pnt] > dx+1: dist[pnt] = dx+1; q[qt] = pnt; qt += 1
                if bid[L] >= 0:
                    z = 1792 + bid[L]
                    if dist[z] > dx+1: dist[z] = dx+1; q[qt] = z; qt += 1
            else:
                for j in range(5):
                    z = 512 + pl[x, j]
                    if dist[z] > dx+1: dist[z] = dx+1; q[qt] = z; qt += 1
        cand = 0
        for b2 in range(nb_):
            if b2 != b0 and cat_of_b[b2] == c and dist[1792+b2] >= 14: cand += 1
        if cand == 0:
            # option pairs: non-forced mate edges in category c
            for u in range(128):
                if forced_of[c, u] >= 0: continue
                for w in range(u+1, 128):
                    if (madj[c, u, w >> 6] >> np.uint64(w & 63)) & np.uint64(1):
                        lu = cat_lines(c, u, cel); lw = cat_lines(c, w, cel)
                        dmin = dist[512+lu] if dist[512+lu] < dist[512+lw] else dist[512+lw]
                        if dmin + 1 >= 14: cand += 1
                if cand > 0: break
        if cand == 0: return 0
    for c in range(10):
        if not tested[c]: return 2
    return 1
@nb.njit(cache=True, inline='always')
def popc64(x):
    c = 0
    while x:
        x &= x - np.uint64(1); c += 1
    return c
@nb.njit(cache=True, inline='always')
def cat_lines(c, u, cel):
    return 32*cel[c, u // 32] + (u % 32)
