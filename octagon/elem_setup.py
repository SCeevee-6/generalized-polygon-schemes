"""Setup for the elementary 32-cover search: constraint lists, column order, arrays, helpers."""
import sys, json, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from base import *
import elem_search as es
def walks_upto(u, L=5):
    out = []
    def rec(v, prev, x, l):
        out.append((v, x))
        if l == L: return
        for w in NB[v]:
            if w == prev: continue
            rec(w, v, x ^ chord_coords(1 << EI[tuple(sorted((v, w)))]), l+1)
    rec(u, None, 0, 0); return out
WALKS = {u: walks_upto(u, 5) for u in VERTS}
def ch(k): return chord_coords(1 << k)
CYC = simple_cycles(7)
GIRTH = {chord_coords(m) for m in CYC}
def mate_cons(k, k2):
    """chord-coordinate vectors z with: lift of k at sheet 0 and lift of k2 at sheet t have delta<=5 iff t = q(z) for some z."""
    a, b = EDGES[k]; c, d = EDGES[k2]; cons = set()
    for src, add in ((a, 0), (b, ch(k))):
        for w, x in WALKS[src]:
            if w == c: cons.add(x ^ add)
            elif w == d: cons.add(x ^ add ^ ch(k2))
    return cons
CATS = sorted(set(CAT)); CATEDGES = {c: [k for k in range(40) if CAT[k] == c] for c in CATS}
PAIRS = []
for ci, c in enumerate(CATS):
    es_ = CATEDGES[c]
    for i in range(4):
        for j in range(i, 4): PAIRS.append((ci, es_[i], es_[j], 16 if i == j else 32))
PAIR_EDGES = [(k, k2) for _, k, k2, _ in PAIRS]
PAIR_CONS = [mate_cons(k, k2) for _, k, k2, _ in PAIRS]
def greedy_order(cons, ncol=25):
    order = []; chosen = 0; rem = set(range(ncol))
    while rem:
        best = None
        for i in sorted(rem):
            cm = chosen | (1 << i); sc = sum(1 for s in cons if s & ~cm == 0 and s >> i & 1)
            if best is None or sc > best[0]: best = (sc, i)
        order.append(best[1]); chosen |= 1 << best[1]; rem.remove(best[1])
    return order
def setup(fp):
    hard = {(z << 1) | 1 for z in PAIR_CONS[fp] if z != 0}
    allc = set(s << 1 for s in GIRTH) | hard
    for pc in PAIR_CONS: allc |= {z << 1 for z in pc}
    order = [0] + [1+i for i in greedy_order({s >> 1 for s in allc if s >> 1})]
    pos = {b: j for j, b in enumerate(order)}
    rm = lambda s: sum(1 << pos[i] for i in range(26) if s >> i & 1)
    gtr = [[] for _ in range(26)]
    for s in set(s << 1 for s in GIRTH) | hard: gtr[rm(s).bit_length()-1].append(rm(s))
    ptr = [[] for _ in range(26)]; zk = [0]*len(PAIRS)
    for pi, pc in enumerate(PAIR_CONS):
        for z in pc:
            if z == 0: zk[pi] |= 1
            else: ptr[rm(z << 1).bit_length()-1].append((pi, rm(z << 1)))
    return gtr, ptr, zk, order
def build_arrays(fp):
    gtr, ptr, zk, order = setup(fp)
    gstart = np.zeros(27, np.int64); gm = []
    for l in range(26): gm += sorted(gtr[l]); gstart[l+1] = len(gm)
    pstart = np.zeros(27, np.int64); pp = []; pmk = []
    for l in range(26):
        for pi, sm in sorted(ptr[l]): pp.append(pi); pmk.append(sm)
        pstart[l+1] = len(pp)
    ep = np.zeros((10, 4, 4), np.int64)
    for ci, c in enumerate(CATS):
        for e_i, k in enumerate(CATEDGES[c]):
            ep[ci, e_i] = [pi for pi in range(ci*10, ci*10+10) if k in PAIR_EDGES[pi]]
    al0 = np.array([((1 << 32)-1) & ~zk[pi] for pi in range(len(PAIRS))], np.int64)
    return dict(gstart=gstart, gmask=np.array(gm, np.int64), pstart=pstart, ppair=np.array(pp, np.int64),
                pmask=np.array(pmk, np.int64), edgepairs=ep, al0=al0, order=np.array(order))
def canon(seq):
    """GL(5,2) first-occurrence echelon form of a sequence of vectors."""
    piv = {}; nb_ = 0; out = []
    for v in seq:
        vv, cc = v, 0
        while vv:
            h = vv.bit_length()-1
            if h not in piv: break
            pv, pc = piv[h]; vv ^= pv; cc ^= pc
        if vv: piv[vv.bit_length()-1] = (vv, cc ^ (1 << nb_)); out.append(1 << nb_); nb_ += 1
        else: out.append(cc)
    return out
def to_seq(vol, tau, order):
    cols = [tau] + [vol[k] for k in CHORDS]
    return canon([cols[o] for o in order])
def from_seq(seq, order):
    """inverse: chord voltages (gauge: tree edges 0) and tau from a column sequence."""
    cols = [0]*26
    for j, o in enumerate(order): cols[o] = int(seq[j])
    vol = [0]*40
    for i, k in enumerate(CHORDS): vol[k] = cols[1+i]
    return vol, cols[0]
def prefix_state(seq, lev0, A):
    y = np.zeros(26, np.int64); r = 0; al = A['al0'].copy()
    for lev in range(lev0):
        c = int(seq[lev]); y[lev] = c
        for q in range(A['gstart'][lev], A['gstart'][lev+1]):
            if es.sval(y, A['gmask'][q]) == 0: return None, ('girth', lev)
        for q in range(A['pstart'][lev], A['pstart'][lev+1]):
            al[A['ppair'][q]] &= ~(1 << es.sval(y, A['pmask'][q]))
        if not es.cap_ok(al, A['edgepairs']): return None, ('cap', lev)
        if c == 1 << r: r += 1
    return (y, r, al), None
def run_from(seq, lev0, A, maxleaves=100000, stop=26):
    assert lev0 < 26
    st, why = prefix_state(seq, lev0, A)
    if st is None: return why
    y, r, al = st
    leaves = np.zeros((maxleaves, 26+1+len(PAIRS)), np.int64); counts = np.zeros(27, np.int64)
    n = es.subtree(y, r, lev0, al, A['gstart'], A['gmask'], A['pstart'], A['ppair'], A['pmask'], A['edgepairs'], 26, leaves, 0, counts, stop)
    return n, leaves[:min(n, maxleaves)], counts

def setup_multi(forced):
    """forced: list of pair indices; column j<m is tau_j (translate of forced orbit j); columns m.. are chords."""
    m = len(forced); ncol = 25 + m
    hard = set()
    for j, fp in enumerate(forced):
        hard |= {(z << m) | (1 << j) for z in PAIR_CONS[fp] if z != 0}
    allc = set(s << m for s in GIRTH) | hard
    for pc in PAIR_CONS: allc |= {z << m for z in pc}
    chord_order = greedy_order({s >> m for s in allc if s >> m})
    order = list(range(m)) + [m+i for i in chord_order]
    pos = {b: j for j, b in enumerate(order)}
    rm = lambda s: sum(1 << pos[i] for i in range(ncol) if s >> i & 1)
    gtr = [[] for _ in range(ncol)]
    for s in set(s << m for s in GIRTH) | hard:
        if s: gtr[rm(s).bit_length()-1].append(rm(s))
    ptr = [[] for _ in range(ncol)]; zk = [0]*len(PAIRS)
    for pi, pc in enumerate(PAIR_CONS):
        for z in pc:
            if z == 0: zk[pi] |= 1
            else: ptr[rm(z << m).bit_length()-1].append((pi, rm(z << m)))
    gstart = np.zeros(ncol+1, np.int64); gm = []
    for l in range(ncol): gm += sorted(gtr[l]); gstart[l+1] = len(gm)
    pstart = np.zeros(ncol+1, np.int64); pp = []; pmk = []
    for l in range(ncol):
        for pi, sm in sorted(ptr[l]): pp.append(pi); pmk.append(sm)
        pstart[l+1] = len(pp)
    ep = np.zeros((10, 4, 4), np.int64)
    for ci, c in enumerate(CATS):
        for e_i, k in enumerate(CATEDGES[c]):
            ep[ci, e_i] = [pi for pi in range(ci*10, ci*10+10) if k in PAIR_EDGES[pi]]
    al0 = np.array([((1 << 32)-1) & ~zk[pi] for pi in range(len(PAIRS))], np.int64)
    return dict(gstart=gstart, gmask=np.array(gm, np.int64), pstart=pstart, ppair=np.array(pp, np.int64),
                pmask=np.array(pmk, np.int64), edgepairs=ep, al0=al0, order=np.array(order), ncol=ncol, m=m, forced=list(forced))
def knuth_multi(A, trials=100, seed=0):
    import random
    rnd = random.Random(seed); ncol = A['ncol']; prof = np.zeros(ncol+1)
    for _ in range(trials):
        y = np.zeros(ncol, np.int64); r = 0; w = 1.0; prof[0] += 1; al = A['al0'].copy()
        for lev in range(ncol):
            left = ncol - lev
            if 5 - r == left: cands = [1 << r]
            else: cands = list(range(1, 1 << r)) + ([1 << r] if r < 5 else [])
            good = []
            for c in cands:
                y[lev] = c
                if any(es.sval(y, A['gmask'][q]) == 0 for q in range(A['gstart'][lev], A['gstart'][lev+1])): continue
                a2 = al.copy()
                for q in range(A['pstart'][lev], A['pstart'][lev+1]):
                    a2[A['ppair'][q]] &= ~(1 << es.sval(y, A['pmask'][q]))
                if not es.cap_ok(a2, A['edgepairs']): continue
                good.append((c, a2))
            if not good: break
            w *= len(good); prof[lev+1] += w
            c, al = rnd.choice(good); y[lev] = c
            if c == 1 << r: r += 1
    return prof/trials

def setup_multi2(forced, order_mode='hard'):
    m = len(forced); ncol = 25 + m
    hard = set()
    for j, fp in enumerate(forced):
        hard |= {(z << m) | (1 << j) for z in PAIR_CONS[fp] if z != 0}
    base = set(s << m for s in GIRTH) | hard
    if order_mode == 'hard': oc = {s >> m for s in base if s >> m}
    else:
        oc = set(base)
        for pc in PAIR_CONS: oc |= {z << m for z in pc}
        oc = {s >> m for s in oc if s >> m}
    chord_order = greedy_order(oc)
    return _finish(forced, m, ncol, base, chord_order)
def _finish(forced, m, ncol, base, chord_order):
    order = list(range(m)) + [m+i for i in chord_order]
    pos = {b: j for j, b in enumerate(order)}
    rm = lambda s: sum(1 << pos[i] for i in range(ncol) if s >> i & 1)
    gtr = [[] for _ in range(ncol)]
    for s in base:
        if s: gtr[rm(s).bit_length()-1].append(rm(s))
    ptr = [[] for _ in range(ncol)]; zk = [0]*len(PAIRS)
    for pi, pc in enumerate(PAIR_CONS):
        for z in pc:
            if z == 0: zk[pi] |= 1
            else: ptr[rm(z << m).bit_length()-1].append((pi, rm(z << m)))
    gstart = np.zeros(ncol+1, np.int64); gm = []
    for l in range(ncol): gm += sorted(gtr[l]); gstart[l+1] = len(gm)
    pstart = np.zeros(ncol+1, np.int64); pp = []; pmk = []
    for l in range(ncol):
        for pi, sm in sorted(ptr[l]): pp.append(pi); pmk.append(sm)
        pstart[l+1] = len(pp)
    ep = np.zeros((10, 4, 4), np.int64)
    for ci, c in enumerate(CATS):
        for e_i, k in enumerate(CATEDGES[c]):
            ep[ci, e_i] = [pi for pi in range(ci*10, ci*10+10) if k in PAIR_EDGES[pi]]
    al0 = np.array([((1 << 32)-1) & ~zk[pi] for pi in range(len(PAIRS))], np.int64)
    return dict(gstart=gstart, gmask=np.array(gm, np.int64), pstart=pstart, ppair=np.array(pp, np.int64),
                pmask=np.array(pmk, np.int64), edgepairs=ep, al0=al0, order=np.array(order), ncol=ncol, m=m, forced=list(forced))

import elem_fast as ef, itertools
def apply_perm(u, p): return sum(((u >> i) & 1) << p[i] for i in range(5))
AUT = [tuple(apply_perm(u, p) ^ tr for u in range(32)) for p in itertools.permutations(range(5)) for tr in VERTS]
def eimg(g, k):
    a, b = EDGES[k]; return EI[tuple(sorted((g[a], g[b])))]
def ginv(g):
    h = [0]*32
    for u in range(32): h[g[u]] = u
    return tuple(h)
def tree_path(u):
    p = []
    while PAR[u] is not None: p.append(EI[tuple(sorted((u, PAR[u])))]); u = PAR[u]
    return p
FUND = []
for k in CHORDS:
    a, b = EDGES[k]; m_ = 1 << k
    for e in tree_path(a) + tree_path(b): m_ ^= 1 << e
    FUND.append(m_)
def pair_stabiliser(fp):
    k, k2 = PAIR_EDGES[fp]
    return [g for g in AUT if {eimg(g, k), eimg(g, k2)} == {k, k2}]
def mate_shift(fp, h):
    """Chord mask correcting a marked translate under pullback by h and tree switching."""
    mask = 0
    for k in PAIR_EDGES[fp]:
        start = EDGES[k][0]
        for e in tree_path(start):
            mask ^= 1 << eimg(h, e)
        # The marked sheet is measured at the smaller endpoint of an oriented edge.
        image = eimg(h, k)
        if h[start] != EDGES[image][0]:
            mask ^= 1 << image
    return chord_coords(mask)

def sym_arrays(A, group):
    """gsupp[g,p]: support (search positions) of column at position p of g.q; identity removed."""
    m = A['m']; ncol = A['ncol']; order = list(A['order']); pos = {b: j for j, b in enumerate(order)}
    ident = tuple(range(32)); G = [g for g in group if g != ident]
    gs = np.zeros((len(G), ncol), np.int64); gt = np.zeros((len(G), ncol), np.int64)
    for gi, g in enumerate(G):
        h = ginv(g)
        for p, o in enumerate(order):
            if o < m:
                shift = mate_shift(A['forced'][o], h)
                sup = (1 << pos[o]) | sum(1 << pos[m+i] for i in range(25) if shift >> i & 1)
            else:
                mm = FUND[o-m]; img = sum(1 << eimg(h, kk) for kk in range(40) if mm >> kk & 1)
                cc = chord_coords(img); sup = sum(1 << pos[m+i] for i in range(25) if cc >> i & 1)
            gs[gi, p] = sup; gt[gi, p] = sup.bit_length()-1
    return gs, gt
def run_fast(A, seq=None, lev0=1, use_sym=True, group=None, maxleaves=200000, maxnodes=10**15, stop=None):
    ncol = A['ncol']; stop = ncol if stop is None else stop
    if seq is None: seq = [1] + [0]*(ncol-1)
    st, why = prefix_state_m(seq, lev0, A)
    if st is None: return why
    y, r, al = st
    if use_sym: gs, gt = sym_arrays(A, group)
    else: gs = np.zeros((0, ncol), np.int64); gt = np.zeros((0, ncol), np.int64)
    leaves = np.zeros((maxleaves, ncol+len(PAIRS)), np.int64); counts = np.zeros(ncol+1, np.int64)
    code, n = ef.search(y, r, lev0, al, A['gstart'], A['gmask'], A['pstart'], A['ppair'], A['pmask'], A['edgepairs'],
                        ncol, leaves, counts, gs, gt, len(gs), ef.MJ, maxnodes, stop)
    return code, n, leaves[:min(n, maxleaves)], counts
def prefix_state_m(seq, lev0, A):
    ncol = A['ncol']; y = np.zeros(ncol, np.int64); r = 0; al = A['al0'].copy()
    for lev in range(lev0):
        c = int(seq[lev]); y[lev] = c
        for q in range(A['gstart'][lev], A['gstart'][lev+1]):
            if es.sval(y, A['gmask'][q]) == 0: return None, ('girth', lev)
        for q in range(A['pstart'][lev], A['pstart'][lev+1]):
            al[A['ppair'][q]] &= ~(1 << es.sval(y, A['pmask'][q]))
        if not es.cap_ok(al, A['edgepairs']): return None, ('cap', lev)
        if c == 1 << r: r += 1
    return (y, r, al), None

def s2_arrays(A):
    m = A['m']; order = list(A['order'])
    ord_chord = np.array([CHORDS[o-m] if o >= m else -1 for o in order], np.int64)
    eA = np.array([VI[a] for a, b in EDGES], np.int64); eB = np.array([VI[b] for a, b in EDGES], np.int64)
    pk = np.array([k for k, k2 in PAIR_EDGES], np.int64); pk2 = np.array([k2 for k, k2 in PAIR_EDGES], np.int64)
    pc = np.array([ci for ci, _, _, _ in PAIRS], np.int64)
    slot = np.zeros(40, np.int64); cel = np.zeros((10, 4), np.int64)
    for ci, c in enumerate(CATS):
        for j, k in enumerate(CATEDGES[c]): slot[k] = j; cel[ci, j] = k
    return ord_chord, eA, eB, pk, pk2, pc, slot, cel
def run_fast2(A, seq=None, lev0=1, use_sym=True, group=None, maxleaves=100000, maxnodes=10**15, stop=None, s2=True):
    ncol = A['ncol']; stop = ncol if stop is None else stop
    if seq is None: seq = [1] + [0]*(ncol-1)
    st, why = prefix_state_m(seq, lev0, A)
    if st is None: return why
    y, r, al = st
    if use_sym: gs, gt = sym_arrays(A, group)
    else: gs = np.zeros((0, ncol), np.int64); gt = np.zeros((0, ncol), np.int64)
    leaves = np.zeros((maxleaves, ncol+len(PAIRS)), np.int64); counts = np.zeros(ncol+1, np.int64)
    stat = np.zeros(4, np.int64)
    code, n = ef.search(y, r, lev0, al, A['gstart'], A['gmask'], A['pstart'], A['ppair'], A['pmask'], A['edgepairs'],
                        ncol, leaves, counts, gs, gt, len(gs), ef.MJ, maxnodes, stop, s2 and stop == ncol, *s2_arrays(A), stat)
    return code, n, leaves[:min(n, maxleaves)], counts, stat
