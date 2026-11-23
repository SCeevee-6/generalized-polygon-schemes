"""Exhaustive search for regular elementary abelian 32-sheeted trace covers (no lifted symmetry assumed).
Columns: 0 = tau (mate translate, fixed e1), 1..25 = chord voltages in a greedy order. GL(5,2) broken by
first-occurrence echelon form. Constraints:
  girth: every simple base cycle of length 4..7 has nonzero voltage;
  hard: the forced mate orbit (fixed pair, t=tau) has delta>=6 (all walks of length<=5 between endpoints avoid t);
  capacity: per category sum_eps 1/(min(mc_eps,3)+1) <= 2, mc counted over still-possible (pair,t) orbits.
"""
import numpy as np, numba as nb
@nb.njit(cache=True, inline='always')
def popc(x):
    c = 0
    while x:
        x &= x-1; c += 1
    return c
@nb.njit(cache=True)
def sval(y, s):
    acc = 0; i = 0
    while s:
        if s & 1: acc ^= y[i]
        s >>= 1; i += 1
    return acc
@nb.njit(cache=True)
def cap_ok(al, edgepairs):
    # edgepairs[c, e, j] pair index (4 pairs per edge); weights 12/(min(mc,3)+1)
    for c in range(10):
        tot = 0
        for e in range(4):
            mc = 0
            for j in range(4): mc += popc(al[edgepairs[c, e, j]])
            if mc == 0: tot += 12
            elif mc == 1: tot += 6
            elif mc == 2: tot += 4
            else: tot += 3
        if tot > 24: return False
    return True
@nb.njit(cache=True)
def subtree(y0, r0, lev0, al0, gstart, gmask, pstart, ppair, pmask, edgepairs, ncol, leaves, nleaf, counts, stop):
    """DFS below a prefix (levels < lev0 fixed in y0, rank r0, alive masks al0). counts[l] = nodes reached at level l.
    stop: process nodes only down to level 'stop' (stop=ncol for full search; smaller to emit prefixes)."""
    npairs = al0.shape[0]
    y = y0.copy(); r = np.zeros(ncol+1, np.int64); r[lev0] = r0
    al = np.zeros((ncol+1, npairs), np.int64); al[lev0] = al0
    idx = np.zeros(ncol+1, np.int64)
    lev = lev0
    while lev >= lev0:
        rr = r[lev]; left = ncol - lev
        forced = (5 - rr) == left
        nc = 1 if forced else ((1 << rr) - 1 + (1 if rr < 5 else 0))
        i = idx[lev]
        if i >= nc:
            lev -= 1; continue
        idx[lev] = i + 1
        if forced or i == (1 << rr) - 1: c = 1 << rr; new = True
        else: c = i + 1; new = False
        y[lev] = c
        ok = True
        for q in range(gstart[lev], gstart[lev+1]):
            if sval(y, gmask[q]) == 0:
                ok = False; break
        if not ok: continue
        for p in range(npairs): al[lev+1, p] = al[lev, p]
        for q in range(pstart[lev], pstart[lev+1]):
            al[lev+1, ppair[q]] &= ~(1 << sval(y, pmask[q]))
        if not cap_ok(al[lev+1], edgepairs): continue
        counts[lev+1] += 1
        r[lev+1] = rr + (1 if new else 0)
        if lev + 1 == stop:
            if nleaf < leaves.shape[0]:
                for j in range(ncol): leaves[nleaf, j] = y[j]
                leaves[nleaf, ncol] = r[lev+1]
                for p in range(npairs): leaves[nleaf, ncol+1+p] = al[lev+1, p]
            nleaf += 1
            continue
        lev += 1; idx[lev] = 0
    return nleaf
