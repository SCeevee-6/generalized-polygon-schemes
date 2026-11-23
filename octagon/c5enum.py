"""Own reproduction of the C5-stable regular elementary 32-cover count (2,167,056 -> 251)."""
import numpy as np, numba, sys, json
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from base import *
def gf2_rank(vs):
    piv = {}
    for v in vs:
        while v:
            h = v.bit_length()-1
            if h not in piv: piv[h] = v; break
            v ^= piv[h]
    return piv
def kernel(M, n):  # M: list of row-ints acting on column vectors? here: M[j] = image of e_j (column j), find x with sum x_j M[j] = 0
    # Gaussian elimination tracking combos
    rows = [(M[j], 1 << j) for j in range(n)]; piv = {}; ker = []
    for v, c in rows:
        while v:
            h = v.bit_length()-1
            if h not in piv: piv[h] = (v, c); break
            pv, pc = piv[h]; v ^= pv; c ^= pc
        if v == 0: ker.append(c)
    return ker
def rot(u): return ((u << 1) & 31) | (u >> 4)
EROT = [EI[tuple(sorted((rot(a), rot(b))))] for a, b in EDGES]
def rot_mask(m): return sum(1 << EROT[k] for k in range(40) if (m >> k) & 1)
FUND = []  # fundamental cycle masks for chords (via tree paths)
def tree_path(u):
    p = []
    while PAR[u] is not None: p.append(EI[tuple(sorted((u, PAR[u])))]); u = PAR[u]
    return p
for k in CHORDS:
    a, b = EDGES[k]; m = 1 << k
    for e in tree_path(a) + tree_path(b): m ^= 1 << e
    FUND.append(m)
def coords(mask): return chord_coords(mask)
R = [coords(rot_mask(FUND[j])) for j in range(25)]  # column j = rho(e_j)
def apply(Mcols, x):
    y = 0
    for j in range(len(Mcols)):
        if (x >> j) & 1: y ^= Mcols[j]
    return y
def main():
    cyc = simple_cycles(7)
    from collections import Counter
    print('cycle lengths', sorted(Counter(cyc.values()).items()))
    # check rho maps cycle space consistently
    for m in list(cyc)[:50]: assert coords(rot_mask(m)) == apply(R, coords(m))
    I = [1 << j for j in range(25)]
    RF = [R[j] ^ I[j] for j in range(25)]
    def poww(x, k):
        for _ in range(k): x = apply(R, x)
        return x
    eF = [sum((poww(1 << j, k) for k in range(5)), 0) for j in range(25)]
    eF = [0]*25
    for j in range(25):
        acc = 0; x = 1 << j
        for k in range(5): acc ^= x; x = apply(R, x)
        eF[j] = acc
    Fix = list(gf2_rank([eF[j] for j in range(25)]).values())
    Nn = list(gf2_rank([eF[j] ^ (1 << j) for j in range(25)]).values())
    print('dims', len(Fix), len(Nn))
    # coordinate solvers in bases
    def solver(basis):
        piv = {}; 
        for i, b in enumerate(basis):
            v, c = b, 1 << i
            while v:
                h = v.bit_length()-1
                if h not in piv: piv[h] = (v, c); break
                pv, pc = piv[h]; v ^= pv; c ^= pc
        def solve(x):
            c = 0
            while x:
                h = x.bit_length()-1; pv, pc = piv[h]; x ^= pv; c ^= pc
            return c
        return solve
    sF, sN = solver(Fix), solver(Nn)
    cl = list(cyc)
    cF = np.array([sF(apply(eF, coords(m))) for m in cl], np.int64)
    cN = np.array([sN(apply(eF, coords(m)) ^ coords(m)) for m in cl], np.int64)
    # rho on N in N-coordinates: RN[i] = coords of rho(N_i); f o rho^j as functional
    RN = [sN(apply(R, b)) for b in Nn]
    # functional f (20-bit, f.x parity) composed with rho: (f o rho)(x) = f(RN x) -> vector g with g_i = f.RN[i]
    RNT = np.array(RN, np.int64)
    return cl, cF, cN, RNT, Fix, Nn, eF, sF, sN
@numba.njit(cache=True)
def par(x):
    x ^= x >> 32; x ^= x >> 16; x ^= x >> 8; x ^= x >> 4; x ^= x >> 2; x ^= x >> 1
    return x & 1
@numba.njit(cache=True)
def compose(f, RNT):
    g = 0
    for i in range(20):
        if par(f & RNT[i]): g |= 1 << i
    return g
@numba.njit(cache=True)
def count(cF, cN, RNT):
    total = 0; surv = []
    for f in range(1, 1 << 20):
        f1 = compose(f, RNT); f2 = compose(f1, RNT); f3 = compose(f2, RNT)
        okg = np.ones(32, np.bool_); okg[0] = False
        for c in range(len(cN)):
            x = cN[c]
            if par(f & x) == 0 and par(f1 & x) == 0 and par(f2 & x) == 0 and par(f3 & x) == 0:
                for g in range(1, 32):
                    if par(g & cF[c]) == 0: okg[g] = False
        for g in range(1, 32):
            if okg[g]:
                total += 1; surv.append((f, g))
    return total, surv
