"""Multiplicity of the Smith factor 2 in the point and line opposition matrices of O(2).
For an integer matrix A put r0 = rank of A mod 2. On x in ker(A mod 2) take an integer lift, then y = A x / 2 mod 2.
The map x -> y is well defined modulo the image of A mod 2, and its rank there is the number e1 of invariant factors
exactly divisible by 2. So e1 = rank[A mod 2 | Y] - r0 with Y the images of a kernel basis.
Run: python smith_two.py generators.npy      expected  points 26 272   lines 246 350"""
import sys, numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
def octagon(gens, n=1755):
    idx = np.arange(n*n); a, b = idx//n, idx % n
    r = np.concatenate([idx]*len(gens)); c = np.concatenate([g[a]*n+g[b] for g in gens])
    _, lab = connected_components(csr_matrix((np.ones(len(r), np.int8), (r, c)), shape=(n*n,)*2), directed=True, connection='weak')
    L = lab.reshape(n, n); col = [l for l in np.unique(L[0]) if (L[0] == l).sum() == 10][0]
    A = L == col; D = shortest_path(csr_matrix(A), unweighted=True).astype(np.int8)
    lines = sorted({tuple(sorted((int(x), int(y), int(np.flatnonzero(A[x] & A[y])[0])))) for x, y in zip(*np.nonzero(np.triu(A)))})
    return D, np.array(lines)
def to_ints(M):
    packed = np.packbits(M.astype(np.uint8), axis=1, bitorder='little')
    return [int.from_bytes(r.tobytes(), 'little') for r in packed]
def rank2(rows):
    piv = {}; r = 0
    for v in rows:
        while v:
            h = v.bit_length()-1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; r += 1; break
    return r
def kernel2(M):
    n = M.shape[1]; piv = []
    for v in to_ints(M % 2):
        for c, rw in piv:
            if v >> c & 1: v ^= rw
        if v:
            c = (v & -v).bit_length()-1
            piv = [(cc, rr ^ v if rr >> c & 1 else rr) for cc, rr in piv]; piv.append((c, v))
    pc = {c for c, _ in piv}; basis = []
    for f in (j for j in range(n) if j not in pc):
        x = 1 << f
        for c, rw in piv:
            if rw >> f & 1: x |= 1 << c
        basis.append(x)
    return len(piv), basis
def e0_e1(A):
    n = A.shape[0]; r0, ker = kernel2(A)
    X = np.array([[(x >> j) & 1 for j in range(n)] for x in ker], dtype=np.float64).T
    AX = A.astype(np.float64) @ X; assert np.all(AX % 2 == 0)
    Y = ((AX/2) % 2).astype(np.uint8)
    return r0, rank2(to_ints(np.concatenate([A % 2, Y], axis=1).T)) - r0
if __name__ == '__main__':
    D, lines = octagon([g.astype(np.int64) for g in np.load(sys.argv[1])])
    print('points', *e0_e1((D == 4).astype(np.int64)))
    nL = len(lines); I = csr_matrix((np.ones(3*nL), (np.repeat(np.arange(nL), 3), lines.ravel())), shape=(nL, 1755))
    LL = (I @ I.T).toarray(); np.fill_diagonal(LL, 0)
    DL = shortest_path(csr_matrix(LL > 0), unweighted=True).astype(np.int8)
    print('lines', *e0_e1((DL == 4).astype(np.int64)))
