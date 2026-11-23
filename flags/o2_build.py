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
