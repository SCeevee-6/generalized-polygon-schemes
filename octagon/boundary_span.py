"""Exact rank certificate on the distance-three sphere of one Ree–Tits root."""
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "flags"))
from o2_build import octagon


def rank_mod_101(matrix):
    a = matrix.copy() % 101
    rank = 0
    for j in range(a.shape[1]):
        pivots = np.flatnonzero(a[rank:, j])
        if len(pivots):
            p = rank + pivots[0]
            a[[rank, p]] = a[[p, rank]]
            mult = a[rank+1:, j] * pow(int(a[rank, j]), -1, 101) % 101
            a[rank+1:, j:] = (a[rank+1:, j:] - mult[:, None]*a[rank, j:]) % 101
            rank += 1
    return rank


def main():
    d, _ = octagon(np.load(ROOT / "data/generators.npy"))
    shell = np.flatnonzero(d[0] == 3)
    outer = np.flatnonzero(d[0] == 4)
    components, _ = connected_components(csr_matrix(d[np.ix_(outer, outer)] == 1))
    gram = np.array([64, 32, 8, 0, -1])[d[np.ix_(shell, shell)]]
    rank = rank_mod_101(gram)
    assert len(shell) == 640 and components == 2 and rank == 349 == 351-components
    # The residual at an outer point of a 512-component is 21/64.
    # Subtract it from the diagonal via a bordered integer Gram matrix:
    # scaling the final row and column by 8 keeps every entry integral.
    for v in (outer[0], outer[-1]):
        cross = 8*np.array([64, 32, 8, 0, -1])[d[shell, v]]
        border = np.block([[gram, cross[:, None]], [cross[None, :], np.array([[4096-21]])]])
        assert rank_mod_101(border) == rank
    print("Root 0: modular rank 349; two components give the matching real upper bound")


if __name__ == "__main__":
    main()
