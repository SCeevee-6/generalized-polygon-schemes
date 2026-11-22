"""Geometry for the escape problem: pairwise non-opposite sets that lie in no line neighbourhood.

geom = 'points': the Ree-Tits octagon O(2) of order (2,4) on its 1755 points.
geom = 'lines' : its dual, order (4,2), whose points are the 2925 lines of O(2).

build(geom, gens) returns (D, nbhd, gens_on_vertices) where D is the distance matrix of the collinearity graph,
nbhd is a list of frozensets N[L] (points on or collinear with L) over all lines L of that geometry, and
gens_on_vertices are the group generators acting on the vertex set.
gens: 'atlas' (fetch c, d of TF42d2G1-p1755B0, via flags/go24_common) or a path to a (2, 1755) .npy array.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'flags'))
import go24_common as gc  # noqa: E402


def parse_gap_perm(path, n=1755):
    import re
    t = open(path).read()
    perm = np.arange(n)
    for cyc in re.findall(r'\(([\d,\s]+)\)', t[t.find('('):]):
        pts = [int(x) - 1 for x in re.sub(r'\s', '', cyc).split(',')]
        perm[pts] = np.roll(pts, -1)
    return perm


def load_gens(gens):
    """'atlas': c, d of 2F4(2) = T.2 (TF42d2G1-p1755B0, fetched by go24_common).
    'gap:DIR': the ATLAS files TF42G1-p1755B0.g1/.g2 (generators of T = 2F4(2)') read from DIR.
    otherwise: a (2, 1755) .npy array of zero-based images."""
    if gens == 'atlas':
        return [np.asarray(g) for g in gc.generators()[1]]
    if gens.startswith('gap:'):
        d = gens[4:]
        return [parse_gap_perm(os.path.join(d, f'TF42G1-p1755B0.g{i}')) for i in (1, 2)]
    return [np.asarray(g) for g in np.load(gens)]


def build(geom, gens):
    G = load_gens(gens)
    D = gc.point_distances(G)
    assert tuple(np.bincount(D[0])) == (1, 10, 80, 640, 1024)
    lines = gc.lines_from_distances(D)
    assert len(lines) == 2925
    if geom == 'points':
        nb = [frozenset(np.nonzero(D[list(L)].min(axis=0) <= 1)[0].tolist()) for L in lines]
        return D, nb, G
    # dual: vertices are lines; two lines are adjacent iff they meet
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path
    idx = {L: i for i, L in enumerate(lines)}
    inc = np.zeros((len(lines), gc.N), np.int8)
    for i, L in enumerate(lines):
        inc[i, list(L)] = 1
    A = (inc.astype(np.int32) @ inc.T.astype(np.int32)) > 0
    np.fill_diagonal(A, False)
    Dl = shortest_path(csr_matrix(A.astype(np.int8)), unweighted=True, directed=False).astype(np.int8)
    assert tuple(np.bincount(Dl[0])) == (1, 12, 96, 768, 2048)
    # dual lines are points p of O(2); N*[p] = lines meeting the closed neighbourhood of p
    near = D <= 1
    nb = [frozenset(np.nonzero(inc[:, near[p]].any(axis=1))[0].tolist()) for p in range(gc.N)]
    Gl = [np.array([idx[tuple(sorted(g[list(L)]))] for L in lines]) for g in G]
    return Dl, nb, Gl
