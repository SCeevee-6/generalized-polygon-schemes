"""Model-only checks for the 512 inputs and the dual packing. Requires igraph."""
import io
import json
import zipfile
from collections import defaultdict
from pathlib import Path

import igraph as ig
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path

from sibling_obstruction import outer_geometry

ROOT = Path(__file__).resolve().parents[1]


def model():
    archive = ROOT / "octagon-regular-256-exclusion/generalized-octagon-research.zip"
    with zipfile.ZipFile(archive) as z:
        generators = np.load(io.BytesIO(z.read("generators.npy")))
    # Orbit of the collinear pair (0,72) under the two stored ATLAS generators.
    edges = [(0, 72)]
    seen = set(edges)
    for u, v in edges:
        for g in generators:
            e = tuple(sorted((int(g[u]), int(g[v]))))
            if e not in seen:
                seen.add(e)
                edges.append(e)
    u, v = np.array(edges).T
    half = csr_matrix((np.ones(len(edges)), (u, v)), shape=(1755, 1755))
    return half + half.T


if __name__ == "__main__":
    graph = model()
    d = shortest_path(graph, indices=0, unweighted=True)
    assert np.bincount(d.astype(int)).tolist() == [1, 10, 80, 640, 1024]
    outer = np.flatnonzero(d == 4)
    induced = graph[outer][:, outer]
    n, labels = connected_components(induced)
    assert n == 2 and list(np.bincount(labels)) == [512, 512]
    data = json.loads(Path(__file__).with_name("stage4_fp0.json").read_text())
    classes = defaultdict(list)
    for i, record in enumerate(data):
        edges, _ = outer_geometry(record["vol"])
        g = ig.Graph(512, edges.tolist())
        key = tuple(sorted(g.permute_vertices(g.canonical_permutation()).get_edgelist()))
        classes[key].append(i)
    print("first_row copies automorphisms diameter")
    for key, ids in classes.items():
        g = ig.Graph(512, list(key))
        print(ids[0], len(ids), g.count_automorphisms(), g.diameter())
        if 6 in ids:
            for c in range(n):
                k = induced[labels == c][:, labels == c].tocoo()
                edges = [(int(u), int(v)) for u, v in zip(k.row, k.col) if u < v]
                assert g.isomorphic_bliss(ig.Graph(512, edges))
    print("row 6 isomorphic to both model components")
    neighbours = [set(graph[i].indices) for i in range(graph.shape[0])]
    lines = sorted((u, v, w) for u, nbr in enumerate(neighbours) for v in nbr if u < v
                   for w in nbr & neighbours[v] if v < w)
    assert len(lines) == 2925
    inc = csr_matrix((np.ones(3*len(lines)), (np.repeat(np.arange(len(lines)), 3),
                     np.asarray(lines).ravel())), shape=(len(lines), 1755))
    dual = inc @ inc.T
    dual.setdiag(0)
    dual.eliminate_zeros()
    centres = np.loadtxt(ROOT / "packing/dual_packing_133.txt", dtype=int)
    assert len(set(centres)) == 133
    dist = shortest_path(dual, indices=centres, unweighted=True).astype(int)
    pairs = dist[:, centres][np.triu_indices(len(centres), 1)]
    print("dual pair histogram", np.bincount(pairs).tolist(), "ball union", int((dist.min(0) <= 1).sum()))
