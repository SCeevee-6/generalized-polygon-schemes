"""Second look at the covers that survive elem_run (stage 1 + H-sibling filter).
For each survivor: rebuild the cover, count mates, test triangle-freeness of the mate graphs, and run the outer-distance
sibling test in every category. Writes survivors_analysis.json; isomorphism with
the Ree-Tits component is checked separately (iso_check.py, needs python-igraph)."""
import sys, os, json, numpy as np, networkx as nx
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import elem_setup as S
from base import *
def edges_of(vol): return [(32*VI[a]+t, 32*VI[b]+(t ^ vol[k])) for k, (a, b) in enumerate(EDGES) for t in range(32)]
def analyse(vol):
    E = np.array(edges_of(vol)); adj = cover(vol); dist = np.array([bfs(adj, s) for s in range(512)])
    catv = np.repeat([S.CATS.index(CAT[k]) for k in range(40)], 32); rep = []
    for c in range(10):
        ids = np.flatnonzero(catv == c); X = E[ids]; flat = X.ravel()
        dd = dist[np.ix_(flat, flat)].reshape(len(X), 2, len(X), 2).min((1, 3)); np.fill_diagonal(dd, 0)
        mate = dd >= 6; g = nx.from_numpy_array(mate.astype(int))
        tri = sum(nx.triangles(g).values()) > 0
        pairs = list(g.edges()); compat = 0
        if True:   # every boundary point carries exactly two K-edges (boundary count), triangles or not
            for i, (a, b) in enumerate(pairs):
                for (c2, d2) in pairs[i+1:]:
                    if {a, b} & {c2, d2}: continue
                    if min(dd[a, c2], dd[a, d2], dd[b, c2], dd[b, d2]) >= 5: compat += 1
        rep.append(dict(cat=c, triangle_free=not tri, mate_edges=len(pairs), sibling_compatible_pairs=compat,
                        comps=dict(Counter(len(x) for x in nx.connected_components(g)))))
    return rep
def verdict(rep):
    """reject if some category has no two disjoint mate pairs at delta >= 5 across. Valid in every category: by the
    boundary count each boundary point carries exactly two K-edges, they are mates, and its sibling carries two more."""
    return any(r['sibling_compatible_pairs'] == 0 for r in rep)
if __name__ == '__main__':
    fp = int(sys.argv[1]); A = S.setup_multi2([fp]); order = list(A['order']); out = []
    for path in sys.argv[2:]:
        R = json.load(open(path)); R = R.get('results', R)
        for r in R.values():
            if not r: continue
            for row in r['surv']:
                vol, tau = S.from_seq(row[:A['ncol']], order); rep = analyse(vol)
                mc = mate_counts(vol, fibre_dist(cover(vol)))
                out.append(dict(vol=[int(v) for v in vol], mates=[int(m) for m in mc], rejected=verdict(rep), cats=rep))
    json.dump(out, open('stage3_fp%d.json' % fp, 'w'))
    print('fp', fp, 'survivors', len(out), 'rejected by outer sibling test', sum(o['rejected'] for o in out))
