"""Is the SRG(45,16,8,4) from the thin octagon of order (2,1) the triangular graph T(10)?

The thin octagon of order (2,1) is the flag geometry of the generalized quadrangle W(2): its points are the 45
flags of W(2), and two flags are collinear when they share their point or their line. W(2) is built as duads and
synthemes of a 6-set. The fusion {3}|{1,2,4} of the distance scheme is the distance-3 graph. This script builds it,
checks the parameters and tests isomorphism with the line graph of K_10 (networkx VF2).

Usage: python verify_t10.py
"""
import itertools
import networkx as nx


def synthemes():
    pts = list(range(6))
    out = set()
    for perm in itertools.permutations(pts):
        m = frozenset(frozenset(perm[i:i + 2]) for i in (0, 2, 4))
        out.add(m)
    return sorted(out, key=lambda m: sorted(map(sorted, m)))


def main():
    duads = [frozenset(c) for c in itertools.combinations(range(6), 2)]
    syn = synthemes()
    assert len(duads) == 15 and len(syn) == 15
    flags = [(p, l) for l in syn for p in l]
    assert len(flags) == 45
    G = nx.Graph()
    G.add_nodes_from(range(45))
    for i, j in itertools.combinations(range(45), 2):
        if flags[i][0] == flags[j][0] or flags[i][1] == flags[j][1]:
            G.add_edge(i, j)
    dist = dict(nx.all_pairs_shortest_path_length(G))
    assert max(max(d.values()) for d in dist.values()) == 4
    k = [sum(1 for j in range(45) if dist[0][j] == h) for h in range(5)]
    G3 = nx.Graph()
    G3.add_nodes_from(range(45))
    G3.add_edges_from((i, j) for i, j in itertools.combinations(range(45), 2) if dist[i][j] == 3)
    lam = {len(set(G3[i]) & set(G3[j])) for i, j in G3.edges()}
    mu = {len(set(G3[i]) & set(G3[j])) for i, j in itertools.combinations(range(45), 2) if not G3.has_edge(i, j)}
    deg = {d for _, d in G3.degree()}
    T10 = nx.line_graph(nx.complete_graph(10))
    iso = nx.is_isomorphic(G3, T10)
    print('octagon valencies', k)
    print('distance-3 graph: degree', deg, 'lambda', lam, 'mu', mu)
    print('isomorphic to T(10):', iso)
    return iso


if __name__ == '__main__':
    main()
