"""Local model for escaping sets that contain a full line L and have every other point at distance 2 from L.

For p on L and a line M_p != L through p, X(p, M_p) is the set of s^2 t points at distance 2 from L whose
geodesic to L is x - q - p with q on M_p. Two points with the same projection p are non-opposite iff their
middle points lie on one line through p, so each projection class lies in a single X(p, M_p). Points in
different classes are non-opposite iff they are at distance 3 ("compatible").

For every choice of (M_p)_{p in L} this script checks in O(2):
  - |X(p, M_p)| = s^2 t and points in one X are pairwise at distance <= 3;
  - each x in X(p, M_p) is compatible with exactly one point of N' - {q'} for every q' in M_p' - {p'} and
    every line N' != M_p' on q' (so with s t points of X(p', M_p'));
computes the twin classes of X(p, M_p) with respect to X(p', M_p') (points with the same compatible set) and
the graph between twin classes; and finds all maximal sets L + (subsets of the X's, pairwise compatible across classes) with at least two
nonempty classes (exactly the ones that escape every line neighbourhood).

Usage: python triple/escape_local.py [--gens atlas|gap:DIR|file.npy]
"""
import argparse, itertools, json, os
from collections import Counter
import numpy as np
import networkx as nx
from escape_common import build, gc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gens', default='atlas')
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_escape'))
    a = ap.parse_args()
    D, nbhd, G = build('points', a.gens)
    lines = gc.lines_from_distances(D)
    L = lines[0]
    s, t = 2, 4
    dL = D[list(L)].min(axis=0)
    at = {}
    for M in lines:
        for x in M:
            at.setdefault(x, []).append(M)
    thru = {p: [M for M in at[p] if M != L] for p in L}

    def X(p, M):
        qs = [q for q in M if q != p]
        return [x for x in np.nonzero((dL == 2) & (D[p] == 2))[0] if any(D[q, x] == 1 for q in qs)]

    stats = Counter(); best = Counter(); patterns = Counter(); n_choice = 0
    twin = Counter(); cycle = Counter(); partitions = {}
    for Ms in itertools.product(*(thru[p] for p in L)):
        n_choice += 1
        Xs = [[int(x) for x in X(p, M)] for p, M in zip(L, Ms)]
        for p, M in zip(L, Ms):
            partitions.setdefault((p, M), set())
        assert all(len(x) == s * s * t for x in Xs)
        assert all(D[np.ix_(x, x)].max() <= 3 for x in Xs)
        for i, j in itertools.permutations(range(3), 2):
            pj, Mj = L[j], Ms[j]
            for x in Xs[i]:
                for q in Mj:
                    if q == pj:
                        continue
                    for N in at[q]:
                        if N != Mj:
                            c = sum(1 for y in N if y != q and D[x, y] == 3)
                            assert c == 1
                assert sum(1 for y in Xs[j] if D[x, y] == 3) == s * t
                assert all(D[x, y] in (3, 4) for y in Xs[j])
        stats['one point per line'] += 1
        for i, j in itertools.permutations(range(3), 2):
            Cm = D[np.ix_(Xs[i], Xs[j])] == 3
            rows = {}
            for r, x in enumerate(Xs[i]):
                rows.setdefault(Cm[r].tobytes(), []).append(x)
            cols = {}
            for c, y in enumerate(Xs[j]):
                cols.setdefault(Cm[:, c].tobytes(), []).append(c)
            twin[tuple(sorted(len(v) for v in rows.values()))] += 1
            partitions[(L[i], Ms[i])].add(frozenset(frozenset(v) for v in rows.values()))
            Hc = nx.Graph()
            for a_, ra in enumerate(rows.values()):
                for b_, cb in enumerate(cols.values()):
                    if Cm[Xs[i].index(ra[0]), cb[0]]:
                        Hc.add_edge(('a', a_), ('b', b_))
            cycle[nx.is_isomorphic(Hc, nx.cycle_graph(2 * len(rows)))] += 1
        H = nx.Graph()
        V = [(i, int(x)) for i in range(3) for x in Xs[i]]
        H.add_nodes_from(V)
        for (i, x), (j, y) in itertools.combinations(V, 2):
            if i == j or D[x, y] == 3:
                H.add_edge((i, x), (j, y))
        sizes = Counter()
        for C in nx.find_cliques(H):
            cls = Counter(i for i, _ in C)
            if len(cls) >= 2:
                sizes[len(C)] += 1
                if len(C) == 12:
                    patterns[tuple(sorted([cls.get(0, 0), cls.get(1, 0), cls.get(2, 0)]))] += 1
        best[max(sizes)] += 1
        stats['max escaping sizes seen'] = max(stats['max escaping sizes seen'], max(sizes))
    res = {'line': list(map(int, L)), 'choices': n_choice, 'checks_passed': dict(stats),
           'max_extra_points_by_choice': {str(k): v for k, v in best.items()},
           'patterns_of_size_12': {str(k): v for k, v in patterns.items()},
           'twin_class_sizes (X_i w.r.t. X_j)': {str(k): v for k, v in twin.items()},
           'twin class graph is a cycle': {str(k): v for k, v in cycle.items()},
           'distinct twin partitions per X(p,M)': dict(Counter(len(v) for v in partitions.values()))}
    print(json.dumps(res), flush=True)
    with open(os.path.join(a.out, 'escape_local.json'), 'w', newline='\n') as f:
        json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
