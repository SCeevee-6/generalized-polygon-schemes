"""Check of the last step in Chen-Makhnev-Klimin, Proc. Inst. Math. NAS Belarus 32:1 (2024) 57-63, Theorem 1.

In section 2 ({44,36,5;1,9,40}) they fix u, v with d(u,v) = 3 and count the edges y ~ z with
y in {uv;31} (d(u,y) = 3, d(v,y) = 1) and z in {uv;32} (d(u,z) = 3, d(v,z) = 2) in two ways:
  via y (Lemma 2):  f1 <= 4 * 3 = 12;
  via z (Lemma 4):  f1 = sum_z [311](u,v,z) = sum_z (r68 - r66) - 16, where for the triple (u, v, z)
                    r68 = [231] and r66 = [321].
So sum_z ([231] - [321]) = f1 + 16 lies in [16, 28]. The contradiction is then drawn from
"[231] = r68 = [321]* = r*66" and "sum_i (r^i_68 - r^i_66) / 16 = 0", i.e. from
    sum over z in {uv;32} of [231](u,v,z)  =  sum over z in {uv;32} of [321](u,v,z).            (*)
The left side counts edges between {uv;23} and {uv;32}; the right side counts ordered adjacent pairs inside
{uv;32}. Section 3 ({44,36,12;1,3,33}) ends with the same step (r101 = [231], r102 = [321]).
(*) is not an identity for distance-regular graphs of diameter 3: this script computes
S(u,v) = sum_z ([231] - [321]) in several known ones, where (*) would force S = 0.
Usage: python search/drg/cmk_check.py   (writes search/drg/cmk_check_output.txt)
"""
import os, itertools
import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)

def array_of(D):
    n = len(D); u = 0
    b, c = [], []
    for i in range(int(D.max()) + 1):
        vs = np.nonzero(D[u] == i)[0]
        w = vs[0]
        nb = np.nonzero(D[w] == 1)[0]
        c.append(int(np.sum(D[u, nb] == i - 1))); b.append(int(np.sum(D[u, nb] == i + 1)))
    return b[:-1], c[1:]

def check(G, name):
    nodes = list(G)
    D = np.array(nx.floyd_warshall_numpy(G, nodelist=nodes), dtype=int)
    A = (D == 1)
    b, c = array_of(D)
    vals = set()
    for u in range(len(nodes)):
        for v in np.nonzero(D[u] == 3)[0]:
            Z = np.nonzero((D[u] == 3) & (D[v] == 2))[0]
            Y23 = (D[u] == 2) & (D[v] == 3)
            Y32 = (D[u] == 3) & (D[v] == 2)
            lhs = int(A[np.ix_(Z, np.nonzero(Y23)[0])].sum())     # sum_z [231]
            rhs = int(A[np.ix_(Z, np.nonzero(Y32)[0])].sum())     # sum_z [321]
            vals.add(lhs - rhs)
    say(f'{name}: intersection array {{{",".join(map(str, b))};{",".join(map(str, c))}}}, {len(nodes)} vertices; '
        f'S(u,v) over all pairs at distance 3: {sorted(vals)}')
    return vals

def johnson(n, k):
    S = list(itertools.combinations(range(n), k)); G = nx.Graph()
    G.add_edges_from((a, b_) for a, b_ in itertools.combinations(S, 2) if len(set(a) & set(b_)) == k - 1)
    return G
def odd(m):
    S = list(itertools.combinations(range(2 * m - 1), m - 1)); G = nx.Graph()
    G.add_edges_from((a, b_) for a, b_ in itertools.combinations(S, 2) if not set(a) & set(b_))
    return G
def hamming(d, q):
    return nx.convert_node_labels_to_integers(nx.cartesian_product(*( [nx.complete_graph(q)] * d)) if d == 2 else
                                              nx.cartesian_product(nx.cartesian_product(nx.complete_graph(q), nx.complete_graph(q)), nx.complete_graph(q)))
res = {}
res['O4'] = check(odd(4), 'Odd graph O4')
res['J93'] = check(johnson(9, 3), 'Johnson J(9,3)')
res['J103'] = check(johnson(10, 3), 'Johnson J(10,3)')
res['H34'] = check(hamming(3, 4), 'Hamming H(3,4)')
res['H35'] = check(hamming(3, 5), 'Hamming H(3,5)')
res['H33'] = check(hamming(3, 3), 'Hamming H(3,3)')
nonzero = [k for k, v in res.items() if v != {0}]
say(f'graphs where (*) fails for some pair: {nonzero}')
say('Conclusion: (*) is not a consequence of distance-regularity, so the last step of Theorem 1 (both arrays) needs an '
    'argument that is not in the paper as I read it. The two counts of f1 only give sum_z([231] - [321]) in [16, 28] '
    'for {44,36,5;1,9,40}, which is not a contradiction by itself.')
open(os.path.join(HERE, 'cmk_check_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
