"""Where do the extra triple relations come from? Geometric split of the triple intersection numbers of the
Ree-Tits octagon GO(2,4) by the shape of the configuration, for the six types with extra relations.

Incidence graph I (points and lines, girth 16, diameter 8). For a triple (u, v, w) the base configuration B is the
union of the geodesics in I between the non-opposite pairs (unique by (G2)). For a point x within distance 3 of
u, v and w, X(x) is B together with the geodesics x-u, x-v, x-w (all unique). x is
  * 'tree'  if X(x) has the same cycle rank as B (x hangs off B by a path; no new circuit),
  * 'c16'   if X(x) has new circuits and every new cycle in a minimum cycle basis has length 16 (apartments),
  * 'long'  otherwise (a new circuit of length > 16 is needed).
For each T-orbit representative and each count [ijk] with i, j, k <= 3 the count is split into these classes.
Checks:
  (1) per base shape (Weisfeiler-Lehman hash of B with u, v, w marked), is the tree part the same for every orbit?
  (2) is the c16 part the same for every orbit of a base shape?
  (3) does each relation of triple_relations.csv hold on the tree + c16 part and on the long part separately?
Usage: python triple/triple_shapes.py   (writes out_orbits/triple_shapes_output.txt, triple_shapes.csv)
"""
import os, sys, csv, re
from collections import defaultdict, Counter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'flags'))
import numpy as np
import networkx as nx
import go24_common as gc

OUT = os.path.join(HERE, 'out_orbits')
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)

Tg, _ = gc.generators()
D = gc.point_distances(Tg)
lines = gc.lines_from_distances(D)
N = D.shape[0]
adj = [[] for _ in range(N + len(lines))]
for j, Lp in enumerate(lines):
    for x in Lp:
        adj[x].append(N + j); adj[N + j].append(x)

def bfs_parent(r):
    par = [-1] * len(adj); dist = [-1] * len(adj); dist[r] = 0; q = [r]
    for a in q:
        for b in adj[a]:
            if dist[b] < 0:
                dist[b] = dist[a] + 1; par[b] = a; q.append(b)
    return par

def path_edges(par, x):
    E = set()
    while par[x] >= 0:
        E.add((min(x, par[x]), max(x, par[x]))); x = par[x]
    return E

def cycle_lengths(E):
    G = nx.Graph(); G.add_edges_from(E)
    return sorted(len(c) for c in nx.minimum_cycle_basis(G))

def shape_hash(E, u, v, w):
    G = nx.Graph(); G.add_edges_from(E)
    for n in [u, v, w]:
        G.add_node(n)
    for n in G.nodes:
        G.nodes[n]['lab'] = {u: 'u', v: 'v', w: 'w'}.get(n, 'P' if n < N else 'L')
    return nx.weisfeiler_lehman_graph_hash(G, node_attr='lab', iterations=6)

TYPES = [(2, 3, 3), (2, 3, 4), (2, 4, 4), (3, 3, 3), (3, 3, 4), (3, 4, 4)]
orbs = defaultdict(list)
with open(os.path.join(OUT, 'triple_orbit_vectors.csv')) as fh:
    for row in csv.DictReader(fh):
        ty = tuple(int(x) for x in re.findall(r'\d', row['type']))
        if ty in TYPES:
            orbs[ty].append((int(row['T_orbit']), tuple(int(x) for x in re.findall(r'\d+', row['rep']))))
rels = defaultdict(list)
with open(os.path.join(OUT, 'triple_relations_long.csv')) as fh:
    for row in csv.DictReader(fh):
        ty = tuple(int(x) for x in re.findall(r'\d', row['type']))
        rels[(ty, int(row['relation']))].append((row['unknown'], int(row['coefficient'])))

table = []
for ty in TYPES:
    per_orbit = []
    for o, (u, v, w) in orbs[ty]:
        par = {r: bfs_parent(r) for r in (u, v, w)}
        B = set()
        for a, b in [(u, v), (u, w), (v, w)]:
            if D[a, b] < 4:
                B |= path_edges(par[a], b)
        Gb = nx.Graph(); Gb.add_edges_from(B); Gb.add_nodes_from([u, v, w])
        rank0 = Gb.number_of_edges() - Gb.number_of_nodes() + nx.number_connected_components(Gb)
        base_cyc = cycle_lengths(B) if B else []
        sh = (shape_hash(B, u, v, w), tuple(base_cyc))
        cnt = {c: Counter() for c in ('tree', 'c16', 'long')}
        newcyc = Counter()
        sel = np.nonzero((D[u] <= 3) & (D[v] <= 3) & (D[w] <= 3))[0]
        for x in sel.tolist():
            X = set(B)
            for r in (u, v, w):
                X |= path_edges(par[r], x)
            key = f'{D[u, x]}{D[v, x]}{D[w, x]}'
            Gx = nx.Graph(); Gx.add_edges_from(X); Gx.add_nodes_from([u, v, w, x])
            rank = Gx.number_of_edges() - Gx.number_of_nodes() + nx.number_connected_components(Gx)
            if rank == rank0:
                cnt['tree'][key] += 1
            else:
                cl = cycle_lengths(X)
                new = list((Counter(cl) - Counter(base_cyc)).elements())
                newcyc[tuple(new)] += 1
                cnt['c16' if all(L_ == 16 for L_ in new) else 'long'][key] += 1
        per_orbit.append({'orbit': o, 'rep': (u, v, w), 'shape': sh, 'cnt': cnt, 'newcyc': newcyc})
    shapes = defaultdict(list)
    for po in per_orbit:
        shapes[po['shape']].append(po)
    say(f'type {ty}: {len(per_orbit)} T-orbits, {len(shapes)} base shapes '
        f'(base cycle lengths {[s_[1] for s_ in shapes]}, orbits per shape {[len(g) for g in shapes.values()]})')
    for si, (sh, grp) in enumerate(shapes.items()):
        same_tree = len({tuple(sorted(po['cnt']['tree'].items())) for po in grp}) == 1
        same_c16 = len({tuple(sorted(po['cnt']['c16'].items())) for po in grp}) == 1
        same_long = len({tuple(sorted(po['cnt']['long'].items())) for po in grp}) == 1
        nc = Counter()
        for po in grp: nc.update(po['newcyc'].keys())
        say(f'  shape {si} (base cycles {list(sh[1])}): tree part same on all its orbits: {same_tree}; c16 part same: '
            f'{same_c16}; long part same: {same_long}; new-cycle patterns seen: {sorted(nc)}')
    # relations on the parts
    nrel = max(k[1] for k in rels if k[0] == ty)
    hold = {'tree+c16': 0, 'long': 0}
    for li in range(1, nrel + 1):
        terms = [(k, c) for k, c in rels[(ty, li)] if k != 'rhs']
        for part in hold:
            vals = set()
            for po in per_orbit:
                cn = po['cnt']['long'] if part == 'long' else po['cnt']['tree'] + po['cnt']['c16']
                vals.add(sum(c * cn[k] for k, c in terms))
            if len(vals) == 1:
                hold[part] += 1
        table.append((ty, li, ' '.join(f'{c:+d}[{k}]' for k, c in terms)))
    say(f'  relations constant on the tree+c16 part alone: {hold["tree+c16"]} of {nrel}; on the long part alone: '
        f'{hold["long"]} of {nrel}')
    with open(os.path.join(OUT, f'triple_shapes_{"".join(map(str, ty))}.csv'), 'w', newline='\n') as fh:
        fh.write('T_orbit,rep,shape,class,ijk,count\n')
        for po in per_orbit:
            si = list(shapes).index(po['shape'])
            for cl in ('tree', 'c16', 'long'):
                for k, c in sorted(po['cnt'][cl].items()):
                    fh.write(f'{po["orbit"]},"{po["rep"]}",{si},{cl},{k},{c}\n')
open(os.path.join(OUT, 'triple_shapes_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
