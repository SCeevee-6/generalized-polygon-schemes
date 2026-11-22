"""The 25-line families in the dual of O(2) (output of escape_bk.py --geom lines --min 25).

Checks: the families form one orbit of the group; the stabiliser of line 0 is transitive on the families through
it (so the stabiliser of a family is transitive on its 25 lines); and for random elements g of order 5
(seeded), g fixes no line, fixes exactly one family, and acts on it without fixed lines.

Usage: python triple/escape_dual_orbit.py [--gens atlas|gap:DIR|file.npy] [--samples 50]
"""
import argparse, json, os
from collections import Counter
from math import lcm
import numpy as np
from escape_common import build
from escape_bk import schreier_orbits

HERE = os.path.dirname(os.path.abspath(__file__))


def order(g):
    seen = np.zeros(len(g), bool); o = 1
    for i in range(len(g)):
        if not seen[i]:
            l, j = 0, i
            while not seen[j]:
                seen[j] = True; j = g[j]; l += 1
            o = lcm(o, l)
    return o


def power(g, k):
    r = np.arange(len(g)); b = g.copy()
    while k:
        if k & 1:
            r = b[r]
        b = b[b]; k >>= 1
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gens', default='atlas')
    ap.add_argument('--samples', type=int, default=50)
    ap.add_argument('--sets', default=os.path.join(HERE, 'out_escape', 'escape_bk_lines_min25_sets.txt'))
    a = ap.parse_args()
    D, nbhd, G = build('lines', a.gens)
    fam = [tuple(map(int, r.split()[1:])) for r in open(a.sets) if r.startswith('E ') and len(r.split()) == 26]
    act = lambda g, x: tuple(sorted(g[list(x)].tolist()))
    orbit = {fam[0]}; fr = [fam[0]]
    while fr:
        nf = []
        for X in fr:
            for g in G:
                Y = act(g, X)
                if Y not in orbit:
                    orbit.add(Y); nf.append(Y)
        fr = nf
    through = {X for X in orbit if 0 in X}
    nlines, orb = schreier_orbits(fam, (0,), G, act)
    d = D[np.ix_(fam[0], fam[0])]
    res = {'families_through_line_0': len(fam), 'orbit_size': len(orbit), 'orbit_through_0_equals_list': through == set(fam),
           'stab_line0_orbits_on_them': orb,
           'pairwise_distances_in_a_family': dict(Counter(d[np.triu_indices(25, 1)].tolist()))}
    rng = np.random.default_rng(2026)
    stats = Counter(); n = 0
    while n < a.samples:
        w = np.arange(D.shape[0])
        for _ in range(40):
            w = G[rng.integers(len(G))][w]
        o = order(w)
        if o % 5:
            continue
        g = power(w, o // 5); n += 1
        fixed = [X for X in orbit if act(g, X) == X]
        stats[(len(fixed), int((g == np.arange(len(g))).sum()),
               bool(fixed) and all(g[x] != x for x in fixed[0]))] += 1
    res['order5_samples'] = n
    res['order5 (families fixed, lines fixed, fixed-point-free on its family)'] = {str(k): v for k, v in stats.items()}
    print(json.dumps(res), flush=True)
    with open(os.path.join(HERE, 'out_escape', 'escape_dual_orbit.json'), 'w', newline='\n') as f:
        json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
