"""Which triple intersection numbers [ijk] are 0 in EVERY generalized octagon of order (s, t)? A mechanical test.

Setting. u, v, w points with collinearity distances (U, V, W) = (d(u,v), d(u,w), d(v,w)); x a point with
d(x,u) = i, d(x,v) = j, d(x,w) = k, all <= 3. Work in the incidence graph (girth 16, diameter 8, point degree t+1,
line degree s+1).

Lemma. The geodesics from x to u, v, w are unique and their union T_x is a tree (two of them that separate and meet
again would close a circuit of length <= 12 < 16). T_x is fixed by the three Gromov products
g_ab = length of the common initial segment of the geodesics x->a and x->b; for a tree, the two smallest of
g_uv, g_uw, g_vw are equal. Let d_T(a,b) = 2d(x,a) + 2d(x,b) - 2 g_ab be the length of the path from a to b in T_x.
That path is non-backtracking, so
  * if d_T(a,b) <= 8 it is a geodesic and 2 d(a,b) = d_T(a,b);
  * if d_T(a,b) > 8, the path and a geodesic a-b close a walk containing a circuit, so 2 d(a,b) + d_T(a,b) >= 16.
At a branch node of T_x the number of outgoing edges plus the incoming edge is at most the node's degree
(s+1 for a line, t+1 for a point; the root x has no incoming edge).
If no choice of (g_uv, g_uw, g_vw) meets all of this, then [ijk] = 0 for the type (U, V, W) in every GO(s, t).
This is a necessary-condition test only: a feasible choice does not mean the count is nonzero.
Usage: python triple/triple_zero_proofs.py [s t]   (default 2 4; writes out_orbits/triple_zero_proofs_output.txt)
"""
import os, sys, csv, re, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out_orbits')
s_, t_ = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (2, 4)
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)

def deg(depth):
    return t_ + 1 if depth % 2 == 0 else s_ + 1

def feasible(U, V, W, i, j, k):
    L = {'u': 2 * i, 'v': 2 * j, 'w': 2 * k}
    base = {('u', 'v'): U, ('u', 'w'): V, ('v', 'w'): W}
    pairs = list(base)
    for g in itertools.product(*[range(min(L[a], L[b]) + 1) for a, b in pairs]):
        G = dict(zip(pairs, g))
        srt = sorted(g)
        if srt[0] != srt[1]:
            continue
        # leaves cannot coincide (base distances >= 1) and a leaf strictly inside another path is allowed
        ok = True
        for (a, b), gab in G.items():
            if gab == L[a] == L[b]:
                ok = False
        if not ok:
            continue
        # branching / degree check
        m = srt[0]; M = srt[2]
        top = [p for p in pairs if G[p] == M]
        nodes = {}
        def outgoing(depth, group):
            return sum(1 for a in group if L[a] > depth)
        if M == m:                       # all three share the first m steps, then separate
            nodes[m] = outgoing(m, ['u', 'v', 'w'])
        else:
            a, b = top[0]
            c = ({'u', 'v', 'w'} - {a, b}).pop()
            nodes[m] = (1 if (L[a] > m or L[b] > m) else 0) + (1 if L[c] > m else 0)
            nodes[M] = outgoing(M, [a, b])
        for dpt, out in nodes.items():
            if out + (1 if dpt > 0 else 0) > deg(dpt):
                ok = False
        if not ok:
            continue
        for (a, b), gab in G.items():
            dT = L[a] + L[b] - 2 * gab
            if dT <= 8:
                if dT != 2 * base[(a, b)]:
                    ok = False
            elif 2 * base[(a, b)] + dT < 16:
                ok = False
        if ok:
            return True
    return False

def main():
    # observed values on the Ree-Tits octagon
    obs = {}
    with open(os.path.join(OUT, 'triple_orbit_vectors.csv')) as fh:
        for row in csv.DictReader(fh):
            ty = tuple(int(x) for x in re.findall(r'\d', row['type']))
            for mm in re.finditer(r'\[(\d\d\d)\]=(\d+)', row['vector']):
                obs.setdefault(ty, {}).setdefault(mm.group(1), set()).add(int(mm.group(2)))

    say(f'order (s, t) = ({s_}, {t_})')
    for ty in sorted(obs):
        proved, zero_unproved, contradictions = [], [], []
        for key, vals in sorted(obs[ty].items()):
            i, j, k = map(int, key)
            if max(i, j, k) > 3 or min(i, j, k) == 0:
                continue
            f = feasible(*ty, i, j, k)
            if not f:
                proved.append(key)
                if vals != {0}:
                    contradictions.append(key)
            elif vals == {0}:
                zero_unproved.append(key)
        say(f'type {ty}: zero in every GO({s_},{t_}) by the tree test: {proved}; zero on all Ree-Tits orbits but not '
            f'excluded by the test: {zero_unproved}' + (f'; CONTRADICTION with data at {contradictions}' if contradictions else ''))
        assert not contradictions
    open(os.path.join(OUT, f'triple_zero_proofs_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')


if __name__ == '__main__':
    main()
