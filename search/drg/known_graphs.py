"""Two graphs with arrays that no primitive group carries (see drg_note.md), built directly.

Br(q), the Brouwer graph: vertices F_q^3 x F_q^3, (a, b) ~ (c, d) iff a != c and b - d = a x c (cross product).
Coset graph of the binary cyclic code of length 2^m - 1 with zeros alpha and alpha^(2^j + 1): the Cayley graph
on F_2^m x F_2^m with connection set {(x, x^(2^j+1)) : x in F_2^m, x != 0}.

For each graph the distance partition from one vertex is computed and checked to be equitable with the stated
array at every vertex (the graphs are vertex-transitive, but every vertex is checked anyway). Adjacency lists
are written to OUT_DIR/<name>.g for GAP/GRAPE (see aut_check.g).

Usage: python known_graphs.py [OUT_DIR]
"""
import os, sys, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def intersection_array(nb):
    """nb: list of neighbour lists. Returns (b, c) if distance-regular, else None."""
    n = len(nb); arr = None
    for s in range(n):
        dist = [-1] * n; dist[s] = 0; fr = [s]; h = 0
        while fr:
            nxt = []
            for x in fr:
                for y in nb[x]:
                    if dist[y] < 0:
                        dist[y] = h + 1; nxt.append(y)
            fr = nxt; h += 1
        d = max(dist)
        if min(dist) < 0:
            return None
        b = [None] * d; c = [None] * d
        for x in range(n):
            i = dist[x]
            ci = sum(1 for y in nb[x] if dist[y] == i - 1); bi = sum(1 for y in nb[x] if dist[y] == i + 1)
            if i < d:
                if b[i] is None: b[i] = bi
                elif b[i] != bi: return None
            if i > 0:
                if c[i - 1] is None: c[i - 1] = ci
                elif c[i - 1] != ci: return None
        if arr is None: arr = (b, c)
        elif arr != (b, c): return None
    return arr


def brouwer(q):
    assert q in (2, 3, 5, 7)
    V = list(itertools.product(itertools.product(range(q), repeat=3), repeat=2))
    idx = {v: i for i, v in enumerate(V)}
    def cross(a, c):
        return ((a[1] * c[2] - a[2] * c[1]) % q, (a[2] * c[0] - a[0] * c[2]) % q, (a[0] * c[1] - a[1] * c[0]) % q)
    nb = []
    for (a, b) in V:
        L = []
        for c in itertools.product(range(q), repeat=3):
            if c == a: continue
            x = cross(a, c)
            d = tuple((b[t] - x[t]) % q for t in range(3))
            L.append(idx[(c, d)])
        nb.append(L)
    return nb


def gf2m(m):
    poly = {5: 0b100101}[m]        # x^5 + x^2 + 1, primitive
    exp = [0] * (2 ** m - 1); x = 1
    for i in range(2 ** m - 1):
        exp[i] = x; x <<= 1
        if x >> m: x ^= poly
    log = {e: i for i, e in enumerate(exp)}
    def mul(a, b):
        if a == 0 or b == 0: return 0
        return exp[(log[a] + log[b]) % (2 ** m - 1)]
    def pw(a, e):
        r = 1
        for _ in range(e): r = mul(r, a)
        return r
    return mul, pw


def coset_graph(m, j):
    mul, pw = gf2m(m)
    Q = 2 ** m
    S = [(x, pw(x, 2 ** j + 1)) for x in range(1, Q)]
    nb = []
    for u in range(Q):
        for w in range(Q):
            nb.append([((u ^ s0) * Q + (w ^ s1)) for s0, s1 in S])
    return nb


def write_gap(nb, path, name):
    with open(path, 'w', newline='\n') as f:
        f.write(f'# {name}: adjacency lists, vertices 1..{len(nb)}\n')
        f.write('ADJ := [\n' + ',\n'.join('[' + ','.join(str(y + 1) for y in L) + ']' for L in nb) + '];;\n')


def main(out=None):
    out = out or os.path.join(HERE, 'known')
    os.makedirs(out, exist_ok=True)
    for name, nb in (('brouwer3', brouwer(3)), ('coset_m5_j1', coset_graph(5, 1)), ('coset_m5_j2', coset_graph(5, 2))):
        arr = intersection_array(nb)
        print(name, len(nb), 'array', arr)
        write_gap(nb, os.path.join(out, name + '.g'), name)


if __name__ == '__main__':
    main(*sys.argv[1:2])
