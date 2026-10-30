"""Explicit check: incidence graph of W(q)=GQ(q,q), q=2^m, vs array {q+1,q,q,q;1,1,1,q+1}. Usage: python verify_thin_explicit.py"""
import itertools, sys, json
from collections import deque

def gf(m):
    poly = {1: 0b11, 3: 0b1011}[m]; q = 1 << m
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a & q: a ^= poly
        return r
    return q, mul

def run(m):
    q, mul = gf(m)
    inv = {a: next(b for b in range(1, q) if mul(a, b) == 1) for a in range(1, q)}
    def norm(v):
        i = next(i for i, x in enumerate(v) if x); c = inv[v[i]]
        return tuple(mul(c, x) for x in v)
    pts = sorted({norm(v) for v in itertools.product(range(q), repeat=4) if any(v)})
    idx = {p: i for i, p in enumerate(pts)}
    form = lambda u, v: mul(u[0], v[1]) ^ mul(u[1], v[0]) ^ mul(u[2], v[3]) ^ mul(u[3], v[2])
    lines = set()
    for i, u in enumerate(pts):
        for v in pts[i+1:]:
            if form(u, v) == 0:
                L = frozenset(idx[norm(tuple(mul(a, x) ^ mul(b, y) for x, y in zip(u, v)))]
                              for a, b in itertools.product(range(q), repeat=2) if a or b)
                lines.add(L)
    lines = sorted(lines, key=sorted); n = len(pts); N = n + len(lines)
    adj = [[] for _ in range(N)]
    for j, L in enumerate(lines):
        for p in L: adj[p].append(n + j); adj[n + j].append(p)
    arr = {}
    for s in range(N):
        d = [-1] * N; d[s] = 0; Q = deque([s])
        while Q:
            x = Q.popleft()
            for y in adj[x]:
                if d[y] < 0: d[y] = d[x] + 1; Q.append(y)
        for x in range(N):
            c = sum(d[y] == d[x] - 1 for y in adj[x]); b = sum(d[y] == d[x] + 1 for y in adj[x])
            arr.setdefault(d[x], set()).add((b, c))
    return dict(q=q, vertices=N, points=n, lines=len(lines), array={k: sorted(v) for k, v in sorted(arr.items())})

for m in (1, 3):
    print(run(m))
