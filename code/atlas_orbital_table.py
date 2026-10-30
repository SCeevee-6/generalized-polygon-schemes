"""ATLAS-based orbital algebra of T = 2F4(2)' on 1755 points.
Recovered from the producing session (lineage of intersection_numbers_independent_verification.csv);
assertions added: transitivity and suborbit lengths 1,10,80,640,1024.
Requires internet access to brauer.maths.qmul.ac.uk (ATLAS v3 generators). Writes the CSV to the working directory.
"""
import urllib.request
import re

base = "https://brauer.maths.qmul.ac.uk/Atlas/exc/TF42/gap/TF42G1-p1755B0.g"
texts = {}
for i in [1,2]:
    url = base + str(i)
    with urllib.request.urlopen(url, timeout=30) as r:
        texts[i] = r.read().decode()

def parse_perm(text, n=1755):
    body = text.split(':=')[1]
    body = body.strip().rstrip(';').strip()
    body = re.sub(r'\s+', '', body)
    cycles = re.findall(r'\(([^)]+)\)', body)
    perm = list(range(n))
    for c in cycles:
        idxs = [int(t)-1 for t in c.split(',')]
        L = len(idxs)
        for k in range(L):
            perm[idxs[k]] = idxs[(k+1)%L]
    return perm

g1 = parse_perm(texts[1])
g2 = parse_perm(texts[2])

from collections import deque

n = 1755
gens = {'g1': g1, 'g2': g2}

ident = list(range(n))
visited = [False]*n
transporter = [None]*n
transporter[0] = ident
visited[0] = True
q = deque([0])

def compose(f, g):
    return [f[g[x]] for x in range(len(g))]

while q:
    p = q.popleft()
    cur_elem = transporter[p]
    for gname, garr in gens.items():
        newp = garr[p]
        if not visited[newp]:
            visited[newp] = True
            transporter[newp] = compose(garr, cur_elem)
            q.append(newp)

def inv_perm(p):
    n = len(p)
    inv = [0]*n
    for x in range(n):
        inv[p[x]] = x
    return inv

inv_transporter = [inv_perm(t) for t in transporter]

schreier_gens = []
for p_pt in range(n):
    tp = transporter[p_pt]
    for gname, garr in gens.items():
        gp = garr[p_pt]
        tmp = compose(garr, tp)
        s = compose(inv_transporter[gp], tmp)
        schreier_gens.append(s)

parent = list(range(n))
def find(x):
    while parent[x]!=x:
        parent[x]=parent[parent[x]]
        x=parent[x]
    return x
def union(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb:
        parent[ra]=rb

for s in schreier_gens:
    for x in range(n):
        union(x, s[x])

from collections import defaultdict
orbits = defaultdict(list)
for x in range(n):
    orbits[find(x)].append(x)

assert all(visited), 'action not transitive'
orbit_list = sorted(orbits.values(), key=len)
assert [len(o) for o in orbit_list] == [1, 10, 80, 640, 1024], 'unexpected suborbit lengths'

R = [set(o) for o in orbit_list]

reps = [next(iter(R[k])) for k in range(5)]

p = {}
for k in range(5):
    yk = reps[k]
    g = transporter[yk]
    for j in range(5):
        Rj_yk = set(g[x] for x in R[j])
        for i in range(5):
            cnt = len(R[i] & Rj_yk)
            p[(i,j,k)] = cnt

import csv
with open('intersection_numbers_independent_verification.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['i','j','k','p^k_ij'])
    for i in range(5):
        for j in range(5):
            for kk in range(5):
                w.writerow([i,j,kk,p[(i,j,kk)]])