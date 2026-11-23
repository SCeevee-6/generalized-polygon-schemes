"""Ree-Tits octagon O(2) from ATLAS generators of T = 2F4(2)' on 1755 points (independent rebuild)."""
import re, numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
import os
GG = os.environ.get('ATLAS_DIR', '.')  # folder with TF42G1-p1755B0.g1 and .g2 (ATLAS v3)
def parse_perm(path, n=1755):
    t = open(path).read(); t = t[t.index(':=')+2:].replace('\n','').replace(' ','').rstrip(';')
    p = np.arange(n)
    for cyc in re.findall(r'\(([\d,]+)\)', t):
        c = [int(x)-1 for x in cyc.split(',')]
        for i,a in enumerate(c): p[a] = c[(i+1)%len(c)]
    return p
def gens():
    return [parse_perm(os.path.join(GG, 'TF42G1-p1755B0.g%d' % i)) for i in (1,2)]
def orbitals(G, n=1755):
    idx = np.arange(n*n); a, b = idx//n, idx%n
    r = np.concatenate([idx]*len(G)); c = np.concatenate([g[a]*n+g[b] for g in G])
    _, lab = connected_components(csr_matrix((np.ones(len(r),np.int8),(r,c)),shape=(n*n,)*2), directed=True, connection='weak')
    return lab.reshape(n,n)
def octagon():
    G = gens(); L = orbitals(G); n = 1755
    vals = {l: int((L[0]==l).sum()) for l in np.unique(L[0])}
    col = [l for l,v in vals.items() if v==10]; assert len(col)==1
    A = (L==col[0]); assert (A==A.T).all() and (A.sum(1)==10).all()
    D = shortest_path(csr_matrix(A), unweighted=True).astype(np.int8)
    lines = set()
    for a,b in zip(*np.nonzero(np.triu(A))):
        c = np.flatnonzero(A[a]&A[b]); assert len(c)==1
        lines.add(tuple(sorted((int(a),int(b),int(c[0])))))
    lines = np.array(sorted(lines)); assert len(lines)==2925
    return G, A, D, lines
