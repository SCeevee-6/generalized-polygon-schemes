"""Orbits of T = 2F4(2)' and of 2F4(2) on ordered triples of points of GO(2,4), their triple intersection
numbers, and a comparison with the LP ranges of go24_triple_lp.py. Also the dimension of the Terwilliger
algebra T(x) of the collinearity graph at a base point x, against the centraliser algebra of the point
stabiliser.

Orbits on triples (u, v, w) are orbits of the stabiliser T_u on ordered pairs (v, w), u = point 0. The
stabiliser generators are Schreier generators whose group order is checked (10240 for T, 20480 for T.2).
For every orbit the vector [i j k] = #{z : d(u,z) = i, d(v,z) = j, d(w,z) = k} is computed exactly.
Usage: python triple/triple_orbits.py     (writes triple/out_orbits/*)
"""
import os, sys, json, itertools, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'flags')); sys.path.insert(0, HERE)
import go24_common as gc
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components

OUT = os.path.join(HERE, 'out_orbits'); os.makedirs(OUT, exist_ok=True)
N = gc.N
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)

Tg, T2g = gc.generators()
D = gc.point_distances(Tg)
S_T = gc.stabiliser_gens(Tg, N, order=10240)
S_T2 = gc.stabiliser_gens(T2g, N, order=20480)
say(f'T_u: {len(S_T)} generators, order 10240 checked; (T.2)_u: {len(S_T2)} generators, order 20480 checked')

def pair_orbits(S):
    idx = np.arange(N * N, dtype=np.int64)
    v, w = idx // N, idx % N
    rows = np.concatenate([idx] * len(S)).astype(np.int32)
    cols = np.concatenate([g[v] * N + g[w] for g in S]).astype(np.int32)
    A = coo_matrix((np.ones(len(rows), np.bool_), (rows, cols)), shape=(N * N, N * N))
    nc, lab = connected_components(A, directed=True, connection='weak')
    return nc, lab.astype(np.int32)

t0 = time.time()
nT, labT = pair_orbits(S_T)
nT2, labT2 = pair_orbits(S_T2)
say(f'T-orbits on ordered triples: {nT}; 2F4(2)-orbits: {nT2} ({time.time() - t0:.0f} s)')

reps = np.full(nT, -1, np.int64)
order = np.argsort(labT, kind='stable')
first = np.searchsorted(labT[order], np.arange(nT))
reps = order[first]
sizes = np.bincount(labT, minlength=nT)
t2_of = labT2[reps]                       # 2F4(2)-orbit of each T-orbit

def triple_vector(u, v, w):
    return np.bincount(25 * D[u].astype(np.int64) + 5 * D[v] + D[w], minlength=125).reshape(5, 5, 5)

orbs = []
for o in range(nT):
    v, w = divmod(int(reps[o]), N)
    ty = (int(D[0, v]), int(D[0, w]), int(D[v, w]))
    orbs.append({'orbit': o, 'rep': (0, v, w), 'type': ty, 'size': int(sizes[o]), 'orbit_2F42': int(t2_of[o]),
                 'vec': triple_vector(0, v, w)})

# ---- compare with the LP of go24_triple_lp.py (types with 1 <= U <= V <= W)
import go24_triple_lp as lp
types = [(x, y, z) for x in range(1, 5) for y in range(x, 5) for z in range(y, 5) if lp.P(x, y, z) > 0]
rows_out = []
summary = []
for ty in types:
    keys, A, bv = lp.system(*ty)
    lo, hi = lp.ranges(A, bv, len(keys))
    free = len(keys) - np.linalg.matrix_rank(A)
    sel = [o for o in orbs if o['type'] == ty]
    vecs = [tuple(int(o['vec'][k]) for k in keys) for o in sel]
    distinct = sorted(set(vecs))
    inside = all(lo[i] - 1e-6 <= vv[i] <= hi[i] + 1e-6 for vv in vecs for i in range(len(keys)))
    # equations satisfied exactly?
    eq_ok = all(np.allclose(A @ np.array(vv, float), bv, atol=1e-6) for vv in distinct)
    span = np.linalg.matrix_rank(np.array([np.array(vv) - np.array(distinct[0]) for vv in distinct], float)) if len(distinct) > 1 else 0
    ntwo = len({o['orbit_2F42'] for o in sel})
    # which coordinates take more than one value across the orbits
    varying = [keys[i] for i in range(len(keys)) if len({vv[i] for vv in vecs}) > 1]
    # coordinates where all orbits sit at an LP bound
    at_bound = [keys[i] for i in range(len(keys)) if hi[i] - lo[i] > 1e-6 and all(abs(vv[i] - lo[i]) < 1e-6 or abs(vv[i] - hi[i]) < 1e-6 for vv in vecs)]
    summary.append({'type': ty, 'unknowns': len(keys), 'free_dim': int(free), 'T_orbits': len(sel), '2F42_orbits': ntwo,
                    'distinct_vectors': len(distinct), 'affine_span_of_vectors': int(span), 'all_in_LP_ranges': bool(inside),
                    'equations_hold': bool(eq_ok), 'coords_varying': [list(k) for k in varying],
                    'coords_at_LP_bound_nontrivial': [list(k) for k in at_bound],
                    'orbit_sizes': [o['size'] for o in sel]})
    for o, vv in zip(sel, vecs):
        rows_out.append((ty, o['orbit'], o['orbit_2F42'], o['size'], o['rep'], vv, keys))
    say(f'type {ty}: {len(keys)} unknowns, LP free dim {free}; T-orbits {len(sel)} (sizes {[o["size"] for o in sel]}), '
        f'2F4(2)-orbits {ntwo}; distinct triple vectors {len(distinct)} spanning affine dim {span}; '
        f'in LP ranges {inside}; equations hold {eq_ok}; varying coords {len(varying)}')

# all orbits, including degenerate triples, by ordered type
by_type = {}
for o in orbs:
    by_type.setdefault(o['type'], []).append(o['size'])
say('T-orbits on ordered triples (u,v,w) by ordered distance type (d(u,v), d(u,w), d(v,w)):')
for ty in sorted(by_type):
    say(f'  {ty}: {len(by_type[ty])} orbits, sizes {sorted(by_type[ty])}')

# ---- Terwilliger algebra at the base point 0, in orbit coordinates
Arow = np.array([np.nonzero(D[x] == 1)[0] for x in range(N)])      # 10 neighbours each
A = csr_matrix((np.ones(N * 10), (np.repeat(np.arange(N), 10), Arow.ravel())), shape=(N, N))
d0 = D[0].astype(np.int64)

def terwilliger_exact(lab, nor, reps_, prime=2147483629):
    """T(x) = <A, E*_0..E*_4> inside the centraliser algebra of T_x, in orbit coordinates.

    An element is a vector c over the orbits of T_x on pairs; (A M)[v,w] = sum_{z ~ v} M[z,w] and
    (E*_i M)[v,w] = [d(x,v) = i] M[v,w], so both act on c through the orbit labels of (z, w), z ~ v, at the
    orbit representatives. Words in the generators are explored from the identity; a word is kept when it
    is independent mod `prime`, and its exact integer vector is stored. Then, over Q (python-flint), the kept
    vectors are checked to be independent and closed under left multiplication by every generator, which
    makes their span exactly T(x)."""
    import flint
    rv, rw = reps_ // N, reps_ % N
    nbo = [[int(lab[int(z) * N + int(rw[o])]) for z in Arow[rv[o]]] for o in range(nor)]
    dist0 = [int(D[0, rv[o]]) for o in range(nor)]
    gens = ['A'] + [f'E{i}' for i in range(5)]
    def left(g, c):
        if g == 'A':
            return [sum(c[j] for j in nbo[o]) for o in range(nor)]
        i = int(g[1:])
        return [c[o] if dist0[o] == i else 0 for o in range(nor)]
    basis, piv, exact = [], [], []
    def add(vec):
        v = np.array([x % prime for x in vec], dtype=np.int64)
        for b, pc in zip(basis, piv):
            if v[pc]:
                v = (v - v[pc] * b) % prime
        nz = np.nonzero(v)[0]
        if len(nz) == 0:
            return False
        pc = int(nz[0]); v = (v * pow(int(v[pc]), prime - 2, prime)) % prime
        for k in range(len(basis)):
            if basis[k][pc]:
                basis[k] = (basis[k] - basis[k][pc] * v) % prime
        basis.append(v); piv.append(pc); exact.append(list(vec))
        return True
    I = [1 if rv[o] == rw[o] else 0 for o in range(nor)]
    add(I); queue = [I]
    while queue:
        c = queue.pop()
        for g in gens:
            nc_ = left(g, c)
            if add(nc_):
                queue.append(nc_)
    dim = len(exact)
    Bm = flint.fmpz_mat(exact)
    assert Bm.rank() == dim
    closed = True
    for c in exact:
        for g in gens:
            if flint.fmpz_mat(exact + [left(g, c)]).rank() != dim:
                closed = False
    return dim, closed


t0 = time.time()
dimT, closed = terwilliger_exact(labT.astype(np.int64), nT, reps)
say(f'Terwilliger algebra T(x) at a point: dimension {dimT} (exact over Q, closure verified: {closed}); '
    f'centraliser algebra of T_x has dimension {nT}, of 2F4(2)_x {nT2} ({time.time() - t0:.0f} s)')

json.dump({'n_T_orbits': int(nT), 'n_2F42_orbits': int(nT2), 'terwilliger_dim': int(dimT), 'terwilliger_closure_verified': bool(closed), 'types': summary,
           'by_ordered_type': {str(k): v for k, v in by_type.items()}},
          open(os.path.join(OUT, 'triple_orbits.json'), 'w'), indent=1)
with open(os.path.join(OUT, 'triple_orbit_vectors.csv'), 'w', newline='\n') as f:
    f.write('type,T_orbit,2F42_orbit,size,rep,vector\n')
    for ty, o, o2, sz, rep, vv, keys in rows_out:
        f.write(f'"{ty}",{o},{o2},{sz},"{rep}","' + ' '.join(f'[{a}{b}{c}]={x}' for (a, b, c), x in zip(keys, vv)) + '"\n')
open(os.path.join(OUT, 'triple_orbits_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
