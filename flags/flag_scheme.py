"""Orbital configuration of T = 2F4(2)' and of 2F4(2) on the 8775 flags of the Ree-Tits octagon GO(2,4).

Steps
 1. Points, lines (triangles of the collinearity graph) and flags (p, L) from the ATLAS 1755-point action.
 2. Flag action of T and T.2, stabiliser of flag 0 (order checked), T-orbitals and their exact p^k_ij.
 3. Weyl distance of every orbital: the 16 relations A_w = A_{s1}...A_{sk} along reduced words in
    D_16 = <sigma, tau>, sigma = change the line through the point, tau = change the point on the line.
 4. T.2-orbitals as unions of T-orbitals.
 5. The 16 Weyl-distance relations against the Iwahori-Hecke algebra of D_16 with q_sigma = t = 4,
    q_tau = s = 2 (structure constants compared exactly).
 6. Centre and Wedderburn decomposition (exact, sympy), degrees of the corresponding characters of T.
 7. All fusions (algebraic WL closure of every union of classes, then closure under meets).
 8. How each orbital projects to the point scheme (1755) and to the line scheme (2925).
Usage: python flags/flag_scheme.py      (writes flags/out/flag_scheme_output.txt and CSV/JSON files)
"""
import os, sys, json, itertools, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'tits1755'))
import go24_common as gc
from coset_schemes import orbitals, intersection_numbers

OUT = os.path.join(HERE, 'out'); os.makedirs(OUT, exist_ok=True)
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)


# ---------------------------------------------------------------- 1. geometry
Tg, T2g = gc.generators()
D = gc.point_distances(Tg)
lines = gc.lines_from_distances(D)
flags = [(p, li) for li, L in enumerate(lines) for p in L]
fidx = {f: i for i, f in enumerate(flags)}
lidx = {L: i for i, L in enumerate(lines)}
NF = len(flags)
say(f'points {gc.N}, lines {len(lines)}, flags {NF}')

def line_perm(g):
    return np.array([lidx[tuple(sorted(int(g[x]) for x in L))] for L in lines])

def flag_perm(g):
    lp = line_perm(g)
    return np.array([fidx[(int(g[p]), int(lp[li]))] for (p, li) in flags])

FT = [flag_perm(g) for g in Tg]
FT2 = [flag_perm(g) for g in T2g]

# ---------------------------------------------------------------- 2. T-orbitals
H = gc.stabiliser_gens(FT, NF, order=17971200 // NF)
say(f'flag stabiliser in T: order {17971200 // NF} certified with {len(H)} Schreier generators')
t0 = time.time()
O, val, reps, _ = orbitals(FT, NF)          # O[x,y] = T-orbital of (x,y), int8
r = len(reps)
p = intersection_numbers(O, reps)          # p[k,i,j] = p^k_ij
tr = [int(O[reps[i], 0]) for i in range(r)]   # transpose class
say(f'T on flags: rank {r}, suborbit lengths {val}, {time.time() - t0:.0f} s')

# ---------------------------------------------------------------- 3. Weyl distance
from scipy.sparse import csr_matrix
pt_of = np.array([f[0] for f in flags]); ln_of = np.array([f[1] for f in flags])
def adj(same):
    rows, cols = [], []
    groups = {}
    for i, key in enumerate(same):
        groups.setdefault(int(key), []).append(i)
    for g in groups.values():
        for a in g:
            for b in g:
                if a != b:
                    rows.append(a); cols.append(b)
    return csr_matrix((np.ones(len(rows), np.int64), (rows, cols)), shape=(NF, NF))
A_sig, A_tau = adj(pt_of), adj(ln_of)      # sigma: same point other line; tau: same line other point
QS = {'s': 4, 't': 2}                       # q_sigma = t = 4, q_tau = s = 2
other = {'s': 't', 't': 's'}

def W_elements():
    els = [('e',)]
    for f in 'st':
        for l in range(1, 8):
            els.append((f, l))
    els.append(('w0',))
    return els
WEL = W_elements()
def wname(w):
    if w == ('e',): return 'e'
    if w == ('w0',): return 'w0'
    f, l = w
    return ''.join(f if i % 2 == 0 else other[f] for i in range(l))
def rmul(w, x):
    """(w, x) -> (wx, length went up?)"""
    if w == ('e',): return (x, 1), True
    if w == ('w0',): return (other[x], 7), False
    f, l = w
    last = f if l % 2 == 1 else other[f]
    if last == x:
        return (('e',) if l == 1 else (f, l - 1)), False
    return (('w0',) if l == 7 else (f, l + 1)), True
def word_of(w):
    return 'stststst' if w == ('w0',) else ('' if w == ('e',) else wname(w))

weyl = np.full(NF, -1)
for wi, w in enumerate(WEL):
    v = np.zeros(NF, np.int64); v[0] = 1
    for ch in word_of(w):
        v = (A_sig if ch == 's' else A_tau).T @ v
    sup = np.nonzero(v)[0]
    assert set(np.unique(v[sup]).tolist()) == {1}, 'reduced gallery not unique'
    assert (weyl[sup] == -1).all()
    weyl[sup] = wi
assert (weyl >= 0).all()
cls_w = []
for k in range(r):
    ws = set(weyl[O[0] == k].tolist())
    assert len(ws) == 1
    cls_w.append(ws.pop())
wsize = {wi: int((weyl == wi).sum()) for wi in range(16)}
split = {wname(WEL[wi]): [k for k in range(r) if cls_w[k] == wi] for wi in range(16)}
say('Weyl distance of the T-orbitals (cell size -> T-suborbit lengths):')
for wi, w in enumerate(WEL):
    ks = [k for k in range(r) if cls_w[k] == wi]
    say(f'  {wname(w):9s} len {0 if w == ("e",) else 8 if w == ("w0",) else w[1]}  |cell| {wsize[wi]:5d}  '
        f'T-classes {ks}  lengths {[val[k] for k in ks]}')

# ---------------------------------------------------------------- 4. T.2 orbitals
comp = list(range(r))
def find(a):
    while comp[a] != a: a = comp[a]
    return a
for g in FT2:
    for k in range(r):
        kk = int(O[g[0], g[reps[k]]])
        a, b = find(k), find(kk)
        if a != b: comp[max(a, b)] = min(a, b)
t2cls = sorted({tuple(k for k in range(r) if find(k) == find(j)) for j in range(r)})
say(f'T.2 on flags: rank {len(t2cls)}; classes {t2cls}')
same_as_weyl = sorted(tuple(v) for v in split.values()) == t2cls
say(f'T.2-orbitals = Weyl-distance relations: {same_as_weyl}')

# ---------------------------------------------------------------- 5. Hecke algebra check
def fuse(p, blocks):
    b = len(blocks)
    q = np.zeros((b, b, b), np.int64)
    for c, C in enumerate(blocks):
        h = C[0]
        for a, A in enumerate(blocks):
            for bb, B in enumerate(blocks):
                q[c, a, bb] = sum(p[h, i, j] for i in A for j in B)
    return q
wblocks = [split[wname(w)] for w in WEL]
pw = fuse(p, wblocks)              # pw[k,i,j], indices = WEL order

def hecke_product(u, v):
    """T_u T_v in the basis WEL, as integer vector."""
    vec = {u: 1}
    for x in word_of(v):
        new = {}
        for w, cf in vec.items():
            wx, up = rmul(w, x)
            if up:
                new[wx] = new.get(wx, 0) + cf
            else:
                new[wx] = new.get(wx, 0) + QS[x] * cf
                new[w] = new.get(w, 0) + (QS[x] - 1) * cf
        vec = new
    out = np.zeros(16, np.int64)
    for w, cf in vec.items():
        out[WEL.index(w)] += cf
    return out
hk_ok = all(np.array_equal(hecke_product(u, v), pw[:, i, j]) for i, u in enumerate(WEL) for j, v in enumerate(WEL))
say(f'Weyl-distance scheme = Iwahori-Hecke algebra H(D16; q_s=4, q_t=2), all 4096 structure constants equal: {hk_ok}')

# ---------------------------------------------------------------- 6. centre, Wedderburn
import sympy
def regular_matrices(p):
    rr = p.shape[0]
    # left multiplication by A_i on basis A_j: A_i A_j = sum_k p[k,i,j] A_k
    return [sympy.Matrix(rr, rr, lambda k, j: int(p[k, i, j])) for i in range(rr)]
def wedderburn(p, nverts, label):
    rr = p.shape[0]
    Ls = regular_matrices(p)
    rows = []
    for j in range(rr):
        for k in range(rr):
            rows.append([int(p[k, i, j] - p[k, j, i]) for i in range(rr)])
    Z = sympy.Matrix(rows).nullspace()     # centre: coefficient vectors c with sum c_i A_i central
    zdim = len(Z)
    # generic central element and its action on the algebra
    rng = np.random.default_rng(7)
    z = sum((int(rng.integers(1, 50)) * v for v in Z), sympy.zeros(rr, 1))
    Lz = sum((z[i] * Ls[i] for i in range(rr)), sympy.zeros(rr, rr))
    lam = sympy.symbols('x')
    cp = sympy.factor_list(Lz.charpoly(lam).as_expr())
    comps, roots = [], []
    for fac, e in cp[1]:
        comps.append((str(fac), int(sympy.degree(fac, lam)), int(e)))
        for rt in sympy.roots(sympy.Poly(fac, lam)):
            roots.append((rt, int(e)))
    assert sum(e for _, e in roots) == rr
    res = []
    for rt, e in roots:
        m = int(round(e ** 0.5)); assert m * m == e
        v = sympy.zeros(rr, 1); v[0] = 1
        for mu, _ in roots:
            if mu != rt:
                v = ((Lz * v - mu * v) / (rt - mu)).applyfunc(sympy.expand)
        chi1 = sympy.nsimplify(sympy.radsimp(sympy.simplify(nverts * v[0] / m)))
        assert chi1.is_integer, chi1
        res.append((m, int(chi1), rt))
    res.sort(key=lambda t: (t[1], t[0]))
    say(f'{label}: centre dimension {zdim}; components (m_i, chi(1)): {[(m, c) for m, c, _ in res]}; '
        f'sum m^2 = {sum(m * m for m, _, _ in res)}; sum m chi(1) = {sum(m * c for m, c, _ in res)}')
    say(f'  char poly of a generic central element factors over Q as {[(f, d, e) for f, d, e in comps]}')
    return zdim, [(m, c) for m, c, _ in res]
wT = wedderburn(p, NF, 'T flag algebra (rank %d)' % r)
wG = wedderburn(pw, NF, '2F4(2) flag algebra = Hecke algebra (rank 16)')
comm = all(p[k, i, j] == p[k, j, i] for k in range(r) for i in range(r) for j in range(r))
say(f'T flag algebra commutative: {comm}; symmetric classes: {sum(1 for k in range(r) if tr[k] == k)} of {r}')

# ---------------------------------------------------------------- 7. fusions
def closure(p, tr, part):
    """Coarsest fusion refining the partition `part` (list of block ids per class, class 0 alone)."""
    rr = p.shape[0]
    part = list(part)
    while True:
        ids = sorted(set(part)); m = {b: i for i, b in enumerate(ids)}
        M = np.zeros((rr, len(ids)), np.int64)
        for k in range(rr): M[k, m[part[k]]] = 1
        S = np.einsum('ia,kij,jb->kab', M, p, M)
        sig = [(part[k], part[tr[k]], S[k].tobytes()) for k in range(rr)]
        lab = {}
        new = [lab.setdefault(s, len(lab)) for s in sig]
        if len(set(new)) == len(set(part)):
            return tuple(canon(new))
        part = new
def canon(part):
    m = {}
    return [m.setdefault(b, len(m)) for b in part]
def fusions(p, tr):
    rr = p.shape[0]
    found = set()
    nontriv = list(range(1, rr))
    for mask in range(1 << (rr - 2)):          # S always contains class 1 (S and its complement give the same closure)
        sel = [k for i, k in enumerate(nontriv[1:]) if (mask >> i) & 1]
        part = [0] + [1 if (k == 1 or k in sel) else 2 for k in range(1, rr)]
        found.add(closure(p, tr, part))
    found.add(tuple([0] + [1] * (rr - 1)))
    # meets
    while True:
        new = set()
        F = list(found)
        for a, b in itertools.combinations(F, 2):
            c = closure(p, tr, canon(list(zip(a, b))))
            if c not in found: new.add(c)
        if not new: break
        found |= new
    return found
def check_fusion(p, tr, part):
    blocks = {}
    for k, b in enumerate(part): blocks.setdefault(b, []).append(k)
    bl = list(blocks.values())
    for C in bl:
        if sorted(set(part[tr[k]] for k in C)).__len__() != 1: return False
    q = fuse(p, bl)
    for c, C in enumerate(bl):
        for h in C:
            for a, A in enumerate(bl):
                for b, B in enumerate(bl):
                    if sum(p[h, i, j] for i in A for j in B) != q[c, a, b]: return False
    return True
def describe(p, tr, part, valv):
    blocks = {}
    for k, b in enumerate(part): blocks.setdefault(b, []).append(k)
    bl = sorted(blocks.values())
    q = fuse(p, bl)
    rk = len(bl)
    commf = bool(np.array_equal(q, q.transpose(0, 2, 1)))
    symf = all(part[tr[k]] == part[k] for k in range(len(part)))
    return {'rank': rk, 'blocks': bl, 'valencies': [int(sum(valv[k] for k in B)) for B in bl],
            'commutative': commf, 'symmetric': symf}
t0 = time.time()
FTf = sorted(fusions(p, tr), key=lambda f: -len(set(f)))
assert all(check_fusion(p, tr, f) for f in FTf)
trw = [int(np.nonzero(pw[0, i, :])[0][0]) for i in range(16)]
FGf = sorted(fusions(pw, trw), key=lambda f: -len(set(f)))
assert all(check_fusion(pw, trw, f) for f in FGf)
valw = [wsize[i] for i in range(16)]
descT = [describe(p, tr, f, val) for f in FTf]
descG = [describe(pw, trw, f, valw) for f in FGf]
say(f'fusions of the T flag configuration: {len(FTf)} (including the configuration itself and the trivial one), {time.time() - t0:.0f} s')
for d in descT:
    say(f'  rank {d["rank"]:2d} comm {d["commutative"]!s:5s} symm {d["symmetric"]!s:5s} valencies {d["valencies"]}')
say(f'fusions of the 2F4(2) flag scheme (Hecke algebra): {len(FGf)}')
for d in descG:
    bl = [[wname(WEL[i]) for i in B] for B in d['blocks']]
    say(f'  rank {d["rank"]:2d} comm {d["commutative"]!s:5s} symm {d["symmetric"]!s:5s} blocks {bl}')

# ---------------------------------------------------------------- 8. projections to points and lines
def line_dist(a, b):
    if a == b: return 0
    return 1 + int(min(D[x, y] for x in lines[a] for y in lines[b]))
proj = []
for k in range(r):
    y = reps[k]
    proj.append((wname(WEL[cls_w[k]]), int(D[flags[0][0], flags[y][0]]), line_dist(flags[0][1], flags[y][1])))
say('orbital -> (Weyl distance, point distance, line distance in the dual octagon):')
for k in range(r):
    say(f'  class {k:2d} length {val[k]:5d} transpose {tr[k]:2d}  {proj[k]}')

json.dump({'rank': r, 'suborbits': [int(v) for v in val], 'transpose': tr, 'p': p.tolist(),
           'weyl': [wname(WEL[c]) for c in cls_w], 'projection': proj,
           'T2_classes': [list(c) for c in t2cls],
           'fusions_T': descT, 'fusions_2F42': [{**d, 'blocks': [[wname(WEL[i]) for i in B] for B in d['blocks']]} for d in descG],
           'wedderburn_T': wT, 'wedderburn_2F42': wG},
          open(os.path.join(OUT, 'flag_scheme.json'), 'w'), indent=1, default=int)
with open(os.path.join(OUT, 'flag_intersection_numbers.csv'), 'w', newline='\n') as f:
    f.write('k,i,j,p_k_ij\n')
    for k, i, j in itertools.product(range(r), repeat=3):
        if p[k, i, j]:
            f.write(f'{k},{i},{j},{p[k, i, j]}\n')
open(os.path.join(OUT, 'flag_scheme_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
