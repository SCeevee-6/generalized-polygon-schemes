"""Wedderburn decomposition of the Terwilliger algebra T(x) of the collinearity graph of GO(2,4) at a point x,
its irreducible modules (dimension, endpoint, dual endpoint, diameters, thinness, multiplicity in the standard
module C^1755), and the decomposition of the centraliser algebra C = End_{H}(C^1755), H = 2F4(2)_x, with the
restriction of each simple C-module to T(x).

Everything is exact over Q (python-flint):
  * elements of C are vectors over the 258 H-orbits on ordered pairs; (a b)_o = sum_z a[lab(v_o, z)] b[lab(z, w_o)]
    at the orbit representative (v_o, w_o);
  * T(x) = <A, E*_0..E*_4> is spanned by words (independence mod p, then exact RREF over Q, closure checked);
  * the centre Z of each algebra is a nullspace; a random integer central element z has char. polynomial
    prod_W (X - lambda_W)^{d_W^2} on the left regular module, so each Q-irreducible factor f with exponent e gives
    deg(f) Galois-conjugate simple components of size d = sqrt(e) (separation checked: sum deg f = dim Z);
  * for a component with block space B_f = ker f(L_z): dim E*_i B_f = deg(f) d r_i, where r_i = dim E*_i W;
    likewise r*_j = dim E_j W with E_j the primitive idempotents of A; trace_V(e_f) = deg(f) d m with m the
    multiplicity of W in C^1755.
Usage: python triple/terwilliger_wedderburn.py   (writes triple/out_orbits/terwilliger_wedderburn_output.txt, .csv)
"""
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'flags'))
import numpy as np
import flint
import go24_common as gc
from fractions import Fraction
from collections import defaultdict
from math import isqrt, gcd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

OUT = os.path.join(HERE, 'out_orbits'); os.makedirs(OUT, exist_ok=True)
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)
random.seed(20261122)
N = gc.N
t0 = time.time()
Tg, T2g = gc.generators()
D = gc.point_distances(Tg)
S = gc.stabiliser_gens(T2g, N, order=20480)

def pair_orbits(S):
    idx = np.arange(N * N, dtype=np.int64)
    v, w = idx // N, idx % N
    rows = np.concatenate([idx] * len(S)).astype(np.int32)
    cols = np.concatenate([g[v] * N + g[w] for g in S]).astype(np.int32)
    A = coo_matrix((np.ones(len(rows), np.bool_), (rows, cols)), shape=(N * N, N * N))
    nc, lab = connected_components(A, directed=True, connection='weak')
    return nc, lab.astype(np.int32)

nor, lab = pair_orbits(S)
order = np.argsort(lab, kind='stable')
reps = order[np.searchsorted(lab[order], np.arange(nor))]
sizes = np.bincount(lab, minlength=nor)
rv, rw = reps // N, reps % N
del order
L = np.stack([lab[int(rv[o]) * N + np.arange(N)] for o in range(nor)])      # lab(v_o, z)
R = np.stack([lab[np.arange(N) * N + int(rw[o])] for o in range(nor)])      # lab(z, w_o)
del lab
diag = [o for o in range(nor) if rv[o] == rw[o]]
d0v = D[0, rv].astype(int); d0w = D[0, rw].astype(int); dvw = D[rv, rw].astype(int)
say(f'H = 2F4(2)_x: {nor} orbits on ordered pairs ({len(diag)} diagonal) ({time.time() - t0:.0f} s)')

def trace_prod(a, b):
    """trace on C^1755 of the product a b (only diagonal orbits are needed)."""
    return sum(int(sizes[o]) * sum(a[l] * b[r_] for l, r_ in zip(L[o].tolist(), R[o].tolist())) for o in diag)

def mul(a, b):
    """exact product of two elements of C (lists of python ints or Fractions)."""
    ao = np.array(a, dtype=object); bo = np.array(b, dtype=object)
    return list((ao[L] * bo[R]).sum(axis=1))

Avec = [1 if dvw[o] == 1 else 0 for o in range(nor)]
Ivec = [1 if rv[o] == rw[o] else 0 for o in range(nor)]
def Estar(i, c, side='L'):
    dd = d0v if side == 'L' else d0w
    return [c[o] if dd[o] == i else 0 for o in range(nor)]

# ---------------- T(x): basis by words
p = 2147483629
basis_mod, piv_mod, words = [], [], []
def add(vec):
    v = np.array([int(x) % p for x in vec], dtype=np.int64)
    for b, pc in zip(basis_mod, piv_mod):
        if v[pc]:
            v = (v - v[pc] * b) % p
    nz = np.nonzero(v)[0]
    if len(nz) == 0:
        return False
    pc = int(nz[0]); v = (v * pow(int(v[pc]), p - 2, p)) % p
    for k in range(len(basis_mod)):
        if basis_mod[k][pc]:
            basis_mod[k] = (basis_mod[k] - basis_mod[k][pc] * v) % p
    basis_mod.append(v); piv_mod.append(pc); words.append([int(x) for x in vec])
    return True
def gens_left(c):
    return [mul(Avec, c)] + [Estar(i, c) for i in range(5)]
add(Ivec); queue = [Ivec]
while queue:
    c = queue.pop()
    for nc_ in gens_left(c):
        if add(nc_):
            queue.append(nc_)
dimT = len(words)
Wm = flint.fmpq_mat(flint.fmpz_mat(words))
Rr, rk = Wm.rref()
assert rk == dimT
Bq = [[Rr[i, j] for j in range(nor)] for i in range(dimT)]
piv = [next(j for j in range(nor) if Bq[i][j] != 0) for i in range(dimT)]
def to_int(row):
    den = 1
    for x in row:
        den = den * int(x.q) // gcd(den, int(x.q))
    return [int(x.p) * (den // int(x.q)) for x in row], den
Bint, Bden = zip(*[to_int(r) for r in Bq])      # b_j = Bint[j] / Bden[j], with b_j[piv_k] = delta_jk
def coords(x):
    return [Fraction(int(x[piv[k]])) if not isinstance(x[piv[k]], Fraction) else x[piv[k]] for k in range(dimT)]
BqM = flint.fmpq_mat([[Rr[i, j] for j in range(nor)] for i in range(dimT)])
def fq(x):
    x = Fraction(int(x)) if not isinstance(x, Fraction) else x; return flint.fmpq(int(x.numerator), int(x.denominator))
def in_span(x):
    """exact membership test; returns coordinates in the RREF-row basis or None."""
    cf = flint.fmpq_mat([[fq(x[piv[k]]) for k in range(dimT)]])
    if cf * BqM != flint.fmpq_mat([[fq(y) for y in x]]):
        return None
    return cf
# closure of T(x) under left multiplication by A and E*_i, and under transpose
closed = True
for j in range(dimT):
    for y in gens_left(list(Bint[j])):
        if in_span(y) is None:
            closed = False
say(f'T(x): dimension {dimT} inside the {nor}-dimensional centraliser algebra; closed under left mult. by '
    f'A, E*_i: {closed} ({time.time() - t0:.0f} s)')

def Lmat(f):
    """matrix (fmpq) of x -> f(x) on T(x) in the basis b_j, for f linear on C mapping T(x) to T(x)."""
    cols = []
    for j in range(dimT):
        y = f(list(Bint[j]))
        cf = coords(y)
        cols.append([c_ / Bden[j] for c_ in cf])
    return flint.fmpq_mat([[flint.fmpq(cols[j][i].numerator, cols[j][i].denominator) for j in range(dimT)] for i in range(dimT)])

LA = Lmat(lambda c: mul(Avec, c))
RA = Lmat(lambda c: mul(c, Avec))
LE = [Lmat(lambda c, i=i: Estar(i, c, 'L')) for i in range(5)]
RE = [Lmat(lambda c, i=i: Estar(i, c, 'R')) for i in range(5)]

def nullspace_q(M):
    """basis (list of fmpq column vectors as lists) of the right nullspace of an fmpq_mat."""
    Rr_, rk_ = M.rref()
    n = M.ncols()
    pivs, r = [], 0
    for j in range(n):
        if r < rk_ and Rr_[r, j] != 0:
            pivs.append(j); r += 1
    free = [j for j in range(n) if j not in pivs]
    out = []
    for fj in free:
        vec = [flint.fmpq(0)] * n; vec[fj] = flint.fmpq(1)
        for r_, pj in enumerate(pivs):
            vec[pj] = -Rr_[r_, fj]
        out.append(vec)
    return out

def stack(ms):
    rows = []
    for M in ms:
        rows += [[M[i, j] for j in range(M.ncols())] for i in range(M.nrows())]
    return flint.fmpq_mat(rows)

Zb = nullspace_q(stack([LA - RA] + [LE[i] - RE[i] for i in range(5)]))
say(f'centre of T(x): dimension {len(Zb)}')

PRIME = 67108859            # < 2^26, so sums of 1755 products of residues stay inside int64
def modp(x):
    x = Fraction(x); return (x.numerator % PRIME) * pow(x.denominator % PRIME, PRIME - 2, PRIME) % PRIME

def mul_modp(a, b):
    ao = np.array(a, dtype=np.int64); bo = np.array(b, dtype=np.int64)
    return (ao[L] * bo[R] % PRIME).sum(axis=1) % PRIME

def rank_modp(M):
    M = np.array(M, dtype=np.int64) % PRIME
    r = 0; rows, cols = M.shape
    for c in range(cols):
        piv_ = next((i for i in range(r, rows) if M[i, c]), None)
        if piv_ is None:
            continue
        M[[r, piv_]] = M[[piv_, r]]
        M[r] = M[r] * pow(int(M[r, c]), PRIME - 2, PRIME) % PRIME
        nz = np.nonzero(M[:, c])[0]
        for i in nz:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % PRIME
        r += 1
        if r == rows:
            break
    return r

def primitive_int(vec):
    den = 1
    for x in vec:
        den = den * Fraction(x).denominator // gcd(den, Fraction(x).denominator)
    return [int(Fraction(x) * den) for x in vec]

def centre_structure(Zvecs, name):
    """Zvecs: integer orbit vectors spanning the centre. Returns the Q-irreducible factors and, for each, the
    central idempotent e_f as an exact orbit vector (Fractions)."""
    c = len(Zvecs)
    Zm = flint.fmpq_mat(flint.fmpz_mat(Zvecs))
    Rz, rz = Zm.rref()
    assert rz == c
    zrows = [[Fraction(int(Rz[i, j].p), int(Rz[i, j].q)) for j in range(nor)] for i in range(c)]
    zpiv = [next(j for j in range(nor) if zrows[i][j] != 0) for i in range(c)]
    zint = [primitive_int(r) for r in zrows]
    zscale = [Fraction(zint[i][zpiv[i]]) for i in range(c)]          # zint[i] = zscale[i] * zrows[i]
    one = [Fraction(Ivec[zpiv[k]]) for k in range(c)]
    for attempt in range(20):
        cf = [random.randint(-3, 3) for _ in range(c)]
        z = [sum(cf[k] * zint[k][o] for k in range(c)) for o in range(nor)]
        cols = []
        for k in range(c):
            pr = mul(z, zint[k])
            cols.append([Fraction(pr[zpiv[i]]) / zscale[k] for i in range(c)])
        Mz = flint.fmpq_mat([[flint.fmpq(cols[k][i].numerator, cols[k][i].denominator) for k in range(c)] for i in range(c)])
        cp = Mz.charpoly()
        fac = cp.factor()[1]
        if sum(f.degree() for f, e in fac) == c and all(e == 1 for f, e in fac):
            break
    else:
        raise RuntimeError('no separating central element')
    say(f'{name}: centre of dimension {c}; a central element with squarefree char. polynomial on the centre found '
        f'(attempt {attempt + 1}); Q-irreducible factors of degrees {sorted(f.degree() for f, e in fac)}')
    out = []
    onev = flint.fmpq_mat([[flint.fmpq(x.numerator, x.denominator)] for x in one])
    for f, _ in fac:
        g = flint.fmpq_poly(cp) // f
        assert g * f == cp
        # h = g * (g^-1 mod f), so h = 1 mod f and h = 0 mod the other factors
        G, S_, T_ = g.xgcd(f)
        assert G.degree() == 0
        h = (g * S_ * (1 / G[0])) % cp
        H = poly_at(h, Mz)
        ec = H * onev
        e = [sum(Fraction(int(ec[k, 0].p), int(ec[k, 0].q)) * zrows[k][o] for k in range(c)) for o in range(nor)]
        out.append({'f': f, 'deg': f.degree(), 'e': e})
    return out

def poly_at(fp, M):
    cf = [fp[i] for i in range(fp.degree() + 1)]
    n = M.nrows()
    I_ = flint.fmpq_mat([[1 if i == j else 0 for j in range(n)] for i in range(n)])
    X = I_ * cf[-1]
    for cc in reversed(cf[:-1]):
        X = X * M + I_ * cc
    return X

def trV(e, rows_dist=None):
    return sum(e[o] * int(sizes[o]) for o in diag if rows_dist is None or d0v[o] == rows_dist)

def trace_prod(a, b):
    return sum(int(sizes[o]) * sum(a[l] * b[r_] for l, r_ in zip(L[o].tolist(), R[o].tolist())) for o in diag)

# ---------------- T(x)
ZT = [primitive_int([sum(Fraction(int(Zb[k][j].p), int(Zb[k][j].q)) * Fraction(Bint[j][o], Bden[j]) for j in range(dimT))
                     for o in range(nor)]) for k in range(len(Zb))]
blocksT = centre_structure(ZT, 'T(x)')
th = [10, 5, 1, -3, -5]
Bmod = [[x % PRIME for x in Bint[j]] for j in range(dimT)]
inv_den = [pow(Bden[j] % PRIME, PRIME - 2, PRIME) for j in range(dimT)]
rowsT = []
for b in blocksT:
    e = b['e']
    em = [modp(x) for x in e]
    # rank of left multiplication by e on T(x), mod p, in the RREF-row basis (coordinates = entries at pivots)
    prods = [mul_modp(em, Bmod[j]) for j in range(dimT)]
    Lm = np.array([[int(prods[j][piv[k]]) * inv_den[j] % PRIME for j in range(dimT)] for k in range(dimT)])
    b['rank'] = rank_modp(Lm)
    # A^k e for the dual idempotents
    Ak = [e]
    for _ in range(4):
        Ak.append(mul(Avec, Ak[-1]))
    trE = []
    for j in range(5):
        pcoef = [Fraction(1)]
        for l in range(5):
            if l != j:
                # multiply polynomial by (X - th_l)/(th_j - th_l)
                q = [Fraction(0)] * (len(pcoef) + 1)
                for i_, cc in enumerate(pcoef):
                    q[i_ + 1] += cc / (th[j] - th[l]); q[i_] -= cc * th[l] / (th[j] - th[l])
                pcoef = q
        Eje = [sum(pcoef[k] * Ak[k][o] for k in range(5)) for o in range(nor)]
        trE.append(trV(Eje))
    b['trV'] = trV(e); b['trEs'] = [trV(e, i) for i in range(5)]; b['trE'] = trE
tot = sum(b['rank'] for b in blocksT)
say(f'T(x): mod-p ranks of the block ideals sum to {tot} (dimension {dimT}); each mod-p rank is at most the rank over '
    f'Q, so equality of the sums makes every block dimension exact')
assert tot == dimT
for b in blocksT:
    d2 = Fraction(b['rank'], b['deg']); d = isqrt(int(d2)); assert d2.denominator == 1 and d * d == d2
    m_ = b['trV'] / (b['deg'] * d); assert m_.denominator == 1 and m_ > 0
    r = [x / (b['deg'] * m_) for x in b['trEs']]; rd = [x / (b['deg'] * m_) for x in b['trE']]
    assert all(x.denominator == 1 for x in r + rd) and sum(r) == d == sum(rd)
    r = [int(x) for x in r]; rd = [int(x) for x in rd]
    supp = [i for i in range(5) if r[i]]; dsupp = [j for j in range(5) if rd[j]]
    b.update({'d': d, 'm': int(m_), 'r': r, 'rd': rd})
    rowsT.append((b['deg'], d, int(m_), supp[0], dsupp[0], supp[-1] - supp[0], dsupp[-1] - dsupp[0],
                  max(r) <= 1, max(rd) <= 1, r, rd, str(b['f'])))
rowsT.sort(key=lambda x: (x[3], x[4], x[1]))
say('T(x)-irreducible modules W (deg = number of Galois-conjugate components; d = dim W; m = multiplicity in C^1755):')
say('  deg  d   m  endpt dual_endpt diam dual_diam thin dual_thin  dim E*_iW          dim E_jW')
for row in rowsT:
    say(f'  {row[0]:3d} {row[1]:2d} {row[2]:3d}  {row[3]:5d} {row[4]:10d} {row[5]:4d} {row[6]:9d} {str(row[7]):5s} {str(row[8]):9s} {row[9]} {row[10]}')
say(f'checks: sum deg*d^2 = {sum(r_[0] * r_[1] ** 2 for r_ in rowsT)}, sum deg*d*m = {sum(r_[0] * r_[1] * r_[2] for r_ in rowsT)}; '
    f'all thin: {all(r_[7] for r_ in rowsT)}; all dual thin: {all(r_[8] for r_ in rowsT)}')

# ---------------- the centraliser algebra C itself
G = 300
Gm = np.random.default_rng(20261122).integers(-3, 4, size=(G, nor, nor), dtype=np.int64)
Rc = np.zeros((G, nor), dtype=np.int64)
for o in range(nor):
    lo, ro = L[o], R[o]
    np.add.at(Rc, (slice(None), lo), Gm[:, o, ro])
    np.add.at(Rc, (slice(None), ro), -Gm[:, o, lo])
del Gm
nsp, nul = flint.fmpz_mat(Rc.tolist()).nullspace()
Zc = [[int(nsp[i, k]) for i in range(nor)] for k in range(nul)]
okc = True
for c in Zc:
    for o in range(nor):
        X = defaultdict(int)
        for l, r_ in zip(L[o].tolist(), R[o].tolist()):
            X[r_] += c[l]; X[l] -= c[r_]
        if any(X.values()):
            okc = False; break
    if not okc: break
say(f'centre of C: dimension {nul} (each basis element commutes with all {nor} orbit matrices, exact: {okc})')
assert okc
blocksC = centre_structure(Zc, 'C')
for b in blocksC:
    em = [modp(x) for x in b['e']]
    Le = np.zeros((nor, nor), dtype=np.int64)
    for o in range(nor):
        np.add.at(Le[o], R[o], np.array(em, dtype=np.int64)[L[o]])
    b['rank'] = rank_modp(Le % PRIME)
tot = sum(b['rank'] for b in blocksC)
say(f'C: mod-p block ranks sum to {tot} (dimension {nor})')
assert tot == nor
for b in blocksC:
    d2 = Fraction(b['rank'], b['deg']); d = isqrt(int(d2)); assert d2.denominator == 1 and d * d == d2
    n_ = trV(b['e']) / (b['deg'] * d); assert n_.denominator == 1
    b['d'] = d; b['n'] = int(n_)
say('C = End_H(C^1755): simple components (deg, m_rho = dim of simple C-module, n_rho = dim of the H-irreducible):')
say('  ' + ', '.join(f'({b["deg"]},{b["d"]},{b["n"]})' for b in sorted(blocksC, key=lambda b: (b['d'], b['n']))))
say(f'checks: sum deg*m^2 = {sum(b["deg"] * b["d"] ** 2 for b in blocksC)}, sum deg*m*n = {sum(b["deg"] * b["d"] * b["n"] for b in blocksC)}')
say('restriction of simple C-modules to T(x): K[rho, W] = trace_V(e_rho e_W) / (n_rho d_W), summed over Galois '
    'conjugates; columns are the T(x)-modules in the order of the table above:')
orderT = sorted(blocksT, key=lambda b: ([i for i in range(5) if b['r'][i]][0], [j for j in range(5) if b['rd'][j]][0], b['d']))
for bc in sorted(blocksC, key=lambda b: (b['d'], b['n'])):
    row = [trace_prod(bc['e'], bt['e']) / (bc['n'] * bt['d']) for bt in orderT]
    say(f'  rho (deg {bc["deg"]}, m={bc["d"]}, n={bc["n"]}): ' + ' '.join(str(x) for x in row))
with open(os.path.join(OUT, 'terwilliger_modules.csv'), 'w', newline='\n') as fh:
    fh.write('galois_deg,dim,multiplicity,endpoint,dual_endpoint,diameter,dual_diameter,thin,dual_thin,dims_EstarW,dims_EW,min_poly_factor\n')
    for row in rowsT:
        fh.write(','.join(str(x) for x in row[:9]) + f',"{row[9]}","{row[10]}","{row[11]}"\n')
say(f'total {time.time() - t0:.0f} s')
open(os.path.join(OUT, 'terwilliger_wedderburn_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
