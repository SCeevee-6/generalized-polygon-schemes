"""Exact integer coefficient matrices for the extra linear relations among triple intersection numbers of the
Ree-Tits octagon GO(2,4) (the relations counted in triple_note.md: LP free dimension minus affine span of the
observed vectors).

For each distance type (U, V, W):
  * the LP system of go24_triple_lp.py (marginals = intersection numbers, plus the Krein equation for q^4_44 = 0)
    is rebuilt over Q and put in reduced row echelon form; the non-pivot unknowns are the free variables F;
  * the observed vectors (one per T-orbit, from out_orbits/triple_orbit_vectors.csv) are restricted to F;
  * the relations are the integer row vectors lam on F with lam . (v - v_1) = 0 for all observed v, i.e. the
    nullspace of the differences, in RREF and scaled to primitive integer rows; the right-hand side is lam . v_1.
Since every LP solution is an affine function of its F-coordinates, these rows are a basis of the extra
relations modulo the LP equations. Each row is checked exactly on every T-orbit vector.
Columns with an index 4 are ordered first, so the free variables are counts [ijk] with i, j, k <= 3 wherever
possible (points within distance 3 of all of u, v, w); the others are affine in these and the intersection numbers.
Usage: python triple/triple_relations.py   (writes out_orbits/triple_relations.csv, triple_relations_long.csv,
                                            triple_relations_output.txt)
"""
import os, sys, csv, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'drgfusion'))
sys.path.insert(0, os.path.join(HERE, '..', 'tits1755'))
import flint
import sympy as sp
from math import gcd
from fractions import Fraction
import go24_triple_lp as lp

OUT = os.path.join(HERE, 'out_orbits')
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)

TYPES = [(2, 3, 3), (2, 3, 4), (2, 4, 4), (3, 3, 3), (3, 3, 4), (3, 4, 4)]
EXPECT = {(2, 3, 3): 15, (2, 3, 4): 12, (2, 4, 4): 9, (3, 3, 3): 14, (3, 3, 4): 7, (3, 4, 4): 3}

def exact_system(Uu, Vv, Ww):
    import itertools
    P = lp.P
    idx = {}
    for i, j, k in itertools.product(range(5), repeat=3):
        if P(Uu, i, j) > 0 and P(Vv, i, k) > 0 and P(Ww, j, k) > 0:
            idx[(i, j, k)] = len(idx)
    rows, rhs = [], []
    for (A, R) in ((lambda i, j, k: (i, j), Uu), (lambda i, j, k: (i, k), Vv), (lambda i, j, k: (j, k), Ww)):
        for x, y in itertools.product(range(5), repeat=2):
            r = [0] * len(idx)
            for key in idx:
                if A(*key) == (x, y):
                    r[idx[key]] = 1
            if any(r) or P(R, x, y):
                rows.append([sp.Integer(v) for v in r]); rhs.append(sp.Integer(P(R, x, y)))
    rows.append([sp.nsimplify(lp.Q[i, 4] * lp.Q[j, 4] * lp.Q[k, 4]) for (i, j, k) in idx]); rhs.append(sp.Integer(0))
    return list(idx), rows, rhs

def fq(x):
    x = sp.Rational(x); return flint.fmpq(int(x.p), int(x.q))

def rref_info(rows, n):
    M = flint.fmpq_mat([[fq(v) for v in r] for r in rows])
    Rr, rk = M.rref()
    piv, r = [], 0
    for j in range(n):
        if r < rk and Rr[r, j] != 0:
            piv.append(j); r += 1
    return piv, rk

def primitive(row):
    den = 1
    for x in row:
        den = den * x.denominator // gcd(den, x.denominator)
    ints = [int(x * den) for x in row]
    g = 0
    for v in ints: g = gcd(g, abs(v))
    ints = [v // g for v in ints]
    first = next(v for v in ints if v)
    return [-v for v in ints] if first < 0 else ints

vecs = {}
with open(os.path.join(OUT, 'triple_orbit_vectors.csv')) as fh:
    for row in csv.DictReader(fh):
        ty = tuple(int(x) for x in re.findall(r'\d', row['type']))
        vv = {tuple(int(c) for c in m.group(1)): int(m.group(2)) for m in re.finditer(r'\[(\d\d\d)\]=(\d+)', row['vector'])}
        vecs.setdefault(ty, []).append((int(row['T_orbit']), vv))

out_rows, long_rows = [], []
CONST = {}
PROVEN_DROP = {}
import triple_zero_proofs as tz
HAND = {(2, 3, 4): {'313': 4}, (2, 4, 4): {'133': 4, '313': 4}}
CONSTKEYS = lambda keys, idx: [''.join(map(str, keys[j])) for j in idx]
total = 0
for ty in TYPES:
    keys, rows, rhs = exact_system(*ty)
    n = len(keys)
    # columns with an index 4 first, so that RREF pivots on them and the free variables are, as far as
    # possible, counts of points within distance 3 of all of u, v, w (unique geodesics)
    perm = sorted(range(n), key=lambda j: (max(keys[j]) <= 3, keys[j]))
    keys = [keys[j] for j in perm]
    rows = [[r[j] for j in perm] for r in rows]
    piv, rk = rref_info(rows, n)
    free = [j for j in range(n) if j not in piv]
    obs = vecs[ty]
    # every observed vector satisfies the LP equations exactly
    for _, vv in obs:
        x = [vv[k] for k in keys]
        assert all(sum(sp.Rational(r[j]) * x[j] for j in range(n)) == b for r, b in zip(rows, rhs))
    distinct = sorted({tuple(vv[k] for k in keys) for _, vv in obs})
    base = distinct[0]
    diffs = [[Fraction(d[j] - base[j]) for j in free] for d in distinct[1:]]
    f = len(free)
    if diffs:
        Dm = flint.fmpq_mat([[flint.fmpq(int(x.numerator), int(x.denominator)) for x in r] for r in diffs])
        Rr, rD = Dm.rref()
        pv, r = [], 0
        for j in range(f):
            if r < rD and Rr[r, j] != 0:
                pv.append(j); r += 1
        fr = [j for j in range(f) if j not in pv]
        null = []
        for fj in fr:
            vec = [Fraction(0)] * f; vec[fj] = Fraction(1)
            for r_, pj in enumerate(pv):
                e = Rr[r_, fj]; vec[pj] = -Fraction(int(e.p), int(e.q))
            null.append(vec)
    else:
        rD = 0
        null = [[Fraction(int(i == j)) for j in range(f)] for i in range(f)]
    # canonical basis: RREF of the relation space, primitive integer rows
    if null:
        Nm = flint.fmpq_mat([[flint.fmpq(int(x.numerator), int(x.denominator)) for x in r] for r in null])
        Rn, rn = Nm.rref()
        rels = [primitive([Fraction(int(Rn[i, j].p), int(Rn[i, j].q)) for j in range(f)]) for i in range(rn)]
    else:
        rels = []
    assert len(rels) == f - rD == EXPECT[ty], (ty, len(rels), f, rD)
    # check every relation exactly on every T-orbit vector
    for li, lam in enumerate(rels):
        vals = {sum(lam[a] * vv[keys[free[a]]] for a in range(f)) for _, vv in obs}
        assert len(vals) == 1
        c = vals.pop()
        terms = ' '.join(f'{"+" if lam[a] > 0 else "-"}{abs(lam[a]) if abs(lam[a]) != 1 else ""}[{"".join(map(str, keys[free[a]]))}]'
                         for a in range(f) if lam[a])
        out_rows.append((ty, li + 1, c, sum(1 for v in lam if v), terms))
        for a in range(f):
            if lam[a]:
                long_rows.append((ty, li + 1, ''.join(map(str, keys[free[a]])), lam[a]))
        long_rows.append((ty, li + 1, 'rhs', c))
    le3 = [j for j in free if max(keys[j]) <= 3]
    total += len(rels)
    # coordinates (indices <= 3, not forced by the LP) that are constant over all T-orbits, and how many of the
    # relations they account for: add [ijk] = value to the LP and measure the drop in free dimension
    const = [j for j in range(n) if max(keys[j]) <= 3 and len({d[j] for d in distinct}) == 1]
    rows_c = rows + [[sp.Integer(int(i == j)) for i in range(n)] for j in const]
    _, rk_c = rref_info(rows_c, n)
    lp_forced = [j for j in const if rref_info(rows + [[sp.Integer(int(i == j)) for i in range(n)]], n)[1] == rk]
    newc = [j for j in const if j not in lp_forced]
    CONST[ty] = {''.join(map(str, keys[j])): distinct[0][j] for j in newc}
    # constants proved for every GO(2,4): zeros from the tree test (triple_zero_proofs.py) and the hand counts
    # [313] = t (types (2,3,4), (2,4,4)) and [133] = t (type (2,4,4)) of triple_relations_note.md
    proven = {kk: 0 for kk in [''.join(map(str, keys[j])) for j in range(n)]
              if max(map(int, kk)) <= 3 and min(map(int, kk)) >= 1 and not tz.feasible(*ty, *map(int, kk))}
    proven.update(HAND.get(ty, {}))
    pj = [keys.index(tuple(map(int, kk))) for kk in proven]
    for kk, val in proven.items():
        assert {d[keys.index(tuple(map(int, kk)))] for d in distinct} == {val}
    _, rk_p = rref_info(rows + [[sp.Integer(int(i == j)) for i in range(n)] for j in pj], n)
    implied = [kk for kk in CONSTKEYS(keys, newc) if kk not in proven and
               rref_info(rows + [[sp.Integer(int(i == j)) for i in range(n)] for j in pj] +
                         [[sp.Integer(int(i == keys.index(tuple(map(int, kk))))) for i in range(n)]], n)[1] == rk_p]
    PROVEN_DROP[ty] = rk_p - rk
    say(f'   proved for every GO(2,4): {proven}; with these added the LP free dimension drops by {rk_p - rk}; '
        f'other observed constants implied by them: {implied}')
    say(f'   constant over all orbits and not fixed by the LP: '
        f'{", ".join("[" + "".join(map(str, keys[j])) + "]=" + str(distinct[0][j]) for j in newc)}; '
        f'they account for {rk_c - rk} of the {len(rels)} relations')
    say(f'type {ty}: {n} unknowns, LP rank {rk}, free dim {f} (free vars {["".join(map(str, keys[j])) for j in free]}); '
        f'{len(obs)} T-orbits, {len(distinct)} distinct vectors, affine span {rD}; relations {len(rels)} '
        f'(checked on all orbits); max |coeff| {max(abs(v) for r in rels for v in r)}; '
        f'free vars with all indices <= 3: {len(le3)} of {f}')
say(f'total relations: {total}; accounted for by constants proved for every GO(2,4): {sum(PROVEN_DROP.values())} {PROVEN_DROP}')
with open(os.path.join(OUT, 'triple_relations.csv'), 'w', newline='\n') as fh:
    fh.write('type,relation,rhs,support,equation\n')
    for ty, li, c, sup, terms in out_rows:
        fh.write(f'"{ty}",{li},{c},{sup},"{terms} = {c}"\n')
with open(os.path.join(OUT, 'triple_relations_long.csv'), 'w', newline='\n') as fh:
    fh.write('type,relation,unknown,coefficient\n')
    for ty, li, k, v in long_rows:
        fh.write(f'"{ty}",{li},{k},{v}\n')
open(os.path.join(OUT, 'triple_relations_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
