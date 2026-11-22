"""LP + Krein + general constraints + free-tree counts, per base shape, at (s, t) = (2, 4) and (8, 64).

Adds to family/triple_family.py the exact counts of family/tree_counts.py (valid in every GO(s, t) for the given base
shape). For each type with a tree base, and each base shape, reports the free dimension left after LP + G + tree
counts + implicit zeros, and at s = 2 checks the result against the Ree-Tits vectors of that shape.
If the free dimension is 0, the triple intersection numbers of that shape are the same in every GO(s, t).
Usage: python family/triple_tree_family.py   (writes family/out/triple_tree_family_output.txt, .csv)
"""
import os, sys, csv, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'triple'))
src = open(os.path.join(HERE, 'triple_family.py'), encoding='utf-8').read()
ns = {'__file__': os.path.join(HERE, 'triple_family.py'), '__name__': 'triple_family_lib'}
exec(compile(src[:src.index('HAND = {')], 'triple_family.py', 'exec'), ns)
F, scheme, system, rank, implicit_zeros, determined = (ns[k] for k in ('F', 'scheme', 'system', 'rank', 'implicit_zeros', 'determined'))
import triple_zero_proofs as tz
import tree_counts as tc

OUT = os.path.join(HERE, 'out')
log = []
def say(*x):
    s_ = ' '.join(str(y) for y in x); log.append(s_); print(s_, flush=True)

HAND = {(2, 3, 4): ['313'], (2, 4, 4): ['133', '313']}
vecs = {}
with open(os.path.join(HERE, '..', 'triple', 'out_orbits', 'triple_orbit_vectors.csv')) as fh:
    for row in csv.DictReader(fh):
        ty = tuple(int(x) for x in re.findall(r'\d', row['type']))
        vv = {tuple(int(c) for c in m_.group(1)): int(m_.group(2)) for m_ in re.finditer(r'\[(\d\d\d)\]=(\d+)', row['vector'])}
        vecs.setdefault(ty, []).append(vv)

rows_out = []
for (s, t) in [(2, 4), (8, 64)]:
    p, Q, q444, v = scheme(s, t)
    tz.s_, tz.t_ = s, t
    say(f'(s, t) = ({s}, {t})')
    types = [(x, y, z) for x in range(1, 5) for y in range(x, 5) for z in range(y, 5) if p[z][x][y] > 0]
    for ty in types:
        C = tc.constraints(*ty, s, t)
        if not C:
            say(f'  {ty}: no tree base (two opposite pairs), skipped')
            continue
        keys, rows, rhs = system(p, Q, *ty)
        n = len(keys)
        G = {}
        for kk in keys:
            if 1 <= min(kk) and max(kk) <= 3 and not tz.feasible(*ty, *kk):
                G[kk] = F(0)
        for h in HAND.get(ty, []):
            G[tuple(map(int, h))] = F(t)
        P = lambda h, i, j: p[h][i][j]
        ub = [max(1, min(P(ty[0], i, j), P(ty[1], i, k), P(ty[2], j, k))) for (i, j, k) in keys]
        for legs, cons in C.items():
            fix = dict(G)
            clash = [kk for kk, val in cons.items() if kk in fix and fix[kk] != val]
            for kk, val in cons.items():
                if kk in keys or val != 0:
                    fix[kk] = F(val)
            missing = [kk for kk, val in fix.items() if kk not in keys and val != 0]
            if clash or missing:
                say(f'  {ty} legs {legs}: tree counts clash with G at {clash} or name impossible unknowns {missing}: shape impossible')
                continue
            fix = {kk: val for kk, val in fix.items() if kk in keys}
            R2 = rows + [[F(int(keys[c] == kk)) for c in range(n)] for kk in fix]
            b2 = rhs + list(fix.values())
            rk = rank(R2)
            consistent = rank([r_ + [b_] for r_, b_ in zip(R2, b2)]) == rk
            if not consistent:
                say(f'  {ty} legs {legs}: LP + G + tree counts inconsistent: shape impossible at ({s},{t})')
                rows_out.append((s, t, ty, legs, len(cons), 'inconsistent', '', '', ''))
                continue
            z, margin = implicit_zeros(keys, R2, b2, ub)
            if z is None:
                say(f'  {ty} legs {legs}: LP infeasible (nonnegativity): shape impossible at ({s},{t})')
                rows_out.append((s, t, ty, legs, len(cons), 'infeasible', '', '', ''))
                continue
            R3 = R2 + [[F(int(keys[c] == kk)) for c in range(n)] for kk in z]
            free = n - rank(R3)
            det = determined(keys, R3) if free > 0 else keys
            rt = ''
            if s == 2:
                match = [vv for vv in vecs.get(ty, []) if all(vv.get(kk, 0) == val for kk, val in cons.items())]
                if free == 0 and match:
                    # the unique solution: solve and compare
                    sol = {}
                    import sympy as sp_
                    M = sp_.Matrix([[sp_.Rational(x.numerator, x.denominator) for x in r_] for r_ in R3])
                    bb = sp_.Matrix([sp_.Rational(x.numerator, x.denominator) for x in (b2 + [F(0)] * len(z))])
                    X = M.gauss_jordan_solve(bb)[0]
                    sol = {keys[c]: X[c] for c in range(n)}
                    rt = f'Ree-Tits orbits of this shape {len(match)}; unique solution equals them: ' + str(
                        all(all(vv.get(kk, 0) == sol[kk] for kk in keys) for vv in match))
                else:
                    rt = f'Ree-Tits orbits of this shape {len(match)}, distinct vectors {len({tuple(sorted(vv.items())) for vv in match})}'
            say(f'  {ty} legs {legs}: {len(cons)} tree counts, free dim after LP+G+tree+implicit zeros = {free}, '
                f'determined unknowns {len(det)}/{n}; {rt}')
            rows_out.append((s, t, ty, legs, len(cons), free, len(det), n, rt))
with open(os.path.join(OUT, 'triple_tree_family.csv'), 'w', newline='\n') as fh:
    fh.write('s,t,type,legs,tree_counts,free_dim,determined,unknowns,ree_tits_check\n')
    for r_ in rows_out:
        fh.write(','.join(f'"{x}"' for x in r_) + '\n')
open(os.path.join(OUT, 'triple_tree_family_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
