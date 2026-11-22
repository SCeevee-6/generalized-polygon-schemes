"""Flag Hecke algebra of a generalized octagon of order (s, s^2) and its generic degrees, as functions of s.

Chambers (flags) of a GO(s, t): sigma changes the line through the point (t other choices), tau changes the point on
the line (s other choices). For the Ree-Tits octagon the group 2F4(q), q = s, has a BN-pair with Weyl group
D16 = I2(8), so its orbitals on flags are the 16 Weyl distances and End_G(C[flags]) is the Iwahori-Hecke algebra
H(I2(8); q_sigma = t, q_tau = s), for every q = 2^(2e+1). (T = 2F4(2)' has index 2 in 2F4(2) and is not a
BN-pair group; it has rank 18 on flags, see flags/flag_note.md. For q > 2, 2F4(q) is simple.)

Irreducible representations of H(I2(8); a, b), a = q_sigma, b = q_tau:
  * four linear ones: T_sigma -> a or -1, T_tau -> b or -1;
  * three of dimension 2: T_sigma -> [[a, 0], [x, -1]], T_tau -> [[-1, 1], [0, b]] with
    x = a + b + sqrt(2ab), a + b, a + b - sqrt(2ab); the braid relation of length 8 is checked symbolically.
Generic degree of the G-character attached to chi: D_chi = P_W chi(1) / sum_w chi(T_w) chi(T_{w^-1}) / ind(w),
with P_W = sum_w ind(w) = number of flags. That G-character occurs chi(1) times in C[flags].
With a = s^2, b = s: sqrt(2ab) = sqrt(2 s^3) =: r, rational exactly when 2s^3 is a square (s = 2^(2e+1)).
Check: at s = 2 the degrees must be 1, 78, 1248, 4096 (linear) and 351, 650, 675 (dimension 2), the 2F4(2) flag
decomposition of flags/out and gap_runs/2026-11-22.
Usage: python family/hecke_family.py   (writes family/out/hecke_family_output.txt)
"""
import os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out'); os.makedirs(OUT, exist_ok=True)
log = []
def say(*x):
    s_ = ' '.join(str(y) for y in x); log.append(s_); print(s_, flush=True)

s = sp.symbols('s', positive=True)
r = sp.symbols('r', positive=True)          # stands for sqrt(2 s^3)
a, b = s**2, s                               # q_sigma = t, q_tau = s

def red(p):
    p = sp.Poly(sp.expand(p), r)
    out = 0
    for (k,), c in p.terms():
        out += c * (2 * s**3)**(k // 2) * r**(k % 2)
    return sp.expand(out)

def simp(e):
    """normal form d0 + d1 r with d0, d1 rational functions of s, using r^2 = 2 s^3."""
    num, den = sp.fraction(sp.together(sp.expand(e)))
    num, den = red(num), red(den)
    d0, d1 = den.coeff(r, 0), den.coeff(r, 1)
    if d1 != 0:
        num = red(num * (d0 - d1 * r)); den = sp.expand(d0**2 - d1**2 * 2 * s**3)
    return sp.factor(num.coeff(r, 0) / den) + sp.factor(num.coeff(r, 1) / den) * r

W = ['']
for L in range(1, 8):
    W += [('ST' * 4)[:L], ('TS' * 4)[:L]]
W.append('STSTSTST')
assert len(W) == 16 and len(set(W)) == 16
ind = {w: sp.Mul(*[a if c == 'S' else b for c in w]) for w in W}
PW = sp.factor(sum(ind.values()))
nflags = sp.factor((s + 1) * (1 + s * s**2) * (1 + s**2 * s**4) * (s**2 + 1))
say(f'P_W = {PW}; number of flags (s+1)(1+st)(1+s^2t^2)(t+1) = {nflags}; equal: {sp.expand(PW - nflags) == 0}')
assert sp.expand(PW - nflags) == 0

def rep_mats(kind):
    if kind[0] == 'lin':
        return sp.Matrix([[kind[1]]]), sp.Matrix([[kind[2]]])
    x = {1: a + b + r, 2: a + b, 3: a + b - r}[kind[1]]
    return sp.Matrix([[a, 0], [x, -1]]), sp.Matrix([[-1, 1], [0, b]])

reps = [('lin', a, b), ('lin', a, -1), ('lin', -1, b), ('lin', -1, -1), ('two', 1), ('two', 2), ('two', 3)]
rows = []
for kind in reps:
    Ms, Mt = rep_mats(kind)
    n = Ms.shape[0]
    I = sp.eye(n)
    ok_quad = all(simp(e) == 0 for e in (Ms - a * I) * (Ms + I)) and all(simp(e) == 0 for e in (Mt - b * I) * (Mt + I))
    ok_braid = all(simp(e) == 0 for e in ((Ms * Mt) ** 4 - (Mt * Ms) ** 4))
    def T(w):
        M = sp.eye(n)
        for c in w:
            M = M * (Ms if c == 'S' else Mt)
        return M
    tot = simp(sum(T(w).trace() * T(w[::-1]).trace() / ind[w] for w in W))
    D = simp(PW * n / tot)
    rows.append((kind, n, D))
    say(f'{kind}: dim {n}, quadratic relations {ok_quad}, braid relation {ok_braid}, generic degree D = {D}')
    assert ok_quad and ok_braid

def at(D, sv):
    return sp.nsimplify(D.subs(r, sp.sqrt(2 * sv**3)).subs(s, sv))
for sv in (2, 8, 32):
    degs = [(n, at(D, sv)) for kind, n, D in rows]
    total = sum(n * d for n, d in degs)
    say(f's = {sv}: (multiplicity, degree) = {degs}; sum mult*degree = {total}; flags = {nflags.subs(s, sv)}; '
        f'rank = sum mult^2 = {sum(n * n for n, d in degs)}')
    assert total == nflags.subs(s, sv) and all(d.is_Integer and d > 0 for n, d in degs)
d2 = (sorted(int(at(D, 2)) for kind, n, D in rows if n == 1), sorted(int(at(D, 2)) for kind, n, D in rows if n == 2))
ok = d2 == ([1, 78, 1248, 4096], [351, 650, 675])
say(f's = 2 check against the 2F4(2) flag decomposition: linear {d2[0]}, dimension 2 {d2[1]}: {ok}')
assert ok
open(os.path.join(OUT, 'hecke_family_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
