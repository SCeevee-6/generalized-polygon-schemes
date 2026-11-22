"""Symbolic checks (sympy) of every identity used in max_nonopposite_proof.md for a generalized octagon of
order (s, s^2), s >= 2. Each check is an exact polynomial identity in s, or an inequality reduced to an
explicit polynomial that is shown to be positive for s >= 2. The s = 2 values are printed for comparison
with max_set_checks.py.
Usage: python triple/max_set_general_checks.py     (writes triple/out_orbits/max_set_general_checks_output.txt)
"""
import os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out_orbits'); os.makedirs(OUT, exist_ok=True)
log = []
def say(*a):
    s_ = ' '.join(str(x) for x in a); log.append(s_); print(s_, flush=True)
def check(name, expr):
    ok = sp.simplify(sp.expand(expr)) == 0
    say(f'[{"ok" if ok else "FAIL"}] {name}')
    assert ok, name

s = sp.symbols('s', positive=True)
t = s**2
r = sp.sqrt(2 * s**3)
b = [s * (t + 1), s * t, s * t, s * t, 0]
c = [0, 1, 1, 1, t + 1]
k = b[0]
a = [sp.expand(k - b[i] - c[i]) for i in range(5)]
kv = [sp.Integer(1)]
for i in range(4):
    kv.append(sp.factor(kv[-1] * b[i] / c[i + 1]))
v = sp.factor(sum(kv))
say(f'k_i = {kv}; v = {v}; a_i = {a}')
check('v = (1+s)(1+s^3)(1+s^6)', v - (1 + s) * (1 + s**3) * (1 + s**6))

thetas = {'k': k, 'th1': s - 1 + r, 'th2': s - 1, 'th3': s - 1 - r, 'th4': -t - 1}
P4 = {}
P = {}
for nm, th in thetas.items():
    u = [sp.Integer(1), th / k]
    for i in range(1, 4):
        u.append(sp.simplify(((th - a[i]) * u[i] - c[i] * u[i - 1]) / b[i]))
    # theta is an eigenvalue: the recurrence closes at i = 4
    check(f'{nm} is an eigenvalue of A_1', sp.simplify((th - a[4]) * u[4] - c[4] * u[3]))
    P[nm] = [sp.simplify(kv[i] * u[i]) for i in range(5)]
    P4[nm] = P[nm][4]
say('P_4 =', {nm: sp.factor(x) for nm, x in P4.items()})
check('P_4(k) = s^10', P4['k'] - s**10)
check('P_4(th1) = -s^4', P4['th1'] + s**4)
check('P_4(th2) = s^4', P4['th2'] - s**4)
check('P_4(th3) = -s^4', P4['th3'] + s**4)
check('P_4(th4) = s^6', P4['th4'] - s**6)
K = sp.factor(v * s**4 / (s**10 + s**4))
say(f'ratio bound K = {K}')
check('K = (s+1)(s^3+1)', K - (s + 1) * (s**3 + 1))
check('line neighbourhood size (s+1)(1 + s t) = K', (s + 1) * (1 + s * t) - K)

# (A1 - th1)(A1 - th3) = A1^2 - 2(s-1)A1 + (s-1)^2 - 2 s^3, and A1^2 = k I + a1 A1 + A2
const = sp.expand((s - 1)**2 - 2 * s**3)
coefA1 = -2 * (s - 1) + a[1]
coefI = sp.expand(k + const)
check('operator = A2 - (s-1)A1 - (s-1)(s^2+1) I: A1 coefficient', coefA1 + (s - 1))
check('operator = A2 - (s-1)A1 - (s-1)(s^2+1) I: I coefficient', coefI + (s - 1) * (s**2 + 1))
check('operator on the constant part (K/v)(k2 - (s-1)k - (s-1)(s^2+1)) = 1',
      K / v * (kv[2] - (s - 1) * kv[1] - (s - 1) * (s**2 + 1)) - 1)

# profiles
aa, m = sp.symbols('a m')
n2 = (s - 1) * aa + (s - 1) * (s**2 + 1) + 1
n3 = K - 1 - aa - n2
m2 = (s - 1) * m + 1
m3 = K - s**4 - m - m2
check('n2 = (s-1)a + s^3 - s^2 + s', n2 - ((s - 1) * aa + s**3 - s**2 + s))
check('n3 = s^4 + s^2 - s a', n3 - (s**4 + s**2 - s * aa))
check('m3 = s^3 + s - s m', m3 - (s**3 + s - s * m))
# Lemma A: s^4 = n3 - c_M
cM = sp.expand(n3 - s**4)
check('Lemma A: c_M = s(s - a)', cM - s * (s - aa))
# Lemma B, f = 0: (t+1) c_M = n3  <=>  s^3 a = 0
check('Lemma B (f = 0): (t+1) c_M - n3 = -s^3 a', (t + 1) * cM - n3 + s**3 * aa)
# Case I
check('Case I: n2(k) = s^4', n2.subs(aa, k) - s**4)
check('Case I: n3(k) = 0', n3.subs(aa, k))
check('Case I: s^4 / (k - s) = s', s**4 / (k - s) - s)
# Case II: m2(y) = s + (s^3 - d) and m1 = 1 + d with m2 = (s-1) m1 + 1  =>  d = s^2
d = sp.symbols('d')
sol = sp.solve(sp.Eq(s + s**3 - d, (s - 1) * (1 + d) + 1), d)
say(f'Case II: d_y(x) = {sol}')
assert sol == [s**2]
check('Case II: m3(y) = 0 for m1 = s^2 + 1', m3.subs(m, s**2 + 1))
# Lemma D: c_w = s + s d_w - d_N, from m3(w) = (n2 - d_N) + (c_N - c_w) with a = 0, c_N = s^2
dw, dN, cw = sp.symbols('d_w d_N c_w')
lemD = sp.solve(sp.Eq(n2.subs(aa, 0) - dN + s**2 - cw, m3.subs(m, 1 + dw)), cw)[0]
check('Lemma D: c_w = s + s d_w - d_N', lemD - (s + s * dw - dN))
lemD2 = sp.solve(sp.Eq(dN - dw + cw, m2.subs(m, 1 + dw)), cw)[0]
check('Lemma D from m2 gives the same c_w', lemD2 - lemD)
cy = lemD.subs({dw: s**2, dN: n2.subs(aa, 0)})
check('Case II: c_y(z) = s^2', cy - s**2)
gap2 = sp.factor(s + (s - 1) * s**2 - cy)
say(f'Case II: lower bound minus value = {gap2}  (positive for s >= 2)')
check('Case II gap = s(s-1)^2', gap2 - s * (s - 1)**2)
# Case III: k in [s^2 - s + 2, s^2 + 1], g(k) = k(s^3 - s^2 + s + 1 - k) - ((s-1)k + 1)
kk = sp.symbols('kk')
g = kk * (s**3 - s**2 + s + 1 - kk) - ((s - 1) * kk + 1)
check('Case III: g = k(s^3 - s^2 + 2 - k) - 1', g - (kk * (s**3 - s**2 + 2 - kk) - 1))
lo, hi = s**2 - s + 2, s**2 + 1
gl, gh = sp.factor(g.subs(kk, lo)), sp.factor(g.subs(kk, hi))
say(f'Case III: g(s^2-s+2) = {sp.expand(gl)}, g(s^2+1) = {sp.expand(gh)}; g is concave in k')
for nm, e in [('g(lo)', gl), ('g(hi)', gh)]:
    e2 = sp.expand(e.subs(s, s + 2))          # shift so positivity of coefficients shows e > 0 for s >= 2
    pos = all(cf > 0 for cf in sp.Poly(e2, s).coeffs())
    say(f'[{"ok" if pos else "FAIL"}] Case III: {nm} at s -> s+2 has positive coefficients: {e2}')
    assert pos
# s = 2 values
say('s = 2:', {'K': K.subs(s, 2), 'n2': n2.subs(s, 2), 'n3': n3.subs(s, 2), 'm3': m3.subs(s, 2),
               'c_M': cM.subs(s, 2), 'P4': [P4[nm].subs(s, 2) for nm in P4]})
open(os.path.join(OUT, 'max_set_general_checks_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
