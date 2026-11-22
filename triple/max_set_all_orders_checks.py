"""Symbolic checks (sympy) for max_nonopposite_all_orders.md: generalized octagons of any order (s, t), s, t >= 2.
Every identity is checked as an exact rational-function identity in s and t (with r = sqrt(2st)); inequalities
are reduced to an explicit polynomial and its sign argument is checked. The multiplicity argument for
|s - t| = 1 is checked symbolically and also numerically for s <= 10^4.
Usage: python triple/max_set_all_orders_checks.py   (writes triple/out_orbits/max_set_all_orders_checks_output.txt)
"""
import os
from fractions import Fraction
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

s, t, a, m, h, dy = sp.symbols('s t a m h d_y', positive=True)
r = sp.sqrt(2 * s * t)
b = [s * (t + 1), s * t, s * t, s * t, 0]
c = [0, 1, 1, 1, t + 1]
k = b[0]
ai = [sp.expand(k - b[i] - c[i]) for i in range(5)]
kv = [sp.Integer(1)]
for i in range(4):
    kv.append(sp.factor(kv[-1] * b[i] / c[i + 1]))
v = sp.factor(sum(kv))
K = (s + 1) * (s * t + 1)
say(f'k_i = {kv}; v = {v}; a_i = {ai}')
check('v = (s+1)(st+1)(s^2t^2+1)', v - K * (s**2 * t**2 + 1))

# eigenvalues, standard sequences, P_4
thetas = {'k': k, 'th1': s - 1 + r, 'th2': s - 1, 'th3': s - 1 - r, 'th4': -t - 1}
P4, U = {}, {}
for nm, th in thetas.items():
    u = [sp.Integer(1), th / k]
    for i in range(1, 4):
        u.append(sp.simplify(((th - ai[i]) * u[i] - c[i] * u[i - 1]) / b[i]))
    check(f'{nm} is an eigenvalue of A_1', sp.simplify((th - ai[4]) * u[4] - c[4] * u[3]))
    U[nm] = u
    P4[nm] = sp.simplify(kv[4] * u[4])
say('P_4 =', {nm: sp.factor(x) for nm, x in P4.items()})
check('P_4(k) = s^4 t^3', P4['k'] - s**4 * t**3)
check('P_4(th1) = -s^2 t', P4['th1'] + s**2 * t)
check('P_4(th2) = s^2 t', P4['th2'] - s**2 * t)
check('P_4(th3) = -s^2 t', P4['th3'] + s**2 * t)
check('P_4(th4) = t^3', P4['th4'] - t**3)
say('least eigenvalue of A_4 is -s^2 t: t^3 > 0 and s^2 t > 0')
check('Hoffman bound v/(1 + s^4t^3/(s^2t)) = K', sp.factor(v / (1 + s**2 * t**2)) - K)
check('|N[L]| = (s+1)(1+st) = K', (s + 1) + (s + 1) * s * t - K)

# operator (A1-th1)(A1-th3) = A2 - (s-1)A1 + (s^2-s+1-st)I, value on 1 is s^2t^2+1, so it maps 1_C to 1
coefI = sp.expand(k + (s - 1)**2 - 2 * s * t)
check('constant term k+(s-1)^2-2st = s^2-s+1-st', coefI - (s**2 - s + 1 - s * t))
check('(k-th1)(k-th3) = s^2t^2+1', sp.expand((k - thetas['th1']) * (k - thetas['th3'])) - (s**2 * t**2 + 1))
n2 = (s - 1) * a + s * (t - s + 1)
n3 = s**2 * t + s**2 - s * a
m2 = (s - 1) * m + 1
m3 = s * t + s - s * m
check('member: n2 - (s-1)a + coefI = 1', n2 - (s - 1) * a + coefI - 1)
check('member: 1 + a + n2 + n3 + 0 = K', 1 + a + n2 + n3 - K)
check('nonmember: m2 - (s-1)m + 0 = 1', m2 - (s - 1) * m - 1)
check('nonmember: m + m2 + m3 + s^2t = K', m + m2 + m3 + s**2 * t - K)

# local structure
cM = n3 - s**2 * t
check('c_M = s(s-a)', cM - s * (s - a))
check('all lines incomplete: (t+1)s(s-a) - n3 = -sta, so a = 0', (t + 1) * s * (s - a) - n3 + s * t * a)
check('k-point: n2(k) = s^2 t', n2.subs(a, k) - s**2 * t)
check('k-point: n3(k) = 0', n3.subs(a, k))
check('s-point: n2(s) = st', n2.subs(a, s) - s * t)
# s-point, y on M: m = 1 + d_y, m2 = s + st - d_y must equal (s-1)m + 1
sol = sp.solve(sp.Eq(s + s * t - dy, m2.subs(m, 1 + dy)), dy)
say('s-point case: d_y =', sol)
assert sol == [t]
# isolated z: m3(w) = n2(z) - d_N + c_N - c_w with m = 1 + d_w
dw, dN, cw = sp.symbols('d_w d_N c_w')
sol = sp.solve(sp.Eq(m3.subs(m, 1 + dw), n2.subs(a, 0) - dN + s**2 - cw), cw)
say('isolated z: c_w =', sol)
check('c_w = s + s d_w - d_N', sol[0] - (s + s * dw - dN))
cy = (s + s * dw - dN).subs({dw: t, dN: s * (t - s + 1)})
check('c_y = s^2 (d_y = t, d_N = s(t-s+1))', cy - s**2)
gap = sp.expand(s + (s - 1) * t - s**2)
check('s + (s-1)t - s^2 = (s-1)(t-s)', gap - (s - 1) * (t - s))
say('  so c_y >= s + (s-1)t contradicts c_y = s^2 iff t > s; n2(0) = s(t-s+1) < 0 iff t <= s-2;'
    ' t = s-1 gives n2(z) = 0 but x is at distance 2 from z')
# all points isolated, t >= s+2
lhs = (s - 1) * h + 1
rhs = h * (s * (t - s + 1) + 1 - h)
check('rhs - lhs = h(st - s^2 + 2 - h) - 1', rhs - lhs - (h * (s * t - s**2 + 2 - h) - 1))
check('st - s^2 + 2 - (t+1) = (s-1)(t-s-1)', s * t - s**2 + 2 - (t + 1) - (s - 1) * (t - s - 1))
say('  h >= t-s+2 >= 4, s >= 2, t-s-1 >= 1 give rhs - lhs >= 4*1*1 - 1 > 0')

# multiplicity of th4 = -(t+1)
check('standard sequence of -(t+1) is (-1/s)^i', sum(sp.simplify(U['th4'][i] - (-1 / s)**i)**2 for i in range(5)))
mult = sp.factor(v / sum(kv[i] * U['th4'][i]**2 for i in range(5)))
say('mult(-(t+1)) =', mult)
check('mult = s^4(st+1)(s^2t^2+1)/((s+t)(s^2+t^2))',
      mult - s**4 * (s * t + 1) * (s**2 * t**2 + 1) / ((s + t) * (s**2 + t**2)))
for tt in (s + 1, s - 1):
    d = sp.expand(s**2 + tt**2)
    check(f't = {tt}: s^2+t^2 = 2st+1', d - (2 * s * tt + 1))
    check(f't = {tt}: 4((st)^2+1) - 5 is divisible by 2st+1',
          sp.rem(sp.expand(4 * ((s * tt)**2 + 1) - 5), sp.expand(2 * s * tt + 1), s))
bad = 0
for S in range(2, 10001):
    for T in (S - 1, S + 1):
        if T < 2:
            continue
        num = S**4 * (S * T + 1) * (S * S * T * T + 1)
        den = (S + T) * (S * S + T * T)
        if num % den == 0:
            bad += 1
say(f'numeric: integral multiplicity for some |s-t|=1, 2 <= s <= 10000: {bad} cases')
assert bad == 0
say('all checks passed')
with open(os.path.join(OUT, 'max_set_all_orders_checks_output.txt'), 'w', newline='\n') as f:
    f.write('\n'.join(log) + '\n')
