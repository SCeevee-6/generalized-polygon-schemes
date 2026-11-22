"""Which eigenspace carries the vanishing Krein parameter of a generalized octagon of order (s, t)?

Bamberg-Lansdown (Algebr. Comb. 6 (2023) 197-212, Thm 8.1) write the zero as q^2_22 in their own ordering,
where E_2 is the eigenspace for which a line is a {2}-design, i.e. (aQ)_2 = 0 for the inner distribution
a = (1, s, 0, 0, 0). This script checks, exactly:
  1. (aQ)_j = m_j (1 + s theta_j / k) vanishes only for theta_j = -(t+1);
  2. for theta = -(t+1), sum_l k_l u_l^3 (proportional to the Krein parameter q^j_jj) as a rational
     function of s, t, and its factors;
  3. for (s,t) = (2,4), all Krein parameters q^k_ij in the ordering 10, 5, 1, -3, -5 used in this repo, and
     which of them vanish.
Usage: python triple/krein_convention.py
"""
import sympy as sp

s, t = sp.symbols('s t', positive=True)

def array(s, t):
    b = [s * (t + 1), s * t, s * t, s * t]
    c = [1, 1, 1, t + 1]
    return b, c

def seq(theta, b, c):
    k = b[0]; bb = b + [0]; cc = [0] + c
    a = [k - bb[i] - cc[i] for i in range(5)]
    u = [sp.Integer(1), theta / k]
    for i in range(1, 4):
        u.append(sp.simplify(((theta - a[i]) * u[i] - cc[i] * u[i - 1]) / bb[i]))
    return u

def valencies(b, c):
    kv = [sp.Integer(1)]
    for i in range(4):
        kv.append(sp.simplify(kv[-1] * b[i] / c[i]))
    return kv

b, c = array(s, t)
kv = valencies(b, c)
theta = -(t + 1)
u = seq(theta, b, c)
k = b[0]
print('1. line inner distribution: 1 + s*theta/k = 0  <=>  theta =', sp.solve(sp.Eq(1 + s * sp.Symbol('th') / k, 0), sp.Symbol('th')))
S3 = sp.factor(sp.together(sum(kv[l] * u[l] ** 3 for l in range(5))))
print('2. theta = -(t+1): sum_l k_l u_l^3 =', S3)
n = sum(kv)
m = sp.factor(sp.simplify(n / sum(kv[l] * u[l] ** 2 for l in range(5))))
print('   multiplicity of -(t+1):', m, '; at (2,4):', m.subs({s: 2, t: 4}))

# 3. (2,4), exact
b4, c4 = array(sp.Integer(2), sp.Integer(4))
kv4 = valencies(b4, c4)
n4 = sum(kv4)
th = [10, 5, 1, -3, -5]
U = [seq(sp.Integer(x), b4, c4) for x in th]
M = [sp.nsimplify(n4 / sum(kv4[l] * U[j][l] ** 2 for l in range(5))) for j in range(5)]
print('3. (2,4): eigenvalues', th, 'multiplicities', M)
zeros = []
for i in range(1, 5):
    for j in range(i, 5):
        for h in range(1, 5):
            q = sp.Rational(1, 1) * M[i] * M[j] / n4 * sum(kv4[l] * U[i][l] * U[j][l] * U[h][l] for l in range(5))
            if q == 0:
                zeros.append((i, j, h))
print('   vanishing q^h_ij (i<=j, i,j,h != 0) in this ordering:', zeros)
print('   so Bamberg-Lansdown q^2_22 (E_2 = eigenspace of -(t+1)) is q^4_44 here (eigenvalue -5, multiplicity', M[4], ')')
