"""Tightness test (Jurisic-Koolen-Terwilliger, J. Algebraic Combin. 12 (2000) 163-197, doi:10.1023/A:1026544111089)
for the point graph of GO(2,4), array {10,8,8,8;1,1,1,5}.

Fundamental bound for a distance-regular graph of diameter d >= 3 with eigenvalues k = th_0 > th_1 > ... > th_d:
    (th_1 + k/(a_1+1)) (th_d + k/(a_1+1)) >= -k a_1 b_1 / (a_1+1)^2,
with equality iff the graph is tight. Tight graphs also have a_d = 0.
Run: python tightness_check.py   (stdlib only, exact Fractions)
"""
from fractions import Fraction as F

b = [10, 8, 8, 8]
c = [None, 1, 1, 1, 5]
k = b[0]
a = [0] + [k - b[i] - c[i] for i in range(1, 4)] + [k - c[4]]
# eigenvalues of the tridiagonal intersection matrix, known in closed form (s,t) = (2,4), m = 4
theta = sorted([10, 5, 1, -3, -5], reverse=True)
L = [[0] * 5 for _ in range(5)]
for i in range(5):
    L[i][i] = a[i]
    if i < 4: L[i][i + 1] = b[i]
    if i > 0: L[i][i - 1] = c[i]
# check each theta is a root of det(L - x I) via exact Gaussian elimination
def det(M):
    M = [[F(x) for x in r] for r in M]; n = len(M); d = F(1)
    for j in range(n):
        p = next((i for i in range(j, n) if M[i][j] != 0), None)
        if p is None: return F(0)
        if p != j: M[j], M[p] = M[p], M[j]; d = -d
        d *= M[j][j]
        for i in range(j + 1, n):
            f = M[i][j] / M[j][j]
            M[i] = [x - f * y for x, y in zip(M[i], M[j])]
    return d
assert all(det([[L[i][j] - (x if i == j else 0) for j in range(5)] for i in range(5)]) == 0 for x in theta)
s = F(k, a[1] + 1)
lhs = (theta[1] + s) * (theta[-1] + s)
rhs = F(-k * a[1] * b[1], (a[1] + 1) ** 2)
print("a_i =", a, " theta =", theta)
print("LHS = (th_1 + k/(a_1+1))(th_d + k/(a_1+1)) =", lhs)
print("RHS = -k a_1 b_1/(a_1+1)^2 =", rhs)
print("bound holds:", lhs >= rhs, " equality (tight):", lhs == rhs, " a_d =", a[4])
