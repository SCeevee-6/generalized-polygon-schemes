"""Verify |T| (Schreier-Sims, sympy) and spectrum of A_1 / irreducibility of the 5 eigenspaces for T=2F4(2)' on 1755 points.
Generators: ATLAS v3 TF42G1-p1755B0.g1/.g2 (needs network). Output: stdout (see verify_group_order_output.txt).
Irreducibility argument: rank (= #orbitals of the point stabiliser) = 5 = <pi,pi>; pi splits into 5 nonzero invariant
eigenspaces of A_1, so each is irreducible and pi is multiplicity-free. Degrees are then the eigenspace dimensions.
"""
import re, urllib.request, numpy as np
from sympy.combinatorics import Permutation, PermutationGroup
B = "https://brauer.maths.qmul.ac.uk/Atlas/exc/TF42/gap/TF42G1-p1755B0.g"
n = 1755
def parse(t):
    body = re.sub(r"\s+", "", t.split(":=")[1].strip().rstrip(";"))
    p = list(range(n))
    for c in re.findall(r"\(([^)]+)\)", body):
        ix = [int(x) - 1 for x in c.split(",")]
        for k in range(len(ix)): p[ix[k]] = ix[(k + 1) % len(ix)]
    return p
g = [parse(urllib.request.urlopen(B + str(i), timeout=30).read().decode()) for i in (1, 2)]
G = PermutationGroup([Permutation(x) for x in g])
order = G.order(); print("|T| =", order, "| equals 2^11.3^3.5^2.13 = 17971200:", order == 2**11 * 3**3 * 5**2 * 13)
orb = G.stabilizer(0).orbits(); sizes = sorted(len(o) for o in orb); print("rank =", len(orb), "suborbits", sizes)
o10 = [o for o in orb if len(o) == 10][0]
P = [np.array(x) for x in g]; tr = {0: np.arange(n)}; q = [0]
while q:
    p = q.pop(0)
    for Q in P:
        y = Q[p]
        if y not in tr: tr[y] = Q[tr[p]]; q.append(y)
A = np.zeros((n, n))
for x in range(n):
    for y in o10: A[x, tr[x][y]] = 1
assert np.allclose(A, A.T)
u, c = np.unique(np.round(np.linalg.eigvalsh(A), 6), return_counts=True)
print("A_1 spectrum {eig: mult} =", {float(a): int(b) for a, b in zip(u, c)})
print("sum mult =", int(c.sum()), "; #distinct eigenspaces == rank:", len(u) == len(orb), "=> pi multiplicity-free, 5 irreducible constituents")
