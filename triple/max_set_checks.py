"""Numerical checks of the statements used in max_nonopposite_proof.md, on the Ree-Tits octagon GO(2,4).

These do not prove anything about other GO(2,4); they check that each step of the hand proof holds on the
2925 line neighbourhoods of the known octagon.
  1. eigenvalues of A_4 = P_4(A_1) for A_1 in {10, 5, 1, -3, -5}, exact (Fractions): -16 exactly at 5 and -3.
  2. (A_2 - A_1 - 5I) 1_C = 1 and A_4 1_C = 16(1 - 1_C) for every line neighbourhood C.
  3. inside C: (n1, n2, n3, n4) = (a, a + 6, 20 - 2a, 0) with a in {2, 10};
     outside C: (m1, m2, m3, m4) = (m, m + 1, 10 - 2m, 16).
  4. Lemma 2: for x in C and a line M on x that is not contained in C, the number c_M(x) of points of C at
     distance 3 from x whose geodesic from x starts along M equals 4 - 2a(x).
Usage: python triple/max_set_checks.py     (writes triple/out_orbits/max_set_checks_output.txt)
"""
import os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'flags'))
import numpy as np
import go24_common as gc

OUT = os.path.join(HERE, 'out_orbits'); os.makedirs(OUT, exist_ok=True)
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); log.append(s); print(s, flush=True)

# 1. exact eigenvalues
b, c = [10, 8, 8, 8, 0], [0, 1, 1, 1, 5]
a_ = [10 - b[i] - c[i] for i in range(5)]
kv = [1]
for i in range(4): kv.append(kv[-1] * b[i] // c[i + 1])
P = {}
for th in [10, 5, 1, -3, -5]:
    u = [Fraction(1), Fraction(th, 10)]
    for i in range(1, 4):
        u.append(((th - a_[i]) * u[i] - c[i] * u[i - 1]) / b[i])
    P[th] = [kv[i] * u[i] for i in range(5)]
say('P_i(theta), i = 0..4:', {th: [str(x) for x in v] for th, v in P.items()})
assert [th for th in P if P[th][4] == -16] == [5, -3] and min(P[th][4] for th in P) == -16

# 2-4 on all line neighbourhoods
Tg, _ = gc.generators()
D = gc.point_distances(Tg)
lines = gc.lines_from_distances(D)
n = D.shape[0]
A = [(D == i).astype(np.int64) for i in range(5)]
line_of = {}
for j, L in enumerate(lines):
    for x in L:
        line_of.setdefault(x, []).append(j)
Cm = np.zeros((n, len(lines)), dtype=np.int64)
for j, L in enumerate(lines):
    S = set(L)
    for x in L: S |= set(np.nonzero(D[x] == 1)[0].tolist())
    Cm[list(S), j] = 1
N = [Ai @ Cm for Ai in A]            # N[i][x, j] = |C_j cap Gamma_i(x)|
ok2 = np.array_equal(N[2] - N[1] - 5 * Cm, np.ones_like(Cm)) and np.array_equal(N[4], 16 * (1 - Cm))
say(f'(A2 - A1 - 5I)1_C = 1 and A4 1_C = 16(1 - 1_C) on all {len(lines)} line neighbourhoods: {ok2}')
ins = Cm == 1
a = N[1]
ok3 = (set(np.unique(a[ins]).tolist()) == {2, 10} and np.array_equal(N[2][ins], a[ins] + 6)
       and np.array_equal(N[3][ins], 20 - 2 * a[ins]) and (N[4][ins] == 0).all()
       and np.array_equal(N[2][~ins], a[~ins] + 1) and np.array_equal(N[3][~ins], 10 - 2 * a[~ins])
       and (N[4][~ins] == 16).all())
say(f'inside (a, a+6, 20-2a, 0) with a in {{2,10}}, outside (m, m+1, 10-2m, 16): {ok3}; '
    f'outside m1 values {sorted(np.unique(a[~ins]).tolist())}')

# 4. Lemma 2 on line neighbourhoods 0..2924, every point x of C, every line M on x not inside C
bad, tested = 0, 0
for j in range(len(lines)):
    C = set(np.nonzero(Cm[:, j])[0].tolist())
    for x in C:
        ax = int(a[x, j])
        for mi in line_of[x]:
            M = lines[mi]
            if all(y in C for y in M):
                continue
            others = [y for y in M if y != x]
            # z at distance 3 from x behind M  <=>  d(x,z) = 3 and d(y,z) = 2 for some y in M \ {x}
            cM = sum(1 for z in C if D[x, z] == 3 and any(D[y, z] == 2 for y in others))
            tested += 1
            bad += cM != 4 - 2 * ax
say(f'Lemma 2 (c_M(x) = 4 - 2a(x)) tested on {tested} (C, x, M) triples: {tested - bad} hold, {bad} fail')
open(os.path.join(OUT, 'max_set_checks_output.txt'), 'w', newline='\n').write('\n'.join(log) + '\n')
