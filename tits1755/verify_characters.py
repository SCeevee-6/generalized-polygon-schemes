"""Character values of the 5 eigenspace representations of T=2F4(2)' on 1755 points, from the ATLAS generators.
chi_i(g) = tr(E_i P_g) = (1/n) sum_j Q[j,i] N_j(g),  N_j(g) = #{x : d(x, x^g) = j}  (d = distance in the point graph).
Checks: (a) every chi_i(g) is a rational integer; (b) chi_i(1) = 1,351,650,675,78; (c) sample mean of chi_i(g)^2 over
random g (product replacement) is close to <chi_i,chi_i> = 1; (d) sum of squares of the full degree list = |T|.
Output: verify_characters_output.txt, character_values_by_order.csv.
"""
import re, csv, random, urllib.request
from fractions import Fraction as F
from collections import defaultdict, deque
import numpy as np
n = 1755
B = "https://brauer.maths.qmul.ac.uk/Atlas/exc/TF42/gap/TF42G1-p1755B0.g"
def parse(t):
    body = re.sub(r"\s+", "", t.split(":=")[1].strip().rstrip(";"))
    p = list(range(n))
    for c in re.findall(r"\(([^)]+)\)", body):
        ix = [int(x) - 1 for x in c.split(",")]
        for k in range(len(ix)): p[ix[k]] = ix[(k + 1) % len(ix)]
    return np.array(p)
gens = [parse(urllib.request.urlopen(B + str(i), timeout=30).read().decode()) for i in (1, 2)]
def orbital(y):
    seen = {(0, y)}; q = deque([(0, y)])
    while q:
        a, b = q.popleft()
        for g in gens:
            e = (int(g[a]), int(g[b]))
            if e not in seen: seen.add(e); q.append(e)
    return seen
adj = None; done = set()
for y in range(1, n):
    if y in done: continue
    E = orbital(y); done |= {b for a, b in E if a == 0}
    if len(E) == 10 * n: adj = E; break
nb = defaultdict(list)
for a, b in adj: nb[a].append(b)
D = np.full((n, n), -1, dtype=np.int8)
for s in range(n):
    D[s, s] = 0; q = deque([s])
    while q:
        a = q.popleft()
        for b in nb[a]:
            if D[s, b] < 0: D[s, b] = D[s, a] + 1; q.append(b)
assert D.max() == 4 and (D >= 0).all()
Q = [[F(x) for x in r[1:]] for r in list(csv.reader(open("out_array/Q_matrix.csv")))[1:]]
def chars(g):
    N = np.bincount(D[np.arange(n), g], minlength=5)
    return [sum(Q[j][i] * int(N[j]) for j in range(5)) / n for i in range(5)]
def order(g):
    o, h = 1, g.copy()
    while not (h == np.arange(n)).all(): h = g[h]; o += 1
    return o
random.seed(1)
X = [g.copy() for g in gens] * 2 + [gens[0][gens[1]]] * 2
def rnd():
    i, j = random.sample(range(len(X)), 2)
    X[i] = X[i][X[j]] if random.random() < .5 else X[j][X[i]]; return X[i]
for _ in range(500): rnd()
out = ["chi(1) = %s" % [int(c) for c in chars(np.arange(n))]]
NS = 4000; sq = [0] * 5; byo = defaultdict(set); nonint = 0
for _ in range(NS):
    g = rnd(); c = chars(g)
    nonint += sum(x.denominator != 1 for x in c)
    for i in range(5): sq[i] += c[i] ** 2
    byo[order(g)].add(tuple(int(x) for x in c))
out.append("samples=%d  non-integer values=%d" % (NS, nonint))
out.append("mean chi_i^2 (expect ~1): %s" % [round(float(s) / NS, 3) for s in sq])
out.append("element orders seen: %s" % sorted(byo))
deg = [1,26,26,27,27,78,300,325,351,351,351,624,624,650,675,702,702,1300,1300,1728,2048,2048]
out.append("degree list: %d chars, sum of squares = %d, |T| = 17971200, equal: %s" % (len(deg), sum(d*d for d in deg), sum(d*d for d in deg) == 17971200))
with open("verify_characters_output.txt", "w") as f: f.write("\n".join(out) + "\n")
with open("character_values_by_order.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["element_order", "chi_1", "chi_351", "chi_650", "chi_675", "chi_78"])
    for o in sorted(byo):
        for v in sorted(byo[o]): w.writerow([o, *v])
print("\n".join(out))
