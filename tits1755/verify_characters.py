"""Exact class-by-class values of the 5 eigenspace characters of T = 2F4(2)' on 1755 points.

For each of the 22 conjugacy classes a representative g is built from the ATLAS class words
(atlas_words.py) in the standard generators a = b11, b = b21 of TF42G1-p1755B0.
chi_e(g) = tr(E_e P_g) = (1/n) sum_j Q[j][e] N_j(g),  N_j(g) = #{x : d(x, x^g) = j},
with Q from out_array/Q_matrix.csv (QP = nI) and d the distance in the valency-10 orbital graph.
pi(g) = N_0(g) is the permutation character.

If gap/tits_chartable.csv (written by gap/tits_chars.g from the GAP character table library)
exists, every value is compared with the library table: element orders, pi, and each chi_e
against a library irreducible; <chi_e, chi_e> = 1 and <pi, pi> = 5 are recomputed with the
library centralizer orders.

Run from tits1755/:  python verify_characters.py
Output: verify_characters_output.txt, eigen_characters_by_class.csv
"""
import csv, os, re, urllib.request
from fractions import Fraction as F
from collections import defaultdict, deque
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_words import ATLAS_CLASS_WORDS, parse_word, evaluate

n = 1755
URL = "https://brauer.maths.qmul.ac.uk/Atlas/exc/TF42/gap/TF42G1-p1755B0.g"


def load(i):
    path = os.path.join("cache", "TF42G1-p1755B0.g%d" % i)
    if os.path.exists(path):
        return open(path).read()
    t = urllib.request.urlopen(URL + str(i), timeout=60).read().decode()
    os.makedirs("cache", exist_ok=True)
    open(path, "w").write(t)
    return t


def parse(t):
    body = re.sub(r"\s+", "", t.split(":=")[1].strip().rstrip(";"))
    p = list(range(n))
    for c in re.findall(r"\(([^)]+)\)", body):
        ix = [int(x) - 1 for x in c.split(",")]
        for k in range(len(ix)):
            p[ix[k]] = ix[(k + 1) % len(ix)]
    return np.array(p)


a, b = parse(load(1)), parse(load(2))
ident = np.arange(n)
mul = lambda x, y: y[x]            # x then y: i -> y[x[i]]  (GAP convention i^(xy) = (i^x)^y)
inv = lambda x: np.argsort(x)


def order(g):
    o, h = 1, g.copy()
    while not (h == ident).all():
        h = g[h]; o += 1
    return o


# distance matrix of the valency-10 orbital graph (BFS from every point)
def orbital(y):
    seen = {(0, y)}; q = deque([(0, y)])
    while q:
        u, v = q.popleft()
        for g in (a, b):
            e = (int(g[u]), int(g[v]))
            if e not in seen:
                seen.add(e); q.append(e)
    return seen

adj = None; done = set()
for y in range(1, n):
    if y in done:
        continue
    E = orbital(y); done |= {v for u, v in E if u == 0}
    if len(E) == 10 * n:
        adj = E; break
nb = defaultdict(list)
for u, v in adj:
    nb[u].append(v)
D = np.full((n, n), -1, dtype=np.int8)
for s0 in range(n):
    D[s0, s0] = 0; q = deque([s0])
    while q:
        u = q.popleft()
        for v in nb[u]:
            if D[s0, v] < 0:
                D[s0, v] = D[s0, u] + 1; q.append(v)
assert D.max() == 4 and (D >= 0).all()
val = np.bincount(D[0], minlength=5).tolist()
assert val == [1, 10, 80, 640, 1024], val

Q = [[F(x) for x in r[1:]] for r in list(csv.reader(open(os.path.join("out_array", "Q_matrix.csv"))))[1:]]
mults = [int(Q[0][e]) for e in range(5)]
assert mults == [1, 351, 650, 675, 78], mults

out = ["ATLAS generators TF42G1-p1755B0, |Omega| = %d, valencies %s, multiplicities %s" % (n, val, mults)]
rows = []
nonint = 0
for name, w in ATLAS_CLASS_WORDS:
    g = evaluate(parse_word(w), a, b, mul, inv, ident)
    o = order(g)
    N = np.bincount(D[ident, g], minlength=5).tolist()
    chi = [sum(Q[j][e] * N[j] for j in range(5)) / n for e in range(5)]
    nonint += sum(c.denominator != 1 for c in chi)
    rows.append(dict(cls=name, word=w, order=o, N=N, pi=N[0], chi=[int(c) for c in chi]))
    assert int(re.match(r"\d+", name).group()) == o, (name, o)
out.append("22 class words evaluated; element order equals the order in the class name for every class: True")
out.append("non-integer character values: %d" % nonint)
for r in rows:
    assert sum(r["chi"]) == r["pi"]
out.append("sum_e chi_e(g) = pi(g) = fix(g) on every class: True")

gapcsv = os.path.join("gap", "tits_chartable.csv")
if os.path.exists(gapcsv):
    tab = list(csv.reader(open(gapcsv)))
    hdr = tab[0]
    cls = hdr[2:]
    meta = {r[0]: r[2:] for r in tab[1:] if r[0].startswith("#")}
    irr = {r[0]: r[2:] for r in tab[1:] if not r[0].startswith("#")}
    cent = [int(x) for x in meta["#centralizer"]]
    ords = [int(x) for x in meta["#order"]]
    pos = {c: i for i, c in enumerate(cls)}
    assert sorted(pos) == sorted(r["cls"] for r in rows), (cls, [r["cls"] for r in rows])
    out.append("GAP library table found: %d classes, class names agree with ATLAS words: True" % len(cls))
    ok_ord = all(ords[pos[r["cls"]]] == r["order"] for r in rows)
    out.append("element orders agree with library: %s" % ok_ord)
    order_T = cent[pos["1A"]]

    def ip(u, v):
        return sum(F(u[pos[r["cls"]]] * v[pos[r["cls"]]], cent[pos[r["cls"]]]) for r in rows)

    for e in range(5):
        vec = [0] * len(cls)
        for r in rows:
            vec[pos[r["cls"]]] = r["chi"][e]
        hits = [k for k, v in irr.items() if all(x == str(y) for x, y in zip(v, vec))]
        out.append("chi_%d (deg %d): <chi,chi> = %s; equals library irreducible(s) %s" % (e, mults[e], ip(vec, vec), hits or "NONE"))
    piv = [0] * len(cls)
    for r in rows:
        piv[pos[r["cls"]]] = r["pi"]
    out.append("<pi,pi> = %s (rank 5 expected)" % ip(piv, piv))
    if "#pi_library" in meta:
        out.append("pi equals library induced character 1_H^T: %s" % ([int(x) for x in meta["#pi_library"]] == piv))
    out.append("|T| from library centralizer of 1A: %d" % order_T)
else:
    out.append("gap/tits_chartable.csv not found: library comparison skipped (run gap/tits_chars.g)")

with open("verify_characters_output.txt", "w", newline="\n") as f:
    f.write("\n".join(out) + "\n")
with open("eigen_characters_by_class.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["class", "atlas_word", "element_order", "N0", "N1", "N2", "N3", "N4",
                 "pi", "chi_1", "chi_351", "chi_650", "chi_675", "chi_78"])
    for r in rows:
        wr.writerow([r["cls"], r["word"], r["order"]] + r["N"] + [r["pi"]] + r["chi"])
print("\n".join(out))
