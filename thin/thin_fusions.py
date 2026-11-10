"""Thin generalized octagons of order (1,t) and (s,1). Order (1,t): point graph = incidence graph of GQ(t,t),
array {t+1,t,t,t;1,1,1,t+1}. Order (s,1): array {2s,s,s,s;1,1,1,2}. Lists every closing fusion of the 4 classes
and, for each block B of a proper fusion, its valency and the numbers sum_{i,j in B} p^h_ij (h = 1..4), i.e. the
number of common B-neighbours of two points at distance h.
Exact Fractions via scheme_from_array. Usage (from the repo root): python thin/thin_fusions.py [--out DIR]; default DIR is thin/out_thin"""
import sys, os, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "tits1755"))
import scheme_from_array as S

def run(s, t):
    """Point scheme of a generalized octagon of order (s, t), array {s(t+1),st,st,st; 1,1,1,t+1}."""
    S.S, S.T, S.D = s, t, 4
    S.B = [Fr(s*(t+1)), Fr(s*t), Fr(s*t), Fr(s*t), Fr(0)]
    S.C = [Fr(0), Fr(1), Fr(1), Fr(1), Fr(t+1)]
    S.R = range(5)
    kk = [Fr(1)]
    for i in range(4): kk.append(kk[i]*S.B[i]/S.C[i+1])
    S.N = sum(kk)
    a, k = S.parameters()
    p = S.intersection_numbers(a); S.check_p(p, k)
    th, P, m, Q = S.spectrum(a, k)
    q = S.krein(P, m, k)
    qp = S.q_polynomial_orderings(q)
    fus = [(pt, S.fusion_violations(p, pt)) for pt in S.set_partitions([1,2,3,4])]
    closing = [pt for pt, v in fus if v == 0]
    proper = [pt for pt in closing if 1 < len(pt) < 4]
    return dict(s=s, t=t, n=int(S.N), k=[int(x) for x in k], theta=th, m=[int(x) for x in m],
                qpoly=len(qp), closing=closing, proper=proper,
                graphs=[describe(p, k, pt) for pt in proper], p=p, a=a)


def describe(p, k, pt):
    """For each block of a closing fusion, the graph it defines: valency, connectivity, and SRG parameters when connected
    with exactly 2 nontrivial classes, or the parameters of its components when it is a disjoint union of SRGs."""
    out = []
    for blk in pt:
        val = sum(k[i] for i in blk)
        lam = {h: sum(p[i][j][h] for i in blk for j in blk) for h in range(1, 5)}
        out.append(dict(block=blk, valency=int(val), common_neighbours={h: int(v) for h, v in lam.items()}))
    return out

if __name__ == "__main__":
    out = sys.argv[sys.argv.index("--out")+1] if "--out" in sys.argv else os.path.join(HERE, "out_thin")
    os.makedirs(out, exist_ok=True)
    res = []
    cases = [(1, t) for t in (2, 8, 18, 32, 50)] + [(s, 1) for s in (2, 8)]
    for s, t in cases:
        try:
            r = run(s, t)
        except AssertionError as e:
            r = dict(s=s, t=t, error=str(e))
        res.append({x: y for x, y in r.items() if x not in ("p", "a")})
        print({x: y for x, y in res[-1].items()})
    json.dump(res, open(f"{out}/thin_fusions.json", "w"), indent=1, default=str)
