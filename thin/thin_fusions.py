"""Thin generalized octagon GO(1,t): point graph = incidence graph of GQ(t,t), array {t+1,t,t,t;1,1,1,t+1}.
Exact Fractions via scheme_from_array. Usage: python thin_fusions.py [--out DIR]"""
import sys, os, json
from fractions import Fraction as Fr
import scheme_from_array as S

def run(t):
    S.S, S.T, S.D = 1, t, 4
    S.B = [Fr(t+1), Fr(t), Fr(t), Fr(t), Fr(0)]
    S.C = [Fr(0), Fr(1), Fr(1), Fr(1), Fr(t+1)]
    S.R = range(5)
    a, k = None, None
    S.N = 1
    kk = [Fr(1)]
    for i in range(4): kk.append(kk[i]*S.B[i]/S.C[i+1])
    S.N = sum(kk)
    a, k = S.parameters()
    p = S.intersection_numbers(a); S.check_p(p, k)
    th, P, m, Q = S.spectrum(a, k)
    q = S.krein(P, m, k)
    qp = S.q_polynomial_orderings(q)
    fus = [(pt, S.fusion_violations(p, pt)) for pt in S.set_partitions([1,2,3,4])]
    return dict(t=t, n=int(S.N), k=[int(x) for x in k], theta=th, m=[int(x) for x in m],
                qpoly=len(qp), closing=[pt for pt, v in fus if v == 0], p=p, a=a)

if __name__ == "__main__":
    out = sys.argv[sys.argv.index("--out")+1] if "--out" in sys.argv else "out_thin"
    os.makedirs(out, exist_ok=True)
    res = []
    for t in (2, 8, 18, 32, 50):
        try:
            r = run(t)
        except AssertionError as e:
            r = dict(t=t, error=str(e))
        res.append({x: y for x, y in r.items() if x not in ("p", "a")})
        print({x: y for x, y in res[-1].items()})
    json.dump(res, open(f"{out}/thin_fusions.json", "w"), indent=1, default=str)
