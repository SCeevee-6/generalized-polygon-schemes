
"""verify_fusions.py -- regenerates every computed claim in the fusion note.
Requires: python>=3.9, sympy.  Run:  python verify_fusions.py
"""
from fractions import Fraction as Fr
import itertools, sympy as sp

s, t, m = sp.symbols('s t m', positive=True)
red = lambda e: sp.expand(sp.rem(sp.expand(e), m**2 - 2*s*t, m))      # octagon: m^2 = 2st

# ---------------------------------------------------------------- eigenmatrices (closed forms)
P8 = sp.Matrix([
 [1, s*(t+1), s**2*t*(t+1), s**3*t**2*(t+1), s**4*t**3],
 [1, s-1+m, s*t+m*(s-1)-s, s*t*(s-1)-m*s, -s**2*t],
 [1, s-1,  -s*(t+1), -s*t*(s-1), s**2*t],
 [1, s-1-m, s*t-m*(s-1)-s, s*t*(s-1)+m*s, -s**2*t],
 [1, -t-1, t*(t+1), -t**2*(t+1), t**3]])
P6 = sp.Matrix([                                                          # hexagon: m^2 = s t
 [1, s*(t+1), s**2*t*(t+1), s**3*t**2],
 [1, s-1+m, m*s-m-s, -m*s],
 [1, s-1-m, -m*s+m-s, m*s],
 [1, -t-1, t*(t+1), -t**2]])

def check_recursion(P, d, msq):
    k = s*(t+1); b = [k] + [s*t]*(d-1); c = [None] + [1]*(d-1) + [t+1]
    a = [0] + [sp.simplify(k-b[i]-c[i]) for i in range(1, d)] + [sp.simplify(k-c[d])]
    ok = True
    for e in range(d+1):
        th = P[e, 1]
        for j in range(d+1):
            lhs = th*P[e, j]
            rhs = (b[j-1]*P[e, j-1] if j >= 1 else 0) + a[j]*P[e, j] + (c[j+1]*P[e, j+1] if j < d else 0)
            r = sp.expand(sp.rem(sp.expand(lhs-rhs), m**2-msq, m))
            ok &= (r == 0)
    sums = [sp.expand(sp.rem(sp.expand(sum(P[e, j] for j in range(d+1))), m**2-msq, m)) for e in range(1, d+1)]
    return ok, sums

def partitions(lst):
    if not lst: yield []; return
    f, rest = lst[0], lst[1:]
    for p in partitions(rest):
        yield [[f]] + p
        for i in range(len(p)): yield p[:i] + [[f] + p[i]] + p[i+1:]

def inter_nums(n, sv, tv):
    d = n//2; D = d+1; k = Fr(sv)*(tv+1)
    b = [k] + [Fr(sv)*tv]*(d-1); c = [None] + [Fr(1)]*(d-1) + [Fr(tv+1)]
    a = [Fr(0)] + [k-b[i]-c[i] for i in range(1, d)] + [k-c[d]]
    L1 = [[Fr(0)]*D for _ in range(D)]
    for j in range(D):
        if j >= 1: L1[j][j-1] = b[j-1]
        L1[j][j] = a[j]
        if j < d: L1[j][j+1] = c[j+1]
    I = [[Fr(int(i == j)) for j in range(D)] for i in range(D)]
    mm = lambda X, Y: [[sum(X[i][l]*Y[l][j] for l in range(D)) for j in range(D)] for i in range(D)]
    L = [I, L1]
    for i in range(1, d):
        M = mm(L1, L[i])
        L.append([[(M[x][y]-b[i-1]*L[i-1][x][y]-a[i]*L[i][x][y])/c[i+1] for y in range(D)] for x in range(D)])
    return L            # L[i][j][k] = p^k_{ij}

def has_fusion(L, p):
    return all(len({sum(L[i][j][kk] for i in Ba for j in Bb) for kk in Bc}) == 1
               for Ba in p for Bb in p for Bc in p)

def grid(n, lo=2, hi=30):
    D = n//2 + 1
    PP = [p for p in partitions(list(range(1, D))) if 1 < len(p) < D-1]
    out = {}
    for sv in range(lo, hi+1):
        for tv in range(lo, hi+1):
            L = inter_nums(n, sv, tv)
            f = [p for p in PP if has_fusion(L, p)]
            if f: out[(sv, tv)] = f
    return out

if __name__ == "__main__":
    print("== 1. closed-form eigenmatrices satisfy the 3-term recursion and sum rule (exact) ==")
    print("   octagon:", check_recursion(P8, 4, 2*s*t))
    print("   hexagon:", check_recursion(P6, 3, s*t))
    x = sp.symbols('x')
    print("== 2. grid test, all integer (s,t) in [2,30]^2 (direct closure of intersection numbers) ==")
    for n in (6, 8):
        g = grid(n)
        print(f"   n={n}: #(s,t) with a proper fusion = {len(g)}; all have s=t: {all(a==b for a,b in g)}; "
              f"every s=t present: {all((q,q) in g for q in range(2,31))}; fusions seen: {sorted({str(v) for v in g.values()})}")
    print("== 3. hexagon cubic: t s^3-(2t+1)s^2-t(2t+1)s-t^2(t+1)^2 ==")
    f = lambda sv, tv: tv*sv**3-(2*tv+1)*sv**2-tv*(2*tv+1)*sv-tv**2*(tv+1)**2
    # squaring m(s-1)=s+t^2+t with m^2=st gives s t (s-1)^2-(s+t^2+t)^2 = f(s,t):
    print("   identity  s*t*(s-1)^2-(s+t^2+t)^2 - f(s,t) =", sp.expand(s*t*(s-1)**2-(s+t**2+t)**2-f(s,t)))
    print("   f(t+2) =", sp.expand(f(t+2, t)), "; f(t+3) =", sp.expand(f(t+3, t)))
    print("   f(t+2)<0 and f(t+3)>0 for all integers t>=2:",
          all(f(q+2, q) < 0 and f(q+3, q) > 0 for q in range(2, 2000)))
    print("== 4. Tits group action on 1755 points = octagon (s,t)=(2,4) ==")
    Pn = P8.subs({s: 2, t: 4, m: 4}); kv = [1, 10, 80, 640, 1024]; v = sum(kv)
    mult = [sp.Rational(v)/sum(sp.Rational(Pn[e, j]**2, kv[j]) for j in range(5)) for e in range(5)]
    print("   v =", v, "= (1+s)(1+st)(1+s^2t^2) =", 3*9*65, "; |T| = 2^11*3^3*5^2*13 =", 2**11*27*25*13,
          "; |H| = 2*2^8*5*4 =", 2*2**8*5*4, "; |T|/|H| =", (2**11*27*25*13)//(2*2**8*5*4))
    print("   P(2,4) =", Pn.tolist()); print("   multiplicities =", mult)
    kr = lambda e, f_, g: sp.Rational(mult[e]*mult[f_], v)*sum(sp.Rational(Pn[e,j]*Pn[f_,j]*Pn[g,j], kv[j]**2) for j in range(5))
    print("   vanishing Krein parameters (e<=f, g, all >=1):",
          [(e, f_, g) for e in range(1,5) for f_ in range(e,5) for g in range(1,5) if kr(e,f_,g) == 0],
          "; min over all =", min(kr(e,f_,g) for e in range(5) for f_ in range(5) for g in range(5)))
    L = inter_nums(8, 2, 4)
    print("   B_i = (p^k_ij)_{j,k}:")
    for i in range(1, 5): print("    B_%d ="%i, [[int(L[i][j][k]) for k in range(5)] for j in range(5)])
    print("   symmetric p^k_ij=p^k_ji:", all(L[i][j][k]==L[j][i][k] for i in range(5) for j in range(5) for k in range(5)))
    print("   valency identity sum_k p^k_ij k_k = k_i k_j:", all(sum(L[i][j][k]*kv[k] for k in range(5))==kv[i]*kv[j] for i in range(5) for j in range(5)))
    print("== 5. SRG from s=t ==")
    q = sp.symbols('q', positive=True)
    ev = [sp.factor((P6[e,1]+P6[e,2]).subs({s:q, t:q, m:q})) for e in range(4)]
    print("   hexagon s=t=q: eigenvalues of A_1+A_2:", ev)
    K = ev[0]; r_, s_ = ev[1], ev[2] if ev[2] != ev[1] else ev[3]
    mu = sp.factor(K + r_*s_); lam = sp.factor(mu + r_ + s_)
    print("   => SRG(v,k,lam,mu) with v=(q^6-1)/(q-1), k=", sp.factor(K), ", lam=", sp.expand(lam), ", mu=", mu)
    for qq in (2, 3, 4, 5):
        Lh = inter_nums(6, qq, qq); B = [1, 2]
        lam_ = sum(Lh[i][j][1] for i in B for j in B); mu_ = sum(Lh[i][j][3] for i in B for j in B)
        k1 = int(Lh[1][1][0]); k2 = int(Lh[2][2][0]); k3 = int(Lh[3][3][0]); vh = 1+k1+k2+k3   # k_i = p^0_{ii}
        print(f"   direct, GH({qq},{qq}): SRG({vh},{k1+k2},{lam_},{mu_})")
    for qq in (2, 3):
        Lo = inter_nums(8, qq, qq)
        vo = 1 + sum(int(Lo[i][i][0]) for i in range(1, 5))                        # v = 1 + sum k_i , k_i = p^0_{ii}
        print(f"   direct, octagon formal (s,t)=({qq},{qq}): v={vo}, Gamma_4: valency k_4=p^0_44={Lo[4][4][0]}, "
              f"lambda=p^4_44={Lo[4][4][4]}, mu=p^1_44={Lo[4][4][1]}   [expected SRG((s^8-1)/(s-1), s^7, s^7-s^6, s^7-s^6)]")
