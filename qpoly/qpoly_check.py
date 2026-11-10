"""Exact Q-polynomial classification for point graphs of generalized 2d-gons of order (s,t), d=2,3,4.
No sympy.solve: Groebner bases over QQ + univariate factorisation only.
Params: d=2 -> (s,t); d=3 -> (s,m), t=m^2/s (m=sqrt(st)); d=4 -> (s,m), t=m^2/(2s) (m=sqrt(2st)).
Usage: python qpoly_check.py > qpoly_check_output.txt
"""
import itertools, sys
import sympy as sp

s, m, t_ = sp.symbols('s m t', positive=True)

def setup(d):
    if d == 2:
        t = t_; th = [s*(t+1), s-1, -t-1]; vars_ = (s, t_)
    elif d == 3:
        t = m**2/s; th = [s*(t+1), s-1+m, s-1-m, -t-1]; vars_ = (s, m)
    else:
        t = m**2/(2*s); th = [s*(t+1), s-1+m, s-1, s-1-m, -t-1]; vars_ = (s, m)
    k = s*(t+1)
    b = [k] + [s*t]*(d-1)
    c = [None] + [1]*(d-1) + [t+1]
    a = [0] + [k - b[i] - c[i] for i in range(1, d)] + [k - c[d]]
    x = sp.Symbol('x')
    v = [sp.Integer(1), x]
    for l in range(1, d):
        v.append(sp.expand(((x - a[l])*v[l] - b[l-1]*v[l-1]) / c[l+1]))
    # exact check: claimed eigenvalues are roots of (x-a_d)v_d - b_{d-1}v_{d-1}
    pchar = (x - a[d])*v[d] - b[d-1]*v[d-1]
    for e in th:
        assert sp.simplify(pchar.subs(x, e)) == 0, (d, e)
    kl = [sp.simplify(v[l].subs(x, k)) for l in range(d+1)]
    n = sp.simplify(sum(kl))
    V = [[sp.simplify(v[l].subs(x, e)) for l in range(d+1)] for e in th]
    mult = [sp.simplify(n / sum(V[i][l]**2 / kl[l] for l in range(d+1))) for i in range(d+1)]
    def q(i, j, kk):
        return sp.factor(sp.together(mult[i]*mult[j]/n * sum(V[i][l]*V[j][l]*V[kk][l]/kl[l]**2 for l in range(d+1))))
    return dict(t=t, th=th, n=n, kl=kl, mult=mult, q=q, vars=vars_)

def _live(f, vars_):
    """False if irreducible f provably has no zero with all vars>0 (monomial or single-sign coeffs)."""
    P = sp.Poly(f, *vars_)
    if P.is_monomial or P.total_degree() == 0: return False
    cs = P.coeffs()
    return not (all(c > 0 for c in cs) or all(c < 0 for c in cs))

def variety(polys, vars_):
    """Branches of {all polys = 0} on vars>0. Returns list of (GB exprs, positive real points or 'curve')."""
    fac = []
    for z in polys:
        fs = sorted({sp.factor(f) for f, _ in sp.factor_list(z, *vars_)[1] if _live(f, vars_)}, key=str)
        if not fs: return []          # some required-zero poly has no positive zero
        fac.append(fs)
    out, seen = [], set()
    for combo in itertools.product(*fac):
        key = frozenset(map(str, combo))
        if key in seen: continue
        seen.add(key)
        G = sp.groebner(list(set(combo)), *vars_, order='lex')
        ex = list(G.exprs)
        if ex == [1]: continue
        uni = [g for g in ex if g.free_symbols <= {vars_[1]}]
        if not uni:
            out.append((ex, 'curve')); continue
        pts = []
        for r in sp.Poly(uni[0], vars_[1]).real_roots():
            if r <= 0: continue
            sub = [sp.expand(g.subs(vars_[1], r)) for g in ex]
            sub = [g for g in sub if g != 0]
            if not sub: pts.append((None, r)); continue
            g0 = sp.Poly(sub[0], vars_[0])
            for r0 in sp.Poly(sp.nsimplify(g0.as_expr()), vars_[0]).nroots(n=30):
                if abs(sp.im(r0)) < 1e-20 and sp.re(r0) > 0 and all(abs(sp.N(g.subs(vars_[0], r0), 30)) < 1e-15 for g in sub):
                    pts.append((sp.nsimplify(sp.re(r0), rational=False, tolerance=1e-20), r))
        if pts: out.append((ex, pts))
    return out

def classify(d, out):
    S = setup(d); q = S['q']; vars_ = S['vars']
    out.write(f"\n=== d={d} (generalized {2*d}-gon), t={S['t']}, n={sp.factor(S['n'])} ===\n")
    out.write("eigenvalues: " + ", ".join(str(sp.simplify(e)) for e in S['th']) + "\n")
    out.write("multiplicities: " + ", ".join(str(sp.factor(mm)) for mm in S['mult']) + "\n")
    found = []
    for perm in itertools.permutations(range(1, d+1)):
        order = (0,) + perm
        pos = {e: p for p, e in enumerate(order)}
        E1 = order[1]
        zeros = []
        for i in range(1, d+1):
            for j in range(1, d+1):
                if abs(pos[i]-pos[j]) >= 2 and i < j:
                    zeros.append(sp.numer(sp.together(q(E1, i, j))))  # q^j_{1i} ~ q^i_{1j} up to m_j/m_i
        nonzero = [q(E1, order[p], order[p+1]) for p in range(1, d)]
        zeros = [sp.factor(z) for z in zeros if z != 0]
        if not zeros:
            out.write(f"order {order}: no required zeros; need nonzero: {[sp.factor(x) for x in nonzero]}\n"); found.append((order, 'generic')); continue
        br = variety(zeros, vars_)
        if not br:
            out.write(f"order {order}: IMPOSSIBLE (no positive solution)\n"); continue
        for ex, pts in br:
            if pts == 'curve':
                nzc = [sp.factor(sp.reduced(sp.numer(sp.together(x)), ex, *vars_, order='lex')[1]) for x in nonzero]
                out.write(f"order {order}: CURVE {[sp.factor(e) for e in ex]}; nonzero-conds mod curve: {nzc}\n")
                found.append((order, ex, 'curve'))
            else:
                for (r0, r1) in pts:
                    sub = {vars_[1]: r1} if r0 is None else {vars_[0]: r0, vars_[1]: r1}
                    tval = sp.simplify(S['t'].subs(sub))
                    nz_ok = all(abs(sp.N(x.subs(sub), 30)) > 1e-12 for x in nonzero) if r0 is not None else '?'
                    out.write(f"order {order}: POINT {sub} t={tval} nonzero-ok={nz_ok}\n")
                    found.append((order, sub, tval, nz_ok))
    return S, found


# ---- independent float cross-check (no Groebner): Krein-graph path test on integer (s,t) grid ----
import math, numpy as np
def scheme(d,s,t):
    k=s*(t+1); b=[k]+[s*t]*(d-1); c=[0]+[1]*(d-1)+[t+1]; a=[0]+[k-b[i]-c[i] for i in range(1,d)]+[k-c[d]]
    L=np.zeros((d+1,d+1))
    for i in range(d+1):
        L[i,i]=a[i]
        if i<d: L[i+1,i]=b[i]
        if i>0: L[i-1,i]=c[i]
    th=np.sort(np.linalg.eigvals(L).real)[::-1]
    def V(x):
        v=[1.0,x]
        for l in range(1,d): v.append(((x-a[l])*v[l]-b[l-1]*v[l-1])/c[l+1])
        return np.array(v)
    P=np.array([V(x) for x in th]); kl=P[0]; n=kl.sum()
    mu=n/np.array([(P[i]**2/kl).sum() for i in range(d+1)])
    q=lambda i,j,h: mu[i]*mu[j]/n*(P[i]*P[j]*P[h]/kl**2).sum()
    return th,q
def qpoly_orders(d,s,t,tol=1e-7):
    th,q=scheme(d,s,t); res=[]
    for perm in itertools.permutations(range(1,d+1)):
        o=(0,)+perm; E1=o[1]; ok=True
        for p in range(d+1):
            for r in range(d+1):
                val=abs(q(E1,o[p],o[r]))
                if abs(p-r)>=2 and val>tol*max(1,q(E1,E1,E1) if False else 1)*1: ok=False
                if abs(p-r)==1 and val<tol: ok=False
        if ok: res.append(o)
    return th,res

def grid_check(N=40):
    hits = []
    for d, cond in [(3, lambda a, b: math.isqrt(a*b)**2 == a*b), (4, lambda a, b: math.isqrt(2*a*b)**2 == 2*a*b)]:
        for a in range(1, N+1):
            for b in range(1, N+1):
                if cond(a, b):
                    _, res = qpoly_orders(d, a, b)
                    if res: hits.append((d, a, b, res))
    return hits

if __name__ == '__main__':
    out = sys.stdout
    for d in (2, 3, 4):
        classify(d, out)
    print('\n=== float grid check, s,t<=40, Feit-Higman admissible: Q-poly hits (d,s,t,orderings) ===')
    for h in grid_check(): print(h)
    print('GO(2,4) Q-poly orderings:', qpoly_orders(4, 2, 4)[1])
