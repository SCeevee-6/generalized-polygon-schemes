"""Gröbner-basis fusion cross-check (Bannai-Muzychuk column criterion) for the octagon point scheme, m^2=2st, t=m^2/(2s).
Last column (feasibility) is now computed: each solution is tested for s,t >= 2 integer with 2st a perfect square
(Feit-Higman). Isolated points are tested exactly (rounded floats are only a pre-filter, then exact check on s,t);
curves are tested by substituting M^2=2*S*t and scanning integer s in [2,SMAX], t from the curve equation.
Runtime about 7 minutes (sympy Groebner bases). Writes bannai_muzychuk_crosscheck.csv to the working directory.
"""
import sympy as sp
import itertools
import numpy as np
import pandas as pd

s, m = sp.symbols('s m', positive=True)
S, M = sp.symbols('S M')
z = sp.symbols('z')

t = m**2/(2*s)
k = s*(t+1)
b = [s*(t+1), s*t, s*t, s*t]
c = [None, 1, 1, 1, t+1]
a = [0]+[sp.simplify(k-b[i]-c[i]) for i in range(1,4)] + [sp.simplify(k-c[4])]
thetas = [k, s-1+m, s-1, s-1-m, -t-1]

def P_col(th):
    A=[sp.Integer(1), th]
    for i in range(1,4):
        A.append(sp.simplify((th*A[i]-b[i-1]*A[i-1]-a[i]*A[i])/c[i+1]))
    return A

P = sp.Matrix([P_col(th) for th in thetas])

def partitions(lst):
    if not lst: yield []; return
    f, rest = lst[0], lst[1:]
    for p in partitions(rest):
        yield [[f]]+p
        for i in range(len(p)):
            yield p[:i]+[[f]+p[i]]+p[i+1:]

parts=[p for p in partitions([1,2,3,4]) if 1<len(p)<4]

def numer(x): return sp.numer(sp.together(sp.simplify(x))).subs({s:S,m:M})

def pos_definite(f):
    """True if irreducible factor f has no zero in S>0,M>0 (checked by dense scan + root test)."""
    for Mv in np.linspace(0.02,80,500):
        if f.has(S):
            rts=sp.Poly(f.subs(M,sp.Rational(Mv).limit_denominator(1000)),S).nroots() if sp.Poly(f.subs(M,1),S).degree()>0 else []
            if any(abs(sp.im(r))<1e-9 and sp.re(r)>1e-9 for r in rts): return False
        else:
            if abs(f.subs(M,Mv))<1e-6: return False
    if not f.has(S):
        return not any(abs(sp.im(r))<1e-9 and sp.re(r)>1e-9 for r in sp.Poly(f,M).nroots())
    return True

def solve_branch(basis):
    """basis: list of polys (no z). returns set of real positive solution descriptors."""
    out=set()
    facs=[]
    for g in basis:
        fl=[f for f,_ in sp.factor_list(g)[1] if not f.is_number and not pos_definite(f)]
        if not fl: return set()
        facs.append(fl)
    for combo in itertools.product(*facs):
        combo=list(set(combo))
        G=sp.groebner(combo,S,M,order='lex',domain='QQ')
        if list(G.exprs)==[1]: continue
        gs=list(G.exprs)
        uni=[g for g in gs if not g.has(S)]
        if uni:
            for f,_ in sp.factor_list(uni[0])[1]:
                for mr in sp.Poly(f,M).nroots(n=30):
                    if abs(sp.im(mr))>1e-9 or sp.re(mr)<=1e-9: continue
                    mv=float(sp.re(mr))
                    ps=[sp.Poly(g.subs(M,mv),S) for g in gs if g.has(S)]
                    if not ps: out.add(('M=',round(mv,5))); continue
                    for sr in ps[0].nroots(n=30):
                        if abs(sp.im(sr))<1e-9 and sp.re(sr)>1e-9:
                            sv=float(sp.re(sr))
                            if all(abs(float(sp.N(g.subs({S:sv,M:mv}))))<1e-5*(1+sv**4+mv**6) for g in gs):
                                out.add((round(sv,5),round(mv,5),'t=%.5f'%(mv**2/(2*sv))))
        else:
            out.add(('CURVE',tuple(str(sp.factor(g)) for g in gs)))
    return out


SMAX=5000
def _is_sq(n): r=sp.integer_nthroot(int(n),2); return r[1]

def int_ok(sv,tv):
    """sv,tv floats from a numeric root; True if (s,t) is a thick-octagon candidate: integers >=2, 2st a perfect square."""
    s_,t_=round(sv),round(tv)
    if abs(sv-s_)>1e-4 or abs(tv-t_)>1e-4 or s_<2 or t_<2: return False
    return _is_sq(2*s_*t_)

def feasible(pt):
    """pt: descriptor from solve_branch or 'IDENTICAL'. Returns list of integer (s,t) candidates, or a family string."""
    if pt=='IDENTICAL' or pt[0]=='M=': return []
    if pt[0]=='CURVE':
        out=[]
        for expr in pt[1]:
            f=sp.sympify(expr,locals={'S':S,'M':M})
            if not f.has(M):      # s fixed, t free
                for c in sp.solve(f,S):
                    if c.is_integer and c>=2: out.append('s=%s, any t>=2 with 2st square'%c)
                continue
            for Mexpr in sp.solve(f,M):      # M=M(S) on the curve; t=M^2/(2S)
                T=sp.lambdify(S,sp.simplify(Mexpr**2/(2*S)),'sympy')
                for sv in range(2,SMAX+1):
                    tv=sp.nsimplify(T(sp.Integer(sv)))
                    if tv.is_integer and tv>=2 and _is_sq(2*sv*tv): out.append((sv,int(tv)))
        return out
    sv,mv=pt[0],pt[1]
    return [(round(sv),round(mv**2/(2*sv)))] if int_ok(sv,mv**2/(2*sv)) else []

final2={}
for p in parts:
    r=len(p)
    colv={e:[sp.together(sum(P[e,j] for j in blk)) for blk in p] for e in range(5)}
    pts=set()
    for q in partitions([1,2,3,4]):
        if len(q)!=r: continue
        E=[]
        for blk in q:
            for x in blk[1:]:
                E+=[numer(u-v) for u,v in zip(colv[blk[0]],colv[x])]
        E=[sp.expand(e) for e in E if sp.expand(e)!=0]
        if not E: pts.add('IDENTICAL'); continue
        G=sp.groebner(E+[1-z*S*M],z,S,M,order='lex',domain='QQ')
        if list(G.exprs)==[1]: continue
        pts|=solve_branch([g for g in G.exprs if not g.has(z)])
    final2[str(p)]=pts
    print(p,'->',sorted(pts,key=str) if pts else 'NO real positive (s,m)',flush=True)

rows=[]
for p,pts in final2.items():
    cand=sorted({c for pt in pts for c in feasible(pt)},key=str)
    rows.append((p, '; '.join(map(str,sorted(pts,key=str))) if pts else 'none', str(cand) if cand else 'none'))
df=pd.DataFrame(rows,columns=['class_partition','real_positive_(s,m)_solutions_[m=sqrt(2st),t=m^2/(2s)]','feasible_thick_octagon_(s,t>=2_integer,2st_square)'])
df.to_csv('bannai_muzychuk_crosscheck.csv',index=False)
print(df.to_string())