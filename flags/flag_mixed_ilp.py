"""Search for a maximum flag family that is neither point- nor line-centred, in the potential form x(p,L) = a_p + b_L.
Run: python flag_mixed_ilp.py SECONDS generators.npy   (HiGHS; 3000 s ended without an answer)"""
import numpy as np, scipy.sparse as sps, time, sys, json
from scipy.optimize import milp, LinearConstraint, Bounds
from o2_build import octagon
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
D, lines = octagon([g.astype(np.int64) for g in np.load(sys.argv[2] if len(sys.argv) > 2 else 'generators.npy')])
_I = csr_matrix((np.ones(3*len(lines)), (np.repeat(np.arange(len(lines)), 3), lines.ravel())), shape=(len(lines), 1755))
_LL = (_I @ _I.T).toarray(); np.fill_diagonal(_LL, 0); DL = shortest_path(csr_matrix(_LL > 0), unweighted=True).astype(np.int8)
nP, nL = 1755, len(lines); TL=float(sys.argv[1])
p0=int(lines[0][0]); L0=0; L1=[L for L in range(nL) if p0 in lines[L] and L!=0][0]
A1=(D==1); A2=(D==2); B1=(DL==1); B2=(DL==2)
ia=lambda p:p; ib=lambda L:nP+L; iz=lambda p:nP+nL+p; iw=lambda p:2*nP+nL+p; ika=3*nP+nL; ikb=ika+1; nv=ikb+1
rows=[];cols=[];vals=[];lo=[];hi=[]; r=[0]
def add(t,l,h):
    for c,v in t: rows.append(r[0]); cols.append(c); vals.append(v)
    lo.append(l); hi.append(h); r[0]+=1
for L in range(nL):
    for p in lines[L]: add([(ia(p),1),(ib(L),1)],0,1)
for p in range(nP):
    t=[(ia(q),1) for q in np.flatnonzero(A2[p])]+[(ia(q),-1) for q in np.flatnonzero(A1[p])]+[(ia(p),-5),(ika,-1)]; add(t,0,0)
for L in range(nL):
    t=[(ib(M),1) for M in np.flatnonzero(B2[L])]+[(ib(M),-3) for M in np.flatnonzero(B1[L])]+[(ib(L),5),(ikb,-1)]; add(t,0,0)
add([(ika,1),(ikb,1)],1,1)
add([(ia(p0),1)],0,0); add([(ib(L0),1)],1,1); add([(ib(L1),1)],0,0)
for p in range(nP):
    add([(ia(p),1),(iz(p),-5)],-4,np.inf); add([(ia(p),1),(iw(p),5)],-np.inf,4)
add([(iz(p),1) for p in range(nP)]+[(iw(p),1) for p in range(nP)],1,np.inf)
A=sps.csr_matrix((vals,(rows,cols)),shape=(r[0],nv))
lb=np.full(nv,0.0); ub=np.full(nv,1.0)
lb[:nP]=-4; ub[:nP]=4; lb[nP:nP+nL]=-4; ub[nP:nP+nL]=5; lb[ika]=lb[ikb]=-300; ub[ika]=ub[ikb]=300
t=time.time()
res=milp(np.zeros(nv),constraints=LinearConstraint(A,lo,hi),integrality=np.ones(nv),bounds=Bounds(lb,ub),options=dict(time_limit=TL,disp=False))
out=dict(status=int(res.status),msg=res.message,sec=time.time()-t,p0=p0,L0=L0,L1=L1)
if res.x is not None:
    x=np.round(res.x).astype(int); out['a']=x[:nP].tolist(); out['b']=x[nP:nP+nL].tolist(); out['ka']=int(x[ika])
json.dump(out,open('flag_mixed_result.json','w')); print(out['status'],out['msg'],out['sec'])
