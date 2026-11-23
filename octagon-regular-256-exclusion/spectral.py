from base import *
from math import comb
from functools import lru_cache
def inertia(co,t):
 p=list(reversed(co));z=[sum(p[k]*comb(k,j)*t**(k-j) for k in range(j,len(p))) for j in range(len(p))]
 nz=[a for a in z if a];return sum((a<0)!=(b<0) for a,b in zip(nz,nz[1:])),next(i for i,a in enumerate(z) if a)
@lru_cache(None)
def allowed(co):
 n,_=inertia(co,-2);m,z=inertia(co,-1);l,e=inertia(co,-4)
 # Scale eigenvalues by 20 for a rational upper endpoint.
 u,_=inertia([x*20**k for k,x in enumerate(co)],69)
 return n==m+z and l+e==16 and not u
def matrix(rs):
 M=np.zeros((len(rs),16,16),dtype=np.int64)
 for e in edges:
  i,j=map(vs.index,e);M[:,i,j]=M[:,j,i]=1-2*signs(rs,e)
 return M
def loads(M):
 P=M@M@M-8*M
 return np.array([P[:,i]+M[:,i,j,None]*P[:,j] for i,j in (tuple(map(vs.index,e)) for e in edges)]).transpose(1,0,2)
def polynomial(R,I):
 E=R;F=I;tr=[];co=[np.ones(len(R),np.int64)]
 for k in range(1,17):
  tr.append(np.trace(E,axis1=1,axis2=2));s=-sum(a*b for a,b in zip(tr,co[::-1]));assert np.all(s%k==0);co.append(s//k)
  if k<16:
   if I is None:E=E@R
   else:E,F=E@R-F@I,E@I+F@R
 return np.array(co).T.tolist()
if __name__=='__main__':
 rs=np.loadtxt('reps.txt',dtype=np.int32)[:,0];B=loads(matrix(np.array([0])))[0];B.astype(np.int8).tofile('B.bin')
 passed=[];stage=0
 for start in range(0,len(rs),128):
  M=matrix(rs[start:start+128]);L=loads(M)
  for j,p in enumerate(polynomial(M,None)):
   if allowed(tuple(p)):
    stage+=1
    if (B+np.abs(L[j])).max()<=16:passed.append(start+j)
 ids=np.memmap('ids.bin',np.int32,'r');good=np.flatnonzero(np.isin(ids,passed)).astype(np.int32);good.tofile('good.bin')
 np.savetxt('first.txt',rs[passed],fmt='%d')
 with open('loads.bin','wb') as f:
  for g in np.array_split(good,150):loads(matrix(g)).astype(np.int8).tofile(f)
 print(len(rs),stage,len(passed),len(good),flush=True)
