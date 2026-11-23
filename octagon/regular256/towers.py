from base import *
from functools import reduce
from operator import xor
from collections import defaultdict
def lift(edges,n,volts):
 return [(2*a+j,2*b+(j^v)) for (a,b),v in zip(edges,volts) for j in (0,1)]
def setup(E,n):
 A=[[] for _ in range(n)]
 for i,(a,b) in enumerate(E):A[a].append((b,i));A[b].append((a,i))
 seen={0};q=[0];tree=set()
 for a in q:
  for b,i in A[a]:
   if b not in seen:seen.add(b);q.append(b);tree.add(i)
 assert len(q)==n
 ch=[i for i in range(len(E)) if i not in tree];marks=[0]*len(E)
 for j,i in enumerate(ch):marks[i]=1<<j
 return A,marks,len(ch)
def conditions(E,n,A,marks):
 m=256//n;rows=set()
 for i,(a,b) in enumerate(E):
  hits=defaultdict(list)
  def visit(v,p,r,mask):
   hits[v].append(mask)
   if r:
    for w,j in A[v]:
     if w!=p:visit(w,v,r-1,mask^marks[j])
  visit(a,b,3,0);visit(b,a,3,marks[i])
  for arr in hits.values():
   if len(arr)>m:return None
   if len(arr)>m//2:rows.add(tuple(sorted(x^arr[0] for x in arr)))
 return rows
def affine(rows,n):
 piv={}
 for v,t in rows:
  while v and v.bit_length() in piv:
   r,u=piv[v.bit_length()];v^=r;t^=u
  if v:piv[v.bit_length()]=(v,t)
  elif t:return
 def back(v,h):
  for j,(r,t) in sorted(piv.items()):
   if ((r&v).bit_count()%2)^(0 if h else t):v|=1<<(j-1)
  return v
 return back(0,False),[back(1<<j,True) for j in range(n) if j+1 not in piv]
base=[tuple(map(vs.index,e)) for e in edges]
def abcover(a,b,cyclic):
 E=[]
 for e,(u,v) in zip(edges,base):
  g=int(signs(a,e))+2*int(signs(b,e))
  E.extend((4*u+j,4*v+((j+g)%4 if cyclic else (j^g))) for j in range(4))
 return E
def lifts(E,n):
 m=256//n;A,marks,d=setup(E,n);cond=conditions(E,n,A,marks)
 if cond is None:return
 sol=affine([(reduce(xor,r),m//2%2) for r in cond if len(r)==m],d)
 if sol is None:return
 v,bs=sol
 for mask in range(1<<len(bs)):
  t=reduce(xor,(x for j,x in enumerate(bs) if mask>>j&1),v)
  if all(len(r)-m//2<=sum((x&t).bit_count()%2 for x in r)<=m//2 for r in cond):
   yield lift(E,n,[(x&t).bit_count()%2 for x in marks])
def bipartite(E,n):
 A,_,_=setup(E,n);part={0:0};q=[0]
 for a in q:
  for b,_ in A[a]:
   if b not in part:q.append(b);part[b]=part[a]^1
   elif part[a]==part[b]:return False
 return True
if __name__=='__main__':
 for a in np.loadtxt('first.txt',dtype=int):
  E=lift(base,16,[int(signs(a,e)) for e in edges]);A,_,_=setup(E,32)
  marks=[v for e in edges for v in ((1<<ix[e] if e in ix else 0),(1<<ix[e] if e in ix else 0)|(int(signs(a,e))<<25))]
  cond=conditions(E,32,A,marks);assert cond is not None
  eq=[reduce(xor,r) for r in cond if len(r)==8]
  sol=affine([(x%2**25,x>>25) for x in eq],25)
  if sol is None:continue
  v,bs=sol
  with open(f'cyc_{a}.dat','w') as f:
   print(a,len(bs),v,*bs,file=f)
   for r in sorted(cond):
    print(len(r),*[sum(((x&b).bit_count()%2)<<j for j,b in enumerate(bs))|(((x&(v|(1<<25))).bit_count()%2)<<len(bs)) for x in r],file=f)
