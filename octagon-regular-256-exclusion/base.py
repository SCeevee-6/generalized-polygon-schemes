from itertools import permutations, combinations
from pathlib import Path
import numpy as np
vs=[v for v in range(32) if v.bit_count()%2==0]
edges=[(a,b) for a,b in combinations(vs,2) if (a^b).bit_count()==4]
adj={v:[w for w in vs if (v^w).bit_count()==4] for v in vs}
order=[0];parent={0:0}
for a in order:
 for b in adj[a]:
  if b not in parent:parent[b]=a;order.append(b)
tree={tuple(sorted((b,a))) for b,a in parent.items() if b}
ix={e:i for i,e in enumerate(e for e in edges if e not in tree)}
def act(m,p):
 marks={tuple(sorted((p[a],p[b]))) for (a,b),i in ix.items() if m>>i&1}
 t={0:0}
 for v in order[1:]:t[v]=t[parent[v]]^(tuple(sorted((v,parent[v]))) in marks)
 return sum((t[a]^t[b]^((a,b) in marks))<<i for (a,b),i in ix.items())
def signs(m,e):return (m>>ix[e])&1 if e in ix else 0
if __name__=='__main__':
 ps=[{v:v^3 for v in vs}]+[{v:sum(1<<s[i] for i in range(5) if v>>i&1) for v in vs} for s in [(1,2,3,4,0),(1,0,2,3,4)]]
 np.array([[act(1<<i,p) for i in range(25)] for p in ps],np.int32).tofile('act.bin')
