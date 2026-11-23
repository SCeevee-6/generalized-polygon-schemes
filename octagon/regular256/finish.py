from towers import *
from spectral import polynomial,inertia
from collections import Counter
stats=Counter()
for cyclic,a,b in np.loadtxt('four.txt',dtype=int):
 stats[int(cyclic),'four']+=1
 if cyclic:
  R=np.zeros((1,16,16),np.int64);I=R.copy()
  for e,(u,v) in zip(edges,base):
   g=int(signs(a,e))+2*int(signs(b,e));R[0,u,v]=R[0,v,u]=(1,0,-1,0)[g];I[0,u,v]=(0,1,0,-1)[g];I[0,v,u]=-I[0,u,v]
  co=polynomial(R,I)[0];l,_=inertia(co,-2);r,z=inertia(co,-1)
  if l>r+z:stats[1,'gap']+=1;continue
 E=abcover(a,b,cyclic)
 for F in lifts(E,64):
  stats[int(cyclic),'eight']+=1
  for H in lifts(F,128):
   stats[int(cyclic),'sixteen']+=1
   if bipartite(H,256):stats[int(cyclic),'bipartite']+=1
   else:print('nonbipartite survivor',cyclic,a,b,flush=True)
print(sorted(stats.items()))
