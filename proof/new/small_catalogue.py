"""Exhaustive connected Hermitian Eisenstein determinant bounds through 7 vertices.
Tree edges are fixed modulo sixth-root switching; every chord phase is enumerated.
All integer arithmetic: matrices are 5I+E (the actual Gram is three times this).
"""
from bounds import *
import numpy as np, networkx as nx, time
from numba import njit
from itertools import product,permutations

@njit(cache=True)
def det_hermitian(a,b):
 a=a.copy();b=b.copy();n=len(a);prev=1
 for k in range(n-1):
  pivot=a[k,k]
  if b[k,k]!=0 or pivot<=0:raise ValueError('nonpositive/nonreal pivot')
  for i in range(k+1,n):
   for j in range(k+1,n):
    x=a[i,k]*a[k,j]-b[i,k]*b[k,j]
    y=a[i,k]*b[k,j]+b[i,k]*a[k,j]-b[i,k]*b[k,j]
    xx=pivot*a[i,j]-x;yy=pivot*b[i,j]-y
    if xx%prev or yy%prev:raise ValueError('noninteger division')
    a[i,j]=xx//prev;b[i,j]=yy//prev
  prev=pivot
 if b[-1,-1]!=0:raise ValueError('nonreal determinant')
 return a[-1,-1]

@njit(cache=True)
def phase_maximum(n,edges,base,ischord,units):
 chord_count=0
 for v in ischord:chord_count+=v
 best=-1;count=6**chord_count
 for serial in range(count):
  a=np.zeros((n,n),dtype=np.int64);b=np.zeros((n,n),dtype=np.int64)
  for i in range(n):a[i,i]=5
  s=serial
  for j in range(len(edges)):
   x=base[j,0];y=base[j,1]
   if ischord[j]:
    k=s%6;s//=6;z=units[k,0];w=units[k,1]
    x,y=x*z-y*w,x*w+y*z-y*w
   v=edges[j,0];w=edges[j,1]
   a[v,w]=x;b[v,w]=y;a[w,v]=x-y;b[w,v]=-y
  d=det_hermitian(a,b)
  if d>best:best=d
 return best,count

def unique_permutations(values):
 """All multiset permutations, without generating repeated factorial copies."""
 from collections import Counter
 counts=Counter(values); keys=tuple(sorted(counts)); n=len(values); prefix=[]
 def visit():
  if len(prefix)==n:
   yield tuple(prefix); return
  for value in keys:
   if counts[value]:
    counts[value]-=1; prefix.append(value)
    yield from visit()
    prefix.pop(); counts[value]+=1
 yield from visit()

EUNITS=np.array(UNITS,dtype=np.int64)
@lru_cache(None)
def norm_bases(w):
 available={(a,b) for a in range(-10,11) for b in range(-10,11) if norm((a,b))==w}
 out=[]
 while available:
  z=min(available);out.append(z);available-={mul(z,x) for x in UNITS}
 return tuple(out)

@lru_cache(None)
def weight_vectors(m,total,minimum=1):
 if m==0:return ((),) if total==0 else ()
 return tuple((w,)+rest for w in UNITS_ALLOWED if minimum<=w<=total for rest in weight_vectors(m-1,total-w,w))

def run(max_units=10):
 start=time.monotonic();records={};tot=0;states=0
 # Warm compilation.
 phase_maximum(2,np.array([[0,1]],dtype=np.int64),np.array([[1,0]],dtype=np.int64),np.array([0],dtype=np.int64),EUNITS)
 for n in range(2,8):
  graphs=[g for g in nx.graph_atlas_g() if len(g)==n and nx.is_connected(g) and g.number_of_edges()<=max_units]
  for u in range(n-1,max_units+1):
   maximum=0;cnt=0;ng=0;arg=None;triangle_max=F(0)
   for g in graphs:
    edges=sorted(tuple(sorted(e)) for e in g.edges());m=len(edges)
    if m>u:continue
    tree={tuple(sorted(e)) for e in nx.minimum_spanning_tree(g).edges()}
    ischord=np.array([int(e not in tree) for e in edges],dtype=np.int64)
    edgearr=np.array(edges,dtype=np.int64)
    for wp in weight_vectors(m,u):
     for ws in unique_permutations(wp):
      ng+=1
      for bs in product(*(norm_bases(w) for w in ws)):
       d,c=phase_maximum(n,edgearr,np.array(bs,dtype=np.int64),ischord,EUNITS)
       cnt+=c
       if d>maximum:maximum=d;arg={'edges':edges,'weights':ws,'bases':bs}
   records[f'{n},{u}']={'max_determinant':int(maximum)*3**n,'scaled_max':int(maximum),'phase_states':cnt,'weighted_supports':ng,'max_support':arg}
   tot+=cnt;states+=ng
  print('vertices',n,'total phase states',tot,'seconds',round(time.monotonic()-start,2),flush=True)
 out={'method':'integer Bareiss; graph atlas; all norm classes; sixth-root chord phases','max_units':max_units,'phase_states':tot,'weighted_supports':states,'records':records,'seconds':time.monotonic()-start}
 Path(__file__).with_name('small_catalogue.json').write_text(json.dumps(out,indent=2)+'\n')
 return out
if __name__=='__main__':run()
