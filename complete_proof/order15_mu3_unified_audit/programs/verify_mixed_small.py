"""Reconstruction of mixed Q=6,15,...,78 using a cited code bound.

External input: B_3(15,10)=12, Todorov--Bogdanova (2020), Proposition14
and Table1. Orthogonal mu3 rows give ternary equidistant words at distance10.
The published classification is an external theorem, not rerun here.

At Q60,69,78 both row/column color types are (14,1), and there are >=6
internal isolates. The opposite hub has zero coordinate in the 15-eigenspace.
Intertwining thus forces the cross entries at these isolates into one mu3
orbit. The energy budget forces their common norm to be3. Normalizing all
cross entries modulo mu3 leaves finite active-support and cross-state lists.
For k internal isolates and an opposite 15-eigenspace supported on <=s
coordinates, the residual Gram PSD inequality requires k*(15-s)<=15.
"""
from pathlib import Path
from itertools import product,combinations
from functools import lru_cache
from collections import Counter
import json,time
import networkx as nx
from exact_core import *
ROOT=Path(__file__).resolve().parents[1]
def alpha(g):
 n=len(g)
 for k in range(n,-1,-1):
  for I in combinations(range(n),k):
   if not any(g.has_edge(i,j) for i,j in combinations(I,2)):return k

def allocations(options,length,budget):
 if length==0:
  if budget==0:yield ()
  return
 for cost,values in options.items():
  if cost<=budget:
   for v in values:
    for tail in allocations(options,length-1,budget-cost):yield (v,)+tail

def canonical_orbits(values):
 left=set(values);out=[]
 while left:
  a=min(left);out.append(a);left-={mul(a,z) for z in ROOTS}
 return tuple(out)

def spectrum(C,p,u,k,iso_norm):
 n=len(C);I=[[(int(i==j),0) for j in range(n)] for i in range(n)];M=I;r=[]
 for j in range(n):
  v=esum(mul(conj(u[a]),esum(mul(M[a][b],u[b]) for b in range(n))) for a in range(n));assert v[1]==0;r.append(v[0])
  if j<n-1:
   M=mm(C,M)
   for a in range(n):M[a][a]=add(M[a][a],(p[j+1],0))
 out=list(convolution((1,0,-k*iso_norm),p))
 for j,v in enumerate(r):out[j+2]-=v
 return tuple(out)+(0,)*(k-1)

def run():
 start=time.monotonic()
 # A same-color class of size a contains an independent set of size at least
 # a-m when its support has m edges. The external code bound forces m>=a-12.
 partitions=[(15,0,0),(14,1,0),(13,2,0),(13,1,1)]
 minimum={p:3*(p[0]*p[1]+p[0]*p[2]+p[1]*p[2])+9*sum(max(0,a-12) for a in p) for p in partitions}
 for q in (6,15,24,33,42,51):
  assert not any(q>=minimum[p] and q%9==(3*(p[0]*p[1]+p[0]*p[2]+p[1]*p[2]))%9 for p in partitions)
 for q in (60,69,78):
  assert [p for p in partitions if q>=minimum[p] and q%9==(3*(p[0]*p[1]+p[0]*p[2]+p[1]*p[2]))%9]==[(14,1,0)]
 cat=json.loads((ROOT/'results/gram_catalogue.json').read_text());inside={w:tuple({(x['a'],x['b']) for x in cat if x['norm']==9*w and x['delta']==0}) for w in (1,3,4)}
 cross={b:canonical_orbits({(x['a'],x['b']) for x in cat if x['delta']==1 and x['norm']==3+9*b}) for b in range(3)}
 assert len(cross[0])==1
 graphs=[g for g in nx.graph_atlas_g() if 1<=g.number_of_edges()<=4 and not any(d==0 for _,d in g.degree())]
 graphs.append(nx.disjoint_union_all([nx.path_graph(2) for _ in range(4)]))
 graphs=[g for g in graphs if len(g)-alpha(g)>=2]
 rows=[]
 for q in (60,69,78):
  states={};tested=norm_pass=0
  for units in range(2,(q-42)//9+1):
   excess=(q-42)//9-units
   for g in graphs:
    n=len(g);k=14-n;edges=sorted(g.edges())
    if len(edges)>units:continue
    # Edge choices combine exact norm cost and all permitted Eisenstein values.
    for zs in allocations(inside,len(edges),units):
     C=[[ZERO]*n for _ in range(n)]
     for (i,j),v in zip(edges,zs):C[i][j]=v;C[j][i]=conj(v)
     A=[[add(C[i][j],(15*int(i==j),0)) for j in range(n)] for i in range(n)]
     d=bareiss(A);assert d[1]==0 and d[0]>0;inv=inverse(A);cp=None;t=None
     for u in allocations(cross,n,excess):
      tested+=1
      loss=esum(mul(conj(u[i]),esum(mul(inv[i][j],u[j]) for j in range(n))) for i in range(n));assert loss[1]==0
      D=F(d[0])*15**k*(15-F(k,5)-loss[0]);assert D.denominator==1;D=int(D)
      if D<=B or D%3**14 or not is_norm(D):continue
      norm_pass+=1
      if cp is None:cp=charpoly(C);t=sum(v>0 for v in kernel_projection_diag(C))
      p=spectrum(C,cp,u,k,3);assert (-1)**15*peval(p,-15)==D
      # E=G-15I characteristic polynomial; a smaller coordinate support only
      # strengthens the necessary inequality, so this upper estimate is safe.
      states[(p,k,k+t)]=D
  before=len(states);deletions=[]
  while states:
   largest={}
   for p,k,s in states:largest[p]=max(s,largest.get(p,0))
   new={key:D for key,D in states.items() if key[1]*(15-largest[key[0]])<=15}
   deletions.append(len(states)-len(new))
   if len(new)==len(states):break
   states=new
  row={'Q':q,'tested_phase_cross_states':tested,'norm_pass_before_deduplication':norm_pass,'states_before_PSD':before,'PSD_deletions':deletions,'residual_states':len(states),'status':'VERIFIED_WITH_CITED_CODE_BOUND' if not states else 'UNRESOLVED','residuals':[{'E_characteristic_polynomial':p,'internal_isolates':k,'eigenspace_support_upper':s,'squared_det':D} for (p,k,s),D in states.items()]}
  rows.append(row);print('mixed Q',q,row['status'],'states',tested,'norm',before,'residual',len(states),flush=True)
 out={'external_theorem':{'statement':'B_3(15,10)=12','authors':'Todor Todorov and Galina Bogdanova','year':2020,'doi':'10.28919/jmcs/4964','locations':'Proposition14, p2717; Table1, p2720','independent_classification_replayed':False},'elementary_shells':[6,15,24,33,42,51],'enumerated_shells':rows,'seconds':time.monotonic()-start}
 (ROOT/'results/mixed_small.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':run()
