"""Independent (14,1) no-internal-isolate component/Schur bounds.

Tree Gram determinants are exhaustively evaluated. Each connected support
with v vertices needs at least v-1 energy units. A v-1-unit component is a
unit-weight tree, whose phases disappear under diagonal unitary switching.
No assertion is made here that internally isolated supports are excluded.
"""
from pathlib import Path
from functools import lru_cache
from math import prod
import json,time
import networkx as nx
from exact_core import *
ROOT=Path(__file__).resolve().parents[1]
def matching_counts(T):
 def rec(v,parent):
  children=[rec(w,v) for w in T[v] if w!=parent]
  f=(1,)
  for _,total in children:f=convolution(f,total)
  total=list(f)
  for i,(free,_) in enumerate(children):
   term=(1,)
   for j,(cf,ct) in enumerate(children):term=convolution(term,cf if i==j else ct)
   if len(total)<len(term)+1:total.extend([0]*(len(term)+1-len(total)))
   for j,a in enumerate(term):total[j+1]+=a
  return f,tuple(total)
 return rec(next(iter(T)),-1)[1]
def det_tree(T):return sum((-9)**k*a*15**(len(T)-2*k) for k,a in enumerate(matching_counts(T)))
def tree_catalogue():
 out={}
 for n in range(2,15):
  best=-1;count=ties=0;best_tree=None
  for T in nx.nonisomorphic_trees(n):
   count+=1;d=det_tree(T)
   if d>best:best=d;ties=1;best_tree=T.copy()
   elif d==best:ties+=1
  assert ties==1 and nx.is_isomorphic(best_tree,nx.path_graph(n))
  out[n]={'count':count,'max_determinant':best,'maximizer':'path','matching_counts':matching_counts(best_tree)}
 assert out[14]['count']==3159 and out[14]['max_determinant']==16802420983158456
 return out
@lru_cache(None)
def configs(vertices,units,minimum=(2,1)):
 if vertices==0:return ((),) if units==0 else ()
 out=[]
 for v in range(2,vertices+1):
  for u in range(v-1,units+1):
   if (v,u)<minimum:continue
   for rest in configs(vertices-v,units-u,(v,u)):out.append(((v,u),)+rest)
 return tuple(out)
def component_cap(v,u):
 # Trees are bipartite: rho(A)^2 <= sum of all edge norms =9u.
 return 15+sqrt_hi(F(9*u) if u==v-1 else F(18*u*(v-1),v))
def bound(cfg,cross,trees):
 caps=[component_cap(v,u) for v,u in cfg]
 ds=[F(trees[v]['max_determinant']) if u==v-1 else stationary(v,9*u) for v,u in cfg]
 corr=sum(F(3*v)/L for (v,u),L in zip(cfg,caps))+F(cross-42)/max(caps)
 return prod(ds)*max(F(0),15-corr)
def run():
 t=time.monotonic();trees=tree_catalogue();rows=[]
 for q in range(105,169,9):
  for e in range(63,q-41,9):
   cs=configs(14,e//9);values=[(cfg,bound(cfg,q-e,trees)) for cfg in cs]
   upper=max((u for _,u in values),default=F(0));res=[{'configuration':cfg,'upper':rat(u),'over_B':float(u/B)} for cfg,u in values if u>B]
   rows.append({'Q':q,'internal':e,'cross':q-e,'configurations':len(cs),'max_upper':rat(upper),'max_over_B':float(upper/B),'status':'VERIFIED_NO_ISOLATE_BRANCH' if upper<=B else 'UNRESOLVED_BOUND','residuals':res})
 # Recompute all six numeric size-14 expressions in the recovered old Q150 script.
 old=[];env=stationary(15,153)
 for e,c,L in [(63,87,F(25)),(72,78,F(127,5)),(81,69,F(261,10)),(90,60,F(27)),(99,51,F(138,5)),(108,42,F(563,20))]:
  u=stationary(14,e)*(15-F(c)/L);assert u<env
  old.append({'internal':e,'cross':c,'lambda_cap':str(L),'upper':rat(u),'over_B':float(u/B),'below_old_Q153':True,'at_or_below_record':u<=B})
 assert sum(not r['at_or_below_record'] for r in old)==5
 out={'tree_catalogue':trees,'component_rows':rows,'Q150_old_comparator_audit':old,'scope':'only (14,1) with no internal support isolates; other color types and isolate exclusion remain separate obligations','seconds':time.monotonic()-t}
 (ROOT/'results/component_bounds.json').write_text(json.dumps(out,indent=2)+'\n')
 for q in range(105,169,9):
  rr=[r for r in rows if r['Q']==q];print('Q',q,'(14,1) no-isolate:',sum(r['status']=='VERIFIED_NO_ISOLATE_BRANCH' for r in rr),'/',len(rr),'allocations; max/B',max(r['max_over_B'] for r in rr))
 print('Old Q150 comparator: five of six upper bounds exceed B.');return out
if __name__=='__main__':run()
