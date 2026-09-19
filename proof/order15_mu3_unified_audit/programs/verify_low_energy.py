"""Exact reconstruction for Q=0,9,...,72.

For these strict-counterexample shells the color reduction forces both
Grams to be all-same. Both supports have isolates (the sole no-isolate
Q72 possibility P3+6K2 fails the norm test). EH=HF then rules out leaves.
Active components have minimum degree 2 and <=8 edges, and are exhausted
by the graph atlas on <=7 vertices plus C8. Edge values are considered
modulo sixth-root switching, legitimate for the spectral necessary tests.
We apply the two-sided necessary bound 15*diag(P_kernel(E)) >= the number
of support isolates in F, together with common spectrum and norm sieves.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import defaultdict,Counter
import json,time
import networkx as nx
from exact_core import *
ROOT=Path(__file__).resolve().parents[1];BASE={1:(1,0),3:(1,-1),4:(2,0)}
def weight_vectors(m,budget):
 if m==0:
  if budget==0:yield ()
  return
 for w in (1,3,4,7):
  if budget-w>=m-1:
   for tail in weight_vectors(m-1,budget-w):yield (w,)+tail

def extra_graphs_8_9():
 """All min-degree-2 connected 8-vertex,9-edge cores, and cycles C8,C9.
 The degree excess is two: figure-eight, barbell, or theta.
 """
 out=[nx.cycle_graph(8),nx.cycle_graph(9)]
 for a,b in ((3,6),(4,5)):
  g=nx.cycle_graph(a);path=[0]+list(range(a,a+b-1))+[0];g.add_edges_from(zip(path,path[1:]));out.append(g)
 for a,b,l in ((3,3,3),(3,4,2),(3,5,1),(4,4,1)):
  g=nx.disjoint_union(nx.cycle_graph(a),nx.cycle_graph(b));path=[0]+list(range(a+b,a+b+l-1))+[a];g.add_edges_from(zip(path,path[1:]));out.append(g)
 for lengths in ((1,2,6),(1,3,5),(1,4,4),(2,2,5),(2,3,4),(3,3,3)):
  g=nx.Graph();next_id=2
  for length in lengths:
   path=[0]+list(range(next_id,next_id+length-1))+[1];next_id+=length-1;g.add_edges_from(zip(path,path[1:]))
  out.append(g)
 for i,g in enumerate(out):
  assert nx.is_connected(g) and min(dict(g.degree()).values())>=2
  if i>=2:assert len(g)==8 and g.number_of_edges()==9
  assert not any(nx.is_isomorphic(g,h) for h in out[:i])
 return out

def components(max_units=8):
 cat=json.loads((ROOT/'results/gram_catalogue.json').read_text())
 bases={}
 for w in (1,3,4,7):
  allowed={(z['a']//3,z['b']//3) for z in cat if z['norm']==9*w and z['delta']==0};left=set(allowed);reps=[]
  while left:
   a=min(left);reps.append(a);orbit={mul(a,u) for u in UNITS};assert orbit<=allowed;left-=orbit
  bases[w]=reps
 graphs=[g for g in nx.graph_atlas_g() if len(g)>=3 and nx.is_connected(g) and min(dict(g.degree()).values())>=2 and g.number_of_edges()<=max_units]
 graphs += [g for g in extra_graphs_8_9() if g.number_of_edges()<=max_units]
 bins=defaultdict(dict);tested=Counter();graph_counts=Counter()
 for g in graphs:
  n=len(g);edges=sorted(tuple(sorted(e)) for e in g.edges());m=len(edges);graph_counts[(n,m)]+=1
  tree={tuple(sorted(e)) for e in nx.minimum_spanning_tree(g).edges()};chords=[i for i,e in enumerate(edges) if e not in tree]
  for budget in range(m,max_units+1):
   for weights in weight_vectors(m,budget):
    for bs in product(*(bases[w] for w in weights)):
     for phase in product(UNITS,repeat=len(chords)):
      ph=dict(zip(chords,phase));A=[[ZERO]*n for _ in range(n)]
      for i,((a,b),w) in enumerate(zip(edges,weights)):
       z=mul(bs[i],ph.get(i,ONE));A[a][b]=z;A[b][a]=conj(z)
      p=charpoly(A);tested[budget]+=1
      if p[-1]!=0:continue
      d=min(kernel_projection_diag(A))
      if d==0:continue
      key=(n,p);old=bins[budget].get(key,F(-1))
      if d>old:bins[budget][key]=d
 return bins,tested,graph_counts

def assemblies(bins,budget):
 kinds=[(cost,n,p,d) for cost in sorted(bins) for (n,p),d in bins[cost].items()]
 def visit(index,left,n,p,d):
  if left==0:
   if n<=15:yield p+(0,)*(15-n),15-n,d
   return
  for i in range(index,len(kinds)):
   cost,v,c,cd=kinds[i]
   if cost<=left and n+v<=15:yield from visit(i,left-cost,n+v,convolution(p,c),min(d,cd))
 yield from visit(0,budget,0,(1,),F(1))

def run():
 start=time.monotonic();bins,tested,gcounts=components();rows=[]
 no_isolate_det=3105*216**6;assert not is_norm(no_isolate_det)
 for k in range(9):
  initial=0;states={};sieve=Counter()
  for p,r,d in assemblies(bins,k):
   initial+=1;D=(-3)**15*peval(p,-5)
   if D<=B:sieve['at_or_below_B']+=1;continue
   if D%3**14 or not is_norm(D):sieve['arithmetic']+=1;continue
   key=(p,r);states[key]=max(states.get(key,F(-1)),d)
  before=len(states);rounds=[]
  while states:
   minima={}
   for p,r in states:minima[p]=min(r,minima.get(p,15))
   new={key:d for key,d in states.items() if 15*d>=minima[key[0]]}
   rounds.append(len(states)-len(new))
   if len(new)==len(states):break
   states=new
  row={'Q':9*k,'component_phase_assignments':tested[k],'assemblies_with_nonzero_kernel_on_each_coordinate':initial,'after_norm_sieve':before,'projection_deletions':rounds,'remaining_states':len(states),'status':'VERIFIED_CLOSED' if not states else 'UNRESOLVED','sieve_counts':dict(sieve),'residuals':[{'characteristic_polynomial':p,'isolates':r,'projection_min':str(d)} for (p,r),d in states.items()]}
  rows.append(row);print('Q',9*k,row['status'],'assemblies',initial,'norm survivors',before,'final',len(states),flush=True)
 out={'method':'exhaustive no-leaf support, exact norm and two-sided kernel projection','graph_counts':[{'vertices':n,'edges':e,'count':c} for (n,e),c in sorted(gcounts.items())],'phase_assignments':dict(tested),'Q72_no_isolate_determinant':no_isolate_det,'rows':rows,'seconds':time.monotonic()-start,'external_code_bound_used':False}
 (ROOT/'results/low_energy.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':run()
