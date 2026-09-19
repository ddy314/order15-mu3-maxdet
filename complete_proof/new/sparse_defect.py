"""Finite arithmetic part of the 12+1, one-norm9-defect obstruction.
The full rank-three frame and spectral-radius-two classification proof is
included in PROOF_zh.md; these are arithmetic checks, not a substituted axiom.
"""
from bounds import *
import networkx as nx
from math import lcm,gcd
from functools import reduce

def affine_candidates():
 for n in range(3,10):yield 'cycle_'+str(n),nx.cycle_graph(n),[1]*n
 G=nx.star_graph(4);yield 'four_arm_star',G,[2,1,1,1,1]
 for n in range(6,10):
  core=n-4;G=nx.path_graph(core)
  G.add_edges_from([(0,core),(0,core+1),(core-1,core+2),(core-1,core+3)])
  yield 'two_forks_'+str(n),G,[2]*core+[1]*4
 for arms in ((2,2,2),(1,3,3),(1,2,5)):
  c=lcm(*(a+1 for a in arms));weights=[c];G=nx.Graph();G.add_node(0)
  for length in arms:
   prev=0
   for k in range(1,length+1):
    i=len(weights);weights.append(c*(length+1-k)//(length+1));G.add_edge(prev,i);prev=i
  gg=reduce(gcd,weights);weights=[w//gg for w in weights]
  yield 'three_arms_'+','.join(map(str,arms)),G,weights

def certificate():
 candidates=[];allowed=[]
 for name,G,m in affine_candidates():
  assert len(G)==len(m)<=9 and all(2*m[i]==sum(m[j]for j in G.neighbors(i))for i in G)
  assert reduce(gcd,m)==1
  base=sum(m)
  for scale0 in range(1,15//base+1):
   row={'family':name,'vertices':len(G),'primitive_weights':m,'primitive_mass':base,'scale':scale0,'mass':base*scale0,'passes_integral_tight_frame':base%3==0}
   candidates.append(row)
   if row['passes_integral_tight_frame']:allowed.append(row)
 masses=sorted(set(x['mass']for x in allowed));assert masses==[3,6,9,12,15]
 final=[x for x in allowed if x['mass']%5==0]
 assert len(final)==1 and final[0]['family']=='cycle_3' and final[0]['scale']==5
 signs=[]
 for s in (1,-1):
  rs=[sub(scale(z,5),(s,0))for z in ROOTS]
  assert all(norm(z)>0 for z in rs)
  for i in range(3):
   for j in range(i+1,3):assert div(rs[i],rs[j]) not in UNITS
  ns=[norm(z)for z in rs];assert ns[1]==ns[2]
  n0=F(15*(24-ns[1]),ns[0]-ns[1])
  assert n0.denominator==1 and n0 in (3,7) and n0%5!=0
  signs.append({'inner_product_sign':s,'scaled_residual_norms':ns,'required_symbol_one_count':int(n0),'forced_count_multiple':5,'contradiction':True})
 out={'excluded_defect':'13 mu3 rows of length15 with a single norm9 off-diagonal Gram edge','candidate_families':candidates,'admissible_component_masses':masses,'only_mass_divisible_by5':final,'sign_checks':signs,'closed':True}
 Path(__file__).with_name('sparse_defect.json').write_text(json.dumps(out,indent=2)+'\n');print('One norm9 defect excluded for both signs; allowed component masses',masses,flush=True)
 return out
if __name__=='__main__':certificate()
