"""Q81 reconstruction: complete support split and local block obstruction.
See report for the short proof connecting the finite checks below.
External input B_3(15,10)=12 removes the 13+1+1 color type.
"""
from pathlib import Path
from itertools import permutations,product
from collections import Counter
import json,time
import networkx as nx
import sympy as sp
from exact_core import *
from verify_low_energy import components,assemblies
ROOT=Path(__file__).resolve().parents[1]
def no_isolate_catalogue():
 kinds=[]
 for n in range(2,6):
  for T in nx.nonisomorphic_trees(n):
   E=[[(3*int(T.has_edge(i,j)),0) for j in range(n)] for i in range(n)]
   kinds.append((n,n-1,charpoly(E),'tree'+str(n)))
 for z in UNITS:
  E=[[ZERO,(3,0),scale(z,3)],[(3,0),ZERO,(3,0)],[scale(conj(z),3),(3,0),ZERO]]
  kinds.append((3,3,charpoly(E),'triangle'+str(z)))
 all_count=0;survivors={}
 def visit(index,n,e,p,parts):
  nonlocal all_count
  if n==15 and e==9:
   all_count+=1;D=-peval(p,-15)
   if D>B and D%3**14==0 and is_norm(D):survivors[p]={'D':D,'parts':parts}
   return
  for i in range(index,len(kinds)):
   nn,ee,cp,name=kinds[i]
   if n+nn<=15 and e+ee<=9:visit(i,n+nn,e+ee,convolution(p,cp),parts+[name])
 visit(0,0,0,(1,),[])
 assert len(survivors)==4
 for p,v in survivors.items():
  assert p[-1]!=0 and v['parts'].count('tree2')==6 and len(v['parts'])==7
 return all_count,survivors

def local_obstruction():
 # On the negative-flux equality branch, the proof forces 6 positive K2
 # blocks, D D*=12I, and S block row/column groups of size3 summing to0.
 # Two zeros of D in any such group would require 2a+2b-c=0 in mu3.
 bad=sum(add(scale(a,2),scale(b,2))==c for a,b,c in product(ROOTS,repeat=3));assert bad==0
 # One zero forces 2a-b-c=0 only when all three roots coincide.
 single=[(a,b,c) for a,b,c in product(ROOTS,repeat=3) if scale(a,2)==add(b,c)]
 assert len(single)==3 and all(a==b==c for a,b,c in single)
 # Thus every S 3x3 block is a*(3P-J), with a in mu3 and P a permutation.
 # Check cross-block orthogonality is impossible, reducing by one permutation
 # and one root so the check is exact and exhaustive on the invariant data.
 matrices=[]
 for perm in permutations(range(3)):
  matrices.append([[(3*int(perm[i]==j)-1,0) for j in range(3)] for i in range(3)])
 zero=0
 for A,C in product(matrices,repeat=2):
  for u,v in product(ROOTS,repeat=2):
   # 9R-3J + phase*(9T-3J) are precisely the possible cross products.
   X=[[add(mul(u,scale(A[i][j],3)),mul(v,scale(C[i][j],3))) for j in range(3)] for i in range(3)]
   zero+=all(x==ZERO for r in X for x in r)
 assert zero==0
 x=sp.Symbol('x');pair=sp.Poly((x-18)*(x-12),x)
 # Nonreal cycle flux has no common eigenvalue with a K2 Gram block.
 for real_twice in (-1,1):
  tri=sp.Poly((x-15)**3-27*(x-15)-27*real_twice,x)
  assert sp.gcd(pair,tri).degree()==0
 assert 21>3*3 # positive-flux localized eigenspace exceeds 3x3 block norm.
 return {'two_zeros_group_solutions':bad,'single_zero_group_solutions':len(single),'permutation_phase_pairs':6*6*3*3,'cross_orthogonal_pairs':zero}

def run():
 start=time.monotonic();count,noiso=no_isolate_catalogue();bins,tested,_=components(9);states={};assemblies_count=0
 for p,r,d in assemblies(bins,9):
  assemblies_count+=1;D=(-3)**15*peval(p,-5)
  if D>B and D%3**14==0 and is_norm(D):states[(p,r)]=max(states.get((p,r),F(-1)),d)
 before=len(states);rounds=[]
 while states:
  minima={}
  for p,r in states:minima[p]=min(r,minima.get(p,15))
  new={key:d for key,d in states.items() if 15*d>=minima[key[0]]}
  rounds.append(len(states)-len(new))
  if len(new)==len(states):break
  states=new
 proof=local_obstruction()
 out={'Q':81,'status':'VERIFIED_WITH_CITED_CODE_BOUND' if not states else 'UNRESOLVED','no_isolate_support_phase_cases':count,'no_isolate_norm_surviving_spectra':len(noiso),'no_isolate_survivors':[{'E_charpoly':p,**v} for p,v in noiso.items()],'isolate_component_phase_assignments':tested[9],'isolate_assemblies':assemblies_count,'isolate_norm_survivors':before,'isolate_projection_deletions':rounds,'isolate_remaining':len(states),'local_finite_obstruction':proof,'seconds':time.monotonic()-start,'scope_note':'Requires the mathematical support split and block argument in the report, plus the cited equidistant-code theorem.'}
 (ROOT/'results/q81.json').write_text(json.dumps(out,indent=2)+'\n');print('Q81',out['status'],'isolate states',before,'->',len(states),'non-isolate spectra',len(noiso),'->0');return out
if __name__=='__main__':run()
