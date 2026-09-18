"""Deterministic cross-checks against a separate real-integer representation."""
from pathlib import Path
import json,random
import sympy as sp
import networkx as nx
from exact_core import *
from verify_component_bounds import det_tree
ROOT=Path(__file__).resolve().parents[1]
def run():
 rng=random.Random(20260917);count=0
 for n in range(1,7):
  for _ in range(10):
   A=[[(rng.randrange(-3,4),rng.randrange(-3,4)) for j in range(n)] for i in range(n)]
   d=bareiss(A);assert d==gaussian_det(A)
   R=sp.zeros(2*n)
   for i,row in enumerate(A):
    for j,(a,b) in enumerate(row):R[i,j]=a;R[i,n+j]=-b;R[n+i,j]=b;R[n+i,n+j]=a-b
   assert int(R.det(method='domain-ge'))==norm(d);count+=1
   if d!=ZERO:assert mm(A,inverse(A))==[[(int(i==j),0) for j in range(n)] for i in range(n)]
 treechecks=0
 for n in range(2,10):
  for T in nx.nonisomorphic_trees(n):
   A=[[(15 if i==j else 3*int(T.has_edge(i,j)),0) for j in range(n)] for i in range(n)]
   assert bareiss(A)==(det_tree(T),0);treechecks+=1
 for n in range(1,7):
  A=[[ZERO]*n for _ in range(n)]
  for i in range(n):
   A[i][i]=(rng.randrange(-2,3),0)
   for j in range(i):A[i][j]=(rng.randrange(-2,3),rng.randrange(-2,3));A[j][i]=conj(A[i][j])
  p=charpoly(A)
  for t in (-3,0,2):
   Z=[[sub((t*int(i==j),0),A[i][j]) for j in range(n)] for i in range(n)]
   assert bareiss(Z)==(peval(p,t),0)
 out={'seed':20260917,'random_matrices_tested':count,'tests':['Bareiss vs rational elimination','Norm(det A) vs determinant of 2n real-integer representation','Exact inverse identity','Newton characteristic polynomial vs direct determinants'],'tree_DP_vs_Bareiss_checks':treechecks}
 (ROOT/'results/arithmetic_tests.json').write_text(json.dumps(out,indent=2)+'\n');print('Arithmetic cross-checks PASS:',count,'matrices;',treechecks,'trees');return out
if __name__=='__main__':run()
