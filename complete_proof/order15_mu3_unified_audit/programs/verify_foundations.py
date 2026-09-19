"""Independent benchmark, congruence, trace and color-partition checks."""
from pathlib import Path
import json
import networkx as nx
import sympy as sp
from exact_core import *
ROOT=Path(__file__).resolve().parents[1]
def run():
 M=json.loads((ROOT/'data/benchmark.json').read_text())['matrix'];H=[[ROOTS[x] for x in r] for r in M]
 G=gram(H);K=gram(adj(H));d=bareiss(H)
 assert d==gaussian_det(H) and norm(d)==B and bareiss(G)==(B,0)
 qr=sum(norm(G[i][j]) for i in range(15) for j in range(i));qc=sum(norm(K[i][j]) for i in range(15) for j in range(i));assert qr==qc==105
 deg=sorted(sum(G[i][j]!=ZERO for j in range(15) if i!=j) for i in range(15));assert deg==[2]*14+[14]
 cat=[]
 for n0 in range(16):
  for n1 in range(16-n0):
   n2=15-n0-n1;z=(n0-n2,n1-n2);N=norm(z);delta=(n1+2*n2)%3
   assert N%3==0 and (N//3-delta*delta)%3==0
   cat.append({'counts':[n0,n1,n2],'a':z[0],'b':z[1],'norm':N,'delta':delta})
 residues={(3*(x*y+x*(15-x-y)+y*(15-x-y)))%9 for x in range(16) for y in range(16-x)};assert residues=={0,6} and len(cat)==136
 ortho={a:F(15**a)*F(75-a,5)**(15-a) for a in range(3,13)};assert max(ortho.values())<B
 parts=[]
 for a in range(15,-1,-1):
  for b in range(min(a,15-a),-1,-1):
   c=15-a-b
   if 0<=c<=b:
    p=(a,b,c);cross=3*(a*b+a*c+b*c)
    if cross+9*sum(3<=x<=12 for x in p)<=168:parts.append(p)
 assert parts==[(15,0,0),(14,1,0),(13,2,0),(13,1,1),(12,3,0),(12,2,1),(11,4,0),(11,3,1),(11,2,2),(10,5,0)]
 graphs={}
 for e,cap in ((3,F(2)),(4,F(7,3))):
  count=0
  for g in nx.graph_atlas_g():
   if g.number_of_edges()!=e or any(v==0 for _,v in g.degree()):continue
   A=sp.Matrix(nx.to_numpy_array(g,dtype=int).tolist())
   assert (sp.eye(A.rows)*sp.Rational(cap.numerator,cap.denominator)-A).is_positive_semidefinite is True
   count+=1
  graphs[str(e)]={'atlas_supports':count,'adjacency_cap':str(cap),'additional_4K2_radius':1 if e==4 else None}
 cases=[(10,9,F(18),150)]+[(11,e,L,132) for e,L in [(9,F(20)),(18,F(21)),(27,F(23))]]+[(12,e,L,108) for e,L in [(9,F(18)),(18,F(77,4)),(27,F(21)),(36,F(22)),(45,F(241,10)),(54,F(25))]]
 bounds=[]
 for size,e,L,cross in cases:
  if size==11 or e>=45:assert (L-15)**2>F(2*e*(size-1),size)
  r=15-size;u=stationary(size,e)*(15-F(cross)/(r*L))**r;assert u<B
  bounds.append({'size':size,'internal':e,'cross_minimum':cross,'lambda_cap':str(L),'upper':rat(u)})
 tail=stationary(15,171);assert tail<B
 friendship=F(216**7)*F(38,3);assert friendship==B
 out={'benchmark':{'det_a':d[0],'det_b':d[1],'squared_det':B,'Q_row':qr,'Q_column':qc,'support_degrees':deg,'row_colors':[sum(r)%3 for r in M],'absolute_hadamard_ratio':(B/15**15)**0.5,'algorithms':['Bareiss H','rational Gaussian H','Bareiss HH*']},'entry_catalogue_count':len(cat),'Q_residues_mod9':sorted(residues),'trace_tail':{'from_Q':171,'upper':rat(tail),'upper_over_B':float(tail/B)},'color_reduction':{'initial':parts,'survivors':parts[:4],'small_graph_caps':graphs,'schur_checks':bounds,'orthogonal_class_checks':{str(a):rat(u) for a,u in ortho.items()}},'friendship_support':{'squared_determinant_max':int(friendship),'scope':'7K2 internal support with minimal cross energy; not all Q105'}}
 (ROOT/'results/gram_catalogue.json').write_text(json.dumps(cat,indent=2)+'\n');(ROOT/'results/foundations.json').write_text(json.dumps(out,indent=2)+'\n')
 print('Foundations PASS: det pair',d,'B',B,'Q',qr,'tail/B',float(tail/B));return out
if __name__=='__main__':run()
