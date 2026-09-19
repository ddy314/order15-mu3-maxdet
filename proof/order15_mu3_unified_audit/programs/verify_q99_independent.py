"""Independently derive all local Q99 intertwiners. Not a global Q99 proof."""
from pathlib import Path
from itertools import product
import importlib.util,json,sys
import sympy as sp
from exact_core import *
ROOT=Path(__file__).resolve().parents[1]
def run():
 p=ROOT/'original/programs/verify_q99_final_orbits_exact.py';spec=importlib.util.spec_from_file_location('original_q99',p);mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
 original_output=mod.certificate();out={}
 for name,L0,R0,fn,test,expected in [('AB',mod.C4,mod.T6,mod.ab_local_blocks,lambda b:mod.is_c4_residual_shape(mod.subtract(mod.C4,mod.gram(b)),9),144),('AA',mod.C4,mod.C4,mod.aa_local_blocks,lambda b:mod.is_c4_residual_shape(mod.subtract(mod.C4,mod.gram(b)),11),225),('BB',mod.T6,mod.T6,mod.bb_local_blocks,lambda b:mod.is_t6_residual_shape(mod.subtract(mod.T6,mod.gram(b)),9),729)]:
  r=len(L0);c=len(R0);N=r*c;T=sp.zeros(N,N)
  for i in range(r):
   for j in range(c):
    for h in range(r):T[i*c+j,h*c+j]+=L0[i][h].a
    for h in range(c):T[i*c+j,i*c+h]-=R0[h][j].a
  C,piv=T.rref();free=[j for j in range(N) if j not in piv];forms=[(col,[(t,F(-C[row,j])) for t,j in enumerate(free) if C[row,j]]) for row,col in enumerate(piv)]
  states=set()
  for vals in product(ROOTS,repeat=len(free)):
   entries=[ZERO]*N
   for j,z in zip(free,vals):entries[j]=z
   good=True
   for col,form in forms:
    v=esum(scale(vals[t],a) for t,a in form)
    if v not in ROOTS:good=False;break
    entries[col]=v
   if good:states.add(tuple(entries))
  original={tuple((v.a,v.b) for row in b for v in row) for b in fn()};assert states==original and len(states)==expected
  survivors=sum(test([[mod.E(*v) for v in flat[i*c:(i+1)*c]] for i in range(r)]) for flat in states);assert survivors==0
  out[name]={'operator_size':N,'rank':len(piv),'nullity':len(free),'free_coordinates':free,'assignments':3**len(free),'mu3_solutions':len(states),'residual_survivors':survivors,'parameterization_complete':True}
 # Cross-check full A and B spectra, energy and the support of the 15-eigenspace.
 x=sp.Symbol('x');expected_cp=tuple(map(int,sp.Poly((x-21)**2*(x-18)**4*(x-15)**2*(x-12)**6*(x-9),x).all_coeffs()))
 for name,sizes in [('A',[2]*4+[3,4]),('B',[2]*3+[3,6])]:
  G=[[(15*int(i==j),0) for j in range(15)] for i in range(15)];o=0
  for n in sizes:
   edges={2:[(0,1,-3)],3:[(0,1,3),(1,2,3),(0,2,3)],4:[(0,1,3),(1,2,3),(2,3,3),(0,3,3)],6:[(0,1,3),(0,2,-3),(0,3,-3),(1,4,-3),(1,5,-3)]}[n]
   for i,j,v in edges:G[o+i][o+j]=G[o+j][o+i]=(v,0)
   o+=n
  assert charpoly(G)==expected_cp and sum(norm(G[i][j]) for i in range(15) for j in range(i))==99
  P=kernel_projection_diag([[sub(G[i][j],(15*int(i==j),0)) for j in range(15)] for i in range(15)])
  assert [i for i,v in enumerate(P) if v]==[11,12,13,14]
 out={'status':'VERIFIED_LOCAL_ONLY','pairs':out,'prerequisites_not_verified':['Completeness of the two-orbit catalogue for the all-same case','All mixed-color Q99 cases, including size13 internal energies9 and18']}
 (ROOT/'results/q99_independent.json').write_text(json.dumps(out,indent=2)+'\n');print('Q99 local rederivation PASS: 144,225,729 solutions; all residual survivors zero.');return out
if __name__=='__main__':run()
