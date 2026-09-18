"""Audit HH* versus H*H orientation in the historical color-orbit argument."""
from pathlib import Path
import sys,json,random
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'independent_audit/programs'))
from exact_core import gram,bareiss,norm,conj,ROOTS
root=Path(__file__).resolve().parents[1]
rng=random.Random(150055)
while True:
 rows=[]
 for _ in range(12):
  r=[0]*5+[1]*5+[2]*5;rng.shuffle(r);rows.append(r)
 last=[-sum(r[j] for r in rows)%3 for j in range(15)]
 if [last.count(i) for i in range(3)]==[5,5,5]:break
E=[[0]*15]+rows+[last]+[[1]+[0]*14]
H=[[ROOTS[e] for e in row] for row in E]
G=gram(H)
# transpose-conjugate H for the independent column Gram.
Ht=[[conj(H[i][j]) for i in range(15)] for j in range(15)]
K=gram(Ht)
assert all(G[0][i]==(0,0) for i in range(1,14))
rc=[sum(r)%3 for r in E];cc=[sum(r[j] for r in E)%3 for j in range(15)]
assert rc==[0]*14+[1] and cc==[1]+[0]*14
D=bareiss(H);assert norm(D)>0
# Row-majority to singleton entry and column-majority to singleton entry.
g=G[0][14];k=K[1][0]
roots=ROOTS
from exact_core import mul,neg
delta=(1,-1)
assert {mul(delta,z) for z in ROOTS}=={mul(neg(conj(delta)),z) for z in ROOTS}
roworbit={mul(g,z) for z in roots};colorbit={mul(k,z) for z in roots}
assert roworbit!=colorbit
row_residue=tuple(x%3 for x in g);column_residue=tuple(x%3 for x in k)
assert row_residue!=column_residue
assert row_residue==tuple(x%3 for x in conj(k))
out={'seed':150055,'exponents':E,'row_colors':rc,'column_colors':cc,'internally_isolated_row':0,'determinant':D,'squared_absolute_determinant':norm(D),'Q_row':sum(norm(G[i][j]) for i in range(15) for j in range(i)),'Q_column':sum(norm(K[i][j]) for i in range(15) for j in range(i)),'row_majority_to_singleton_gram_entry':g,'column_majority_to_singleton_gram_entry':k,'row_cross_residue_mod_3':row_residue,'column_cross_residue_mod_3':column_residue,'same_oriented_cross_residue_claim_false':True,'corrected_norm3_orbits_coincide':True,'scope':'Refutes unrestricted same-oriented cross-orbit/isolate inference; not a counterexample to the low-energy maximal-determinant theorem.'}
(root/'provenance/remote_orbit_orientation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='exponents'},indent=2))
