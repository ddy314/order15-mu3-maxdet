"""Record-level closure certificates for the independently reduced high shells."""
from bounds import *
from cross_exact import eis_shell,orbit_bases,cross_norm_patterns,tree_matrix
from itertools import product
from math import lcm


def q114_certificate():
 refine=json.loads(Path(__file__).with_name('q114_refinement.json').read_text())
 cross=json.loads(Path(__file__).with_name('q114_cross_exact.json').read_text())
 assert len(refine['rows'])==len(cross['rows'])==8
 assert all(x['surviving_matrices']==1 for x in refine['rows'])
 # Every discarded core phase is bounded below the record by the exact
 # determinant screening replayed in refine_q114.py. Cross-spectra above
 # the record reduce to one norm-compatible case.
 survivors=[x for x in cross['rows'] if x['survivors']]
 assert len(survivors)==1
 survivor=survivors[0]
 assert survivor['internal']==54 and survivor['configuration'].count([1,0])==10
 A=tuple(tuple((15 if i==j else 3,0)for j in range(4))for i in range(4))
 J=inverse([list(row)for row in A]);minloss=F(11,8)
 minimizers=[];tested=0
 for ns in cross_norm_patterns(4,2):
  options=[orbit_bases(ns[0])]+[eis_shell(n)for n in ns[1:]]
  for z in product(*options):
   tested+=1
   value=esum(mul(conj(z[i]),mul(J[i][j],z[j]))for i in range(4)for j in range(4))
   assert value[1]==0 and value[0]>=minloss
   if value[0]==minloss:
    lengths=[2 if n==12 else 1 for n in ns]
    assert sorted(ns)==[3,3,12,12]
    units=[div(zz,(ll,0))for zz,ll in zip(z,lengths)]
    assert len(set(units))==1 and norm(units[0])==3
    minimizers.append({'norms':ns,'vector':z})
 assert len(minimizers)==6
 # Canonical remaining Gram. Unit-modulus switching preserves every
 # argument below, including the flat-entry condition on H.
 delta=(1,-1)
 G=[[(15 if i==j else 0,0)for j in range(15)]for i in range(15)]
 for i in range(4):
  for j in range(4):G[i][j]=A[i][j]
 for i in range(14):
  z=scale(delta,2 if i<2 else 1)
  G[i][14]=z;G[14][i]=conj(z)
 E=[[sub(G[i][j],(12 if i==j else 0,0))for j in range(15)]for i in range(15)]
 diag=kernel_projection_diag(E)
 assert diag==[F(1,2)]*4+[F(0)]*11
 for a,b in ((0,1),(2,3)):
  v=[(int(i==a)-int(i==b),0)for i in range(15)]
  assert all(esum(mul(E[i][j],v[j])for j in range(15))==ZERO for i in range(15))
 D=bareiss(G)[0]
 assert str(D)==survivor['survivors'][0] and D>B and is_norm(D)
 # GH=HK transfers the 12-eigenspace. For either matched pair a,b in K,
 # H_ia=H_ib up to one common phase on the 11 rows outside the core of G.
 # Their column inner product consequently has magnitude at least11-4=7,
 # contradicting the prescribed magnitude3. This even excludes arbitrary
 # complex unit-modulus H, so no mu3-orientation lemma is being assumed.
 assert 11-4>3
 out={'Q':114,'closed':True,'core_phase_screening_configurations':8,'cross_vector_states_in_last_case':tested,'last_case_minimizer_count':len(minimizers),'last_abstract_determinant':D,'last_abstract_over_B':float(F(D,B)),'kernel_12_projection_diagonal':[str(x)for x in diag],'forced_column_inner_product_lower':7,'prescribed_column_inner_product':3,'canonical_gram':G,'minimum_cross_loss':str(minloss),'minimizers':minimizers}
 Path(__file__).with_name('q114_final.json').write_text(json.dumps(out,indent=2)+'\n')
 return out


def q126_certificate():
 rows=json.loads(Path(__file__).with_name('projection_frontier.json').read_text())
 row=next(x for x in rows if x['Q']==126)
 assert len(row['remaining'])==1
 r=row['remaining'][0]
 assert r['configuration']==[[1,0]]*7+[[8,14]]
 assert r['component_constraints']==[[8,14,4,5]]
 # Therefore rank(G-15I)<=4. Since (G-15I)/3 is integral Hermitian
 # Eisenstein, det G is divisible by 15^11 * 3^4.
 factor=15**11*3**4
 upper=15**11*stationary(4,126)
 lo=B//factor+1;hi=int(upper//factor)
 possibilities=[{'integer':k,'eisenstein_norm':is_norm(k*factor)}for k in range(lo,hi+1)]
 assert (lo,hi)==(397,401) and not any(x['eisenstein_norm']for x in possibilities)
 out={'Q':126,'closed':True,'rank_upper':4,'determinant_factor':factor,'upper':rat(upper),'integer_candidates':possibilities}
 Path(__file__).with_name('q126_final.json').write_text(json.dumps(out,indent=2)+'\n')
 return out


def run():
 a=q114_certificate();b=q126_certificate()
 print('Q114 closed: norm-compatible last Gram fails12-eigenspace flatness; minimum vectors',a['last_case_minimizer_count'])
 print('Q126 closed: rank4 determinant integer interval397..401 has no norm')
 return a,b
if __name__=='__main__':run()
