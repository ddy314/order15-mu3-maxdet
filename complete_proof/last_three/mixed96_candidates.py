"""Exact Q96 double-(14,1) residuals, with all internal isolates retained."""
from column_test import *
from itertools import combinations

SPECS={
 'P3_C4':{'v':7,'r':7,'u':6,'components':[(0,1,2),(3,4,5,6)],'edges':[(0,1),(1,2),(3,4),(4,5),(5,6),(6,3)]},
 'C4':{'v':4,'r':10,'u':4,'components':[(0,1,2,3)],'edges':[(0,1),(1,2),(2,3),(3,0)]},
 'K23':{'v':5,'r':9,'u':6,'components':[(0,1,2,3,4)],'edges':[(i,j)for i in (0,1)for j in (2,3,4)]},
}

def cross_reps(k):
 if k==0:return ((1,-1),)
 left={z for z in eis_shell(3+9*k)if z[0]%3==1 and z[1]%3==2};out=[]
 while left:
  z=min(left);out.append(z);left-={mul(z,t)for t in ROOTS}
 return tuple(out)

def core_unscaled(spec):
 v=spec['v'];A=[[(15*int(i==j),0)for j in range(v)]for i in range(v)]
 for i,j in spec['edges']:A[i][j]=A[j][i]=(3,0)
 return A

def full_mixed_gram(row):
 spec=SPECS[row['family']];v=spec['v'];C=core_unscaled(spec);d=[UNITS[k]for k in row['phases']]
 G=[[(15*int(i==j),0)for j in range(15)]for i in range(15)]
 for i in range(v):
  for j in range(v):G[i][j]=mul(mul(d[i],C[i][j]),conj(d[j]))
 z=list(map(tuple,row['core_cross']))+list(map(tuple,row['isolate_cross']))
 assert len(z)==14
 for i,w in enumerate(z):G[i][14]=w;G[14][i]=conj(w)
 return G

def run():
 out=[];stats=[];start=time.monotonic()
 for name,spec in SPECS.items():
  v=spec['v'];r=spec['r'];u=spec['u'];maxk=(96-9*u-42)//9
  C=core_unscaled(spec);J=inverse(C);dc=bareiss(C)[0]*15**r
  # Representatives for cross magnitude allocations modulo isolated-row permutations.
  zs=[]
  for ks in product(range(maxk+1),repeat=v):
   left=maxk-sum(ks)
   if left<0:continue
   icases=[[(1,-1)]*r]if left==0 else []
   if left==1:icases=[[z]+[(1,-1)]*(r-1)for z in cross_reps(1)]
   if left==2:
    icases=[[z]+[(1,-1)]*(r-1)for z in cross_reps(2)]
    icases+=[[a,b]+[(1,-1)]*(r-2)for a in cross_reps(1)for b in cross_reps(1)]
   for zz in product(*(cross_reps(k)for k in ks)):
    for iz in icases:zs.append((zz,iz))
  roots=[c[0]for c in spec['components']];free=[i for i in range(v)if i not in roots]
  checked=0;n0=len(out);spectra=set()
  for ds in product(range(6),repeat=len(free)):
   d=[0]*v
   for i,k in zip(free,ds):d[i]=k
   for z,iz in zs:
    zz=[mul(conj(UNITS[k]),a)for k,a in zip(d,z)]
    Jz=[esum(mul(a,b)for a,b in zip(row,zz))for row in J]
    loss=esum(mul(conj(a),b)for a,b in zip(zz,Jz));assert loss[1]==0
    D=F(dc)*(15-loss[0]-F(sum(norm(a)for a in iz),15));checked+=1
    if D<=B or D.denominator!=1 or D.numerator%3**14 or not is_norm(D.numerator):continue
    out.append({'Q':96,'family':name,'phases':d,'core_cross':z,'isolate_cross':iz,'D':int(D)});spectra.add(int(D))
  rr={'family':name,'normalized_phase_cross_states':checked,'candidate_grams':len(out)-n0,'candidate_determinants':sorted(spectra)};stats.append(rr);print(rr,'sec',round(time.monotonic()-start,2),flush=True)
 result={'stats':stats,'rows':out};Path(__file__).with_name('mixed96_candidates.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':run()
