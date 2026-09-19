"""Exact two-hub Schur enumeration for the unresolved energy18 size13 blocks.
This is exploratory until all residual spectra/realizations have been excluded.
Within one Gram, the ratio of the two minimum-norm cross entries has a fixed
mu3 coset; one hub switch makes every ratio a third root of unity.
"""
from cross_exact import *
from itertools import combinations_with_replacement

@lru_cache(None)
def component_states(v):
 A=tree_matrix(v);J=inverse([list(r)for r in A]);delta=(1,-1);out=set()
 for phases in product(UNITS,repeat=v-1):
  u=(delta,)+tuple(mul(delta,z)for z in phases)
  for ratios in product(ROOTS,repeat=v):
   w=tuple(mul(z,t)for z,t in zip(u,ratios))
   Ju=[esum(mul(a,z)for a,z in zip(row,u))for row in J]
   Jw=[esum(mul(a,z)for a,z in zip(row,w))for row in J]
   x=esum(mul(conj(a),b)for a,b in zip(u,Ju));y=esum(mul(conj(a),b)for a,b in zip(w,Jw));z=esum(mul(conj(a),b)for a,b in zip(u,Jw))
   assert x[1]==y[1]==0
   out.add((F(x[0]),F(y[0]),F(z[0]),F(z[1])))
 return tuple(out)

@njit(cache=True)
def all_scores(states,shifts,den,r):
 vals=np.empty(len(states)*len(shifts),dtype=np.int64);i=0
 for x,y,a,b in states:
  t=(15*den-r*(den//5)-x)*(15*den-r*(den//5)-y)
  for c,d in shifts:
   aa=c-a;bb=d-b;vals[i]=t-(aa*aa-aa*bb+bb*bb);i+=1
 return np.unique(vals)

def run():
 out=[]
 for name,v,r in (('P3',3,10),('2K2',4,9)):
  if name=='P3':ss=component_states(3);Dcore=3105
  else:
   raw=component_states(2);ss=tuple(set(tuple(a+b for a,b in zip(x,y))for x,y in combinations_with_replacement(raw,2)));Dcore=216**2
  den=lcm(5,*(x.denominator for row in ss for x in row));arr=np.array([[int(x*den)for x in row]for row in ss],dtype=np.int64)
  sums={esum([scale(ROOTS[0],a),scale(ROOTS[1],b),scale(ROOTS[2],r-a-b)])for a in range(r+1)for b in range(r-a+1)}
  for h in (0,3,9,12):
   shifts={sub(scale(g,den),scale(s,den//5))for g in eis_shell(h)for s in sums}
   scores=all_scores(arr,np.array(sorted(shifts),dtype=np.int64),den,r)
   scaleD=F(Dcore*15**r,den**2);above=[int(x)for x in scores if F(int(x))*scaleD>B]
   norms=[]
   for t in above:
    D=F(t)*scaleD
    if D.denominator==1 and D.numerator%3**14==0 and is_norm(D.numerator):norms.append(int(D))
   row={'family':name,'outside_norm':h,'internal_energy':18,'cross_energy':78,'Q':96+h,'core_states':len(ss),'isolate_sum_states':len(sums),'scalar_shift_states':len(shifts),'score_count':len(scores),'scale':str(scaleD),'max_over_B':float(F(int(scores[-1]))*scaleD/B),'above_B_score_count':len(above),'norm_survivors':norms,'closed_relaxation':not norms}
   print(name,'Q',96+h,'max',row['max_over_B'],'scoreabove',len(above),'norms',norms,flush=True);out.append(row)
 Path(__file__).with_name('size13_energy18.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':run()
