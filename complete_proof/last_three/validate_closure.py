"""Independent direct quadratic enumeration and exact matrix checks."""
from mixed96_columns import *
from collections import Counter
from math import lcm

@njit(cache=True)
def direct_scan(diagonal,cost,n):
 roots=np.zeros(n,dtype=np.int64);hits=[]
 for serial in range(3**(n-1)):
  s=serial
  for i in range(1,n):roots[i]=s%3;s//=3
  q=diagonal
  for i in range(n):
   for j in range(i+1,n):q+=cost[i,j,roots[i],roots[j]]
  if q==0:hits.append(serial)
 return hits

def direct_hits(a,b,L):
 n=len(a);cost=np.zeros((n,n,3,3),dtype=np.int64)
 for i in range(n):
  for j in range(i+1,n):
   z=(int(a[i,j]),int(b[i,j]))
   for k,x in enumerate(ROOTS):
    for l,y in enumerate(ROOTS):
     w=mul(mul(conj(x),z),y);cost[i,j,k,l]=2*w[0]-w[1]
 return direct_scan(int(np.trace(a))-L,cost,n)

def run():
 start=time.monotonic();d1=json.loads(Path(__file__).with_name('column_test.json').read_text());d2=json.loads(Path(__file__).with_name('mixed96_columns.json').read_text());assert d1['complete']and d2['complete']
 cases=[];witnesses=[];fullchecks=[]
 for group,data,builder in [('size13',d1,full_gram),('mixed96',d2,full_mixed_gram)]:
  for row in data['rows']:
   i=row['candidate_index'];G=builder(row);assert len(G)==15 and all(G[j][j]==(15,0)for j in range(15))
   assert sum(norm(G[i][j])for i in range(15)for j in range(i))==row['Q']
   assert bareiss(G)==(row['D'],0)
   colors=[0]*13+[1,2 if row['Q']==99 else 1]if group=='size13'else[0]*14+[1]
   for u in range(15):
    for v in range(15):
     target=scale((1,-1),colors[v]-colors[u]);assert all((G[u][v][k]-target[k])%3==0 for k in (0,1))
   a,b,L=inverse_arrays(G)
   hits=list(map(int,direct_hits(a,b,L)));assert hits==row['column_serials']
   assert hits==sorted(set(hits))
   fullchecks.append({'group':group,'candidate_index':i,'assignments':3**14,'matching_columns':len(hits),'identical_to_incremental':True});print('direct',group,i,'pass',time.monotonic()-start,flush=True)
   cols=[decode(s)for s in row['column_serials']]
   # Independent annihilator certificate: no rank assertion alone is needed.
   W=[[esum(mul(v[j],conj(v[k]))for v in cols)for k in range(15)]for j in range(15)]
   RR,piv=rref(W)
   if len(piv)<15:
    free=next(j for j in range(15)if j not in piv);w=[ZERO]*15;w[free]=ONE
    for ri,c in enumerate(piv):w[c]=neg(RR[ri][free])
    den=lcm(*(F(z).denominator for pair in w for z in pair));w=[tuple(int(F(z)*den)for z in pair)for pair in w]
    assert any(z!=ZERO for z in w)
    for v in cols:assert esum(mul(conj(c),z)for c,z in zip(w,v))==ZERO
    witnesses.append({'group':group,'candidate_index':i,'columns':len(cols),'annihilator':w,'span_rank':len(piv)})
   cases.append({'group':group,'candidate_index':i,'determinant_exact':True,'color_congruence_checked':True,'inverse_identity_checked':True,'span_rank':len(piv)})
 assert len(cases)==194 and len(witnesses)==193 and len(fullchecks)==194
 result={'matrix_checks':cases,'annihilators':witnesses,'direct_full_enumerations':fullchecks,'all_passed':True,'seconds':time.monotonic()-start}
 Path(__file__).with_name('closure_validation.json').write_text(json.dumps(result,indent=2)+'\n');print('VALIDATED194 matrices,193 annihilators,194 full direct scans',result['seconds'],flush=True);return result
if __name__=='__main__':run()
