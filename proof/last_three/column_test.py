"""Exact necessary column test for HH*=G, H in mu3^(15x15).
Tests all 3^14 columns up to global mu3 scaling. Uses an integer-scaled
inverse and incrementally updates its quadratic form, with no floating point.
"""
from gram_candidates import *
from math import lcm,gcd

@njit(cache=True)
def column_patterns(a,b,target):
 n=len(a);za=np.ones(n,dtype=np.int64);zb=np.zeros(n,dtype=np.int64)
 ra=np.zeros(n,dtype=np.int64);rb=np.zeros(n,dtype=np.int64)
 for i in range(n):
  for j in range(n):ra[i]+=a[i,j];rb[i]+=b[i,j]
 value=0
 for i in range(n):value+=ra[i]
 # Since the quadratic form is real, sum rb must vanish.
 if rb.sum()!=0:raise ValueError('non-real quadratic')
 ids=np.zeros(n,dtype=np.int64);roots=np.array([[1,0],[0,1],[-1,-1]],dtype=np.int64)
 hits=[]
 total=3**(n-1)
 for serial in range(total):
  if value==target:hits.append(serial)
  if serial==total-1:break
  for k in range(1,n):
   old=ids[k];new=(old+1)%3;ids[k]=new
   da=roots[new,0]-za[k];db=roots[new,1]-zb[k]
   # q(z+d e_k)=q(z)+2Re(conj(d)(Jz)_k)+J_kk |d|^2.
   ca=da-db;cb=-db;pa=ca*ra[k]-cb*rb[k];pb=ca*rb[k]+cb*ra[k]-cb*rb[k]
   value+=2*pa-pb+a[k,k]*(da*da-da*db+db*db)
   za[k]=roots[new,0];zb[k]=roots[new,1]
   for i in range(n):
    ra[i]+=a[i,k]*da-b[i,k]*db;rb[i]+=a[i,k]*db+b[i,k]*da-b[i,k]*db
   if new!=0:break
 return hits

def inverse_arrays(G):
 J=inverse(G);den=lcm(*(F(z).denominator for row in J for pair in row for z in pair))
 a=np.array([[int(z[0]*den)for z in row]for row in J],dtype=np.int64)
 b=np.array([[int(z[1]*den)for z in row]for row in J],dtype=np.int64)
 M=max(max(abs(int(x))for x in a.flat),max(abs(int(x))for x in b.flat))
 # Conservative bound covering all quadratic and incremental-update intermediates.
 assert den<2**63 and 100*len(G)**2*M<2**63
 # Verify the inverse identity in independent arbitrary-precision arithmetic.
 JJ=[[(int(a[i,j]),int(b[i,j]))for j in range(len(G))]for i in range(len(G))]
 assert mm(G,JJ)==[[(den*int(i==j),0)for j in range(len(G))]for i in range(len(G))]
 return a,b,den

def run():
 data=json.loads(Path(__file__).with_name('gram_candidates.json').read_text());out=[];start=time.monotonic()
 for i,row in enumerate(data['rows']):
  G=full_gram(row);a,b,den=inverse_arrays(G)
  hits=column_patterns(a,b,den)
  rr={**row,'candidate_index':i,'inverse_scale':den,'column_count_up_to_scalar':len(hits),'column_serials':list(map(int,hits))}
  out.append(rr)
  print(i,row['Q'],row['family'],row['D'],'cols',len(hits),'time',round(time.monotonic()-start,2),flush=True)
  Path(__file__).with_name('column_test.json').write_text(json.dumps({'rows':out,'processed':len(out),'total':len(data['rows']),'complete':len(out)==len(data['rows'])},indent=2)+'\n')
 return out
if __name__=='__main__':run()
