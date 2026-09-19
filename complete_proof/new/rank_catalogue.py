"""Exact component nullities and kernel support for low mixed-color energy."""
from small_catalogue import *
from itertools import combinations

@njit(cache=True)
def characteristic(a,b):
 n=len(a);pa=a.copy();pb=b.copy();coef=np.zeros(n+1,dtype=np.int64);coef[0]=1;tr=np.zeros(n+1,dtype=np.int64)
 for k in range(1,n+1):
  ta=0;tb=0
  for i in range(n):ta+=pa[i,i];tb+=pb[i,i]
  assert tb==0;tr[k]=ta
  s=0
  for j in range(1,k+1):s+=coef[k-j]*tr[j]
  assert s%k==0;coef[k]=-s//k
  if k<n:
   qa=np.zeros((n,n),dtype=np.int64);qb=np.zeros((n,n),dtype=np.int64)
   for i in range(n):
    for j in range(n):
     x=0;y=0
     for z in range(n):
      x+=pa[i,z]*a[z,j]-pb[i,z]*b[z,j]
      y+=pa[i,z]*b[z,j]+pb[i,z]*a[z,j]-pb[i,z]*b[z,j]
     qa[i,j]=x;qb[i,j]=y
   pa=qa;pb=qb
 return coef

@njit(cache=True)
def rank_of(a,b):
 p=characteristic(a,b);n=len(a);r=n
 while r>0 and p[r]==0:r-=1
 return r

@njit(cache=True)
def rank_support(a,b):
 n=len(a);r=rank_of(a,b);s=0
 for i in range(n):
  aa=np.empty((n-1,n-1),dtype=np.int64);bb=np.empty((n-1,n-1),dtype=np.int64);u=0
  for j in range(n):
   if j==i:continue
   v=0
   for k in range(n):
    if k==i:continue
    aa[u,v]=a[j,k];bb[u,v]=b[j,k];v+=1
   u+=1
  if rank_of(aa,bb)==r:s+=1
 return n-r,s

@njit(cache=True)
def phase_signatures(n,edges,bases,chords):
 nc=len(chords);count=6**nc;ans=np.zeros((count,3),dtype=np.int64)
 for phase_id in range(count):
  a=np.zeros((n,n),dtype=np.int64);b=np.zeros((n,n),dtype=np.int64);p=phase_id
  for k in range(len(edges)):
   x=bases[k,0];y=bases[k,1]
   if k in chords:
    z=p%6;p//=6;c=EUNITS[z,0];d=EUNITS[z,1];x,y=x*c-y*d,x*d+y*c-y*d
   i=edges[k,0];j=edges[k,1];a[i,j]=x;b[i,j]=y;a[j,i]=x-y;b[j,i]=-y
  null,support=rank_support(a,b)
  for i in range(n):a[i,i]=5
  det=det_hermitian(a,b)
  ans[phase_id,0]=null;ans[phase_id,1]=support;ans[phase_id,2]=det
 return ans

def independence(G):
 n=len(G)
 for k in range(n,0,-1):
  if any(not any(G.has_edge(i,j)for i,j in combinations(ss,2))for ss in combinations(range(n),k)):return k

def run():
 start=time.monotonic();data={};counts=0
 graphs=[g for g in nx.graph_atlas_g()if 2<=len(g)<=7 and nx.is_connected(g)and g.number_of_edges()<=7]
 graphs+=list(nx.nonisomorphic_trees(8))
 for G in graphs:
  G=nx.convert_node_labels_to_integers(G);n=len(G);m=G.number_of_edges();alpha=independence(G)
  ed=sorted(tuple(sorted(e))for e in G.edges());tree=set(tuple(sorted(e))for e in nx.minimum_spanning_tree(G).edges());ch=np.array([i for i,e in enumerate(ed)if e not in tree],dtype=np.int64);edges=np.array(ed,dtype=np.int64)
  for u in range(m,8):
   for ws in weight_vectors(m,u):
    for order in unique_permutations(ws):
     for bs in product(*(norm_bases(w)for w in order)):
      ans=phase_signatures(n,edges,np.array(bs,dtype=np.int64),ch);counts+=len(ans)
      for null,support,det in ans:
       key=f'{n},{u},{int(null)},{int(support)},{alpha}';D=int(det)*3**n
       old=data.get(key)
       if old is None:data[key]={'v':n,'u':u,'nullity':int(null),'kernel_support':int(support),'independence':alpha,'max_determinant':D,'phase_states':1}
       else:old['max_determinant']=max(old['max_determinant'],D);old['phase_states']+=1
 out={'max_energy_units':7,'phase_states':counts,'types':list(data.values()),'seconds':time.monotonic()-start}
 Path(__file__).with_name('rank_catalogue.json').write_text(json.dumps(out,indent=2)+'\n');print('rank catalogue',len(data),'types',counts,'states',out['seconds'],'s',flush=True)
 return out
if __name__=='__main__':run()
