"""Additional exact catalogue for Q108, including both support-isolate regimes."""
from pure_catalogue import *
MAXU=12

def large_extended(maxu=MAXU):
 allg=[];stats=[]
 for n in range(8,maxu+1):
  buckets=defaultdict(list);tested=passed=0;limit_noiso=maxu-(16-n)//2
  for T in nx.nonisomorphic_trees(n):
   T=nx.convert_node_labels_to_integers(T);te=set(tuple(sorted(e))for e in T.edges());missing=[e for e in combinations(range(n),2)if e not in te];base_deg=[T.degree(i)for i in range(n)]
   for k in range(maxu-n+2):
    noleaf=n-1+k>limit_noiso
    if noleaf and sum(d==1 for d in base_deg)>2*k:continue
    for extra in combinations(missing,k):
     tested+=1;dd=base_deg[:]
     for i,j in extra:dd[i]+=1;dd[j]+=1
     if noleaf and min(dd)<2:continue
     passed+=1;G=T.copy();G.add_edges_from(extra)
     key=(G.number_of_edges(),tuple(sorted(dd)),nx.weisfeiler_lehman_graph_hash(G))
     if any(nx.is_isomorphic(G,H)for H in buckets[key]):continue
     buckets[key].append(G)
  gs=[g for vv in buckets.values()for g in vv];allg+=gs
  stats.append({'v':n,'extra_edge_sets':tested,'passing_degree':passed,'graphs':len(gs)})
  print('extended large',stats[-1],flush=True)
 return allg,stats

@njit(cache=True)
def batch_extended(n,edges,bases,chords,minimum_det,allow_zero):
 count=6**len(chords);out=np.zeros((count,n+4),dtype=np.int64)
 for serial in range(count):
  a=np.zeros((n,n),dtype=np.int64);b=np.zeros((n,n),dtype=np.int64);s=serial
  for j in range(len(edges)):
   x=bases[j,0];y=bases[j,1]
   if j in chords:
    k=s%6;s//=6;u=EUNITS[k,0];v=EUNITS[k,1];x,y=x*u-y*v,x*v+y*u-y*v
   i,k=edges[j,0],edges[j,1];a[i,k]=x;b[i,k]=y;a[k,i]=x-y;b[k,i]=-y
  aa=a.copy()
  for i in range(n):aa[i,i]=5
  D=det_hermitian(aa,b)
  if D<=minimum_det:continue
  p=characteristic(a,b);num,den=kernel_min_from_poly(a,b,p)
  if not allow_zero and num==0:continue
  for i in range(n+1):out[serial,i]=p[i]
  out[serial,n+1]=D;out[serial,n+2]=num;out[serial,n+3]=den
 return out

def run():
 start=time.monotonic();graphs=[g for g in nx.graph_atlas_g()if 2<=len(g)<=7 and nx.is_connected(g)and g.number_of_edges()<=MAXU]
 gg,stats=large_extended();graphs+=gg;types={};tested=0
 for G in graphs:
  G=nx.convert_node_labels_to_integers(G);n=len(G);m=G.number_of_edges();no_leaf=min(dict(G.degree()).values())>=2
  umax=MAXU if no_leaf else MAXU-(16-n)//2
  if m>umax:continue
  ed=sorted(tuple(sorted(e))for e in G.edges());tree={tuple(sorted(e))for e in nx.minimum_spanning_tree(G).edges()};ch=np.array([i for i,e in enumerate(ed)if e not in tree],dtype=np.int64);edges=np.array(ed,dtype=np.int64)
  minimum=B//(3**n*15**(15-n))
  for u in range(m,umax+1):
   if cb(n,u)[0]*15**(15-n)<=B:continue
   allow_zero=u+(16-n)//2<=MAXU
   for ws in weight_vectors(m,u):
    for order in unique_permutations(ws):
     for bs in product(*(norm_bases(w)for w in order)):
      arr=batch_extended(n,edges,np.array(bs,dtype=np.int64),ch,minimum,allow_zero);tested+=len(arr)
      for row in arr:
       if row[0]!=1:continue
       p=tuple(map(int,row[:n+1]));D=int(row[n+1])*3**n;val=F(int(row[n+2]),int(row[n+3]));key=(n,u,p)
       old=types.get(key)
       if old is None:types[key]={'v':n,'u':u,'p':p,'determinant':D,'min_kernel':str(val),'no_leaf':no_leaf,'phase_states':1}
       else:
        old['min_kernel']=str(max(F(old['min_kernel']),val));old['phase_states']+=1;old['no_leaf']|=no_leaf
 print('extended catalogue',len(types),'types',tested,'states',time.monotonic()-start,'sec',flush=True)
 out={'max_units':MAXU,'types':list(types.values()),'phase_states':tested,'graph_generation':stats,'seconds':time.monotonic()-start}
 Path(__file__).with_name('pure_extended.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':run()
