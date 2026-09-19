"""Exact low-energy all-same Gram component polynomials and kernel diagonals."""
from rank_catalogue import *
from collections import defaultdict

@njit(cache=True)
def kernel_min_from_poly(a,b,p):
 n=len(a);r=n
 while r>0 and p[r]==0:r-=1
 if r==n:return 0,1
 # q(A)/q(0), q=p/x^(n-r), is the orthogonal zero-eigenvalue projector.
 qa=np.eye(n,dtype=np.int64);qb=np.zeros((n,n),dtype=np.int64)
 for k in range(1,r+1):
  na=np.zeros((n,n),dtype=np.int64);nb=np.zeros((n,n),dtype=np.int64)
  for i in range(n):
   for j in range(n):
    x=0;y=0
    for z in range(n):
     x+=qa[i,z]*a[z,j]-qb[i,z]*b[z,j]
     y+=qa[i,z]*b[z,j]+qb[i,z]*a[z,j]-qb[i,z]*b[z,j]
    na[i,j]=x;nb[i,j]=y
   na[i,i]+=p[k]
  qa=na;qb=nb
 den=p[r];sgn=1 if den>0 else -1;num=10**12
 for i in range(n):
  assert qb[i,i]==0
  num=min(num,sgn*qa[i,i])
 return num,sgn*den

@njit(cache=True)
def phase_records(n,edges,bases,chords,cutoff):
 count=6**len(chords);ans=np.zeros((count,n+4),dtype=np.int64);kept=0
 for serial in range(count):
  a=np.zeros((n,n),dtype=np.int64);b=np.zeros((n,n),dtype=np.int64);s=serial
  for k in range(len(edges)):
   x=bases[k,0];y=bases[k,1]
   if k in chords:
    z=s%6;s//=6;c=EUNITS[z,0];d=EUNITS[z,1];x,y=x*c-y*d,x*d+y*c-y*d
   i=edges[k,0];j=edges[k,1];a[i,j]=x;b[i,j]=y;a[j,i]=x-y;b[j,i]=-y
  p=characteristic(a,b);D=0
  for x in p:D=D*(-5)+x
  if n%2:D=-D
  if D<=cutoff:continue
  num,den=kernel_min_from_poly(a,b,p)
  for j in range(n+1):ans[kept,j]=p[j]
  ans[kept,n+1]=D;ans[kept,n+2]=num;ans[kept,n+3]=den;kept+=1
 return ans[:kept],count

def sparse_large_graphs(max_units=10):
 """Complete via an unlabeled spanning tree plus every extra-edge set."""
 out=[];details=[]
 for n in range(8,max_units+1):
  buckets=defaultdict(list);generated=0;tested=0
  for tree in nx.nonisomorphic_trees(n):
   tree=nx.convert_node_labels_to_integers(tree);degs=[tree.degree(i)for i in range(n)]
   leaves=sum(1<<i for i,d in enumerate(degs)if d==1)
   missing=[e for e in combinations(range(n),2)if not tree.has_edge(*e)]
   for k in range(1,max_units-(n-1)+1):
    if leaves.bit_count()>2*k:continue
    for extra in combinations(missing,k):
     tested+=1;covered=0
     for i,j in extra:covered|=(1<<i)|(1<<j)
     if covered&leaves!=leaves:continue
     G=tree.copy();G.add_edges_from(extra)
     assert min(dict(G.degree()).values())>=2
     generated+=1;key=nx.weisfeiler_lehman_graph_hash(G)
     if any(nx.is_isomorphic(G,h)for h in buckets[key]):continue
     buckets[key].append(G)
  gs=[g for group in buckets.values()for g in group];out+=gs
  details.append({'v':n,'extra_edge_sets':tested,'passing_degree':generated,'nonisomorphic_graphs':len(gs)})
  print('large graphs',details[-1],flush=True)
 return out,details

def run(maxu=10):
 start=time.monotonic();data={};counts=0
 graphs=[g for g in nx.graph_atlas_g()if 2<=len(g)<=7 and nx.is_connected(g)and g.number_of_edges()<=maxu]
 large,details=sparse_large_graphs(maxu);graphs+=large
 for gi,G in enumerate(graphs):
  G=nx.convert_node_labels_to_integers(G);n=len(G);m=G.number_of_edges()
  ed=sorted(tuple(sorted(e))for e in G.edges());tree=set(tuple(sorted(e))for e in nx.minimum_spanning_tree(G).edges());ch=np.array([i for i,e in enumerate(ed)if e not in tree],dtype=np.int64);edges=np.array(ed,dtype=np.int64)
  cut=B//(3**n*15**(15-n));has_leaf=min(dict(G.degree()).values())==1
  for u in range(m,maxu+1):
   for ws in weight_vectors(m,u):
    for order in unique_permutations(ws):
     for bs in product(*(norm_bases(w)for w in order)):
      ans,cnt=phase_records(n,edges,np.array(bs,dtype=np.int64),ch,cut);counts+=cnt
      for rec in ans:
       p=tuple(map(int,rec[:n+1]));D,num,den=map(int,rec[n+1:]);key=(n,u,p)
       if n>=8 and num==0:continue
       d=F(num,den)
       old=data.get(key)
       if old is None:data[key]={'v':n,'u':u,'polynomial':p,'normalized_det':D,'min_kernel':str(d),'phase_states':1}
       else:
        old['min_kernel']=str(max(F(old['min_kernel']),d));old['phase_states']+=1
 out={'max_energy_units':maxu,'phase_states':counts,'large_graph_completeness':details,'types':list(data.values()),'seconds':time.monotonic()-start,'scope':'all connected v<=7; connected minimum-degree2 large supports; sufficient for Q90 isolate/no-isolate split'}
 Path(__file__).with_name('pure_catalogue.json').write_text(json.dumps(out,indent=2)+'\n');print('pure catalogue',len(data),'types',counts,'states',out['seconds'],'s',flush=True)
 return out
if __name__=='__main__':run()
