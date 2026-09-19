from small_catalogue import *

@njit(cache=True)
def phase_candidates(n,edges,base,ischord,units,minimum_det):
 chord_count=0
 for c in ischord:chord_count+=c
 survivors=[];max_discarded=0;count=6**chord_count
 for serial in range(count):
  a=np.zeros((n,n),dtype=np.int64);b=np.zeros((n,n),dtype=np.int64)
  for i in range(n):a[i,i]=5
  s=serial
  for j in range(len(edges)):
   x=base[j,0];y=base[j,1]
   if ischord[j]:
    k=s%6;s//=6;z=units[k,0];w=units[k,1];x,y=x*z-y*w,x*w+y*z-y*w
   v=edges[j,0];w=edges[j,1];a[v,w]=x;b[v,w]=y;a[w,v]=x-y;b[w,v]=-y
  d=det_hermitian(a,b)
  if d>=minimum_det:survivors.append((serial,d))
  elif d>max_discarded:max_discarded=d
 return survivors,max_discarded,count

def make_matrix(n,edges,bs,chords,serial):
 A=[[(15 if i==j else 0,0) for j in range(n)] for i in range(n)]
 for e,z,is_chord in zip(edges,bs,chords):
  if is_chord:z=mul(z,UNITS[serial%6]);serial//=6
  i,j=e;A[i][j]=scale(z,3);A[j][i]=conj(A[i][j])
 return A

@lru_cache(None)
def cross_norm_patterns(v,t):
 if v==0:return ((),) if t==0 else ()
 return tuple((3+9*k,)+tail for k in range(t+1) if 3+9*k in {z['norm'] for z in CAT} for tail in cross_norm_patterns(v-1,t-k))

def inverse_loss_curve(A,L,max_t):
 J=inverse(A);n=len(A)
 diag=[F(J[i][i][0]) for i in range(n)]
 assert all(J[i][i][1]==0 and diag[i]>0 for i in range(n))
 off=[(i,j,norm(J[i][j])) for i in range(n) for j in range(i)]
 out=[]
 for t in range(max_t+1):
  spectral=F(3*n+9*t)/L
  lower=min(sum(diag[i]*p[i] for i in range(n))-2*sum(sqrt_hi(w*p[i]*p[j]) for i,j,w in off) for p in cross_norm_patterns(n,t))
  out.append(max(spectral,lower))
 return out

def weak_curve(v,u,t):return F(3*v+9*t)/cb(v,u)[1]

def allocations(k,t):
 if k==0:
  if t==0:yield ()
  return
 for x in range(t+1):
  for tail in allocations(k-1,t-x):yield (x,)+tail

def run():
 exploration=json.loads(Path(__file__).with_name('exploration.json').read_text())
 todo=exploration['114']['res14'];out=[];start=time.monotonic()
 for e,cfg0,old_ratio in todo:
  cfg=tuple(map(tuple,cfg0));cores=[c for c in cfg if c[0]>=4];assert len(cores)==1
  core=cores[0];v,u=core;other=list(cfg);other.remove(core)
  coreD,L=cb(v,u)
  ds,Ls=zip(*(cb(a,b)for a,b in cfg))
  correction=sum(F(3*a)/lam for (a,b),lam in zip(cfg,Ls))+F(114-e-42)/max(Ls)
  full_raw=prod(ds)*max(F(0),15-correction)
  # Screening must scale the raw Schur product, never the independent
  # rank/trace minimum. Both coincide in these eight present branches.
  assert full_raw==bound14(cfg,114-e)
  unit_upper=full_raw/coreD
  threshold=B/unit_upper/(3**v);minimum=threshold.numerator//threshold.denominator+1
  tmax=(114-e-42)//9
  prodother=prod(cb(a,b)[0] for a,b in other)
  isolate_count=other.count((1,0));groups=[(a,b) for a,b in other if a>1]
  if isolate_count:groups.append((isolate_count,-1))
  others=[]
  for tcore in range(tmax+1):
   best=F(10**9)
   for aa in allocations(len(groups),tmax-tcore):
    s=sum(F(3*a+9*t)/15 if b==-1 else weak_curve(a,b,t) for (a,b),t in zip(groups,aa))
    best=min(best,s)
   if best<F(10**9):others.append((tcore,best))
  best=F(0);candidate_count=total=0;prunedmax=0;worst=None
  graphs=[g for g in nx.graph_atlas_g() if len(g)==v and nx.is_connected(g) and g.number_of_edges()<=u]
  seen=set()
  for g in graphs:
   edges=sorted(tuple(sorted(x)) for x in g.edges());m=len(edges)
   tree={tuple(sorted(x)) for x in nx.minimum_spanning_tree(g).edges()};chords=[int(x not in tree) for x in edges]
   for wp in weight_vectors(m,u):
    for ws in unique_permutations(wp):
     for bs in product(*(norm_bases(w) for w in ws)):
      survivors,discarded,count=phase_candidates(v,np.array(edges,dtype=np.int64),np.array(bs,dtype=np.int64),np.array(chords,dtype=np.int64),EUNITS,minimum)
      total+=count;prunedmax=max(prunedmax,discarded)
      for serial,d in survivors:
       A=make_matrix(v,edges,bs,chords,serial);key=tuple(tuple(row) for row in A)
       if key in seen:continue
       seen.add(key);candidate_count+=1
       assert bareiss(A)==(int(d)*3**v,0)
       losses=inverse_loss_curve(A,L,tmax)
       for tcore,sother in others:
        value=F(int(d)*3**v)*prodother*max(F(0),15-sother-losses[tcore])
        if value>best:best=value;worst={'core_cross_extra_units':tcore,'core_det':int(d)*3**v,'core_loss':str(losses[tcore]),'outside_loss':str(sother),'matrix':A}
  prunedupper=F(prunedmax*3**v)*unit_upper;best=max(best,prunedupper)
  row={'internal':e,'configuration':cfg,'phase_states':total,'surviving_matrices':candidate_count,'upper':rat(best),'over_B':float(best/B),'closed':best<=B,'worst':worst}
  out.append(row);print(e,core,'phases',total,'survivors',candidate_count,'ratio',float(best/B),'sec',round(time.monotonic()-start,2),flush=True)
 Path(__file__).with_name('q114_refinement.json').write_text(json.dumps({'rows':out,'seconds':time.monotonic()-start},indent=2)+'\n')
 return out
if __name__=='__main__':run()
