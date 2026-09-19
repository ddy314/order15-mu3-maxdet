"""Exact cross-vector quadratic forms over finite Eisenstein norm shells."""
from small_catalogue import *
from refine_q114 import cross_norm_patterns,allocations
from math import lcm

@njit(cache=True)
def quadratic_scores(Ja,Jb,opts,counts):
 n=len(counts);total=1
 for c in counts:total*=c
 scores=np.empty(total,dtype=np.int64)
 for serial in range(total):
  za=np.empty(n,dtype=np.int64);zb=np.empty(n,dtype=np.int64);s=serial
  for i in range(n):
   k=s%counts[i];s//=counts[i];za[i]=opts[i,k,0];zb[i]=opts[i,k,1]
  value=0
  for i in range(n):
   value+=Ja[i,i]*(za[i]*za[i]-za[i]*zb[i]+zb[i]*zb[i])
   for j in range(i+1,n):
    # conj(z_i) J_ij z_j; 2 Re(a+b omega)=2a-b.
    a=za[i]-zb[i];b=-zb[i];c=Ja[i,j];d=Jb[i,j]
    x=a*c-b*d;y=a*d+b*c-b*d
    a=x*za[j]-y*zb[j];b=x*zb[j]+y*za[j]-y*zb[j]
    value+=2*a-b
  scores[serial]=value
 return np.unique(scores),total

def eis_shell(n):
 out=[]
 for a in range(-10,11):
  for b in range(-10,11):
   if a*a-a*b+b*b==n:out.append((a,b))
 return tuple(out)

def orbit_bases(n):
 left=set(eis_shell(n));reps=[]
 while left:
  a=min(left);reps.append(a);left-={mul(a,u) for u in UNITS}
 return tuple(reps)

@lru_cache(None)
def cross_spectrum(Akey,t):
 A=[list(row) for row in Akey];J=inverse(A);n=len(A)
 den=lcm(*(x.denominator if isinstance(x,F) else 1 for row in J for z in row for x in z))
 Ja=np.array([[int(z[0]*den) for z in row]for row in J],dtype=np.int64)
 Jb=np.array([[int(z[1]*den) for z in row]for row in J],dtype=np.int64)
 scores=set();count=0
 for ns in cross_norm_patterns(n,t):
  options=[orbit_bases(ns[0])]+[eis_shell(x) for x in ns[1:]]
  counts=np.array([len(o) for o in options],dtype=np.int64)
  opts=np.zeros((n,max(counts),2),dtype=np.int64)
  for i,o in enumerate(options):opts[i,:len(o),:]=o
  values,c=quadratic_scores(Ja,Jb,opts,counts);scores.update(map(int,values));count+=c
 return tuple(sorted(F(x,den)for x in scores)),count

def tree_matrix(v):
 A=[[(15 if i==j else 0,0)for j in range(v)]for i in range(v)]
 for i in range(v-1):A[i][i+1]=A[i+1][i]=(3,0)
 return tuple(tuple(row)for row in A)

def convolution_sets(xs,ys):return {x+y for x in xs for y in ys}

def run_q114():
 r=json.loads(Path(__file__).with_name('q114_refinement.json').read_text());rows=[]
 for row in r['rows']:
  cfg=tuple(map(tuple,row['configuration']));e=row['internal'];core=next(c for c in cfg if c[0]>=4)
  A=tuple(tuple(tuple(z)for z in ar)for ar in row['worst']['matrix']);D=bareiss([list(x)for x in A])[0]
  other=list(cfg);other.remove(core);riso=other.count((1,0));ts=[c for c in other if c[0]>1];assert all(c in ((2,1),(3,2))for c in ts)
  tmax=(114-e-42)//9;values=set();states=0
  corecurves=[]
  for t in range(tmax+1):
   sc,cnt=cross_spectrum(A,t);corecurves.append(sc);states+=cnt
  for aa in allocations(1+len(ts)+1,tmax):
   si=F(3*riso+9*aa[-1],15)
   if not riso and aa[-1]:continue
   losses=set(corecurves[aa[0]])
   for (v,u),t in zip(ts,aa[1:-1]):
    sc,cnt=cross_spectrum(tree_matrix(v),t);states+=cnt;losses=convolution_sets(losses,sc)
   fullD=D*prod(T[str(v)]['max_determinant']for v,u in ts)*15**riso
   values.update(F(fullD)*(15-si-s)for s in losses)
  over=sorted(v for v in values if v>B)
  survivors=[v for v in over if v.denominator==1 and v.numerator%3**14==0 and is_norm(v.numerator)]
  result={'internal':e,'configuration':cfg,'core_cross_minima':[str(min(sc))for sc in corecurves],'quadratic_vector_states':states,'distinct_full_determinants':len(values),'max_over_B':float(max(values)/B),'above_B': [{'value':str(v),'over_B':float(v/B),'integral':v.denominator==1,'eisenstein_norm':is_norm(v.numerator)if v.denominator==1 else False,'divisible_3_14':v.denominator==1 and v.numerator%3**14==0}for v in over], 'survivors':[str(v)for v in survivors]}
  print(e,core,'max',result['max_over_B'],'over',len(over),'normsurvivors',len(survivors), 'mins',result['core_cross_minima'],flush=True);rows.append(result)
 Path(__file__).with_name('q114_cross_exact.json').write_text(json.dumps({'rows':rows},indent=2)+'\n')
if __name__=='__main__':run_q114()
