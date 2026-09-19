"""Exact integer characteristic-polynomial sieve for low-rank pure Grams."""
from bounds import *
import sympy as sp,time
X=sp.Symbol('x')

def roots_real_and_positive_shift(coeffs):
 p=sp.Poly.from_list(coeffs,X,domain=sp.ZZ)
 intervals=sp.polys.polytools.intervals(p)
 real=sum(m for _,m in intervals)
 return real==p.degree(),intervals

def run117():
 out={};start=time.monotonic()
 # Rank <=4: an extra zero root is allowed, so the fourth-degree factor
 # retains integer coefficients also when the actual rank is smaller.
 factor4=15**11*3**4
 U4=15**11*stationary(4,117)
 Ns4=[n for n in range(B//factor4+1,int(U4//factor4)+1)if is_norm(n*factor4)]
 assert Ns4==[405]
 tests4=[]
 for n in Ns4:
  for t in range(-22,23):
   c=n-300-5*t
   coeff=[1,0,-13,-t,c];ok,ints=roots_real_and_positive_shift(coeff)
   tests4.append({'R':n,'e3':t,'e4':c,'all_roots_real':ok})
 assert not any(r['all_roots_real']for r in tests4)
 out['rank4']={'factor':factor4,'norm_candidates':Ns4,'triangle_e3_abs_cap':22,'tested':tests4,'closed':True}
 # Rank <=5, the sole support-size9 case has triangle sum at most6.
 v,u,rmin=9,13,6;d=5
 L=15+sqrt_hi(F(18*u*3,4))
 a=F(u,25)-F(2,125)*triangle_weight_cap(v,u,True)+F(3*u*u,125*d)/L
 U5=15**15*expneg_upper(a)
 factor5=15**10*3**5
 Ns5=[n for n in range(B//factor5+1,int(U5//factor5)+1)if is_norm(n*factor5)]
 tests5=[];surv=[]
 for n in Ns5:
  for t in range(-12,13):
   for e4 in range(-85,51):
    e5=n-1500-25*t-5*e4
    if abs(e5)>62:continue
    coeff=[1,0,-13,-t,e4,-e5]
    ok,ints=roots_real_and_positive_shift(coeff)
    tests5.append({'R':n,'e3':t,'e4':e4,'e5':e5,'all_roots_real':ok})
    if ok:surv.append({'coefficients':coeff,'intervals':str(ints),'R':n})
 out['rank5']={'factor':factor5,'upper':rat(U5),'norm_candidates':Ns5,'tested':tests5,'survivors':surv,'closed':not surv}
 print('rank4 candidates',Ns4,'tests',len(tests4),'rank5 candidates',Ns5,'tests',len(tests5),'survivors',surv,'seconds',time.monotonic()-start,flush=True)
 Path(__file__).with_name('q117_integer_sieve.json').write_text(json.dumps(out,indent=2)+'\n')
 return out
if __name__=='__main__':run117()
