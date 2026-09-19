"""Exact Z[omega]/Q(omega) arithmetic. No float is used for decisions."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
if not __debug__:raise RuntimeError('Verification requires assertions: do not use python -O.')
B=2**22*3**20*19
ZERO=(0,0);ONE=(1,0);ROOTS=(ONE,(0,1),(-1,-1));UNITS=ROOTS+tuple((-a,-b) for a,b in ROOTS)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return x[0]-y[0],x[1]-y[1]
def scale(x,c):return x[0]*c,x[1]*c
def mul(x,y):
 a,b=x;c,d=y;return a*c-b*d,a*d+b*c-b*d
def conj(x):return x[0]-x[1],-x[1]
def norm(x):return x[0]**2-x[0]*x[1]+x[1]**2
def div(x,y):
 z=mul(x,conj(y));n=norm(y)
 if not n:raise ZeroDivisionError
 return F(z[0],n),F(z[1],n)
def ediv(x,y):
 z=div(x,y)
 if any(a.denominator!=1 for a in z):raise ArithmeticError('Nonintegral Bareiss quotient')
 return tuple(map(int,z))
def esum(xs):
 a=b=0
 for x,y in xs:a+=x;b+=y
 return a,b
def mm(A,C):
 cols=list(zip(*C));return [[esum(mul(x,y) for x,y in zip(r,c)) for c in cols] for r in A]
def adj(A):return [[conj(A[i][j]) for i in range(len(A))] for j in range(len(A[0]))]
def gram(A):return mm(A,adj(A))
def bareiss(A):
 n=len(A)
 if n==0:return ONE
 C=[r[:] for r in A];last=ONE;sign=1
 for k in range(n-1):
  p=next((i for i in range(k,n) if C[i][k]!=ZERO),None)
  if p is None:return ZERO
  if p!=k:C[k],C[p]=C[p],C[k];sign=-sign
  pivot=C[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):C[i][j]=ediv(sub(mul(pivot,C[i][j]),mul(C[i][k],C[k][j])),last)
  for i in range(k+1,n):C[i][k]=ZERO
  last=pivot
 return scale(C[-1][-1],sign)
def gaussian_det(A):
 C=[r[:] for r in A];out=ONE;n=len(C)
 for k in range(n):
  p=next((i for i in range(k,n) if C[i][k]!=ZERO),None)
  if p is None:return ZERO
  if p!=k:C[k],C[p]=C[p],C[k];out=neg(out)
  pivot=C[k][k];out=mul(out,pivot)
  for i in range(k+1,n):
   z=div(C[i][k],pivot)
   for j in range(k+1,n):C[i][j]=sub(C[i][j],mul(z,C[k][j]))
   C[i][k]=ZERO
 return out

def charpoly(A):
 n=len(A);p=[1];tr=[0];power=[r[:] for r in A]
 for k in range(1,n+1):
  v=esum(power[i][i] for i in range(n));assert v[1]==0;tr.append(v[0])
  t=sum(p[k-j]*tr[j] for j in range(1,k+1));assert t%k==0;p.append(-t//k)
  if k<n:power=mm(power,A)
 return tuple(p)
def peval(p,x):
 out=0
 for c in p:out=out*x+c
 return out
def convolution(p,q):
 out=[0]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):out[i+j]+=a*b
 return tuple(out)
def rref(A):
 C=[r[:] for r in A];m=len(C);n=len(C[0]);r=0;piv=[]
 for j in range(n):
  p=next((i for i in range(r,m) if C[i][j]!=ZERO),None)
  if p is None:continue
  C[r],C[p]=C[p],C[r];v=C[r][j];C[r]=[div(x,v) for x in C[r]]
  for i in range(m):
   if i!=r and C[i][j]!=ZERO:
    v=C[i][j];C[i]=[sub(a,mul(v,b)) for a,b in zip(C[i],C[r])]
  piv.append(j);r+=1
  if r==m:break
 return C,tuple(piv)
def inverse(A):
 n=len(A);C,p=rref([row+[(int(i==j),0) for j in range(n)] for i,row in enumerate(A)])
 if p[:n]!=tuple(range(n)):raise ArithmeticError('singular inverse')
 return [r[n:] for r in C]
def kernel_projection_diag(A):
 n=len(A);C,p=rref(A);free=[j for j in range(n) if j not in p]
 if not free:return [F(0)]*n
 V=[[ZERO]*len(free) for _ in range(n)]
 for k,j in enumerate(free):
  V[j][k]=ONE
  for i,c in enumerate(p):V[c][k]=neg(C[i][j])
 P=mm(mm(V,inverse(mm(adj(V),V))),adj(V))
 assert all(P[i][i][1]==0 for i in range(n))
 return [F(P[i][i][0]) for i in range(n)]
def sqrt_lo(x,scale=10**15):
 x=F(x)
 if x<0:raise ValueError('negative square root')
 a=F(isqrt(x.numerator*scale*scale//x.denominator),scale);assert a*a<=x;return a
def sqrt_hi(x,scale=10**15):
 a=sqrt_lo(x,scale);return a if a*a==x else a+F(1,scale)
@lru_cache(None)
def stationary(d,q,mean=F(15)):
 q=F(q);mean=F(mean)
 if d==0:return F(1) if q==0 else F(0)
 if mean<=0 or q<0:return F(0)
 if q==0:return mean**d
 if d==1 or 2*q>=d*(d-1)*mean**2:return F(0)
 vals=[]
 for m in range(1,d):
  t=sqrt_lo(2*q/F(d*m*(d-m)));high=mean+(d-m)*t;low=mean-m*t
  if low>0:vals.append(high**m*low**(d-m))
 return max(vals,default=F(0))
@lru_cache(None)
def is_norm(n):
 if n<0:return False
 if n==0:return True
 from sympy import factorint
 return all(p%3!=2 or e%2==0 for p,e in factorint(n).items())
def rat(x):
 x=F(x);return {'numerator':x.numerator,'denominator':x.denominator,'approx':float(x)}
