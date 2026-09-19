"""Independent arbitrary-precision checks of accelerated exact primitives."""
from rank_catalogue import *
from pure_catalogue import kernel_min_from_poly
import random, platform, sympy as sp

def run():
 rng=random.Random(20260918);checks=[]
 for n in range(2,13):
  for trial in range(5):
   a=np.zeros((n,n),dtype=np.int64);b=np.zeros((n,n),dtype=np.int64)
   edges=list(combinations(range(n),2));rng.shuffle(edges)
   for i,j in edges[:min(len(edges),trial+2)]:
    z=rng.choice(UNITS);a[i,j],b[i,j]=z;a[j,i],b[j,i]=conj(z)
   A=[[(int(a[i,j]),int(b[i,j]))for j in range(n)]for i in range(n)]
   p=charpoly(A);fast=tuple(map(int,characteristic(a,b)));assert p==fast
   di=kernel_projection_diag(A);num,den=map(int,kernel_min_from_poly(a,b,np.array(p,dtype=np.int64)))
   assert F(num,den)==min(di)
   null,support=map(int,rank_support(a,b));r=n
   while r>0 and p[r]==0:r-=1
   assert null==n-r and support==sum(x>0 for x in di)
   aa=a+5*np.eye(n,dtype=np.int64)
   M=[[add(A[i][j],(5 if i==j else 0,0))for j in range(n)]for i in range(n)]
   exact=bareiss(M);assert exact==gaussian_det(M)==(int(det_hermitian(aa,b)),0)
   checks.append({'n':n,'trial':trial,'characteristic_equal':True,'determinants_equal':True,'kernel_equal':True})
 for values in ((1,1,1),(1,1,3,4),(1,3,3,4,4),(1,)*10):
  actual=list(unique_permutations(values));assert len(actual)==len(set(actual))
  if len(values)<=5:assert set(actual)==set(permutations(values))
 result={'checks':checks,'matrix_cases':len(checks),'seed':20260918,'all_passed':True,
 'python':platform.python_version(),'numpy':np.__version__,'networkx':nx.__version__,'sympy':sp.__version__}
 Path(__file__).with_name('primitive_validation.json').write_text(json.dumps(result,indent=2)+'\n')
 print('Validated',len(checks),'exact matrices against independent arbitrary-precision algorithms',flush=True)
 return result
if __name__=='__main__':run()
