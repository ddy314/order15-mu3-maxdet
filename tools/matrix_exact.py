"""Small exact Eisenstein-integer toolkit. No floating-point decisions."""
from __future__ import annotations
from itertools import combinations

ROOTS=((1,0),(0,1),(-1,-1))

def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def sub(x,y): return add(x,neg(y))

def mul(x,y):
    a,b=x;c,d=y
    return (a*c-b*d,a*d+b*c-b*d)

def conj(x):return (x[0]-x[1],-x[1])

def norm(x):
    a,b=x
    return a*a-a*b+b*b

def div(x,y):
    z=mul(x,conj(y));n=norm(y)
    if not n or z[0]%n or z[1]%n:
        raise ArithmeticError('non-exact Eisenstein division')
    return (z[0]//n,z[1]//n)

def gram(E):
    return [[tuple(sum(ROOTS[(a-b)%3][s] for a,b in zip(u,v)) for s in (0,1)) for v in E] for u in E]

def det(A):
    A=[[tuple(z) for z in row] for row in A];n=len(A);old=(1,0);sign=1
    for k in range(n-1):
        p=next((i for i in range(k,n) if A[i][k]!=(0,0)),None)
        if p is None:return (0,0)
        if p!=k:A[p],A[k]=A[k],A[p];sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                A[i][j]=div(sub(mul(pivot,A[i][j]),mul(A[i][k],A[k][j])),old)
        for i in range(k+1,n):A[i][k]=(0,0)
        old=pivot
    return tuple(sign*x for x in A[-1][-1])

def orthogonal_sets(E):
    G=gram(E);n=len(E)
    for k in range(n,0,-1):
        best=[list(c) for c in combinations(range(n),k)
              if all(G[i][j]==(0,0) for i,j in combinations(c,2))]
        if best:return k,best
