#!/usr/bin/env python3
"""Exact recovered witness checks; all arithmetic uses Python integers.

The JSON certificate is regenerated, never trusted as a cache. Finite controls
are distinguished from the analytic theorems in the manuscript.
"""
from __future__ import annotations
import hashlib
import json
from itertools import product, permutations
from pathlib import Path
from matrix_exact import ROOTS, add, sub, neg, mul, conj, norm, det, gram

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/revision'
ZERO = (0, 0)
ONE = (1, 0)
B = 277868041444786176
SIX = (ONE, (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1))


def scale(z, k): return (k*z[0], k*z[1])
def total(zs):
    s = ZERO
    for z in zs: s = add(s, z)
    return s

def integer_det(rows):
    """Fraction-free elimination, independent of Eisenstein elimination."""
    a = [row[:] for row in rows]
    n = len(a); old = 1; sign = 1
    for k in range(n-1):
        pivot = next((i for i in range(k,n) if a[i][k]), None)
        if pivot is None: return 0
        if pivot != k: a[k], a[pivot] = a[pivot], a[k]; sign = -sign
        p = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value = p*a[i][j] - a[i][k]*a[k][j]
                q, r = divmod(value, old)
                assert not r
                a[i][j] = q
        for i in range(k+1,n): a[i][k] = 0
        old = p
    return sign*a[-1][-1]


def regular_representation(h):
    n = len(h); a = [[0]*(2*n) for _ in range(2*n)]
    for i,row in enumerate(h):
        for j,(x,y) in enumerate(row):
            a[2*i][2*j] = x; a[2*i][2*j+1] = -y
            a[2*i+1][2*j] = y; a[2*i+1][2*j+1] = x-y
    return a


def multiply_transpose(a, b):
    return [[sum(x*y for x,y in zip(r,s)) for s in b] for r in a]


def transpose(a): return [list(r) for r in zip(*a)]
def load(name): return json.loads((OUT/(name+'.json')).read_text())
def save(name, value): (OUT/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
def energy(g): return sum(norm(g[i][j]) for i in range(len(g)) for j in range(i))


def verify_bw25():
    w = load('bw25_16_witness'); n = len(w)
    assert n == 25 and all(len(r)==n and set(r)<={-1,0,1} for r in w)
    support = [[abs(x) for x in row] for row in w]
    zeros = [[1-x for x in row] for row in support]
    for a,diag,off in ((w,16,0),(support,16,10),(zeros,9,3)):
        for rows in (a,transpose(a)):
            assert all(sum(abs(x) for x in r)==diag for r in rows)
            g = multiply_transpose(rows,rows)
            assert g == [[diag if i==j else off for j in range(n)] for i in range(n)]
    assert integer_det(w) == -2**50
    exps = [[0]*51 for _ in range(51)]
    for i in range(25):
        for j in range(25):
            for a,b in product(range(2),repeat=2):
                exps[2*i+a][2*j+b] = ((1 if a==b else -1)*w[i][j])%3
        exps[2*i][50] = exps[2*i+1][50] = 2
    h = [[ROOTS[e] for e in r] for r in exps]
    g = gram(exps)
    delta = (1,-1)
    for i,j in product(range(51),repeat=2):
        expected = (51,0) if i==j else delta if j==50 else conj(delta) if i==50 else (3,0) if i//2==j//2 else ZERO
        assert g[i][j] == expected, (i,j,g[i][j])
    d = det(h); factor = 2**63*3**49
    assert d == (8*factor,-9*factor)
    squared = 48**25*54**24*2604
    assert norm(d) == squared == 2**126*3**98*217
    assert integer_det(regular_representation(h)) == squared
    assert energy(g) == energy(gram(transpose(exps))) == 375
    save('order51_witness',exps)
    save('bw25_support',zeros)
    return {'determinant':d,'squared_determinant':squared,'energy':375,'integer_matrix_order':102}


def verify_sixth():
    e = load('sixth_root_order15')
    assert len(e)==15 and all(len(r)==15 and set(r)<=set(range(6)) for r in e)
    h = [[SIX[x] for x in r] for r in e]
    d = det(h)
    assert d == (-335544320,-671088640)
    assert norm(d)==337769972052787200 > B
    def gmat(rows):
        return [[total(mul(x,conj(y)) for x,y in zip(r,s)) for s in rows] for r in rows]
    assert energy(gmat(h)) == energy(gmat(transpose(h))) == 53
    assert integer_det(regular_representation(h))==norm(d)
    return {'determinant':d,'squared_determinant':norm(d),'energy':53}


def gf_mul(a,b):
    t=0
    while b:
        if b&1:t^=a
        b>>=1;a<<=1
        if a&8:a^=13
    return t


def gf_trace(a):return a ^ gf_mul(a,a) ^ gf_mul(gf_mul(a,a),gf_mul(a,a))


def fano_factors():
    w = [[0 if gf_trace(gf_mul(a,b))==0 else -1 if gf_mul(a,b)==1 else 1 for b in range(1,8)] for a in range(1,8)]
    n = [[int(x==0) for x in row] for row in w]
    assert multiply_transpose(w,w)==[[4*int(i==j) for j in range(7)] for i in range(7)]
    assert multiply_transpose(n,n)==[[3 if i==j else 1 for j in range(7)] for i in range(7)]
    j0=0;i0=next(i for i in range(7) if n[i][j0]==0)
    u=[n[i][j0] for i in range(7)];v=n[i0]
    records=[]
    for r,s in product(range(3),repeat=2):
        rho,sigma=ROOTS[r],ROOTS[s]
        alpha=[rho if x else ONE for x in u];beta=[sigma if x else ONE for x in v]
        a=[[mul(mul(alpha[i],(3*n[i][j]-1,0)),beta[j]) for j in range(7)] for i in range(7)]
        ww=[[scale(mul(alpha[i],beta[j]),w[i][j]) for j in range(7)] for i in range(7)]
        for i in range(7):a[i][j0]=add(a[i][j0],scale(sub(ONE,rho),2*u[i]))
        for j in range(7):a[i0][j]=add(a[i0][j],scale(sub(ONE,sigma),2*v[j]))
        h=[[ZERO]*15 for _ in range(15)]
        for i,j in product(range(7),repeat=2):
            diff=mul((1,2),ww[i][j])
            for t,q in product(range(2),repeat=2):
                z=add(a[i][j],diff if t==q else neg(diff))
                assert z[0]%2==z[1]%2==0
                h[2*i+t][2*j+q]=(z[0]//2,z[1]//2)
        for i in range(14):h[i][14]=ROOTS[2]
        h[14]=[ONE]*15
        e=[[ROOTS.index(x) for x in row] for row in h]
        g=gram(e)
        assert energy(g)==105 and norm(det(h))==B
        assert all(g[i][j]==((15,0) if i==j else (1,-1) if j==14 else (2,1) if i==14 else (3,0) if i//2==j//2 else ZERO) for i,j in product(range(15),repeat=2))
        records.append({'row_phase':r,'column_phase':s,'exponents':e,'automorphism_order':1008 if r==s==0 else 144 if r*s==0 else 36})
    # Independently enumerate the 168 Fano incidence permutations.
    line_sets=[frozenset(j for j,x in enumerate(row) if x) for row in n]
    actions=[]
    for p in permutations(range(7)):
        image=[frozenset(p[j] for j in line) for line in line_sets]
        if set(image)==set(line_sets):actions.append((p,tuple(line_sets.index(z) for z in image)))
    assert len(actions)==168
    patterns=[(tuple([0]*7),tuple([0]*7))]
    opts=[tuple([0]*7)]+[tuple(t*x for x in row) for row in n for t in (1,2)]
    for aa,bb in product(opts,repeat=2):
        if any(aa) and any(bb):
            point=next(j for j in range(7) if tuple(int(x!=0) for x in aa)==tuple(n[i][j] for i in range(7)))
            line=next(i for i in range(7) if tuple(int(x!=0) for x in bb)==tuple(n[i]))
            if n[line][point]:continue
        patterns.append((aa,bb))
    patterns=set(patterns);assert len(patterns)==141
    def orbit(x):
        aa,bb=x
        return {(tuple(aa[q[i]] for i in range(7)),tuple(bb[p[j]] for j in range(7))) for p,q in actions}
    orbits=[]
    while patterns:
        x=min(patterns);o=orbit(x);assert o<=patterns
        orbits.append({'size':len(o),'stabilizer':168//len(o)})
        patterns-=o
    assert sorted(z['size'] for z in orbits)==[1,7,7,7,7,28,28,28,28]
    save('friendship_nine_factors',records)
    return {'fano_group_order':168,'normalized_factors':141,'orbits':orbits,'classes_per_orientation':9}


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled.')
    OUT.mkdir(exist_ok=True,parents=True)
    result={'bw25':verify_bw25(),'sixth_root':verify_sixth(),'fano_switches':fano_factors(),
            'scope':'Exact explicit witnesses and Fano switch orbits; completeness of arbitrary Gram factorizations is an analytic proof, not inferred from this finite check.'}
    result['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*witness.json')}
    save('continuation_verified',result)
    print('PASS: BW(25,16), order 51 (independent 102x102 determinant), sixth-root witness, all 141 Fano switches and nine orbits.')

if __name__=='__main__':main()
