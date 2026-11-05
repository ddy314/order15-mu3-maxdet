#!/usr/bin/env python3
"""Exhaust all 675 normalized support-potential patterns independently."""
from itertools import product
from math import comb
from pathlib import Path
import json
from verify_continuation import ROOTS, ONE, ZERO, ROOT, add, mul, conj, scale, total, gf_mul, gf_trace

def main():
 if not __debug__:raise RuntimeError('Assertions must be enabled.')
 n=[[int(gf_trace(gf_mul(a,b))==0) for b in range(1,8)] for a in range(1,8)]
 options=[(0,)*7]+[tuple(t*x for x in row) for row in n for t in (1,2)]
 count=0;passing=[];tested=0
 for gamma,alpha,beta in product(range(3),options,options):
  tested+=1;choices=[]
  for i in range(7):
   positions=[j for j in range(7) if n[i][j]];row=[scale(ROOTS[(gamma+alpha[i]+beta[j])%3],-1) if not n[i][j] else ZERO for j in range(7)]
   rr=[]
   for e in product(range(3),repeat=3):
    r=row[:]
    for j,k in zip(positions,e):r[j]=scale(ROOTS[k],2)
    if total(r)==(2,0):rr.append(r)
   choices.append(rr)
  completions=[]
  def visit(rows):
   if len(rows)==7:completions.append(rows[:]);return
   for row in choices[len(rows)]:
    if all(total(mul(a,conj(b)) for a,b in zip(row,r))==(-2,0) for r in rows):visit(rows+[row])
  visit([])
  if completions:
   assert gamma==0 and len(completions)==1
   passing.append({'row':alpha,'column':beta})
   count+=len(completions)
 assert tested==675 and count==141,(tested,count)
 # Independently enumerate the 792 unordered entry-pair multiplicities.
 pair_types=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
 def compositions(n,k):
  if k==1:yield (n,);return
  for a in range(n+1):
   for tail in compositions(n-a,k-1):yield (a,)+tail
 tables={7:set(),11:set()};tested_pairs=0
 for counts in compositions(7,6):
  tested_pairs+=1;cc=[0,0,0];d=0
  for number,(a,b) in zip(counts,pair_types):
   cc[a]+=number;cc[b]+=number;d+=number*(a!=b)
  z=(cc[0]-cc[2]-6,cc[1]-cc[2]+6)
  if z[0]**2-z[0]*z[1]+z[1]**2!=19*(16-3*d):continue
  for r in (7,11):
   if (z[0]+r*z[1])%19==0:tables[r].add(tuple(cc)+(d,))
 assert tested_pairs==792
 assert tables[7]=={(3,2,9,3),(0,0,14,0),(1,7,6,1),(6,4,4,4),(4,9,1,1)}
 assert tables[11]=={(3,3,8,3),(1,0,13,1),(0,8,6,0),(5,6,3,3),(8,1,5,5),(10,4,0,4)}
 # Exhaustive quadrangle control, with columns indexed by K4 intersections.
 neighbors=[[j for j in range(4) if j!=i] for i in range(4)]
 quad=[]
 for p in range(3):
  triples=[e for e in product(range(3),repeat=3) if sum(e)%3==p]
  number=0
  for rows in product(triples,repeat=4):
   z=[{j:ROOTS[e] for j,e in zip(neighbors[i],rows[i])} for i in range(4)]
   ok=True
   for i in range(4):
    for j in range(i):
     left=total(z[i][k] for k in range(4) if k!=i and k!=j)
     left=add(left,total(conj(z[j][k]) for k in range(4) if k!=i and k!=j))
     right=add((2,0),scale(mul(z[i][j],conj(z[j][i])),2))
     if left!=right:ok=False;break
    if not ok:break
   if ok:
    assert p==0 and all(len(set(row))==1 for row in rows)
    assert sum(row[0]!=0 for row in rows)<=1
    number+=1
  quad.append(number)
 assert quad==[9,0,0],quad
 (ROOT/'data/revision/fano_completion_certificate.json').write_text(json.dumps({'patterns':675,'completions':141,'column_pair_multiplicities':792,'column_types':{str(r):sorted(t) for r,t in tables.items()},'quadrangle_counts':quad,'passing_potentials':passing},indent=2)+'\n')
 print('PASS: all 675 patterns, exactly 141 unique completions; exhaustive quadrangle rigidity.')
if __name__=='__main__':main()
