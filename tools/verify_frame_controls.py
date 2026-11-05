#!/usr/bin/env python3
"""Exact positive controls for the analytic general energy gaps."""
from itertools import combinations,product
import json
from pathlib import Path
from matrix_exact import gram,norm
from verify_mixed_inert_gap import gram_energy_and_colors
ROOT=Path(__file__).resolve().parents[1]

def evaluate(e):
 n=len(e);q,rows=gram_energy_and_colors(e);q2,cols=gram_energy_and_colors([list(r) for r in zip(*e)])
 assert q==q2
 if q:
  assert q>=3*(n-1)
  if len(set(rows))==len(set(cols))==1:assert q>=9*((n+1)//2)
  elif n>=6:assert q==3*(n-1) or q>=6*(n-2)
 return q

def main():
 if not __debug__:raise RuntimeError('Assertions must be enabled.')
 e=[list(map(int,s)) for s in ('000000','001122','010221','012012','022101','021210')]
 assert evaluate(e)==0
 single=0
 for pos in range(36):
  i,j=divmod(pos,6)
  for d in (1,2):
   h=[r[:] for r in e];h[i][j]=(h[i][j]+d)%3
   assert evaluate(h)==15
   single+=1
 double=0
 for a,b in combinations(range(36),2):
  for da,db in product((1,2),repeat=2):
   h=[r[:] for r in e]
   for pos,d in ((a,da),(b,db)):
    i,j=divmod(pos,6);h[i][j]=(h[i][j]+d)%3
   evaluate(h);double+=1
 assert double==2520
 pure=0;minimum=10**9
 for x in product(range(3),repeat=4):
  h=[[(a*b+x[2*i+j])%3 for b in range(3) for j in range(2)] for a in range(3) for i in range(2)]
  q=evaluate(h);minimum=min(minimum,q);pure+=1
 assert pure==81 and minimum==27
 (ROOT/'data/revision/frame_controls.json').write_text(json.dumps({'single_entry_repairs':single,'two_entry_controls':double,'pure_tensor_controls':pure,'pure_energy_minimum_at_order6':minimum,'scope':'Finite controls; all-order theorems are proved analytically.'},indent=2)+'\n')
 print('PASS: 72 one-entry repairs, 2520 two-entry changes, 81 pure-color tensor controls.')
if __name__=='__main__':main()
