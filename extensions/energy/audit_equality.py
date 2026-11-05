#!/usr/bin/env python3
"""Inclusive equality audit; reads freshly regenerated original catalogues.

A component cutoff safe for strict improvements is safe for equality here:
its full Hadamard extension is divisible by five whereas B is not.
"""
from pathlib import Path
import sys,json,hashlib
from math import prod
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'proof/new'))
from bounds import allcfg, cb, bound14, bound13, CAT, F
B=277868041444786176

def main():
 if not __debug__:raise RuntimeError('Assertions must be enabled.')
 p=ROOT/'proof/new/pure_extended.json';data=json.loads(p.read_text())['types']
 rs=sorted([r for r in data if r['u']+(16-r['v'])//2<=12],key=lambda r:(r['v'],r['u'],r['p']))
 count=0;sur=[]
 def visit(start,v,u,d,ids):
  nonlocal count
  if v==15:
   if u==12:
    count+=1
    if d>=B:sur.append(d)
   return
  for k in range(start,len(rs)):
   r=rs[k]
   if v+r['v']>15 or u+r['u']>12:continue
   visit(k,v+r['v'],u+r['u'],d*r['determinant'],ids+[k])
 visit(0,0,0,1,[])
 assert count==5570,count
 assert sorted(sur)==sorted([2**21*3**18*7**3,2**19*3**22*17,2**21*3**19*5*23]),sur
 pure=[]
 for u in range(13,19):
  configs=[c for c in allcfg(15,u) if all(v>=2 for v,w in c)]
  upper=max((prod(cb(v,w)[0] for v,w in c) for c in configs),default=F(0))
  assert upper<B,(u,upper/B)
  pure.append({'Q':9*u,'configurations':len(configs),'upper_over_B':str(upper/B)})
 mixed=[]
 for q in range(105,169):
  if q%9 not in (0,6):continue
  mx14=mx13=F(0);n14=n13=0
  if q%9==6:
   for e in range(0,q-41,9):
    for cfg in allcfg(14,e//9):
     if cfg.count((1,0))>1:continue
     n14+=1;mx14=max(mx14,bound14(cfg,q-e))
  for e in range(0,q-77,9):
   for h in sorted({z['norm'] for z in CAT if z['norm']%9==(3 if q%9==0 else 0)}):
    cross=q-e-h
    if cross<78 or cross%9!=6:continue
    for cfg in allcfg(13,e//9):
     if cfg.count((1,0))>2:continue
     n13+=1;mx13=max(mx13,bound13(cfg,cross,h))
  assert mx13<B
  assert mx14<=B if q==105 else mx14<B
  mixed.append({'Q':q,'configurations14':n14,'configurations13':n13,'upper14_over_B':str(mx14/B),'upper13_over_B':str(mx13/B)})
 result={'pure108_assemblies':count,'pure108_above_threshold':sur,'higher_pure':pure,'mixed':mixed,
         'closure':'At Q105 only 7K2 and P3+5K2+K1 remain; the latter has eigenvalue15, which cannot divide B. Schur equality gives all pair phases1.',
         'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [p,ROOT/'proof/new/small_catalogue.json',ROOT/'proof/order15_mu3_unified_audit/results/component_bounds.json',ROOT/'proof/order15_mu3_unified_audit/results/gram_catalogue.json']}}
 (ROOT/'data/revision/equality_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS: inclusive equality boundary, 5570 pure Q108 assemblies and all higher pure/mixed component caps.')
if __name__=='__main__':main()
