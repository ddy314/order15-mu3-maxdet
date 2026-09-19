"""All all-same-color Q99 assemblies; retain all component decompositions."""
from pure99_frontier import *
import sympy as sp
x=sp.Symbol('x')

def run():
 units=11;types=sorted([t for t in CATALOG if t['u']+(15-t['v']+1)//2<=units],key=lambda t:(t['u'],t['v'],t['polynomial']))
 out=[];counts=0
 def rec(start,v,u,p,inds):
  nonlocal counts
  if u==units:
   if v!=15:return
   counts+=1;D=-3**15*peval(p,-5)
   if D<=B or not is_norm(D):return
   out.append({'polynomial':p,'D':D,'types':inds,'components':[{'v':types[i]['v'],'u':types[i]['u'],'p':types[i]['p']}for i in inds]});return
  for i in range(start,len(types)):
   t=types[i];V=v+t['v'];U=u+t['u']
   if V>15 or U>units or U+(15-V+1)//2>units:continue
   rec(i,V,U,pmul(p,tuple(t['p'])),inds+(i,))
 rec(0,0,0,(1,),())
 def gcd_poly(pa,pb):return sp.gcd(sp.Poly.from_list(pa,x),sp.Poly.from_list(pb,x))
 pairs=[];proofs=[]
 for ia,a in enumerate(out):
  for ib,b in enumerate(out):
   if ib<ia or a['polynomial']!=b['polynomial']:continue
   fail=None
   for ca in a['components']:
    for cb0 in b['components']:
     p=gcd_poly(ca['p'],cb0['p']);rank=p.degree();va=ca['v'];vb=cb0['v']
     La=cb(va,ca['u'])[1];Lb=cb(vb,cb0['u'])[1]
     bound=rank*min(La,Lb)
     if va*vb>bound:fail={'left':ca,'right':cb0,'common':list(map(int,p.all_coeffs())),'rectangle_size':[va,vb],'rank_bound':rank,'spectral_upper':str(min(La,Lb)),'squared_frobenius':va*vb};break
    if fail:break
   if fail:proofs.append({'left_id':ia,'right_id':ib,'obstruction':fail})
   else:pairs.append([ia,ib])
 print('assemblies',counts,'normabove',len(out),'pairs surviving',pairs,flush=True)
 for i in set(j for p in pairs for j in p):
  print(i,out[i]['D'],[(c['v'],c['u'],str(sp.factor(sp.Poly.from_list(c['p'],x).as_expr())))for c in out[i]['components']],flush=True)
 data={'assembly_count':counts,'candidates':out,'spectral_obstructions':proofs,'surviving_pairs':pairs}
 Path(__file__).with_name('pure99_pairs.json').write_text(json.dumps(data,indent=2)+'\n');return data
if __name__=='__main__':run()
