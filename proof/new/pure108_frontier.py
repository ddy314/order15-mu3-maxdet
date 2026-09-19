from bounds import *
from itertools import product
CATALOG=json.loads(Path(__file__).with_name('pure_extended.json').read_text())['types']
for t in CATALOG:t['polynomial']=t['p']

def pmul(p,q):
 z=[0]*(len(p)+len(q)-1)
 for i,x in enumerate(p):
  for j,y in enumerate(q):z[i+j]+=x*y
 return tuple(z)
def peval(p,x):
 s=0
 for c in p:s=s*x+c
 return s

def run():
 q=108;units=12;result={}
 for iso in (False,True):
  if iso:
   types=[t for t in CATALOG if F(t['min_kernel'])>0]
  else:
   types=[t for t in CATALOG if t['u']+(15-t['v']+1)//2<=units]
  types=sorted(types,key=lambda t:(t['u'],t['v'],t['polynomial']))
  states={};assemblies=0
  def rec(start,v,u,p,mink,inds):
   nonlocal assemblies
   if u==units:
    r=15-v
    if (iso and r<1)or(not iso and r!=0):return
    assemblies+=1;pol=p+(0,)*r;D=-3**15*peval(pol,-5)
    if D<=B or not is_norm(D):return
    key=(r,pol);row={'isolates':r,'polynomial':pol,'determinant':D,'over_B':float(F(D,B)),'min_kernel':str(mink),'types':inds,'components':[(types[i]['v'],types[i]['u'])for i in inds]}
    if key not in states or mink>F(states[key]['min_kernel']):states[key]=row
    return
   for i in range(start,len(types)):
    t=types[i];V=v+t['v'];U=u+t['u']
    if U>units or V>(14 if iso else 15):continue
    if not iso and U+(15-V+1)//2>units:continue
    rec(i,V,U,pmul(p,tuple(t['polynomial'])),min(mink,F(t['min_kernel'])),inds+(i,))
  rec(0,0,0,(1,),F(1),())
  old=list(states.values());surv=[a for a in old if any(a['polynomial']==b['polynomial'] and(15*F(a['min_kernel'])>=b['isolates'])and(15*F(b['min_kernel'])>=a['isolates'])for b in old)]
  print('pure108','isolate'if iso else'noisolate','types',len(types),'assemblies',assemblies,'normaboveB',len(old),'afterprojection',len(surv),flush=True)
  for x in surv:print(x,flush=True)
  result['isolate'if iso else'no_isolate']={'type_count':len(types),'assembly_count':assemblies,'norm_above_B':old,'survivors':surv}
 Path(__file__).with_name('pure108_frontier.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':run()
