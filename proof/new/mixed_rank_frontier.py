"""Flat submatrix rank obstructions for paired (14,1) Grams, isolates retained."""
from bounds import *
RANKS=json.loads(Path(__file__).with_name('rank_catalogue.json').read_text())['types']
TYPES=sorted(RANKS,key=lambda x:(x['u'],x['v'],x['nullity'],x['kernel_support'],x['independence']))

def assemble(maxu):
 def rec(start,left,vertices,selected):
  if selected:yield tuple(selected)
  for i in range(start,len(TYPES)):
   t=TYPES[i]
   if t['u']<=left and vertices+t['v']<=13:
    yield from rec(i,left-t['u'],vertices+t['v'],selected+[i])
 yield from rec(0,maxu,0,[])

def necessary_pair(a,b):
 r,k,t=a['isolates'],a['nullity'],a['kernel_support'];R,K,T=b['isolates'],b['nullity'],b['kernel_support']
 if r+k!=R+K:return False
 if r*(15-R-T)>15 or R*(15-r-t)>15:return False
 if r*(15-R)>15*(K+1) or R*(15-r)>15*(k+1):return False
 return True

def run():
 out=[]
 for q in (87,96,105):
  states=[];count=0
  for inds in assemble((q-42)//9):
   cc=[TYPES[i]for i in inds];v=sum(x['v']for x in cc);u=sum(x['u']for x in cc);r=14-v
   if r+sum(x['independence']for x in cc)>12:continue
   count+=1;k=sum(x['nullity']for x in cc);t=sum(x['kernel_support']for x in cc)
   cfg=((1,0),)*r+tuple((x['v'],x['u'])for x in cc)
   ds=[F(15)]*r+[F(x['max_determinant'])for x in cc];Ls=[F(15)]*r+[cb(x['v'],x['u'])[1]for x in cc]
   corr=sum(F(3*vv)/L for (vv,uu),L in zip(cfg,Ls))+F(q-9*u-42)/max(Ls)
   D=prod(ds)*max(F(0),15-corr)
   if D<=B:continue
   states.append({'isolates':r,'nullity':k,'kernel_support':t,'internal':9*u,'upper':rat(D),'over_B':float(D/B),'component_types':inds,'configuration':cfg})
  before=len(states);rounds=[]
  while states:
   new=[a for a in states if any(necessary_pair(a,b)for b in states)]
   rounds.append([len(states),len(new)])
   if len(new)==len(states):break
   states=new
  row={'Q':q,'scope':'both row and column color type (14,1); strict-counterexample branches with isolates','assemblies_after_code_bound':count,'above_record_before_pairing':before,'pairing_rounds':rounds,'remaining':states,'closed_double_14_1':not states}
  print(q,'before',before,'rounds',rounds,'remaining',len(states),flush=True)
  for x in states:print('  ',x['configuration'],'rkt',x['isolates'],x['nullity'],x['kernel_support'],'over',x['over_B'],flush=True)
  out.append(row)
 Path(__file__).with_name('mixed_rank_frontier.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':run()
