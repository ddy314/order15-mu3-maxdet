"""Two-sided isolate-projector dimension bounds on all-same shells."""
from bounds import *
from cap_exact import capped_stationary

def clique_no_leaf(v,u):
 ks=[]
 for k in range(2,v+1):
  needed=k*(k-1)//2+(v-k+1 if v>k else 0)
  if needed<=u:ks.append(k)
 return max(ks,default=1)

@lru_cache(None)
def component_projection_cap(v,u,rmin):
 if v==1:return F(15),0,1
 if not any(v<=len(ws)<=v*(v-1)//2 for ws in wparts(u)):return F(0),0,1
 # A unit 3-cycle has nonzero determinant since Re(mu6) never vanishes.
 if (v,u)==(3,3):return F(0),0,3
 nullity=(v*rmin+14)//15
 if nullity>v-2:return F(0),nullity,1
 clique=clique_no_leaf(v,u)
 L=15+sqrt_hi(F(18*u*(clique-1),clique))
 reduced,_=capped_stationary(v-nullity,9*u,L)
 d=v-nullity
 a=F(u,25)-F(2,125)*triangle_weight_cap(v,u,True)+F(3*u*u,125*d)/L
 moment=15**v*expneg_upper(a) if a>0 else F(15**v)
 return min(cb(v,u)[0],15**nullity*reduced,moment),nullity,clique

def run():
 original=json.loads(Path(__file__).with_name('pure_exploration.json').read_text());out=[]
 for row in original:
  q=row['Q']
  if q<117:continue
  states=[tuple(map(tuple,cfg))for cfg,x in row['residuals']];rounds=[]
  while states:
   rmin=min(cfg.count((1,0))for cfg in states)
   assert rmin>0
   new=[];bounds=[]
   for cfg in states:
    D=prod(component_projection_cap(v,u,rmin)[0]for v,u in cfg)
    if D>B:new.append(cfg)
    bounds.append((cfg,D))
   rounds.append({'input_count':len(states),'minimum_isolates':rmin,'output_count':len(new),'max_over_B':float(max(D for cfg,D in bounds)/B)})
   if len(new)==len(states):break
   states=new
  residual=[]
  for cfg in states:
   cc=[(v,u,*component_projection_cap(v,u,rmin)[1:])for v,u in cfg if v>1]
   D=prod(component_projection_cap(v,u,rmin)[0]for v,u in cfg)
   residual.append({'configuration':cfg,'upper':rat(D),'over_B':float(D/B),'component_constraints':cc})
  print(q,rounds,'residual',residual,flush=True)
  out.append({'Q':q,'rounds':rounds,'remaining':residual,'closed_all_same':len(states)==0})
 Path(__file__).with_name('projection_frontier.json').write_text(json.dumps(out,indent=2)+'\n')
 return out
if __name__=='__main__':run()
