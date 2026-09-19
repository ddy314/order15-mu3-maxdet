from bounds import *
from cap_exact import capped_stationary
rows=[]
for u in range(10,19):
 q=9*u;mx=F(0);res=[];count=0
 for cfg in allcfg(15,u):
  count+=1;d=prod(cb(v,e)[0]for v,e in cfg)
  r=cfg.count((1,0))
  if r:d=min(d,15**r*stationary(15-r,q))
  if d>mx:mx=d;arg=cfg
  if d>B:res.append((cfg,float(d/B)))
 print(q, 'max',float(mx/B),'res',len(res),'/',count,flush=True)
 if q>=135:print(res,flush=True)
 rows.append({'Q':q,'count':count,'max':rat(mx),'max_config':arg,'residuals':res})
Path(__file__).with_name('pure_exploration.json').write_text(json.dumps(rows,indent=2)+'\n')
