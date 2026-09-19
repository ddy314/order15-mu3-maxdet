from bounds import *
ALL_Q=(87,90,96,99,105,108,114,117,123,126,132,135,141,144,150,153,159,162,168)
out={}
for q in ALL_Q:
 r13=[];mx13=F(0);arg13=None;count13=0
 for e in range(0,q-77,9):
  for h in sorted({z['norm'] for z in CAT if z['norm']%9==(3 if q%9==0 else 0)}):
   c=q-e-h
   if c<78 or c%9!=6:continue
   for cfg in allcfg(13,e//9):
    u=bound13(cfg,c,h);count13+=1
    if u>mx13:mx13=u;arg13=(e,c,h,cfg)
    if u>B:r13.append((e,c,h,cfg,float(u/B)))
 mx14=F(0);arg14=None;r14=[];count14=0
 if q%9==6:
  for e in range(0,q-41,9):
   for cfg in allcfg(14,e//9):
    u=bound14(cfg,q-e);count14+=1
    if u>mx14:mx14=u;arg14=(e,cfg)
    if u>B:r14.append((e,cfg,float(u/B)))
 print(q,'13:',float(mx13/B),len(r13),'14:',float(mx14/B),len(r14),flush=True)
 if q%9==6 and q>=114:print(r14,flush=True)
 out[q]={'max13':rat(mx13),'arg13':arg13,'count13':count13,'res13':r13,'max14':rat(mx14),'arg14':arg14,'count14':count14,'res14':r14}
Path(__file__).with_name('exploration.json').write_text(json.dumps(out,indent=2))
