"""Check the exact energy18 Schur formula against independent 15-by15 Bareiss determinants."""
from size13_energy18 import *
import random

def run():
 rng=random.Random(2026091802);rows=[];delta=(1,-1)
 for name,r in (('P3',10),('2K2',9)):
  core=tree_matrix(3)if name=='P3'else tuple(tuple((15 if i==j else 3 if i//2==j//2 else 0,0)for j in range(4))for i in range(4))
  v=len(core);J=inverse([list(x)for x in core]);da=bareiss([list(x)for x in core])[0]*15**r
  for h in (0,3,9,12):
   for trial in range(2):
    u=[mul(delta,rng.choice(UNITS))for i in range(v)];w=[mul(z,rng.choice(ROOTS))for z in u]
    ui=[delta]*r;wi=[mul(delta,rng.choice(ROOTS))for _ in range(r)];g=rng.choice(eis_shell(h))
    Ju=[esum(mul(a,z)for a,z in zip(row,u))for row in J];Jw=[esum(mul(a,z)for a,z in zip(row,w))for row in J]
    x=esum(mul(conj(a),b)for a,b in zip(u,Ju))[0]+F(r,5)
    y=esum(mul(conj(a),b)for a,b in zip(w,Jw))[0]+F(r,5)
    z=add(esum(mul(conj(a),b)for a,b in zip(u,Jw)),scale(esum(mul(conj(a),b)for a,b in zip(ui,wi)),F(1,15)))
    schur=F(da)*((15-x)*(15-y)-norm(sub(g,z)))
    G=[[(15 if i==j else 0,0)for j in range(15)]for i in range(15)]
    for i in range(v):
     for j in range(v):G[i][j]=core[i][j]
    for i,(a,b)in enumerate(zip(u+ui,w+wi)):
     G[i][13]=a;G[13][i]=conj(a);G[i][14]=b;G[14][i]=conj(b)
    G[13][14]=g;G[14][13]=conj(g)
    exact=bareiss(G);assert exact==(schur,0)
    rows.append({'family':name,'outside_norm':h,'trial':trial,'match':True})
 out={'seed':2026091802,'cases':len(rows),'all_passed':True,'rows':rows}
 Path(__file__).with_name('size13_validation.json').write_text(json.dumps(out,indent=2)+'\n');print('Validated',len(rows),'full15-by15 energy18 determinant formulas',flush=True);return out
if __name__=='__main__':run()
