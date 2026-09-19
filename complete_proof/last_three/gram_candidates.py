"""Normalized size-13 Gram representatives, with original mu3 row cosets.
Majority-to-first-hub entries are delta=1-omega. For Q99 the second hub
has the other color: cross entries are -delta*t and g lies in delta*mu3.
For Q96/Q105 both hubs have the same color. All t are in mu3.
"""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'new'))
from cross_exact import *

def core_matrix(family,ei):
 v=3 if family=='P3' else 4
 edges=((0,1),(1,2)) if v==3 else ((0,1),(2,3))
 C=[[(15*int(i==j),0)for j in range(v)]for i in range(v)]
 for (i,j),e in zip(edges,ei):C[i][j]=scale(UNITS[e],3);C[j][i]=conj(C[i][j])
 return C

def full_gram(row):
 family=row['family'];v=3 if family=='P3' else 4;r=13-v;eps=-1 if row['Q']==99 else 1
 G=[[(15*int(i==j),0)for j in range(15)]for i in range(15)]
 C=core_matrix(family,row['edges'])
 for i in range(v):
  for j in range(v):G[i][j]=C[i][j]
 ts=list(row['ratios'])+sum(([k]*n for k,n in enumerate(row['counts'])),[])
 for i in range(13):
  G[i][13]=(1,-1);G[13][i]=conj(G[i][13]);G[i][14]=scale(mul((1,-1),ROOTS[ts[i]]),eps);G[14][i]=conj(G[i][14])
 G[13][14]=tuple(row['g']);G[14][13]=conj(G[13][14]);return G

def run():
 start=time.monotonic();out=[];stats=[]
 old=json.loads((Path(__file__).resolve().parents[1]/'new/size13_energy18.json').read_text())
 allowed={(x['family'],x['Q']):set(x['norm_survivors'])for x in old}
 for family in ('P3','2K2'):
  v=3 if family=='P3' else 4;r=13-v;L=345 if v==3 else 360;dc=3105 if v==3 else 216**2;factor=F(dc*15**r,L**2)
  sums=[((a,b,r-a-b),esum([scale(ROOTS[0],a),scale(ROOTS[1],b),scale(ROOTS[2],r-a-b)]))for a in range(r+1)for b in range(r-a+1)]
  for q in (96,99,105):
   eps=-1 if q==99 else 1;h=q-96
   gs=[mul((1,-1),t)for t in ROOTS]if h==3 else eis_shell(h)
   good=allowed[family,q]
   count=0
   for ei in product(range(6),repeat=2):
    C=core_matrix(family,ei);J=inverse(C);J=[[scale(x,L)for x in row]for row in J]
    assert all(all(F(a).denominator==1 for a in z)for row in J for z in row)
    J=[[tuple(map(int,z))for z in row]for row in J]
    u=[(1,-1)]*v;Ju=[esum(mul(a,b)for a,b in zip(row,u))for row in J]
    xx=esum(mul(conj(a),b)for a,b in zip(u,Ju));assert xx[1]==0
    for tt in product(range(3),repeat=v-1):
     t=(0,)+tt;w=[scale(mul((1,-1),ROOTS[k]),eps)for k in t]
     Jw=[esum(mul(a,b)for a,b in zip(row,w))for row in J]
     yy=esum(mul(conj(a),b)for a,b in zip(w,Jw));zz=esum(mul(conj(a),b)for a,b in zip(u,Jw));assert yy[1]==0
     diag=(15*L-xx[0]-r*(L//5))*(15*L-yy[0]-r*(L//5))
     for ns,s in sums:
      for g in gs:
       off=sub(sub(scale(g,L),zz),scale(s,eps*(L//5)));score=diag-norm(off)
       D=F(score)*factor
       if D<=B or D.denominator!=1 or D.numerator%3**14 or not is_norm(D.numerator):continue
       # Consistency check against the older, more permissive Schur relaxation.
       assert D.numerator in good
       out.append({'Q':q,'family':family,'edges':ei,'ratios':t,'counts':ns,'g':g,'D':int(D)});count+=1
   stats.append({'family':family,'Q':q,'candidate_grams':count});print(stats[-1],time.monotonic()-start,flush=True)
 data={'normalization':'majority cross delta; mu3 row scaling; hub2 negative cross for Q99; t0=1','stats':stats,'rows':out}
 Path(__file__).with_name('gram_candidates.json').write_text(json.dumps(data,indent=2)+'\n');return data
if __name__=='__main__':run()
