from bounds import *
def capped_stationary(n,q,L):
 L=F(L);q=F(q);values=[]
 for k in range(n+1):
  d=n-k;S=15*n-k*L;SS=225*n+2*q-k*L*L
  if d==0:
   if S==0 and SS==0:values.append((k,0,L**k))
   continue
  mean=S/d;rvar=SS-d*mean*mean
  if mean<=0 or mean>L or rvar<0:continue
  if rvar==0:
   values.append((k,0,L**k*mean**d));continue
  if d==1:continue
  for m in range(1,d):
   t=sqrt_lo(rvar/F(d*m*(d-m)))
   high=mean+(d-m)*t;low=mean-m*t
   if high>L or low<=0:continue
   values.append((k,m,L**k*high**m*low**(d-m)))
 return max((v[2] for v in values),default=F(0)),values
if __name__=='__main__':
 for q,L in ((153,F(154,5)),(153,F(31)),(162,F(63,2)),(144,F(30)),(135,F(147,5))):
  u,vals=capped_stationary(15,q,L)
  print(q,L,float(u/B),'arg',max(vals,key=lambda x:x[2])[:2],'states',len(vals))
