"""Exploratory rational bounds, with internal isolated vertices retained."""
import sys,json
from pathlib import Path
from functools import lru_cache
from math import prod,isqrt,comb
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'order15_mu3_unified_audit/programs'))
from exact_core import *
from verify_component_bounds import configs,component_cap
BASE=Path(__file__).resolve().parents[1]/'order15_mu3_unified_audit'
T=json.loads((BASE/'results/component_bounds.json').read_text())['tree_catalogue']
CAT=json.loads((BASE/'results/gram_catalogue.json').read_text())
_SMALL_PATH=Path(__file__).with_name('small_catalogue.json')
SMALL=json.loads(_SMALL_PATH.read_text())['records'] if _SMALL_PATH.exists() else {}
UNITS_ALLOWED=sorted({z['norm']//9 for z in CAT if z['norm']%9==0 and z['norm']>0})
@lru_cache(None)
def wparts(u,minimum=1):
 if u==0:return ((),)
 return tuple((w,)+rest for w in UNITS_ALLOWED if minimum<=w<=u for rest in wparts(u-w,w))
def expneg_upper(a,terms=24):
 assert a>=0
 term=F(1);s=F(1)
 for k in range(1,terms+1):term*=a/k;s+=term
 return 1/s
@lru_cache(None)
def cycle_triangle_cap(g):
 if g==0:return 0
 # Triangles and cycle rank add over biconnected blocks. In a cyclic
 # block of rank g, min degree d satisfies d(d-1)/2 <= g. Deleting a
 # min-degree vertex leaves a connected graph of rank g-d+1 and
 # removes at most d(d-1)/2 triangles.
 vals=[cycle_triangle_cap(a)+cycle_triangle_cap(g-a) for a in range(1,g)]
 for d in range(2,(1+isqrt(1+8*g))//2+1):
  vals.append(d*(d-1)//2+cycle_triangle_cap(g-d+1))
 return max(vals)
def no_leaf_triangle_cap(v,m):
 g=m-v+1
 if g<1:return 0
 split=max((cycle_triangle_cap(a)+cycle_triangle_cap(g-a) for a in range(1,g)),default=0)
 dmax=min(2*m//v,(1+isqrt(1+8*g))//2)
 block=max((d*(d-1)//2+cycle_triangle_cap(g-d+1) for d in range(2,dmax+1)),default=0)
 return max(split,block)
@lru_cache(None)
def triangle_weight_cap(v,u,no_leaf=False):
 best=F(0)
 for ws in wparts(u):
  m=len(ws)
  if not (v if no_leaf else v-1)<=m<=v*(v-1)//2:continue
  g=m-v+1
  if not g:continue
  tmax=min(cycle_triangle_cap(g),2**(g-1),g*(g+1)//2,comb(v,3),m*(m-1)//6)
  if no_leaf:tmax=min(tmax,no_leaf_triangle_cap(v,m))
  deltas=[sqrt_hi(F(w))-1 for w in ws if w>1]
  from itertools import combinations
  expanded=F(tmax)+g*sum(deltas)+sum(prod(z) for z in combinations(deltas,2))+sum(prod(z) for z in combinations(deltas,3))
  best=max(best,min(F(tmax)*sqrt_hi(F(prod(ws[-3:]))),expanded))
 return best
@lru_cache(None)
def cb(v,u,cubic=True):
 if v==1:return F(15),F(15)
 d=F(T[str(v)]['max_determinant']) if u==v-1 and v<=14 else stationary(v,9*u)
 k=min(v,(1+isqrt(1+8*u))//2)
 L=min(component_cap(v,u),15+sqrt_hi(F(18*u*(k-1),k)))
 if not any(v-1<=len(ws)<=v*(v-1)//2 for ws in wparts(u)):return F(0),L
 if cubic:
  a=F(u,25)-F(2,125)*triangle_weight_cap(v,u)+F(3*u*u,125*v)/L
  if a>0:d=min(d,F(15**v)*expneg_upper(a))
 if f'{v},{u}' in SMALL:d=min(d,F(SMALL[f'{v},{u}']['max_determinant']))
 return d,L
@lru_cache(None)
def allcfg(n,units):
 return tuple(((1,0),)*r+c for r in range(n+1) for c in configs(n-r,units))
def bound14(cfg,cross,cubic=True):
 ds,Ls=zip(*(cb(v,u,cubic) for v,u in cfg))
 s=sum(F(3*v)/L for (v,u),L in zip(cfg,Ls))+F(cross-42)/max(Ls)
 upper=prod(ds)*max(F(0),15-s);r=cfg.count((1,0))
 if r>1:upper=min(upper,15**(r-1)*stationary(16-r,cross+9*sum(u for v,u in cfg)))
 return upper

def bound13(cfg,cross,h,cubic=True):
 ds,Ls=zip(*(cb(v,u,cubic) for v,u in cfg))
 s=sum(F(6*v)/L for (v,u),L in zip(cfg,Ls))+F(cross-78)/max(Ls)
 if s>=30:return F(0)
 f=(15-s/2)**2-max(F(0),sqrt_lo(F(h))-s/2)**2
 upper=prod(ds)*max(F(0),f);r=cfg.count((1,0))
 if r>2:upper=min(upper,15**(r-2)*stationary(17-r,cross+h+9*sum(u for v,u in cfg)))
 return upper
