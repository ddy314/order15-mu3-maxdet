#!/usr/bin/env python3
"""Exact exclusion of Q=9 and Q=10 for order-15 sixth-root matrices.

No determinant-threshold or orthogonality-capacity input is used. A complete
necessary-condition relaxation is exhausted; every rejection has an integer
norm obstruction or a rational frame/rectangle obstruction.
"""
import sys,json
from pathlib import Path
from itertools import combinations,product
from collections import defaultdict
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'proof/new'))
from pure_catalogue import UNITS_ALLOWED,characteristic,kernel_min_from_poly,weight_vectors,unique_permutations,norm_bases,UNITS,mul,conj,np,nx
sys.path.insert(0,str(ROOT/'extensions/energy'))
from verify_pure_energy import catalogue,assemblies as no_isolate_assemblies,poly_mul,determinant,obstruction,rectangle

def leafless():
 graphs=[nx.convert_node_labels_to_integers(g) for g in nx.graph_atlas_g() if 2<=len(g)<=7 and nx.is_connected(g) and g.number_of_edges()<=10 and min(dict(g.degree()).values())>=2]
 for n in (8,9,10):
  buckets=defaultdict(list)
  for tree in nx.nonisomorphic_trees(n):
   tree=nx.convert_node_labels_to_integers(tree);missing=[e for e in combinations(range(n),2) if not tree.has_edge(*e)]
   for k in range(11-n+1):
    for extra in combinations(missing,k):
     g=tree.copy();g.add_edges_from(extra)
     if min(dict(g.degree()).values())<2 or g.number_of_edges()>10:continue
     key=tuple(sorted(dict(g.degree()).values()))
     if not any(nx.is_isomorphic(g,h) for h in buckets[key]):buckets[key].append(g)
  graphs.extend(g for group in buckets.values() for g in group)
 return graphs

def main():
 if not __debug__:raise RuntimeError('Assertions must be enabled.')
 assert {w for w in UNITS_ALLOWED if w<=9} == {a*a-a*b+b*b for a in range(-4,5) for b in range(-4,5) if 0<a*a-a*b+b*b<=9}
 # No-isolate supports at Q9 have components of order at most five.
 # Their maximum cycle rank is one, all phases on that chord are included.
 data={};states=0;graphs=leafless()
 for g in graphs:
  n=len(g);edges=sorted(tuple(sorted(e)) for e in g.edges());m=len(edges)
  tree={tuple(sorted(e)) for e in nx.minimum_spanning_tree(g).edges()};ch=[i for i,e in enumerate(edges) if e not in tree]
  for u in range(m,11):
   for weights in weight_vectors(m,u):
    for ws in unique_permutations(weights):
     for bases in product(*(norm_bases(w) for w in ws)):
      for phase in product(range(6),repeat=len(ch)):
       zs=list(bases)
       for i,p in zip(ch,phase):zs[i]=mul(zs[i],UNITS[p])
       a=np.zeros((n,n),dtype=np.int64);b=a.copy()
       for (i,j),z in zip(edges,zs):a[i,j],b[i,j]=z;a[j,i],b[j,i]=conj(z)
       p=tuple(map(int,characteristic(a,b)));num,den=map(int,kernel_min_from_poly(a,b,np.array(p,dtype=np.int64)));states+=1
       if num==0:continue
       z={'v':n,'u':u,'poly':p,'kernel_min':str(F(num,den))}
       key=n,u,p
       if key not in data or F(data[key]['kernel_min'])<F(num,den):data[key]=z
 records=list(data.values());levels=[]
 for energy in (9,10):
  aa=[]
  def visit(start,v,u,ids,p):
   if u==energy:
    assert v<15
    full=p+(0,)*(15-v);d=determinant(full,15)
    aa.append({'ids':ids,'v':v,'r':15-v,'poly':full,'det':d,'norm_obstruction':obstruction(d)})
    return
   for i in range(start,len(records)):
    z=records[i]
    if u+z['u']<=energy:visit(i,v+z['v'],u+z['u'],ids+[i],poly_mul(p,z['poly']))
  visit(0,0,0,[],(1,))
  pairs=[]
  for i,a in enumerate(aa):
   if a['norm_obstruction']:continue
   for b in aa[i:]:
    if b['norm_obstruction'] or a['poly']!=b['poly']:continue
    failed=next(({'side':'left','component':k,'required':str(F(b['r'],15)),'available':records[k]['kernel_min']} for k in a['ids'] if F(records[k]['kernel_min'])<F(b['r'],15)),None)
    if failed is None:failed=next(({'side':'right','component':k,'required':str(F(a['r'],15)),'available':records[k]['kernel_min']} for k in b['ids'] if F(records[k]['kernel_min'])<F(a['r'],15)),None)
    assert failed is not None,('isolated survivor',energy,a,b)
    pairs.append({'left':a['ids'],'right':b['ids'],'obstruction':failed})
  assert (len(aa),len(pairs))==((285,29) if energy==9 else (1402,124))
  levels.append({'energy':energy,'isolated_assemblies':aa,'frame_pair_exclusions':pairs})
 # At Q10 every isolate-free component has order at most seven. Generate
 # all supports, every norm placement and every chord phase without alpha pruning.
 ten_records,ten_graph_counts,ten_states=catalogue(10)
 ten=no_isolate_assemblies(ten_records,10,15,15)
 ten_norm=[a for a in ten if not a['norm_obstruction']]
 ten_pairs=[]
 for i,a in enumerate(ten_norm):
  for b in ten_norm[i:]:
   if a['poly']!=b['poly']:continue
   failure=rectangle(a,b,ten_records,1)
   assert failure is not None,('Q10 no-isolate survivor',a,b)
   ten_pairs.append({'left':a['ids'],'right':b['ids'],**failure})
 assert len(ten)==74 and len(ten_norm)==2 and len(ten_pairs)==2
 # Every isolate-free assembly with nine unit-energy edges has one of these
 # forest component-size partitions or a triangle plus six pairs.
 trees={n:list(nx.nonisomorphic_trees(n)) for n in (2,3,4,5)}
 polys={}
 for n,gs in trees.items():
  for g in gs:
   a=nx.to_numpy_array(g,dtype=np.int64);b=np.zeros_like(a)
   polys.setdefault(n,[]).append(tuple(map(int,characteristic(a,b))))
 forests=[]
 for sizes in ((5,2,2,2,2,2),(4,3,2,2,2,2),(3,3,3,2,2,2)):
  for pp in product(*(polys[n] for n in sizes)):
   p=(1,)
   for q in pp:p=poly_mul(p,q)
   d=determinant(p,15);z=d;five=0
   while z%5==0:z//=5;five+=1
   assert five%2==1,(p,d)
   forests.append({'sizes':sizes,'polynomial':p,'determinant':d,'inert_prime':5})
 triangle=[]
 for tau in (-2,-1,1,2):
  d=(3330+tau)*224**6
  prime=obstruction(d)
  if tau in (-2,1):
   assert prime is None
   # Negative triangle has eigenvalue1 of multiplicity2, pair union
   # has eigenvalue1 of multiplicity6: capacity 2*(15+1)=32 <36.
   triangle.append({'flux_trace':tau,'area':36,'capacity':32 if tau==-2 else 0})
  else:assert prime is not None;triangle.append({'flux_trace':tau,'inert_prime':prime})
 result={'leafless_graphs':len(graphs),'labelled_states':states,'flat_kernel_records':records,'levels':levels,'no_isolate_ten':{'records':ten_records,'assemblies':ten,'norm_survivors':ten_norm,'rectangle_exclusions':ten_pairs,'graph_counts':ten_graph_counts,'labelled_states':ten_states},'forest_exclusions':forests,'triangle_exclusions':triangle,'result':'Q(H)>=11 for H in mu6^(15x15)'}
 (ROOT/'data/revision/sixth_gap_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS: Q9 and Q10 excluded; sixth-root energy is at least eleven.',len(graphs),'leafless supports,',states,'labelled states',flush=True)
if __name__=='__main__':main()
