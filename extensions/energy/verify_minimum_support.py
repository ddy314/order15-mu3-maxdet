#!/usr/bin/env python3
"""All unweighted majority supports with 5--7 edges and independence <=9."""
from pathlib import Path
import json
from itertools import combinations
import networkx as nx
ROOT=Path(__file__).resolve().parents[2]
def alpha(g):
 n=len(g)
 for k in range(n,0,-1):
  if any(all(not g.has_edge(a,b) for a,b in combinations(v,2)) for v in combinations(range(n),k)):return k

def main():
 if not __debug__:raise RuntimeError('Assertions must be enabled.')
 graphs=[nx.convert_node_labels_to_integers(g) for g in nx.graph_atlas_g() if 2<=len(g)<=7 and nx.is_connected(g) and g.number_of_edges()<=7]
 graphs += [nx.convert_node_labels_to_integers(g) for g in nx.nonisomorphic_trees(8)]
 graphs.sort(key=lambda g:(g.number_of_edges(),len(g),nx.to_graph6_bytes(g,header=False)))
 records=[(len(g),g.number_of_edges(),alpha(g),g) for g in graphs]
 out=[]
 def visit(start,v,m,a,ids):
  if 5<=m<=7 and 14-v+a<=9:
   g=nx.disjoint_union_all([records[i][3] for i in ids]);g.add_nodes_from(range(14))
   degrees=dict(g.degree());r=sum(d==0 for d in degrees.values());s=len({next(iter(g[i])) for i,d in degrees.items() if d==1})
   out.append({'edges':m,'isolates':r,'leaf_neighbors':s,'component_ids':ids,'support_edges':[list(e) for e in g.edges()]})
  for i in range(start,len(records)):
   w,k,b,_=records[i]
   if v+w<=14 and m+k<=7:visit(i,v+w,m+k,a+b,ids+[i])
 visit(0,0,0,0,[])
 counts=[sum(r['edges']==m for r in out) for m in (5,6,7)]
 assert counts==[1,4,19],counts
 isolated=[a for a in out if a['isolates']]
 pairs=[(a,b) for a in isolated for b in isolated if (a['leaf_neighbors']+1)*b['isolates']<=15 and (b['leaf_neighbors']+1)*a['isolates']<=15]
 assert len(pairs)==1 and pairs[0][0]['isolates']==pairs[0][1]['isolates']==1
 assert len([a for a in out if not a['isolates']])==1
 (ROOT/'data/revision/minimum_support_certificate.json').write_text(json.dumps({'counts':counts,'supports':out,'isolated_pair_survivors':pairs},indent=2)+'\n')
 print('PASS: 24 majority supports; sole two-sided-isolate survivor P3 + 5K2 + K1.')
if __name__=='__main__':main()
