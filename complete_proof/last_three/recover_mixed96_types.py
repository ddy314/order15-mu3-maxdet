"""Re-enumerate the three component signatures needed by the Q96 frontier."""
from mixed96_candidates import *
from rank_catalogue import phase_signatures,independence
from refine_q114 import make_matrix

def run():
 targets={(3,2):(1,2,2,'P3'),(4,4):(2,4,2,'C4'),(5,6):(3,5,3,'K23')};out=[]
 for (v,u),(null,ks,alpha,name)in targets.items():
  total=0;matches=[]
  for G in nx.graph_atlas_g():
   if len(G)!=v or not nx.is_connected(G)or G.number_of_edges()>u:continue
   G=nx.convert_node_labels_to_integers(G);edges=sorted(tuple(sorted(e))for e in G.edges());tree={tuple(sorted(e))for e in nx.minimum_spanning_tree(G).edges()}
   mask=[int(e not in tree)for e in edges];ch=np.array([i for i,z in enumerate(mask)if z],dtype=np.int64);ed=np.array(edges,dtype=np.int64)
   for ws in weight_vectors(len(edges),u):
    for ordered in unique_permutations(ws):
     for bases in product(*(norm_bases(w)for w in ordered)):
      arr=phase_signatures(v,ed,np.array(bases,dtype=np.int64),ch);total+=len(arr)
      for serial,row in enumerate(arr):
       if (int(row[0]),int(row[1]),independence(G))!=(null,ks,alpha):continue
       assert all(w==1 for w in ordered)
       expected=nx.path_graph(3)if name=='P3'else nx.cycle_graph(4)if name=='C4'else nx.complete_bipartite_graph(2,3)
       assert nx.is_isomorphic(G,expected)
       C=make_matrix(v,edges,bases,mask,serial);d={0:(1,0)};queue=[0]
       for i in queue:
        for j in G.neighbors(i):
         if j not in d:d[j]=mul(conj(scale(C[i][j],F(1,3))),d[i]);queue.append(j)
       assert all(d[i]in UNITS for i in d)
       assert all(C[i][j]==scale(mul(d[i],conj(d[j])),3)for i,j in edges)
       matches.append({'edges':edges,'scaled_core_determinant':int(row[2]),'zero_cycle_flux':True})
  assert matches
  result={'name':name,'v':v,'u':u,'signature':[null,ks,alpha],'tested_phase_states':total,'matching_gauge_states':len(matches),'matches':matches};out.append(result);print(name,'states',total,'matches',len(matches),flush=True)
 Path(__file__).with_name('mixed96_type_recovery.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':run()
