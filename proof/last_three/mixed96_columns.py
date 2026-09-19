from mixed96_candidates import *
from column_rank import decode

def run():
 data=json.loads(Path(__file__).with_name('mixed96_candidates.json').read_text());out=[];start=time.monotonic()
 for i,row in enumerate(data['rows']):
  G=full_mixed_gram(row);a,b,L=inverse_arrays(G);hits=column_patterns(a,b,L)
  cols=[decode(s)for s in hits]
  W=[[esum(mul(v[i],conj(v[j]))for v in cols)for j in range(15)]for i in range(15)];R,piv=rref(W);rank=len(piv)
  rr={**row,'candidate_index':i,'inverse_scale':L,'column_count':len(hits),'column_serials':list(map(int,hits)),'span_rank':rank,'excluded':rank<15};out.append(rr)
  print(i,row['family'],row['D'],'columns',len(hits),'rank',rank,'seconds',round(time.monotonic()-start,2),flush=True)
  Path(__file__).with_name('mixed96_columns.json').write_text(json.dumps({'total':len(data['rows']),'processed':len(out),'complete':len(out)==len(data['rows']),'rows':out},indent=2)+'\n')
 return out
if __name__=='__main__':run()
