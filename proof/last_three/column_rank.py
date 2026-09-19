from column_test import *

def decode(s):
 out=[(1,0)]
 for i in range(14):out.append(ROOTS[s%3]);s//=3
 return out

def run():
 data=json.loads(Path(__file__).with_name('column_test.json').read_text());assert data['complete'];out=[]
 for row in data['rows']:
  cols=[decode(s)for s in row['column_serials']]
  W=[[esum(mul(v[i],conj(v[j]))for v in cols)for j in range(15)]for i in range(15)]
  R,piv=rref(W);rank=len(piv)
  out.append({'candidate_index':row['candidate_index'],'Q':row['Q'],'family':row['family'],'D':row['D'],'column_count':len(cols),'span_rank':rank,'excluded':rank<15})
  print(row['candidate_index'],row['Q'],row['family'],len(cols),'rank',rank,flush=True)
 Path(__file__).with_name('column_rank.json').write_text(json.dumps({'rows':out,'excluded':sum(x['excluded']for x in out),'total':len(out)},indent=2)+'\n');return out
if __name__=='__main__':run()
