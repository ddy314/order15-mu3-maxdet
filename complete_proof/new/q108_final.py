"""Complete Q108 closure from exhaustive pure spectra and exact size13 Schur bounds."""
from bounds import *

def certificate():
 e=json.loads(Path(__file__).with_name('exploration.json').read_text())['108']['res13']
 assert all(x[0]in(0,9,18)for x in e)
 sparse=json.loads(Path(__file__).with_name('sparse_defect.json').read_text())
 assert sparse['closed']
 size=json.loads(Path(__file__).with_name('size13_energy18.json').read_text())
 rows=[r for r in size if r['Q']==108]
 assert {r['family']for r in rows}=={'P3','2K2'}
 # The enumeration reports all possible score values. This closure requires
 # no score above B, not merely the absence of norm-compatible examples.
 assert all(r['above_B_score_count']==0 for r in rows)
 p=json.loads(Path(__file__).with_name('pure108_frontier.json').read_text())
 assert not p['isolate']['survivors'] and not p['no_isolate']['survivors']
 out={'Q':108,'closed':True,'size13_energy18_families':rows,'no_isolate_assemblies':p['no_isolate']['assembly_count'],'no_isolate_norm_survivors':len(p['no_isolate']['norm_above_B']),'isolate_assemblies':p['isolate']['assembly_count'],'isolate_norm_survivors_before_projection':len(p['isolate']['norm_above_B']),'isolate_survivors_after_projection':0,'external_dependency':'B3(15,10)=12 for the zero-energy 13-row block'}
 Path(__file__).with_name('q108_final.json').write_text(json.dumps(out,indent=2)+'\n')
 print('Q108 closed: both size13 energy18 families below B; all pure-Gram assemblies excluded.',flush=True)
 return out
if __name__=='__main__':certificate()
