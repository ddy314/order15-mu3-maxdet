"""Fail-closed global coverage ledger. A local PASS never closes a whole Q."""
from pathlib import Path
import json,csv
from exact_core import B,stationary,rat
ROOT=Path(__file__).resolve().parents[1]
def run():
 low=json.loads((ROOT/'results/low_energy.json').read_text());mixed=json.loads((ROOT/'results/mixed_small.json').read_text());q81=json.loads((ROOT/'results/q81.json').read_text());comp=json.loads((ROOT/'results/component_bounds.json').read_text());q99=json.loads((ROOT/'results/q99_independent.json').read_text());found=json.loads((ROOT/'results/foundations.json').read_text())
 lowmap={x['Q']:x for x in low['rows']};mixmap={x['Q']:x for x in mixed['enumerated_shells']}
 rows=[]
 for q in range(169):
  if q%9 not in (0,6):continue
  status='UNVERIFIED_GLOBAL';basis='';missing='Complete record-level proof for every color, support, and phase case.'
  if q in lowmap and lowmap[q]['status']=='VERIFIED_CLOSED':status='VERIFIED_CLOSED';basis='Exact no-leaf catalogue, norm sieve and two-sided kernel projection.';missing=''
  elif q in (6,15,24,33):status='VERIFIED_CLOSED';basis='Color reduction and minimum cross energy 42.';missing=''
  elif q in (42,51):status='VERIFIED_CLOSED_WITH_PUBLISHED_THEOREM';basis='Color/energy bound plus B3(15,10)=12.';missing=''
  elif q in mixmap and mixmap[q]['status']=='VERIFIED_WITH_CITED_CODE_BOUND':status='VERIFIED_CLOSED_WITH_PUBLISHED_THEOREM';basis='Full small mixed support/cross enumeration and eigenprojection PSD; B3(15,10)=12.';missing=''
  elif q==81 and q81['status']=='VERIFIED_WITH_CITED_CODE_BOUND':status='VERIFIED_CLOSED_WITH_PUBLISHED_THEOREM';basis='Isolate support catalogue and exact K3+6K2 block obstruction; B3(15,10)=12.';missing=''
  partial=[]
  if q==99:
   partial.append('All AB/AA/BB local mu3 intertwiners independently recovered; 144/225/729 solutions, zero compatible residuals.')
   missing='Exhaustive reduction to A/B and mixed-color branches, especially size13 internal energy9/18.'
  cr=[r for r in comp['component_rows'] if r['Q']==q]
  if cr and all(r['status']=='VERIFIED_NO_ISOLATE_BRANCH' for r in cr):partial.append('(14,1) no-internal-isolate branch completely bounded by B over all component allocations.')
  if q==105:partial.append('An explicit matrix attains B; the 7K2 friendship support family has exact upper bound B.')
  if q in (87,90,96,108,114,117,123,132,141):missing='Record-level size13 defect/local inverse reductions and other remaining support branches; a lower historical comparison target is insufficient.'
  if q in (126,135):missing='Complete polynomial-majorant/high-triangle orbits, Fourier branches and mixed-color record-level reductions.'
  if q==144:missing='Complete K6 support catalogue and record-level state-pairing certificates.'
  if q==150:missing='Record-level (13,2) local Schur/state-pairing and (14,1) isolate exclusion; the no-isolate branch has now been reconstructed.'
  if q in (153,159,162,168):missing='Full upstream sparse/spectral/support reductions of the tail certificate, including all-same/mixed residuals. The historical closure claim has not been fully replayed here.'
  rows.append({'Q':q,'status':status,'legacy_ledger_claim':'CLOSED','complete_proof_basis':basis,'verified_partial_results':partial,'missing_obligation':missing,'generic_trace_upper':rat(stationary(15,q)),'generic_trace_upper_over_B':float(stationary(15,q)/B)})
 assert len(rows)==38 and [r['Q'] for r in rows]==[q for q in range(169) if q%9 in (0,6)]
 closed=[r['Q'] for r in rows if r['status'].startswith('VERIFIED_CLOSED')];pending=[r['Q'] for r in rows if r['status']=='UNVERIFIED_GLOBAL']
 result={'benchmark_B':B,'global_maximality_verified':not pending,'admissible_shells_checked_for_coverage':38,'closed_shells':closed,'unverified_shells':pending,'Q_congruence_excluded_residues':[1,2,3,4,5,7,8],'all_Q_at_least_171':'VERIFIED_EXCLUDED_BY_TRACE_BOUND','scope_warning':'UNVERIFIED does not mean a counterexample exists; it means no complete replayable proof has been supplied in this audit.','rows':rows}
 (ROOT/'results/all_Q_audit.json').write_text(json.dumps(result,indent=2)+'\n')
 with (ROOT/'results/all_Q_audit.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['Q','status','legacy_claim','complete_proof_basis','verified_partial_results','missing_obligation','trace_upper_over_B'])
  for r in rows:w.writerow([r['Q'],r['status'],r['legacy_ledger_claim'],r['complete_proof_basis'],' | '.join(r['verified_partial_results']),r['missing_obligation'],r['generic_trace_upper_over_B']])
 print('Global coverage:',len(closed),'closed;',len(pending),'unverified. Global maximality:',result['global_maximality_verified']);return result
if __name__=='__main__':run()
