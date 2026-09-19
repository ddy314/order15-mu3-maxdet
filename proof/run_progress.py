#!/usr/bin/env python3
"""Reproduce the new proof certificates and assemble the honest global ledger.

All catalogues are regenerated from source. Exit 2 means that genuine global
proof obligations remain; --audit-only asks only whether the documented
partial proof and its checks completed.
"""
from __future__ import annotations
import argparse,csv,json,os,platform,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
NEW=ROOT/'new';RESULTS=ROOT/'results'
CATALOGUES=['small_catalogue.py','rank_catalogue.py','pure_catalogue.py','pure_extended.py']
DERIVED=['explore.py','explore_pure.py','projection_frontier.py','spectral_integer.py',
         'refine_q114.py','cross_exact.py','finalize_high.py','sparse_defect.py',
         'mixed_rank_frontier.py','pure90_frontier.py','q90_final.py',
         'size13_energy18.py','validate_size13.py','pure108_frontier.py','q108_final.py']

def read(name:str):
 return json.loads((NEW/name).read_text())

def make_ledger()->dict:
 baseline=json.loads((ROOT/'order15_mu3_unified_audit/results/all_Q_audit.json').read_text())
 baseclosed=set(baseline['closed_shells']);expected_old={0,6,9,15,18,24,27,33,36,42,45,51,54,60,63,69,72,78,81}
 assert baseclosed==expected_old
 exp=read('exploration.json');pure={r['Q']:r for r in read('projection_frontier.json')}
 methods={};evidence={}
 def accept(q,method,files):methods[q]=method;evidence[q]=files
 for q in (123,132,141,150,159,168):
  assert not exp[str(q)]['res13'] and not exp[str(q)]['res14']
  accept(q,'Complete component/Schur bound, including all internal isolates.',['exploration.json','small_catalogue.json'])
 for q in (135,144,153,162):
  assert not exp[str(q)]['res13'] and pure[q]['closed_all_same'] and not pure[q]['remaining']
  accept(q,'Complete mixed-color bound and two-sided projection/rank/moment bound.',['exploration.json','projection_frontier.json'])
 si=read('q117_integer_sieve.json')
 assert not exp['117']['res13']
 assert all(x['closed'] for x in si.values())
 for state in pure[117]['remaining']:
  cs=[x for x in state['configuration'] if x[0]>1]
  assert len(cs)==1 and cs[0][1]==13
  v=cs[0][0];constraint=state['component_constraints'][0]
  assert (v==9 and v-constraint[2]<=5)or(v in (7,8)and v-constraint[2]<=4)
 accept(117,'Projection reduction and exact integer real-rooted polynomial sieve.',['projection_frontier.json','q117_integer_sieve.json'])
 assert not exp['126']['res13'] and read('q126_final.json')['closed']
 assert len(pure[126]['remaining'])==1
 accept(126,'Rank-four determinant lattice and exact Eisenstein norm sieve.',['projection_frontier.json','q126_final.json'])
 assert not exp['114']['res13'] and read('q114_final.json')['closed']
 accept(114,'Exhaustive core/cross-vector enumeration, norm sieve, and 12-eigenspace flatness contradiction.',['q114_refinement.json','q114_cross_exact.json','q114_final.json'])
 defect=read('sparse_defect.json');assert defect['closed']
 mixed={x['Q']:x for x in read('mixed_rank_frontier.json')}
 for q in (87,90):
  assert all(x[0] in (0,9)for x in exp[str(q)]['res13'])
 assert mixed[87]['closed_double_14_1'] and not mixed[87]['remaining']
 accept(87,'One-edge defect lemma, published B3 bound, and exhaustive double-(14,1) projection/rank pairing.',['sparse_defect.json','mixed_rank_frontier.json'])
 assert read('q90_final.json')['closed']
 accept(90,'One-edge defect lemma and exact pure-Gram catalogue with spectral rectangle obstruction.',['sparse_defect.json','pure_catalogue.json','pure90_frontier.json','q90_final.json'])
 assert all(x[0]in(0,9,18)for x in exp['108']['res13'])
 assert read('q108_final.json')['closed']
 accept(108,'Complete extended pure-Gram catalogue and minimum-cross two-hub energy18 enumeration.',['pure_extended.json','pure108_frontier.json','size13_energy18.json','q108_final.json'])
 new=set(methods);assert new.isdisjoint(baseclosed)
 allowed={q for q in range(169)if q%9 in (0,6)};closed=baseclosed|new;remaining=sorted(allowed-closed)
 rows=[]
 for q in sorted(allowed):
  row={'Q':q,'status':'UNRESOLVED_GLOBAL','proof_basis':'','evidence_files':'','remaining_obligation':''}
  if q in baseclosed:
   row.update(status='INHERITED_VERIFIED',proof_basis='Unchanged 2026-09-17 independent audit; published code theorem dependency retained.',evidence_files='order15_mu3_unified_audit/results/all_Q_audit.json')
  elif q in new:
   row.update(status='NEWLY_CLOSED',proof_basis=methods[q],evidence_files=';'.join('new/'+p for p in evidence[q]))
  else:
   row['remaining_obligation']={96:'Size13 energy18 (P3 or 2K2), and any surviving mixed-color pairings; double14 partial certificates are not a full shell proof.',
    99:'Size13 energy18, plus complete all-same support reduction; the old A/B local certificate does not by itself prove exhaustive coverage.',
    105:'Size13 energy18 and color-pair coverage; equality at B and the complete double14 branch do not settle the remaining colors.',
    108:'Complete all-same support/phase coverage; size13 energy18 is separately investigated in the optional frontier supplement.'}[q]
  rows.append(row)
 result={'benchmark_squared_determinant':277868041444786176,'total_shells':len(allowed),'baseline_closed':sorted(baseclosed),'newly_closed':sorted(new),'closed_shells':sorted(closed),'remaining_shells':remaining,'closed_count':len(closed),'remaining_count':len(remaining),'global_maximality_verified':not remaining,'proof_interpretation':'Computer-assisted mathematical proof with explicitly stated manual lemmas and inherited B3(15,10)=12 theorem; not a formal proof-assistant derivation.','rows':rows}
 RESULTS.mkdir(exist_ok=True)
 (RESULTS/'progress_ledger.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
 # Derive the remaining candidate values from the executed exact enumeration,
 # rather than copying an informal progress claim.
 energy18=read('size13_energy18.json')
 selected=[{k:r[k] for k in ('family','outside_norm','internal_energy','cross_energy','Q','norm_survivors','closed_relaxation')} for r in energy18 if r['Q'] in remaining]
 q105_values=sorted({d for r in energy18 if r['Q']==105 for d in r['norm_survivors']})
 assert mixed[105]['closed_double_14_1'] and not mixed[105]['remaining']
 assert all(x[0] in (0,9,18) for x in exp['105']['res13'])
 assert all(tuple(x[:3])==(18,78,9) for x in exp['105']['res13'] if x[0]==18)
 frontier={'complete_shells_remaining':remaining,'size13_exact_relaxation':selected,
  'Q105_conditional_determinant_squared_candidates':q105_values,
  'Q105_statement':'Every strict Q105 counterexample has squared determinant in this candidate set. The P3 branch is still open; no realization is claimed.',
  'process_correction':'The earlier informal claim of a unique Q105 determinant was too strong. The exact replay retains three different values.',
  'global_proof_completed':False}
 (RESULTS/'remaining_frontier.json').write_text(json.dumps(frontier,ensure_ascii=False,indent=2)+'\n')
 with (RESULTS/'all_Q_progress.csv').open('w',newline='',encoding='utf-8-sig')as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 return result

def main()->int:
 if not __debug__:raise RuntimeError('Assertions must be enabled. Do not use python -O.')
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit-only',action='store_true');p.add_argument('--include-baseline',action='store_true');args=p.parse_args()
 logs=RESULTS/'replay_logs';logs.mkdir(exist_ok=True,parents=True);steps=[];start=time.monotonic()
 commands=[]
 if args.include_baseline:commands.append(('baseline',[sys.executable,str(ROOT/'order15_mu3_unified_audit/run_all.py'),'--audit-only']))
 if not args.include_baseline:
  # These immutable-source computations regenerate every baseline numeric
  # catalogue actually consumed by the new proofs.
  for name in ('verify_foundations.py','verify_component_bounds.py'):
   commands.append((name,[sys.executable,str(ROOT/'order15_mu3_unified_audit/programs'/name)]))
 commands.append(('validate_primitives.py',[sys.executable,str(NEW/'validate_primitives.py')]))
 commands.extend((name,[sys.executable,str(NEW/name)])for name in CATALOGUES)
 commands.extend((name,[sys.executable,str(NEW/name)])for name in DERIVED)
 for name,command in commands:
  t=time.monotonic();r=subprocess.run(command,cwd=ROOT,text=True,capture_output=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
  (logs/(name+'.stdout.txt')).write_text(r.stdout);(logs/(name+'.stderr.txt')).write_text(r.stderr)
  steps.append({'program':name,'returncode':r.returncode,'seconds':time.monotonic()-t})
  print(name+(': PASS'if r.returncode==0 else': FAILED'),flush=True)
  if r.returncode:
   print(r.stderr,file=sys.stderr);return 1
 ledger=make_ledger()
 summary={'python':sys.version,'platform':platform.platform(),'baseline_replayed_in_this_run':args.include_baseline,'steps':steps,'all_executed_steps_passed':True,'seconds':time.monotonic()-start,'newly_closed':ledger['newly_closed'],'remaining':ledger['remaining_shells'],'global_maximality_verified':ledger['global_maximality_verified']}
 (RESULTS/'replay_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(f"NEWLY CLOSED: {ledger['newly_closed']}\nTOTAL CLOSED: {ledger['closed_count']}/38\nUNRESOLVED: {ledger['remaining_shells']}",flush=True)
 return 0 if args.audit_only or ledger['global_maximality_verified'] else 2
if __name__=='__main__':raise SystemExit(main())