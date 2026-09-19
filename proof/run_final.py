#!/usr/bin/env python3
"""Rebuild the complete computer-assisted 38/38 proof; fails on any missing step.

Rerun the original audit, all intermediate catalogues, and all final
candidates with two independent exhaustive column algorithms.
"""
from pathlib import Path
import json,os,platform,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
SCRIPTS=['gram_candidates.py','column_test.py','column_rank.py','q105_last.py',
 'recover_mixed96_types.py','mixed96_candidates.py','mixed96_columns.py',
 'pure99_frontier.py','pure99_pairs.py','pure99_spectral.py',
 'validate_closure.py','final_certificate.py']

def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled; do not use python -O.')
    dest=ROOT/'results/final_replay_logs';dest.mkdir(exist_ok=True,parents=True)
    # Avoid presenting a stale positive summary after an interrupted/failed run.
    summary_path=ROOT/'results/final_replay_summary.json'
    summary_path.write_text(json.dumps({'completed':False,'all_steps_passed':False})+'\n')
    steps=[];start=time.monotonic()
    prior=[sys.executable,str(ROOT/'run_progress.py'),'--include-baseline','--audit-only']
    commands=[('inherited_35_shells',prior)]
    commands.extend((s,[sys.executable,str(ROOT/'last_three'/s)]) for s in SCRIPTS)
    for name,cmd in commands:
        t=time.monotonic()
        with (dest/(name+'.stdout.txt')).open('w') as out, (dest/(name+'.stderr.txt')).open('w') as err:
            proc=subprocess.run(cmd,cwd=ROOT,stdout=out,stderr=err,
                env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1'})
        steps.append({'program':name,'returncode':proc.returncode,'seconds':time.monotonic()-t})
        print(name+(': PASS' if proc.returncode==0 else ': FAIL'),flush=True)
        if proc.returncode:
            summary_path.write_text(json.dumps({'completed':False,'all_steps_passed':False,
              'failed_program':name,'steps':steps},indent=2)+'\n')
            return 1
    cert=json.loads((ROOT/'last_three/final_certificate.json').read_text())
    if not cert['global_maximality_verified'] or cert['remaining_shells']:
        raise RuntimeError('Certificate has remaining obligations.')
    summary={'completed':True,'all_steps_passed':True,
      'python':sys.version,'platform':platform.platform(),'steps':steps,
      'seconds':time.monotonic()-start,'closed_shells':cert['closed_shells'],
      'remaining_shells':cert['remaining_shells'],'global_maximality_verified':True,
      'independent_full_column_scans':cert['grams_checked_by_two_complete_algorithms']}
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    print('SUCCESS: 38/38 shells, all three final obligations proved; no remaining shells.',flush=True)
    return 0
if __name__=='__main__':raise SystemExit(main())