#!/usr/bin/env python3
"""Replay available certificates; exit 2 unless global proof is complete.
Use --audit-only to request success of the audit, not a global theorem.
"""
import argparse,hashlib,json,subprocess,sys,time,platform,os
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
SCRIPTS=['verify_arithmetic_tests','verify_foundations','verify_low_energy','verify_mixed_small','verify_q81','verify_q99_independent','verify_component_bounds','audit_all_shells']
def main():
 if not __debug__:raise RuntimeError('Do not run verification with -O.')
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit-only',action='store_true');args=p.parse_args()
 source_manifest=ROOT/'manifest_source_sha256.json'
 if source_manifest.exists():
  for path,digest in json.loads(source_manifest.read_text()).items():
   assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,f'Source hash mismatch: {path}'
 logs=ROOT/'results/logs';logs.mkdir(parents=True,exist_ok=True);steps=[];t=time.monotonic()
 for name in SCRIPTS:
  st=time.monotonic();r=subprocess.run([sys.executable,str(ROOT/'programs'/f'{name}.py')],cwd=ROOT,text=True,capture_output=True,env={**os.environ,"PYTHONDONTWRITEBYTECODE":"1"})
  (logs/f'{name}.stdout.txt').write_text(r.stdout);(logs/f'{name}.stderr.txt').write_text(r.stderr)
  print(r.stdout,end='',flush=True)
  steps.append({'program':name,'exit_code':r.returncode,'seconds':time.monotonic()-st})
  if r.returncode:
   print(r.stderr,file=sys.stderr);raise SystemExit(r.returncode)
 report=json.loads((ROOT/'results/all_Q_audit.json').read_text())
 import sympy,networkx
 out={'python':sys.version,'platform':platform.platform(),'sympy':sympy.__version__,'networkx':networkx.__version__,'steps':steps,'seconds':time.monotonic()-t,'available_tests_passed':True,'global_maximality_verified':report['global_maximality_verified'],'closed_Q':report['closed_shells'],'unverified_Q':report['unverified_shells'],'no_mocked_repository_certificates':True}
 (ROOT/'results/replay_summary.json').write_text(json.dumps(out,indent=2)+'\n')
 if not report['global_maximality_verified']:
  print('INCOMPLETE GLOBAL PROOF: local tests passed; unverified shells remain.')
  if not args.audit_only:raise SystemExit(2)
if __name__=='__main__':main()
