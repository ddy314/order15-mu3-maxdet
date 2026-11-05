#!/usr/bin/env python3
"""Run all reconstructed continuation certificates, failing on any missing step."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parent
STEPS=['tools/verify_frame_controls.py','tools/verify_continuation.py','tools/verify_fano_completions.py',
       'tools/verify_mixed_inert_gap.py','extensions/energy/verify_pure_energy.py',
       'extensions/energy/verify_minimum_support.py','extensions/energy/audit_equality.py',
       'extensions/energy/verify_minimum_phases.py','extensions/sixth_roots/verify_gap.py']

def main():
 if not __debug__:raise RuntimeError('Assertions must be enabled.')
 dest=ROOT/'data/revision';dest.mkdir(exist_ok=True,parents=True)
 summary=dest/'replay_summary.json';summary.write_text('{"completed":false}\n')
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
 start=time.monotonic();steps=[]
 with tempfile.TemporaryDirectory(prefix='order15-orthogonality-') as tmp:
  exe=Path(tmp)/'verify'
  subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'extensions/orthogonality/verify.cpp'),'-o',str(exe)],check=True)
  subprocess.run([str(exe),str(dest/'orthogonality_certificate.json')],check=True)
  sys.path.insert(0,str(ROOT/'tools'))
  from matrix_exact import gram
  orth=json.loads((dest/'orthogonality_certificate.json').read_text())
  for case in orth['cases']:
   h=case['witness'];g=gram(h)
   assert g==[[(15,0) if i==j else (0,0) for j in range(len(h))] for i in range(len(h))]
  steps.append('independent orthogonality graph and small-graph controls')
 for path in STEPS:
  print('\n==> '+path,flush=True)
  subprocess.run([sys.executable,str(ROOT/path)],cwd=ROOT,env=env,check=True)
  steps.append(path)
 certs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.rglob('*.json') if p!=summary}
 source_paths=sorted(set(ROOT.rglob('*.py')) | set(ROOT.rglob('*.cpp')) | set((ROOT/'paper').rglob('*.tex')))
 sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths if '.venv' not in p.parts and '.git' not in p.parts}
 summary.write_text(json.dumps({'source_sha256':sources,'completed':True,'all_steps_passed':True,'steps':steps,'certificate_sha256':certs,'seconds':time.monotonic()-start},indent=2)+'\n')
 print('SUCCESS: every recovered continuation obligation and the new sixth-root Q>=11 certificate passed.',flush=True)
if __name__=='__main__':main()
