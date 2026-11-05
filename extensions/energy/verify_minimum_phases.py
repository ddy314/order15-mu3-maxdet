#!/usr/bin/env python3
"""Rebuild the column certificates, then independently verify their coverage."""
import json,subprocess,tempfile,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent

def main():
 if not __debug__:raise RuntimeError('Assertions must be enabled.')
 dest=ROOT/'data/revision/phase_certificates';dest.mkdir(exist_ok=True)
 with tempfile.TemporaryDirectory(prefix='order15-phases-') as tmp:
  tmp=Path(tmp)
  for source,name in [('minimum_phases.cpp','producer'),('check_phase_certificates.cpp','checker')]:
   subprocess.run(['g++','-std=c++17','-O3',str(HERE/source),'-o',str(tmp/name)],check=True)
  subprocess.run([str(tmp/'producer'),str(dest)],check=True)
  results=[]
  for path in sorted(dest.glob('phase_*.json')):
   c=json.loads(path.read_text());lines=[' '.join(map(str,c['multiplicities'])),str(len(c['columns'])),' '.join(map(str,c['columns'])),' '.join(map(str,c['coloring'])),str(len(c['local_colorings']))]
   for local in c['local_colorings']:
    lines.extend([f"{local['vertex']} {len(local['neighbors'])}",' '.join(map(str,local['neighbors'])),' '.join(map(str,local['colors']))])
   cert=tmp/'input.txt';cert.write_text('\n'.join(lines)+'\n')
   proc=subprocess.run([str(tmp/'checker'),str(cert)],check=True,text=True,capture_output=True)
   count,pairs,local=map(int,proc.stdout.split())
   results.append({'file':path.name,'columns':count,'same_color_pairs':pairs,'neighborhoods':local,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
   print('Checked',path.name,proc.stdout.strip(),flush=True)
  assert len(results)==34
  # Three directed negative controls, each independently rejected.
  import copy
  def encode(c):
   lines=[' '.join(map(str,c['multiplicities'])),str(len(c['columns'])),' '.join(map(str,c['columns'])),' '.join(map(str,c['coloring'])),str(len(c['local_colorings']))]
   for local in c['local_colorings']:
    lines += [f"{local['vertex']} {len(local['neighbors'])}",' '.join(map(str,local['neighbors'])),' '.join(map(str,local['colors']))]
   return '\n'.join(lines)+'\n'
  c=next(json.loads(p.read_text()) for p in dest.glob('phase_*.json') if json.loads(p.read_text())['local_colorings'])
  bad=tmp/'bad.txt'
  variants=[]
  omitted=copy.deepcopy(c);omitted['columns']=[];omitted['coloring']=[];omitted['local_colorings']=[];variants.append(('missing_column',omitted))
  improper=copy.deepcopy(c);local=improper['local_colorings'][0];improper['coloring'][local['vertex']]=improper['coloring'][local['neighbors'][0]];variants.append(('improper_global_color',improper))
  neighbor=copy.deepcopy(c);neighbor['local_colorings'][0]['neighbors'].pop();neighbor['local_colorings'][0]['colors'].pop();variants.append(('missing_actual_neighbor',neighbor))
  negatives=[]
  for name,variant in variants:
   bad.write_text(encode(variant))
   assert subprocess.run([str(tmp/'checker'),str(bad)],cwd=tmp,capture_output=True).returncode!=0
   negatives.append(name)
  report={'phase_multisets':120,'norm_survivors':35,'excluded_nonfriendship':34,'independent_columns_scanned':34*3**14,'records':results,'negative_controls_rejected':negatives}
  (ROOT/'data/revision/minimum_phase_audit.json').write_text(json.dumps(report,indent=2)+'\n')
 print('PASS: all 34 non-friendship phase multisets excluded by checked coloring certificates.')
if __name__=='__main__':main()
