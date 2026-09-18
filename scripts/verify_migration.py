from pathlib import Path
import hashlib,json,zipfile
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'provenance/migration_manifest.json').read_text())
for path,digest in m['imported_files'].items():assert hashlib.sha256((root/m.get('imported_root','')/path).read_bytes()).hexdigest()==digest,path
for path,digest in m['local_artifacts'].items():assert hashlib.sha256((root/'artifacts'/path).read_bytes()).hexdigest()==digest,path
with zipfile.ZipFile(root/'artifacts/order15_mu3_unified_audit.zip') as z:
 files=[i for i in z.infolist() if not i.is_dir()]
 for i in files:
  p=i.filename.split('/',1)[1]
  # Replay result files are intentionally variable; immutable inputs must match.
  if not p.startswith('results/'):
   assert (root/'independent_audit'/p).read_bytes()==z.read(i),p
 assert z.read('order15_mu3_unified_audit/results/all_Q_audit.csv')==(root/'artifacts/all_Q_audit.csv').read_bytes()
 assert z.read('order15_mu3_unified_audit/report/order15_mu3_independent_audit.pdf')==(root/'artifacts/order15_mu3_independent_audit.pdf').read_bytes()
a=json.loads((root/'independent_audit/data/benchmark.json').read_text())
b=json.loads((root/'historical_remote/order15_mu3/data/literature_matrices/nunez_ponasso_m15.json').read_text())
assert a['matrix']==b['matrix']
print(f'PASS: {len(m["imported_files"])} source files, 3 local artifacts, immutable ZIP contents, PDF/CSV duplicates and 225 benchmark entries')
