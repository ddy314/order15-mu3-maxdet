#!/usr/bin/env python3
"""All 5,806,080 permitted row transforms; exact matrix entry comparison."""
from pathlib import Path
import json,subprocess,tempfile
ROOT=Path(__file__).resolve().parents[2]
def main():
 if not __debug__:raise RuntimeError('Assertions must be enabled.')
 data=json.loads((ROOT/'data/revision/friendship_nine_factors.json').read_text())
 with tempfile.TemporaryDirectory(prefix='order15-factor-audit-') as tmp:
  exe=Path(tmp)/'direct'
  subprocess.run(['g++','-O3','-std=c++17',str(Path(__file__).with_name('direct.cpp')),'-o',str(exe)],check=True)
  inp=' '.join(str(v) for m in data for row in m['exponents'] for v in row)+'\n'
  result=json.loads(subprocess.run([str(exe)],input=inp,capture_output=True,text=True,check=True).stdout)
 assert result['row_transformations']==5806080
 for i,row in enumerate(result['projective_equivalence_counts']):
  for j,c in enumerate(row):assert 3*c==(data[i]['automorphism_order'] if i==j else 0),(i,j,c)
 (ROOT/'data/revision/direct_factor_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS: 5,806,080 direct row transforms, all 81 equivalence comparisons and stabilizer orders.')
if __name__=='__main__':main()
