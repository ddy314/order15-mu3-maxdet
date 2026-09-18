"""Compare historical Q150 bounds against the trusted record, using exact ratios."""
from pathlib import Path
from fractions import Fraction
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'historical_remote'))
from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
from independent_audit.programs.exact_core import B
root=Path(__file__).resolve().parents[1]
rows=[]
for internal,cross,cap in [(63,87,Fraction(25)),(72,78,Fraction(127,5)),(81,69,Fraction(261,10)),(90,60,Fraction(27)),(99,51,Fraction(138,5)),(108,42,Fraction(563,20))]:
 upper=stationary_upper_dimension(14,internal)*(15-Fraction(cross)/cap)
 rows.append({'internal':internal,'cross':cross,'upper_numerator':upper.numerator,'upper_denominator':upper.denominator,'exceeds_record':upper>B,'ratio_to_record':str(upper/B)})
assert sum(r['exceeds_record'] for r in rows)==5
out={'benchmark':B,'source':'order15_mu3/scripts/verify_q150_boundary.py','comparison':'record B, rather than Q153 trace envelope','rows':rows,'historical_q150_bounds_close_entire_shell_against_record':False,'independent_audit_remains_authority':True}
(root/'provenance/remote_bound_comparison.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS: exact comparison confirms 5 of 6 historical Q150 size14 bounds exceed B.')
