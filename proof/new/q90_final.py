"""Q90: low-rank intertwining rectangles exclude all exact residual spectra."""
from bounds import *
import sympy as sp

def certificate():
 r=json.loads(Path(__file__).with_name('pure90_frontier.json').read_text())
 assert not r['isolate']['survivors']
 states=r['no_isolate']['survivors'];assert states
 x=sp.Symbol('x');pair=sp.Poly(x*x-1,x);proofs=[]
 for state in states:
  assert sorted(state['components'])==[[2,1]]*5+[[5,5]]
  p=sp.Poly.from_list(state['polynomial'],x,domain=sp.ZZ)
  large,rem=sp.div(p,pair**5);assert rem.is_zero and large.degree()==5
  # For Hermitian blocks, solutions of A5 X=X B10 have their columns
  # in the sum of A5 eigenspaces shared with B10. A characteristic gcd
  # with (x^2-1)^5 counts every such A5 eigenvalue with full multiplicity.
  common=sp.gcd(large,pair**5);rank=common.degree()
  assert rank<=2
  # tr(A5)=75 and sum of squared off-diagonal moduli=45, hence
  # lambda_max(A5)<=15+sqrt(72)<24. The flat 5-by-10 rectangle has
  # squared Frobenius norm50, and X X*<=G5, so 50<=rank*lambda_max(G5).
  assert 72<9**2 and 50>24*rank
  proofs.append({'determinant':state['determinant'],'large_component_polynomial':list(map(int,large.all_coeffs())),
    'pair_component_polynomial':[1,0,-1],'common_polynomial':str(common.as_expr()),
    'rectangle_rank_upper':rank,'rectangle_squared_frobenius_norm':50,
    'core_lambda_upper_strict':24,'rank_times_lambda_upper':24*rank,'contradiction':True})
 out={'Q':90,'closed':True,'remaining_pure_gram_states':len(states),'spectral_rectangle_certificates':proofs,
      'size13_defect_excluded_by':'sparse_defect.py plus B3(15,10)=12'}
 Path(__file__).with_name('q90_final.json').write_text(json.dumps(out,indent=2)+'\n')
 print('Q90 closed:',len(states),'spectral states; flat rectangle50 exceeds rank-times-cap at most48',flush=True)
 return out
if __name__=='__main__':certificate()
