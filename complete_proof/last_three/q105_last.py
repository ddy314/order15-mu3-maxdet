"""Analytic exclusion of the final Q105 candidate, also checked by enumeration.

If both core pairs disagree, the exact inverse quadratic form either exceeds
one (different hub entries), or forces an Eisenstein norm to equal 1666.
The latter is impossible modulo four.  Thus both pairs cannot disagree in
one column; their required eight disagreements each cannot fit in 15 columns.
"""
from column_rank import *

def run():
    ranks=json.loads(Path(__file__).with_name('column_rank.json').read_text())
    sur=[r for r in ranks['rows'] if not r['excluded']]
    assert len(sur)==1 and sur[0]['Q']==105 and sur[0]['family']=='2K2'
    data=json.loads(Path(__file__).with_name('column_test.json').read_text())
    assert data['complete']
    r=data['rows'][sur[0]['candidate_index']];G=full_gram(r)
    # Fix the actual residual representative; no generic Gram is substituted.
    assert G[0][1]==G[2][3]==G[13][14]==(3,0)
    assert all(G[i][13]==G[i][14]==(1,-1) for i in range(13))
    A=[row[:13] for row in G[:13]];J=inverse(A)
    for i in range(13):
        for j in range(13):
            expected=(15,0) if i==j else ((3,0) if {i,j} in ({0,1},{2,3}) else ZERO)
            assert A[i][j]==expected
    # Inverse identity as a sum of exact Hermitian quadratic forms:
    # q=x*A^{-1}x + (15/392)|u+v-2 bar(delta) k|^2 + |u-v|^2/24,
    # k=(a+b+c+d)/18 + sum(other nine entries)/15.
    bd=conj((1,-1))
    ell=[scale(bd,-F(1,9))]*4+[scale(bd,-F(2,15))]*9+[ONE,ONE]
    em=[ZERO]*13+[ONE,neg(ONE)]
    recovered=[]
    for i in range(15):
        line=[]
        for j in range(15):
            z=J[i][j] if i<13 and j<13 else ZERO
            z=add(z,scale(mul(conj(ell[i]),ell[j]),F(15,392)))
            z=add(z,scale(mul(conj(em[i]),em[j]),F(1,24)))
            line.append(z)
        recovered.append(line)
    assert recovered==inverse(G)
    # Both core pairs disagree: each pair contributes 11/72, the other
    # nine coordinates contribute 9/15, irrespective of their phases.
    base=2*F(11,72)+F(9,15)
    different_hubs_lower=base+F(3,24)
    assert different_hubs_lower>1
    # Equal hub entries u=v yield
    # N(30 delta u - 5(a+b+c+d) - 6 sum(other nine)) = 1666.
    impossible_norm=(1-base)*F(98,15)*F(90**2,3)
    assert impossible_norm.denominator==1
    residues=sorted({(a*a-a*b+b*b)%4 for a in range(4) for b in range(4)})
    assert int(impossible_norm)%4 not in residues
    counts={}
    for serial in r['column_serials']:
        v=decode(serial);a=v[0]!=v[1];b=v[2]!=v[3]
        key=f'{int(a)},{int(b)}';counts[key]=counts.get(key,0)+1
        assert not(a and b)
    agreements=F(2*G[0][1][0]+15,3)
    assert agreements.denominator==1
    needed=2*(15-int(agreements));available=15
    assert needed>available
    out={'Q':105,'residual_candidate_index':r['candidate_index'],'D':r['D'],
      'columns':len(r['column_serials']),'disagreement_counts':counts,
      'analytic_inverse_identity_checked':True,'both_pairs_different_base':str(base),
      'different_hubs_quadratic_lower':str(different_hubs_lower),
      'equal_hubs_required_Eisenstein_norm':int(impossible_norm),
      'Eisenstein_norm_residues_mod4':residues,'analytic_column_obstruction':True,
      'row_pairs':[[0,1],[2,3]],'agreements_per_pair':int(agreements),
      'required_disagreements':needed,'maximum_disagreements':available,'contradiction':True}
    Path(__file__).with_name('q105_last.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out,flush=True);return out
if __name__=='__main__':run()
