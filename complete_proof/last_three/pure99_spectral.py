"""Exact trace obstruction for EVERY equal-spectrum Q99 component pair.

For normalized components T_A,T_B let p=gcd(char(T_A),char(T_B)).
If X is a corresponding submatrix of a factor H with HH*=15I+3T,
then ||X||_F^2 <= 15 deg(p)-3 p_1.  Every entry of X has modulus one.
All unions of connected components are checked, not just one component.
"""
from pure99_pairs import *
from functools import lru_cache

@lru_cache(None)
def subsets(components):
    # Subsets with identical characteristic polynomial have identical
    # dimension and yield identical inequalities; retain one witness.
    states={(1,):()}
    for idx,p in enumerate(components):
        for a,ids in list(states.items()):
            states.setdefault(pmul(a,p),ids+(idx,))
    return [(p,inds) for p,inds in states.items() if len(p)>1]

@lru_cache(None)
def spectral_cap(pa,pb):
    g=sp.gcd(sp.Poly.from_list(pa,x),sp.Poly.from_list(pb,x))
    co=tuple(map(int,g.all_coeffs()));r=g.degree()
    assert co[0]==1
    return 15*r-(3*co[1] if r else 0),co

def run():
    data=json.loads(Path(__file__).with_name('pure99_pairs.json').read_text())
    candidates=data['candidates'];out=[];rem=[]
    for ia,A in enumerate(candidates):
        for ib in range(ia,len(candidates)):
            BB=candidates[ib]
            if A['polynomial']!=BB['polynomial']:
                continue
            sa=subsets(tuple(tuple(c['p']) for c in A['components']))
            sb=subsets(tuple(tuple(c['p']) for c in BB['components']))
            found=None
            for pa,ida in sa:
                for pb,idb in sb:
                    cap,common=spectral_cap(pa,pb)
                    volume=(len(pa)-1)*(len(pb)-1)
                    if volume>cap:
                        found={'left_components':ida,'right_components':idb,
                          'left_polynomial':pa,'right_polynomial':pb,'gcd':common,
                          'flat_squared_frobenius':volume,'spectral_trace_upper':cap}
                        break
                if found:
                    break
            rr={'left_id':ia,'right_id':ib,'obstruction':found};out.append(rr)
            if found:
                print(rr,flush=True)
            else:
                rem.append([ia,ib])
    assert len(candidates)==14 and len(out)==16
    result={'candidate_component_assemblies':len(candidates),
      'equal_spectrum_unordered_pairs':len(out),
      'spectral_pair_obstructions':out,'remaining_pure_pure_pairs':rem,
      'all_pairs_excluded':bool(out) and not rem}
    Path(__file__).with_name('pure99_spectral.json').write_text(json.dumps(result,indent=2)+'\n')
    print('remaining pure pairs',rem,flush=True)
    return result
if __name__=='__main__':run()
