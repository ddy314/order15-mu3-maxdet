#!/usr/bin/env python3
"""Exact checks for the sharper difference-matrix reduction.

External mathematical input: Lampio--Ostergard (2011), M_3(15)=9.
The external classification is NOT rerun by this program.
Run proof/run_final.py first to regenerate all inherited finite catalogues.
This program checks the new inequalities, explicit independent-set witnesses,
spectral rectangles and determinant-lattice arithmetic, not stored 'closed' labels.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import sys,json
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'proof'
B=277868041444786176
sys.path.insert(0,str(P/'last_three'))
from gram_candidates import full_gram
from mixed96_candidates import full_mixed_gram
from matrix_exact import det, norm

def load(rel):
    return json.loads((P/rel).read_text())

def independent_witness(G,size=10):
    return next(([i+1 for i in s] for s in combinations(range(len(G)),size)
                 if all(tuple(G[i][j])==(0,0) for i,j in combinations(s,2))),None)

def norm_obstruction(n):
    fac=sp.factorint(n)
    return next((int(p) for p,e in fac.items() if p%3==2 and e%2),None)

def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')

    partitions=[(a,b,15-a-b) for a in range(15,-1,-1)
                for b in range(a,-1,-1) if 0<=15-a-b<=b]
    sets={}
    for cap in (12,9):
        candidates=set()
        for a,b,c in partitions:
            cross=3*(a*b+a*c+b*c)
            floor=cross+9*sum(max(0,x-cap) for x in (a,b,c))
            candidates.update(q for q in range(169)
                              if q>=floor and (q-cross)%9==0)
        sets[cap]=sorted(candidates)

    assert len(sets[12])==29 and len(sets[9])==23
    all38=[q for q in range(169) if q%9 in (0,6)]
    t=F(319,250)
    tail=(15-t)**14*(15+14*t)
    assert t*t<F(171,105) and tail<B

    # Construct every one of the legacy 194 Gram matrices and exhibit ten
    # independent support vertices. These exclude them under the sharp theorem.
    records=[]
    for name,constructor in [
        ('last_three/gram_candidates.json',full_gram),
        ('last_three/mixed96_candidates.json',full_mixed_gram),
    ]:
        for i,row in enumerate(load(name)['rows']):
            G=constructor(row)
            w=independent_witness(G)
            assert w is not None,(name,i)
            assert det(G)==(int(row['D']),0)
            records.append({
                'source':name,
                'row_index_0based':i,
                'Q':row['Q'],
                'independent_vertices_1based':w,
            })
    assert len(records)==194

    r114=load('new/q114_final.json')
    w114=independent_witness(r114['canonical_gram'])
    assert w114 and det(r114['canonical_gram'])==(r114['last_abstract_determinant'],0)

    # Check the three surviving mixed-Q96 support signatures directly.
    mixed=load('new/mixed_rank_frontier.json')
    ranktypes=load('new/rank_catalogue.json')['types']
    ranktypes=sorted(ranktypes,key=lambda a:(a['u'],a['v'],a['nullity'],a['kernel_support'],a['independence']))
    signature_bounds=[]
    for row in mixed:
        for state in row['remaining']:
            alpha=state['isolates']+sum(ranktypes[i]['independence'] for i in state['component_types'])
            assert alpha>9
            signature_bounds.append({
                'Q':row['Q'],
                'configuration':state['configuration'],
                'independence_lower_bound':alpha,
            })
    assert len(signature_bounds)==3

    # Re-evaluate every stated spectral-rectangle contradiction with integer polynomials.
    x=sp.Symbol('x')
    rectangles=[]
    for row in load('last_three/pure99_spectral.json')['spectral_pair_obstructions']:
        a=row['obstruction']
        pa=sp.Poly.from_list(a['left_polynomial'],x)
        pb=sp.Poly.from_list(a['right_polynomial'],x)
        g=sp.gcd(pa,pb)
        co=list(map(int,g.all_coeffs()))
        cap=15*g.degree()-(3*co[1] if g.degree() else 0)
        volume=pa.degree()*pb.degree()
        assert volume==a['flat_squared_frobenius']
        assert cap==a['spectral_trace_upper']
        assert volume>cap
        rectangles.append({
            'Q':99,
            'left_id':row['left_id'],
            'right_id':row['right_id'],
            'volume':int(volume),
            'cap':int(cap),
            'gcd':co,
        })
    assert len(rectangles)==16

    q90=load('new/pure90_frontier.json')['no_isolate']['survivors']
    for row in q90:
        poly=sp.Poly.from_list(row['polynomial'],x)
        small=sp.Poly(x*x-1,x)**5
        big,rem=sp.div(poly,small)
        assert rem.is_zero
        g=sp.gcd(big,small)
        co=list(map(int,g.all_coeffs()))
        cap=15*g.degree()-3*co[1]
        assert 50>cap
        rectangles.append({'Q':90,'volume':50,'cap':int(cap),'gcd':co})

    # The norm test is on the FULL determinant, not on R alone (400 is a norm).
    fac=15**11*3**4
    obstruction126=[
        {
            'R':r,
            'full_determinant':r*fac,
            'inert_prime_odd_valuation':norm_obstruction(r*fac),
        }
        for r in range(397,402)
    ]
    assert all(a['inert_prime_odd_valuation'] for a in obstruction126)

    # Exact real-root checks of ALL bounded characteristic polynomials at Q117.
    counts={}
    for label,row in load('new/q117_integer_sieve.json').items():
        for a in row['tested']:
            coeff=[1,0,-13,-a['e3'],a['e4']]
            if label=='rank5':
                coeff.append(-a['e5'])
            poly=sp.Poly.from_list(coeff,x)
            assert sum(m for interval,m in sp.polys.polytools.intervals(poly))<poly.degree()
        counts[label]=len(row['tested'])

    out={
        'verified':True,
        'external_input':{
            'statement':'M_3(15)=9',
            'source':'Lampio and Ostergard, JSPI 141 (2011), 1194-1207, doi:10.1016/j.jspi.2010.09.023',
            'classification_rerun':False,
        },
        'coarse_38':all38,
        'equidistant_code_29':sets[12],
        'sharp_difference_matrix_23':sets[9],
        'removed_by_sharp_orthogonality':[q for q in all38 if q not in sets[9]],
        'tail_rational_upper_over_B':str(tail/B),
        'legacy_194_independent_set_witnesses':records,
        'q114_independent_vertices_1based':w114,
        'mixed_signature_independence':signature_bounds,
        'exact_spectral_rectangles':rectangles,
        'q126_full_determinant_obstructions':obstruction126,
        'q117_polynomials_rechecked':counts,
        'scope':'Checks the new reductions and stated residual witnesses; not a replacement for generation/completeness checks in proof/run_final.py.',
    }
    (ROOT/'data/sharp_reduction_verified.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: 38 -> 29 -> 23 energy candidates; 194 exact Gram independent-set witnesses; Q114; all Q90/Q99 rectangles; Q117 real-root tests; Q126 FULL determinant norm obstructions.')

if __name__=='__main__':
    main()
