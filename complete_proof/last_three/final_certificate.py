"""Fail-closed aggregation of the three new shell proofs and inherited basis.

This checks the coverage joins between certificates.  Source generation,
exhaustive enumeration, and the mathematical lemmas in the report remain
part of the proof; a status flag on its own is not used as a certificate.
"""
from pathlib import Path
from fractions import Fraction as F
import json

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
B=277868041444786176

def load(path):
    return json.loads((ROOT/path).read_text())

def require(condition,message):
    if not condition:
        raise RuntimeError('Unproved obligation: '+message)

def rat(x):
    return F(x['numerator'],x['denominator'])

def run():
    old=load('results/progress_ledger.json')
    require(old['closed_count']==35 and old['remaining_shells']==[96,99,105],
            'the inherited 35-shell basis has changed')
    foundation=load('order15_mu3_unified_audit/results/foundations.json')
    require(foundation['benchmark']['squared_det']==B,'record matrix')
    require(foundation['trace_tail']['from_Q']==171 and
            rat(foundation['trace_tail']['upper'])<B,'full infinite tail')
    require(foundation['Q_residues_mod9']==[0,6],'energy congruences')
    require(foundation['color_reduction']['survivors']==
            [[15,0,0],[14,1,0],[13,2,0],[13,1,1]],'color reduction')
    require(load('new/sparse_defect.json')['closed'],'one-edge defect lemma')
    exp=load('new/exploration.json')
    size13_coverage={}
    for q in (96,99,105):
        residual=exp[str(q)]['res13']
        require(all(row[0] in (0,9,18) for row in residual),'size13 energy coverage')
        e18=[row for row in residual if row[0]==18]
        require(len(e18)==2 and all(row[:3]==[18,78,q-96] for row in e18),
                'size13 remaining cross and outside energies')
        configurations={tuple(map(tuple,row[3])) for row in e18}
        expected={((1,0),)*10+((3,2),),((1,0),)*9+((2,1),(2,1))}
        require(configurations==expected,'P3 and two-edge support coverage')
        size13_coverage[str(q)]={'remaining_internal_energy':18,
          'cross_energy':78,'outside_norm':q-96,'core_families':['P3','2K2'],
          'energy0_reason':'Published B3(15,10)=12',
          'energy9_reason':'Replayed single-edge defect lemma'}
    grams=load('last_three/gram_candidates.json')['rows']
    columns=load('last_three/column_test.json')
    ranks=load('last_three/column_rank.json')['rows']
    require(columns['complete'] and len(grams)==len(columns['rows'])==len(ranks)==142,
            'complete size13 candidate and column lists')
    for i,(g,c,r) in enumerate(zip(grams,columns['rows'],ranks)):
        require(all(c[k]==v for k,v in g.items()),'candidate-column provenance')
        require(c['candidate_index']==r['candidate_index']==i and
                g['Q']==r['Q'] and g['D']==r['D'],'rank provenance')
        require(r['excluded']==(r['span_rank']<15),'rank exclusion predicate')
    last=load('last_three/q105_last.json')
    survivors=[r for r in ranks if not r['excluded']]
    require(len(survivors)==1 and survivors[0]['candidate_index']==last['residual_candidate_index']
            and survivors[0]['Q']==105,'unique remaining column-span case')
    require(last['analytic_inverse_identity_checked'] and last['analytic_column_obstruction'],
            'analytic final column lemma')
    require(last['equal_hubs_required_Eisenstein_norm']==1666 and
            last['Eisenstein_norm_residues_mod4']==[0,1,3] and
            F(last['different_hubs_quadratic_lower'])>1,'norm and positivity contradictions')
    require(last['required_disagreements']==16>last['maximum_disagreements']==15,
            'disjoint disagreement count')
    mixed={row['Q']:row for row in load('new/mixed_rank_frontier.json')}
    require(mixed[105]['closed_double_14_1'] and not mixed[105]['remaining'],
            'Q105 double14 branch with isolates')
    rows=load('order15_mu3_unified_audit/results/component_bounds.json')['component_rows']
    q105=[row for row in rows if row['Q']==105]
    require(len(q105)==1 and q105[0]['internal']==63 and q105[0]['cross']==42 and
            not q105[0]['residuals'] and rat(q105[0]['max_upper'])<=B,
            'Q105 double14 branch without isolates, including equality')
    # For Q96 the 14 vertices cannot be covered without at least seven
    # internal edges: that would already require Q >= 7*9+42=105.
    require(96<7*9+42,'Q96 no-isolate energy exclusion')
    signatures=load('last_three/mixed96_type_recovery.json')
    require({(r['v'],r['u'],tuple(r['signature']),r['name']) for r in signatures}==
      {(3,2,(1,2,2),'P3'),(4,4,(2,4,2),'C4'),(5,6,(3,5,3),'K23')},
      'Q96 recovered rank signatures')
    require(all(r['matches'] and all(m['zero_cycle_flux'] for m in r['matches'])
                for r in signatures),'complete balanced-core recovery')
    ranktypes=sorted(load('new/rank_catalogue.json')['types'],
      key=lambda t:(t['u'],t['v'],t['nullity'],t['kernel_support'],t['independence']))
    observed=[]
    for row in mixed[96]['remaining']:
        tt=[ranktypes[i] for i in row['component_types']]
        observed.append((row['isolates'],tuple((t['v'],t['u'],t['nullity'],
                         t['kernel_support'],t['independence']) for t in tt)))
    require(set(observed)=={(7,((3,2,1,2,2),(4,4,2,4,2))),
                 (10,((4,4,2,4,2),)),(9,((5,6,3,5,3),))},
            'Q96 all paired signatures covered')
    mg=load('last_three/mixed96_candidates.json');mc=load('last_three/mixed96_columns.json')
    require(mc['complete'] and len(mg['rows'])==len(mc['rows'])==52,
            'Q96 full Gram representatives')
    require({r['family'] for r in mg['stats']}=={'P3_C4','C4','K23'},
            'all Q96 support families generated')
    for i,(g,c) in enumerate(zip(mg['rows'],mc['rows'])):
        require(c['candidate_index']==i and all(c[k]==v for k,v in g.items()),
                'Q96 candidate-column provenance')
        require(c['span_rank']<15 and c['excluded'],'Q96 allowed-column span')
    pure=load('last_three/pure99_frontier.json')
    require(not pure['isolate']['survivors'],'Q99 two-sided isolate projection')
    pairs=load('last_three/pure99_pairs.json');spec=load('last_three/pure99_spectral.json')
    require(pure['no_isolate']['assembly_count']==pairs['assembly_count']==632,
            'Q99 full no-isolate assembly coverage')
    require(len(pairs['candidates'])==14 and
       {tuple(c['polynomial']) for c in pairs['candidates']}==
       {tuple(c['polynomial']) for c in pure['no_isolate']['norm_above_B']},
       'Q99 all component decompositions retained')
    expected_pairs={(i,j) for i,a in enumerate(pairs['candidates'])
        for j,b in enumerate(pairs['candidates']) if i<=j and a['polynomial']==b['polynomial']}
    witnesses=spec['spectral_pair_obstructions']
    require(len(witnesses)==len(expected_pairs)==16 and
            {(w['left_id'],w['right_id']) for w in witnesses}==expected_pairs,
            'Q99 all equal-spectrum pairs checked')
    for w in witnesses:
        v=w['obstruction']
        require(v is not None and v['flat_squared_frobenius']>v['spectral_trace_upper'],
                'Q99 spectral trace witness')
    validation=load('last_three/closure_validation.json')
    require(validation['all_passed'],'independent arithmetic verification')
    expected_ids={('size13',i) for i in range(142)}|{('mixed96',i) for i in range(52)}
    require({(r['group'],r['candidate_index']) for r in validation['matrix_checks']}==expected_ids,
            'all 194 determinant, color and inverse checks')
    direct=validation['direct_full_enumerations']
    require(len(direct)==194 and {(r['group'],r['candidate_index']) for r in direct}==expected_ids
            and all(r['assignments']==3**14 and r['identical_to_incremental'] for r in direct),
            'second complete enumeration of every candidate')
    annihilators=validation['annihilators']
    require(len(annihilators)==193 and {(r['group'],r['candidate_index']) for r in annihilators}==
        expected_ids-{('size13',last['residual_candidate_index'])},'193 exact annihilators')
    methods={96:'Normalized size13 columns; exhaustive double14 rank-signature recovery and column annihilators.',
             99:'Normalized size13 columns; pure-Gram projection and all 16 exact spectral trace witnesses.',
             105:'Normalized size13 columns; analytic norm-mod4 and 16>15 contradiction; inherited double14 bounds.'}
    universe={q for q in range(169) if q%9 in (0,6)}
    closed=set(old['closed_shells'])|set(methods)
    require(closed==universe,'all finite shells joined to the infinite tail')
    rows=[]
    for row in old['rows']:
        row=dict(row)
        if row['Q'] in methods:
            row.update(status='CLOSED_IN_FINAL_SUPPLEMENT',proof_basis=methods[row['Q']],
                       evidence_files='last_three/final_certificate.json',remaining_obligation='')
        rows.append(row)
    result={'benchmark_squared_determinant':B,'maximum_modulus_exact':'120932352 * sqrt(19)',
      'closed_shells':sorted(closed),'closed_count':len(closed),'total_shells':len(universe),
      'newly_closed_shells':sorted(methods),'remaining_shells':[],
      'global_maximality_verified':True,'infinite_tail_from_Q':171,
      'size13_coverage':size13_coverage,'size13_candidate_grams':len(grams),
      'extra_Q96_candidate_grams':len(mg['rows']),'exact_annihilators':len(annihilators),
      'assignments_per_gram_per_algorithm':3**14,'grams_checked_by_two_complete_algorithms':len(direct),
      'total_column_assignments_per_algorithm':len(direct)*3**14,
      'Q99_equal_spectrum_pairs':len(expected_pairs),'Q105_final_analytic_obstruction':last,
      'external_theorem_dependency':'Todorov--Bogdanova (2020), B3(15,10)=12; Proposition 14 and Table 1.',
      'scope':'Computer-assisted proof, using the retained mathematical lemmas and published code theorem; not a proof-assistant formalization. No uniqueness classification is claimed.',
      'rows':rows}
    (HERE/'final_certificate.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    (ROOT/'results/final_ledger.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('NEWLY CLOSED:',sorted(methods));print('COMPLETE: 38/38; infinite tail verified; maximum squared determinant =',B)
    return result

if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Do not run with Python -O.')
    run()
