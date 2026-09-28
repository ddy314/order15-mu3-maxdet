#!/usr/bin/env python3
"""Recompute the witness invariants and its COMPLETE mu_3-monomial stabilizer.

The exhaustiveness argument (unique hub, seven pairs, forced common phase)
is in the manuscript. This enumerates all 7! * 2^7 permitted row permutations,
then verifies the returned pairs by direct comparison of all 225 entries.
Requires Python 3.10+, SymPy and a C++17 compiler (default: g++).
"""
from __future__ import annotations
import argparse, itertools, json, math, os, subprocess, tempfile
from collections import Counter
from pathlib import Path
from sympy.combinatorics import Permutation, PermutationGroup
from matrix_exact import ROOTS, gram, det, norm, mul, conj, orthogonal_sets

ROOT=Path(__file__).resolve().parents[1]
B=277868041444786176

def main() -> None:
    if not __debug__:
        raise RuntimeError('Assertions must be enabled.')

    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--compiler',default=os.environ.get('CXX','g++'))
    args=ap.parse_args()

    original=json.loads((ROOT/'data/benchmark_original.json').read_text())
    E=[
        [
            (original[i][j]-original[i][0]-original[0][j]+original[0][0])%3
            for j in range(15)
        ]
        for i in range(15)
    ]
    assert E==json.loads((ROOT/'data/benchmark_normalized.json').read_text())

    G=gram(E)
    K=[[conj(z) for z in row] for row in gram(list(map(list,zip(*E))))]

    h=next(
        i for i in range(15)
        if sum(G[i][j]!=(0,0) for j in range(15) if i!=j)==14
    )
    pairs=[
        (i,j)
        for i,j in itertools.combinations(range(15),2)
        if norm(G[i][j])==9
    ]
    assert len(pairs)==7 and len(set(sum((list(p) for p in pairs),[])))==14

    delta=(1,-1)
    phases=[
        next(k for k in range(3) if G[i][h]==mul(delta,ROOTS[k]))
        if i!=h else 0
        for i in range(15)
    ]
    order=[i for p in pairs for i in p]+[h]
    F=[
        [(E[i][j]-phases[i]-E[h][j])%3 for j in range(15)]
        for i in order
    ]
    assert F==json.loads((ROOT/'data/benchmark_canonical.json').read_text())

    GC=gram(F)
    for i,j in itertools.product(range(15),repeat=2):
        target=(
            (15,0) if i==j
            else delta if j==14
            else conj(delta) if i==14
            else (3,0) if i//2==j//2
            else (0,0)
        )
        assert GC[i][j]==target

    determinants=[
        det([[ROOTS[e] for e in row] for row in A])
        for A in [original,E,F]
    ]
    assert determinants==[
        (604661760,241864704),
        (-362797056,-604661760),
        (241864704,-362797056),
    ]
    assert all(norm(z)==B for z in determinants)
    assert det(G)==(B,0)
    assert 216**7*38==3*B
    assert 12**7*18**6*228==B

    row_alpha,row_sets=orthogonal_sets(E)
    col_alpha,col_sets=orthogonal_sets(list(map(list,zip(*E))))
    assert (row_alpha,col_alpha,len(row_sets),len(col_sets))==(7,7,128,128)

    histogram=Counter(
        norm(G[i][j])
        for i,j in itertools.combinations(range(15),2)
    )
    assert histogram=={0:84,3:14,9:7}
    assert sum(k*v for k,v in histogram.items())==105
    assert G!=K

    with tempfile.TemporaryDirectory(prefix='order15-aut-') as td:
        exe=Path(td)/'enumerate_aut'
        subprocess.run(
            [
                args.compiler,'-std=c++17','-O2',
                str(ROOT/'tools/enumerate_aut.cpp'),
                '-o',str(exe),
            ],
            check=True,
        )
        inp='\n'.join(' '.join(map(str,row)) for row in F)+'\n'
        proc=subprocess.run(
            [str(exe)],
            input=inp,
            text=True,
            capture_output=True,
            check=True,
        )

    records=json.loads(proc.stdout)
    assert len(records)==336 and 'searched=645120' in proc.stderr
    assert len({tuple(a['rows']) for a in records})==336

    for a in records:
        p,c=a['rows'],a['columns']
        assert sorted(p)==sorted(c)==list(range(15))
        assert all(
            F[p[i]][j]==F[i][c[j]]
            for i,j in itertools.product(range(15),repeat=2)
        )

    # A small generating set, found deterministically rather than trusted.
    generators=[]
    group=PermutationGroup([Permutation(list(range(15)))])
    for a in records:
        p=Permutation(a['rows'])
        if not group.contains(p):
            generators.append(a)
            group=PermutationGroup(
                [Permutation(g['rows']) for g in generators]
            )
        if group.order()==336:
            break
    assert group.order()==336

    center=group.center()
    derived=group.derived_subgroup()
    assert center.order()==2 and derived.order()==168

    z=next(
        p for p in center.generate_schreier_sims()
        if not p.is_Identity
    )
    assert not derived.contains(z)
    assert [z(i) for i in range(15)]==[i^1 for i in range(14)]+[14]

    pair_group=PermutationGroup([
        Permutation([g['rows'][2*i]//2 for i in range(7)])
        for g in generators
    ])
    assert pair_group.order()==168

    lines={
        frozenset(t)
        for t in [
            (0,1,2),(0,3,4),(0,5,6),
            (1,3,5),(1,4,6),(2,3,6),(2,4,5),
        ]
    }
    assert all(
        sum({i,j}<=line for line in lines)==1
        for i,j in itertools.combinations(range(7),2)
    )
    assert all(
        {frozenset(p(i) for i in line) for line in lines}==lines
        for p in pair_group.generators
    )

    derived_images={
        tuple(p(2*i)//2 for i in range(7))
        for p in derived.generate_schreier_sims()
    }
    assert len(derived_images)==168

    # Independently count all permutations of the seven points preserving these lines.
    fano_aut=sum(
        {frozenset(p[i] for i in line) for line in lines}==lines
        for p in itertools.permutations(range(7))
    )
    assert fano_aut==168

    out={
        'verified':True,
        'squared_determinant':B,
        'determinants_original_dephased_canonical':determinants,
        'row_energy':105,
        'column_energy':sum(
            norm(K[i][j])
            for i,j in itertools.combinations(range(15),2)
        ),
        'offdiagonal_unordered_norm_histogram':dict(histogram),
        'hub_dephased_1based':h+1,
        'matching_dephased_1based':[[i+1,j+1] for i,j in pairs],
        'largest_orthogonal_row_subset':row_alpha,
        'largest_orthogonal_column_subset':col_alpha,
        'number_of_maximum_row_subsets':len(row_sets),
        'number_of_maximum_column_subsets':len(col_sets),
        'one_orthogonal_row_subset_1based':[i+1 for i in row_sets[0]],
        'dephased_row_gram':G,
        'dephased_column_gram':K,
        'row_switch_subtractions':phases,
        'canonical_row_order_1based':[i+1 for i in order],
        'canonical_column_switch_subtractions':E[h],
        'searched_row_permutations':math.factorial(7)*2**7,
        'projective_automorphism_order':336,
        'monomial_pair_automorphism_order':1008,
        'center_projective_order':2,
        'derived_projective_order':168,
        'pair_action_order':168,
        'derived_pair_action_faithful':True,
        'fano_automorphism_order':168,
        'fano_lines_1based':[
            sorted(i+1 for i in line)
            for line in sorted(lines,key=lambda s:tuple(sorted(s)))
        ],
        'projective_group':'C2 x GL(3,2)',
        'monomial_pair_group':'C6 x GL(3,2)',
        'generators_0based':generators,
        'global_maximizer_equivalence_classification_performed':False,
    }

    (ROOT/'data/automorphisms.json').write_text(
        json.dumps(records,indent=2)+'\n'
    )
    (ROOT/'data/structure_verified.json').write_text(
        json.dumps(out,indent=2)+'\n'
    )

    print(
        'PASS: exact determinant and Gram; row/column orthogonality 7; '
        '645120 row candidates; 336 projective / 1008 monomial '
        'automorphisms; Fano action verified.'
    )

if __name__=='__main__':
    main()
