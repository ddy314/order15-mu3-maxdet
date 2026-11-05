#!/usr/bin/env python3
"""Fresh complete catalogue for isolate-free pure Grams with Q <= 99.

No determinant cutoff is used. Sixth-root switching is a spectral relaxation.
Isolate cases are excluded analytically by the sparse flat-kernel lemmas.
"""
from pathlib import Path
import sys,json
from itertools import product,combinations
from collections import defaultdict
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'proof/new'))
from rank_catalogue import characteristic, independence, norm_bases, UNITS, mul, conj, np, nx
import sympy as sp
X=sp.Symbol('x')


def catalogue(budget=11):
    assert budget in (10,11)
    data={};states=0;counts=[]
    atlas=nx.graph_atlas_g()
    for v in range(2,2*budget-14+2):
        cap=budget-(16-v)//2
        if v<=7:graphs=[g for g in atlas if len(g)==v and nx.is_connected(g) and g.number_of_edges()<=cap]
        else:
            assert cap==v-1, (v,cap)
            graphs=list(nx.nonisomorphic_trees(v))
        counts.append(len(graphs))
        for graph in graphs:
            g=nx.convert_node_labels_to_integers(graph)
            edges=sorted(tuple(sorted(e)) for e in g.edges());m=len(edges)
            tree={tuple(sorted(e)) for e in nx.minimum_spanning_tree(g).edges()}
            chords=[k for k,e in enumerate(edges) if e not in tree]
            weights=[(None,1)] + [(k,w) for k in range(m) for w in (3,4) if m+w-1<=cap]
            # Two heavy edges cannot fit the whole 15-vertex budget.
            assert budget-8<4
            for heavy,w in weights:
                for phases in product(range(6),repeat=len(chords)):
                    labels=[(1,0)]*m
                    if heavy is not None:labels[heavy]=norm_bases(w)[0]
                    for k,t in zip(chords,phases):labels[k]=mul(labels[k],UNITS[t])
                    a=np.zeros((v,v),dtype=np.int64);b=a.copy()
                    for (i,j),z in zip(edges,labels):
                        a[i,j],b[i,j]=z;a[j,i],b[j,i]=conj(z)
                    poly=tuple(map(int,characteristic(a,b)));u=m+w-1 if heavy is not None else m
                    alpha=independence(g);key=(v,u,poly)
                    entry={'v':v,'u':u,'poly':poly,'alpha':alpha}
                    if key not in data or alpha<data[key]['alpha']:data[key]=entry
                    states+=1
    if budget==11:
        assert counts==[1,2,5,13,19,44,23,47],counts
        assert states==676,states
        assert len(data)==373,len(data)
    return sorted(data.values(),key=lambda z:(z['v'],z['u'],z['poly'])),counts,states


def poly_mul(p,q):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]+=a*b
    return tuple(r)


def determinant(p,scale=5):
    z=0
    for a in p:z=-scale*z+a
    return z*(-1)**(len(p)-1)


@lru_cache(None)
def obstruction(d):
    return next((int(p) for p,k in sp.factorint(d).items() if p%3==2 and k%2),None)


@lru_cache(None)
def poly_gcd(p,q):
    z=sp.Poly(sp.gcd(sp.Poly.from_list(p,X),sp.Poly.from_list(q,X)),X)
    return tuple(map(int,z.all_coeffs()))


def assemblies(records,budget,alpha_cap=9,scale=5):
    out=[]
    def visit(start,v,u,alpha,ids,p):
        if v==15:
            if u==budget and alpha<=alpha_cap:
                d=determinant(p,scale)
                out.append({'ids':ids,'poly':p,'det':d,'norm_obstruction':obstruction(d)})
            return
        for k in range(start,len(records)):
            r=records[k]
            if v+r['v']>15 or u+r['u']>budget or alpha+r['alpha']>alpha_cap:continue
            visit(k,v+r['v'],u+r['u'],alpha+r['alpha'],ids+[k],poly_mul(p,r['poly']))
    visit(0,0,0,0,[],(1,))
    return out


def unions(a,records):
    out={}
    for bits in range(1,1<<len(a['ids'])):
        ids=tuple(a['ids'][i] for i in range(len(a['ids'])) if bits>>i&1)
        if ids in out:continue
        p=(1,);v=0
        for k in ids:p=poly_mul(p,records[k]['poly']);v+=records[k]['v']
        out[ids]=(v,p)
    return list(out.values())


def rectangle(a,b,records,scale=3):
    for v,p in unions(a,records):
        for w,q in unions(b,records):
            g=poly_gcd(p,q);degree=len(g)-1
            cap=15*degree-(scale*g[1] if degree else 0)
            if v*w>cap:return {'left':p,'right':q,'gcd':g,'area':v*w,'capacity':cap}
    return None


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled.')
    records,counts,states=catalogue();tables=[];certificates=[];survivors=[]
    for u in range(8,12):
        aa=assemblies(records,u);nn=[a for a in aa if not a['norm_obstruction']];pairs=0;sur=[]
        for i,a in enumerate(nn):
            for b in nn[i:]:
                if a['poly']!=b['poly']:continue
                pairs+=1;r=rectangle(a,b,records)
                if r:certificates.append({'u':u,'a':a['ids'],'b':b['ids'],**r})
                else:sur.append((a,b))
        tables.append([u,len(aa),len(nn),pairs,len(sur)])
        survivors.extend({'u':u,'left':a,'right':b} for a,b in sur)
    expected=[[8,1,0,0,0],[9,10,4,4,1],[10,67,2,2,0],[11,581,53,63,2]]
    assert tables==expected,(tables,expected)
    p=poly_mul(poly_mul((1,0),(1,0,-5)),(1,0,-1))
    # Last two pairs must have x(x^2-1)^6(x^2-5).
    last=(1,0)
    for _ in range(6):last=poly_mul(last,(1,0,-1))
    last=poly_mul(last,(1,0,-5))
    assert all(s['left']['poly']==last for s in survivors if s['u']==11)
    out={'component_graph_counts':counts,'labelled_states':states,'records':records,'tables':tables,'rectangle_certificates':certificates,'survivors':survivors,
         'analytic_closure':'u=9: triangle block argument; u=11: signed adjugate trace five forces an Eisenstein norm of fifteen.'}
    path=ROOT/'data/revision/pure_energy_certificate.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: 676 labelled states, 373 component spectra; pure-energy tables',tables,flush=True)

if __name__=='__main__':main()
