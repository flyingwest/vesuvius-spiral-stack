"""Independent geometry/provenance checks and one exact bounded-search replay."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import p0_t2_recovery as m

ROOT=m.ROOT;OUT=ROOT/'outputs/p0_t2';RUN=OUT/'run_02'

def load_graph(z):
    g=json.loads((ROOT/f'outputs/p0_t1/slice_{z}_graph.json').read_text())
    original=json.loads((RUN/f'slice_{z}_oriented_chain_map.json').read_text())
    ps=[np.array(g['chains'][cid])[::(-1 if i%2 else 1)] for i,cid in enumerate(original)]
    edges=np.load(RUN/f'slice_{z}_virtual_edges.npz')['edges'];adj=[[] for _ in ps]
    for eid,e in enumerate(edges):
        if e[13]:adj[int(e[0])].append((float(e[12]),int(e[1]),eid))
    adj=[sorted(a) for a in adj]
    L=np.array([m.arc(p)[-1] for p in ps]);tout=m.unit(np.array([p[-1]-p[max(0,len(p)-6)] for p in ps]))
    seed=np.array(g['pieces'][0][0]);sid=g['selected_chains'][0];seeds=[]
    for i,cid in enumerate(original):
        if cid==sid:
            p=ps[i];j=int(np.argmin(np.linalg.norm(p-seed,axis=1)))
            if len(p[j:])>1:seeds.append((i,p[j:]))
    return ps,original,L,adj,edges,seeds,tout

def intersections(segments):
    a=np.concatenate([p[:-1] for p in segments]);b=np.concatenate([p[1:] for p in segments]);v=b-a
    pairs=cKDTree((a+b)/2).query_pairs(12.01,output_type='ndarray')
    if not len(pairs):return []
    i,j=pairs.T;keep=(j-i)>1;i=i[keep];j=j[keep]
    def cross(x,y):return x[:,0]*y[:,1]-x[:,1]*y[:,0]
    den=cross(v[i],v[j]);good=np.abs(den)>1e-9;i=i[good];j=j[good];den=den[good]
    t=cross(a[j]-a[i],v[j])/den;u=cross(a[j]-a[i],v[i])/den
    keep=(t>1e-8)&(t<1-1e-8)&(u>1e-8)&(u<1-1e-8)
    return (a[i[keep]]+t[keep,None]*v[i[keep]]).tolist()

def main():
    volume=np.load(ROOT/m.PLAN['source'],mmap_mode='r');metrics=json.loads((RUN/'metrics.json').read_text());checks=[]
    for z in reversed(m.ZS):
        ps,original,L,adj,edges,seeds,tout=load_graph(z)
        top=json.loads((RUN/f'slice_{z}_top_paths.json').read_text());response=np.load(ROOT/f'outputs/p0_t1/slice_{z}_candidates.npz')['response']
        a=np.load(RUN/f'slice_{z}_samples.npz');candidate_checks=[]
        for cand in top:
            path=cand['path'];eids=cand['edge_ids'];first=np.array(cand['seed_piece'])
            assert len(set(original[i] for i in path))==len(path)
            assert float(m.target_dist(ps[path[-1]][-1]))==0
            score=m.arc(first)[-1]+sum(L[i] for i in path[1:])-sum(edges[e,12] for e in eids)+1.5*float(m.target_dist(first[0]))
            assert np.isclose(score,cand['score'],atol=1e-8,rtol=0)
            for i,eid in enumerate(eids):
                e=edges[eid];assert e[13] and e[0]==path[i] and e[1]==path[i+1] and e[2]<=12
            state=(cand['score'],tuple(path),tuple(eids),frozenset(original[i] for i in path),frozenset(),first)
            seg,labs,tr,ss,xy,prov,lens=m.materialize(state,ps,edges,adj,response,z)
            assert all(np.array_equal(p[-1],q[0]) for p,q in zip(seg,seg[1:]))
            crosses=intersections(seg)
            candidate_checks.append({'rank':cand['rank'],'continuous':True,'proper_self_crossings':len(crosses),'crossing_xy':crosses,'score_recomputed':float(score)})
            if cand['rank']==1:
                assert np.array_equal(xy,a['xy']) and np.array_equal(prov,a['provenance'])
                assert np.allclose(np.diff(a['s_um']),63.28,atol=1e-8)
                valid=np.isin(prov,[1,3]);assert np.isfinite(a['ct'][valid]).all() and np.isnan(a['ct'][~valid]).all()
                q=xy[valid];x=np.floor(q[:,0]).astype(int);y=np.floor(q[:,1]).astype(int);dx=q[:,0]-x;dy=q[:,1]-y;im=volume[z-128].astype(float)
                direct=im[y,x]*(1-dx)*(1-dy)+im[y,x+1]*dx*(1-dy)+im[y+1,x]*(1-dx)*dy+im[y+1,x+1]*dx*dy
                assert np.allclose(direct,a['ct'][valid],rtol=0,atol=1e-10)
                # Save auditable panels for the longest/highest penalty bridges, not flattering selections.
                selected=sorted([eid for eid in eids if edges[eid,2]>0],key=lambda eid:-edges[eid,12])[:6]
                fig,axs=plt.subplots(2,3,figsize=(14,8))
                for ax,eid in zip(axs.flat,selected):
                    e=edges[eid];p=ps[int(e[0])];q=ps[int(e[1])];mid=(p[-1]+q[0])/2
                    ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
                    ax.plot(*p.T,'c-',lw=1);ax.plot(*q.T,'c-',lw=1);ax.plot(*np.array([p[-1],q[0]]).T,'m--',lw=2)
                    ax.set(xlim=(mid[0]-12,mid[0]+12),ylim=(mid[1]+12,mid[1]-12),title=f'edge{eid}: gap{e[2]:.2f}px, penalty{e[12]:.1f}')
                fig.tight_layout();fig.savefig(OUT/f'slice_{z}_bridge_audit.png',dpi=140);plt.close(fig)
        for f in RUN.glob(f'slice_{z}_*.png'):
            with Image.open(f) as im:im.verify()
        old=np.load(ROOT/f'outputs/p0_t1/slice_{z}_samples.npz');oldmask=np.isfinite(old['ct'])
        checks.append({'z':z,'candidate_checks':candidate_checks,'T1_raw_valid_samples':int(oldmask.sum()),
                       'T2_valid_samples':int(np.isfinite(a['ct']).sum()),'sample_utilization_gain_vs_T1':float(np.isfinite(a['ct']).sum()/oldmask.sum())})
        print(f'Z{z}: geometry/provenance/sampling PASS; crossings {[c["proper_self_crossings"] for c in candidate_checks]}',flush=True)
    # One independent reverse-order experiment replay; other slices checked against all saved graph edges and samples.
    z=128;ps,original,L,adj,edges,seeds,tout=load_graph(z)
    top,complete,log,counts=m.search(ps,original,L,adj,edges,seeds,tout)
    saved=json.loads((RUN/f'slice_{z}_top_paths.json').read_text())
    assert [list(s[1]) for s in top]==[s['path'] for s in saved]
    assert [s[0] for s in top]==[s['score'] for s in saved]
    assert len(complete)==next(r['complete_count'] for r in metrics if r['z']==z)
    for run in ['run_01','run_02']:
        for name,h in json.loads((OUT/run/'preservation.json').read_text()).items():assert m.sha(ROOT/name)==h,name
    m.dump(OUT/'verification.json',{'slices':checks,'search_replay_Z128':'exact paths/scores/complete count PASS','prior_artifacts':'all preserved hashes PASS',
                                  'note':'proper crossings exclude endpoint touches and collinear overlap; no oracle-based identity check'})
    print('PASS: exact search replay and all preservation hashes',flush=True)

if __name__=='__main__':main()
