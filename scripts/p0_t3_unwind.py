"""Select maximum-recovery T2 alternatives with exact noncrossing geometry; unwind."""
import csv, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import p0_t2_recovery as m
from p0_t2_validate import load_graph, intersections

ROOT=m.ROOT;OUT=ROOT/'outputs/p0_t3/run_02';T2=ROOT/'outputs/p0_t2/run_02'
PLAN={'scope':'Same five independent XY slices; select highest T2-score complete path without proper crossing or repeated nonadjacent vertex',
      'objective':'unchanged T2 maximum-recovery score; hard topology feasibility filter, no percentage gate',
      'sampling_px':.5,'spacing_um':126.56,'preserve':'T2 original runs and report; all rejected candidate evaluations',
      'search':'evaluate existing score-sorted T2 complete archive, stop at first topology-feasible path; if none, disclose and expand search separately',
      'sampling':'CT only for provenance1/3 and response>3; bridge geometry retained, CT NaN; no downstream masking'}

def cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]

class ConflictIndex:
    """Exact bridge-vs-observed line intersection cache, independent of raster occupancy."""
    def __init__(self,ps,original,edges,seeds):
        self.ps=ps;self.original=original;self.edges=edges;self.cache={};self.seed_map=dict(seeds)
        aa=[];bb=[];owner=[]
        for i in range(0,len(ps),2):
            p=ps[i];aa.extend(p[:-1]);bb.extend(p[1:]);owner.extend([original[i]]*(len(p)-1))
        self.a=np.array(aa);self.b=np.array(bb);self.owner=np.array(owner);self.tree=cKDTree((self.a+self.b)/2)
    def bridge_conflicts(self,eid):
        if eid in self.cache:return self.cache[eid]
        e=self.edges[eid];a=self.ps[int(e[0])][-1];b=self.ps[int(e[1])][0];v=b-a
        if not np.linalg.norm(v):self.cache[eid]=set();return set()
        idx=np.array(self.tree.query_ball_point((a+b)/2,np.linalg.norm(v)/2+1),int)
        w=self.b[idx]-self.a[idx];den=cross(v,w);good=np.abs(den)>1e-10;idx=idx[good];w=w[good];den=den[good]
        t=cross(self.a[idx]-a,w)/den;u=cross(self.a[idx]-a,v)/den
        conflict=set(self.owner[idx[(t>1e-8)&(t<1-1e-8)&(u>1e-8)&(u<1-1e-8)]].tolist())
        self.cache[eid]=conflict;return conflict
    def quick_conflict(self,path,eids):
        members={self.original[i] for i in path};sid=self.original[path[0]]
        for eid in eids:
            hits=self.bridge_conflicts(eid)&members
            # Seed chain is truncated; defer that rare case to full exact geometry.
            if hits-{sid}:return True
        es=[self.edges[eid] for eid in eids if self.edges[eid,2]>0]
        if not es:return False
        a=np.array([self.ps[int(e[0])][-1] for e in es]);b=np.array([self.ps[int(e[1])][0] for e in es]);v=b-a
        i,j=np.triu_indices(len(es),1);den=cross(v[i],v[j]);good=np.abs(den)>1e-10;i=i[good];j=j[good];den=den[good]
        t=cross(a[j]-a[i],v[j])/den;u=cross(a[j]-a[i],v[i])/den
        return bool(np.any((t>1e-8)&(t<1-1e-8)&(u>1e-8)&(u<1-1e-8)))

def geometry(path,eids,first,ps,edges):
    seg=[first]
    for k,j in enumerate(path[1:]):
        prev=path[k]
        if edges[eids[k],2]>0:seg.append(np.array([ps[prev][-1],ps[j][0]]))
        seg.append(ps[j])
    return seg

def repeated_vertices(seg):
    p=np.concatenate([seg[0]]+[s[1:] for s in seg[1:]])
    p=p[np.r_[True,np.any(np.diff(p,axis=0)!=0,axis=1)]]
    _,counts=np.unique(p,axis=0,return_counts=True)
    return int(np.sum(counts-1))

def main():
    # Enforce documented T2 entry condition rather than merely trusting its presence.
    assert (ROOT/'docs/p0_t2_maximum_recovery_thread.md').is_file()
    for line in (ROOT/'outputs/p0_t2/SHA256SUMS').read_text().splitlines():
        h,p=line.split('  ',1);assert m.sha(ROOT/p)==h,p
    OUT.mkdir(exist_ok=False)
    prior=[p for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents]
    prior += [ROOT/'README.md',ROOT/'AGENTS.md']
    m.dump(OUT/'preservation.json',{str(p.relative_to(ROOT)):m.sha(p) for p in prior});m.dump(OUT/'plan.json',PLAN)
    volume=np.load(ROOT/m.PLAN['source'],mmap_mode='r');baseline={r['z']:r for r in json.loads((ROOT/'outputs/p0_t2/summary.json').read_text())};metrics=[]
    for z in m.ZS:
        dest=OUT/f'slice_{z}';dest.mkdir()
        ps,original,L,adj,edges,seeds,tout=load_graph(z);seedmap=dict(seeds);conflicts=ConflictIndex(ps,original,edges,seeds)
        archive=json.loads((T2/f'slice_{z}_complete_paths.json').read_text());chosen=None;log=[]
        for rank,c in enumerate(archive,1):
            path=c['path'];eids=c['edge_ids'];first=seedmap[path[0]]
            if conflicts.quick_conflict(path,eids):
                log.append({'rank':rank,'score':c['score'],'result':'proper_crossing_cached'});continue
            seg=geometry(path,eids,first,ps,edges);crossings=intersections(seg);touches=repeated_vertices(seg)
            log.append({'rank':rank,'score':c['score'],'proper_crossings':len(crossings),'repeated_vertices':touches,'result':'accept' if not crossings and not touches else 'topology_reject'})
            if not crossings and not touches:chosen=c;break
        m.dump(dest/'candidate_evaluations.json',log)
        if chosen is None:
            m.dump(dest/'unresolved.json',{'reason':'no feasible complete archived path; no shortened substitute','candidates':len(archive)})
            print(f'Z{z}: no feasible archive alternative',flush=True);continue
        path=chosen['path'];eids=chosen['edge_ids'];first=seedmap[path[0]]
        response=np.load(ROOT/f'outputs/p0_t1/slice_{z}_candidates.npz')['response']
        state=(chosen['score'],tuple(path),tuple(eids),frozenset(original[i] for i in path),frozenset(),first)
        seg,labs,tr,ss,xy,prov,lengths=m.materialize(state,ps,edges,adj,response,z)
        row=np.full(len(ss),np.nan);valid=np.isin(prov,[1,3]);row[valid]=m.sample(volume[z-128],xy[valid]);mask=~valid
        np.save(dest/'row.npy',row);np.save(dest/'mask.npy',mask);np.save(dest/'provenance.npy',prov)
        np.savez_compressed(dest/'path_and_samples.npz',s_px=ss,s_um=ss*126.56,xy=xy,ct=row,mask=mask,provenance=prov,
                            vertex_xy=np.concatenate(seg),segment_offsets=np.r_[0,np.cumsum([len(p) for p in seg])],segment_provenance=labs)
        m.dump(dest/'selected_path.json',{**chosen,'archive_rank':rank,'seed_piece':first.tolist(),'transitions':tr})
        with (dest/'path_coordinates.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['segment','vertex_in_segment','x_px','y_px','s_px','s_um','provenance']);start=0.
            for k,(p,label) in enumerate(zip(seg,labs)):
                a=m.arc(p)
                for j,(q,s) in enumerate(zip(p,a)):w.writerow([k,j,*q,start+s,(start+s)*126.56,label])
                start+=a[-1]
        with (dest/'sample_provenance.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['sample','s_px','s_um','x_px','y_px','provenance','masked','ct'])
            for i in range(len(ss)):w.writerow([i,ss[i],ss[i]*126.56,*xy[i],int(prov[i]),int(mask[i]),'' if mask[i] else row[i]])
        Image.fromarray(np.nan_to_num(row,nan=0).round().astype('uint8')[None,:]).save(dest/'unwrapped_row.png')
        Image.fromarray(mask.astype('uint8')[None,:]*255).save(dest/'unwrapped_mask.png')
        for provenance in [False,True]:
            fig,ax=plt.subplots(figsize=(12,9));ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255)
            ax.add_collection(LineCollection(seg,colors=[m.COLORS[l] for l in labs] if provenance else 'cyan',linewidths=1.1))
            ax.scatter(*seg[0][0],c='lime',marker='*',s=90);ax.scatter(*seg[-1][-1],c='lime',marker='o',s=40)
            ax.add_patch(Rectangle((430,150),30,60,fill=False,edgecolor='lime'))
            ax.set(xlim=(150,475),ylim=(320,25),title=f'Z{z}: topology-feasible T2 alternative rank{rank}; continuous physical s'+('\ncyan observed / magenta inferred / orange low confidence' if provenance else ''))
            fig.tight_layout();fig.savefig(dest/('provenance_overlay.png' if provenance else 'final_thread_overlay.png'),dpi=170);plt.close(fig)
        # Contract checks on actual arrays and continuous geometry.
        assert all(np.array_equal(p[-1],q[0]) for p,q in zip(seg,seg[1:]))
        assert not intersections(seg) and not repeated_vertices(seg)
        assert np.isnan(row[np.isin(prov,[0,2])]).all() and np.isfinite(row[valid]).all()
        assert np.allclose(np.diff(ss*126.56),63.28,atol=1e-8)
        labs=np.array(labs);obs=float(lengths[labs!=2].sum());bridge=float(lengths[labs==2].sum())
        b=baseline[z];r={'z':z,'T2_archive_rank':rank,'candidates_evaluated':len(log),'score':chosen['score'],
            'total_arc_px':float(lengths.sum()),'total_arc_mm':float(lengths.sum()*.12656),'observed_arc_px':obs,'observed_arc_mm':obs*.12656,
            'bridge_arc_px':bridge,'bridge_arc_mm':bridge*.12656,'observed_graph_percent':100*obs/lengths.sum(),
            'observed_retained_vs_T2_percent':100*obs/b['observed_graph_arc_px'],'T2_observed_arc_px':b['observed_graph_arc_px'],
            'crossings_before':b['proper_self_crossings'],'crossings_after':0,'repeated_nonadjacent_vertices':0,
            'total_samples':len(ss),'observed_samples':int(valid.sum()),'masked_bridge_samples':int((prov==2).sum()),
            'unresolved_samples':int((prov==0).sum()),'low_confidence_samples':int((prov==3).sum()),
            'fragments':len(path),'bridges':int((labs==2).sum()),'exterior_reached':bool(m.target_dist(seg[-1][-1])==0),
            'seed_xy':seg[0][0].tolist(),'end_xy':seg[-1][-1].tolist(),'T2_observed_samples':b['usable_ct_sample_count']}
        metrics.append(r);m.dump(OUT/'metrics.json',metrics);print(json.dumps(r),flush=True)

if __name__=='__main__':main()
