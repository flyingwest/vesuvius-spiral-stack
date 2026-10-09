"""Recovery-preserving geometric two-edge uncrossing of the selected T2 paths."""
import csv, json
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import p0_t2_recovery as m
from p0_t3_local_topology import flatten,crossings
from p0_t3_unwind import repeated_vertices

ROOT=m.ROOT;OUT=ROOT/'outputs/p0_t3/run_03'
PLAN={'method':'two-edge uncrossing: replace crossing A->B,C->D by A->C,B->D; reverse intervening edges; preserve their provenance',
      'selection':'at each step choose strict decrease in proper crossing count, then minimum lost observed length, then minimum new bridge length; max20 corrections',
      'scope':'only four original T2 crossings; no detection, ideal spiral, new CT, Z tracking or downstream masking',
      'gap':'both replacement edges labeled2 and unsampled even if they touch bright CT',
      'uncertainty':'observed edges within5px along path on either side of each new bridge relabeled3 locally',
      'arc_sampling_px':.5,'spacing_um':126.56,
      'amendment':'Existing-archive alternatives lose31-48% observed length on three slices; loop diagnosis shows96-682px intervening arc; preserve both findings, apply local two-edge exchange instead of loop deletion'}

def change(a,b,l,origin,i,j):
    aa=np.concatenate([a[:i],a[i:i+1],b[i+1:j][::-1],b[i:i+1],a[j+1:]])
    bb=np.concatenate([b[:i],a[j:j+1],a[i+1:j][::-1],b[j:j+1],b[j+1:]])
    ll=np.r_[l[:i],2,l[i+1:j][::-1],2,l[j+1:]].astype(np.uint8)
    oo=np.r_[origin[:i],-1,origin[i+1:j][::-1],-1,origin[j+1:]]
    return aa,bb,ll,oo

def main():
    OUT.mkdir(exist_ok=False)
    # Whole prior evidence, including the unsatisfactory archive-only choices, remains immutable.
    previous=[p for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents]
    previous += [ROOT/'README.md',ROOT/'AGENTS.md']
    m.dump(OUT/'preservation.json',{str(p.relative_to(ROOT)):m.sha(p) for p in previous});m.dump(OUT/'plan.json',PLAN)
    volume=np.load(ROOT/m.PLAN['source'],mmap_mode='r');baseline={r['z']:r for r in json.loads((ROOT/'outputs/p0_t2/summary.json').read_text())};metrics=[]
    for z in m.ZS:
        dest=OUT/f'slice_{z}';dest.mkdir();source=np.load(ROOT/f'outputs/p0_t2/slice_{z}_samples.npz')
        a,b,l=flatten(source);origin=np.arange(len(a));initial_cross=crossings(a,b);audit=[]
        for step in range(20):
            current=crossings(a,b)
            if not current:break
            choices=[]
            for c in current:
                i,j=c['i'],c['j'];new=change(a,b,l,origin,i,j);aa,bb,ll,oo=new
                after=crossings(aa,bb);lost=sum(np.linalg.norm(b[k]-a[k]) for k in [i,j] if l[k]!=2)
                newbridge=float(np.linalg.norm(aa[i]-bb[i])+np.linalg.norm(aa[j]-bb[j]))
                choice={'crossing':c,'before_count':len(current),'after_count':len(after),'lost_observed_px':float(lost),'new_bridge_length_px':newbridge,
                        'old_edges_xy':[[a[i].tolist(),b[i].tolist()],[a[j].tolist(),b[j].tolist()]],'new_edges_xy':[[aa[i].tolist(),bb[i].tolist()],[aa[j].tolist(),bb[j].tolist()]],
                        'old_labels':[int(l[i]),int(l[j])],'reversed_original_edge_ids':origin[i+1:j].tolist()}
                if len(after)<len(current):choices.append((len(after),lost,newbridge,i,j,new,choice))
            if not choices:raise RuntimeError(f'No local crossing-decreasing exchange for {z}; preserve failure')
            best=min(choices,key=lambda q:q[:5]);a,b,l,origin=best[5];audit.append(best[6])
        assert not crossings(a,b)
        base=l.copy();length=np.linalg.norm(b-a,axis=1);s=np.r_[0,np.cumsum(length)]
        # Mark local topology-changing connection neighborhoods, never the entire suffix.
        for k in np.flatnonzero(origin<0):
            near=((s[:-1]>=s[k]-5)&(s[:-1]<=s[k+1]+5))|((s[1:]>=s[k]-5)&(s[1:]<=s[k+1]+5))
            l[near&np.isin(l,[1,3])]=3
        ss=np.arange(0,s[-1]+1e-8,.5);idx=np.minimum(np.searchsorted(s,ss,side='right')-1,len(a)-1)
        fraction=(ss-s[idx])/length[idx];xy=a[idx]+fraction[:,None]*(b[idx]-a[idx]);prov=l[idx].copy()
        response=np.load(ROOT/f'outputs/p0_t1/slice_{z}_candidates.npz')['response'];support=m.sample(response,xy)>3
        prov[np.isin(prov,[1,3])&~support]=0;mask=~np.isin(prov,[1,3]);row=np.full(len(ss),np.nan)
        row[~mask]=m.sample(volume[z-128],xy[~mask])
        # At exact bridge/observed boundaries, preserve T2's downstream-edge convention.
        np.save(dest/'row.npy',row);np.save(dest/'mask.npy',mask);np.save(dest/'provenance.npy',prov)
        np.savez_compressed(dest/'path_and_samples.npz',edge_start_xy=a,edge_end_xy=b,edge_provenance=l,edge_original_provenance=base,
                            original_T2_edge_id=origin,edge_s_px=s,s_px=ss,s_um=ss*126.56,xy=xy,ct=row,mask=mask,provenance=prov)
        m.dump(dest/'corrections.json',audit)
        with (dest/'path_coordinates.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['edge','start_x_px','start_y_px','end_x_px','end_y_px','start_s_um','end_s_um','provenance','original_T2_edge_id'])
            for k in range(len(a)):w.writerow([k,*a[k],*b[k],s[k]*126.56,s[k+1]*126.56,int(l[k]),int(origin[k])])
        with (dest/'sample_provenance.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['sample','s_um','x_px','y_px','provenance','masked','ct'])
            for k in range(len(ss)):w.writerow([k,ss[k]*126.56,*xy[k],int(prov[k]),int(mask[k]),'' if mask[k] else row[k]])
        Image.fromarray(np.nan_to_num(row,nan=0).round().astype('uint8')[None,:]).save(dest/'unwrapped_row.png')
        Image.fromarray(mask.astype('uint8')[None,:]*255).save(dest/'unwrapped_mask.png')
        for colored in [False,True]:
            fig,ax=plt.subplots(figsize=(12,9));ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255)
            ax.add_collection(LineCollection(np.stack([a,b],axis=1),colors=[m.COLORS[c] for c in l] if colored else 'cyan',linewidths=1))
            ax.scatter(*a[0],c='lime',marker='*',s=90);ax.scatter(*b[-1],c='lime',s=35)
            ax.add_patch(Rectangle((430,150),30,60,fill=False,edgecolor='lime'))
            for c in audit:ax.scatter(*c['crossing']['xy'],s=140,facecolors='none',edgecolors='red',lw=1)
            ax.set(xlim=(150,475),ylim=(320,25),title=f'Z{z}: continuous uncrossed thread; red circles = corrected sites'+('\ncyan observed / magenta inferred / orange local low confidence' if colored else ''))
            fig.tight_layout();fig.savefig(dest/('provenance_overlay.png' if colored else 'final_thread_overlay.png'),dpi=170);plt.close(fig)
        # Before/after panels at every repaired site, including raw CT background.
        if audit:
            fig,axs=plt.subplots(len(audit),2,figsize=(10,5*len(audit)),squeeze=False)
            olda,oldb,oldl=flatten(source)
            for c,axes in zip(audit,axs):
                center=np.array(c['crossing']['xy'])
                for ax,aa,bb,ll,label in zip(axes,[olda,a],[oldb,b],[oldl,l],['T2 before','T3 after']):
                    ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
                    ax.add_collection(LineCollection(np.stack([aa,bb],axis=1),colors=[m.COLORS[v] for v in ll],linewidths=1.3))
                    ax.set(xlim=(center[0]-10,center[0]+10),ylim=(center[1]+10,center[1]-10),title=label)
            fig.tight_layout();fig.savefig(dest/'topology_corrections.png',dpi=160);plt.close(fig)
        bline=baseline[z];observed=float(length[l!=2].sum());bridges=float(length[l==2].sum())
        segments=np.stack([a,b],axis=1);assert np.allclose(b[:-1],a[1:],rtol=0,atol=1e-10)
        assert not repeated_vertices(segments)
        record={'z':z,'selection':'original T2 best with local 2-edge uncrossing' if audit else 'original T2 best unchanged',
                'total_arc_px':float(length.sum()),'total_arc_mm':float(length.sum()*.12656),'observed_arc_px':observed,'observed_arc_mm':observed*.12656,
                'bridge_arc_px':bridges,'bridge_arc_mm':bridges*.12656,'observed_graph_percent':100*observed/length.sum(),
                'T2_observed_arc_px':bline['observed_graph_arc_px'],'observed_retained_vs_T2_percent':100*observed/bline['observed_graph_arc_px'],
                'lost_observed_arc_px':bline['observed_graph_arc_px']-observed,'crossings_before':len(initial_cross),'crossings_after':0,'corrections':len(audit),
                'total_samples':len(ss),'observed_samples':int((~mask).sum()),'masked_bridge_samples':int((prov==2).sum()),'unresolved_samples':int((prov==0).sum()),
                'low_confidence_samples':int((prov==3).sum()),'fragments_incorporated':bline['fragments'],'inferred_edges':int((l==2).sum()),
                'exterior_reached':bool(m.target_dist(b[-1])==0),'seed_xy':a[0].tolist(),'end_xy':b[-1].tolist(),'sampling_um':63.28,
                'nonadjacent_repeated_vertices':0,'T2_observed_samples':bline['usable_ct_sample_count']}
        metrics.append(record);m.dump(OUT/'metrics.json',metrics);print(json.dumps(record),flush=True)

if __name__=='__main__':main()
