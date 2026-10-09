"""Independent T2/T3 reconstruction, immutable row barrier, then T5 evaluation."""
import argparse, csv, hashlib, json, os, platform, sys, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import p0_t1_threads as t1
import p0_t2_recovery as t2
import p0_t2_papyrus_domain as domain
import p0_t3_contacts as contacts
import p0_t3_uncross as uncross
from p0_t3_local_topology import flatten
from p0_t3_verify import exhaustive_topology

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/p0_t5'; ZS=list(range(128,178)); SPACING=126.56; DS=63.28
SOURCE=ROOT/'outputs/orientation/ct_level4_zyx.npy'
def dump(p,x): p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def table(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def initialize():
    OUT.mkdir(exist_ok=False)
    prior=[p for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents]
    prior += [ROOT/'README.md',ROOT/'AGENTS.md',ROOT/'spiral.png',SOURCE]
    dump(OUT/'preservation.json',{str(p.relative_to(ROOT)):sha(p) for p in prior})
    for seal in ['outputs/p0_t3/SHA256SUMS','outputs/p0_t4/SHA256SUMS']:
        for line in (ROOT/seal).read_text().splitlines():
            h,p=line.split('  ',1);assert sha(ROOT/p)==h,p
    dump(OUT/'plan.json',{'slices':ZS,'spacing_um':SPACING,'ds_um':DS,'source':str(SOURCE.relative_to(ROOT)),
        'source_sha256':sha(SOURCE),'T1':t1.PLAN,'T2':t2.PLAN,'domain_polygon':domain.POLYGON,
        'T3':uncross.PLAN,'T3_contact_amendment':'inclusive nonparallel contact exchange, same final run04 method',
        'scope_override':'50 independent slices; no restriction to original five or four defects',
        'parameters':'unchanged graph, seed, target, bridge, beam, depth, scoring, max20 exchanges; historical Z146/Z160 low-confidence zones retained',
        'failure_policy':'preserve capped/unresolvable topology as negative evidence; no retuning; retain row and report defects',
        'independence':'each worker sees only its XY plane and fixed constants; no other reconstructed slice supplied; comparisons only after 50 final receipts and row hashes',
        'coordinate':'unchanged seed-relative s, no shifts/stretch/reversal during stacking',
        'evaluation':'all49 raw/difference Pearson; 256-sample windows stride128 min64 valid; 199 circular-shift controls; lag ±64 diagnostic only',
        'environment':{'python':sys.version,'numpy':np.__version__,'platform':platform.platform()}})

def reconstruct(z):
    dest=OUT/f'slice_{z}';dest.mkdir(exist_ok=False)
    # This compatibility root lets the unmodified T2 builder load this slice's newly
    # extracted T1 graph. It contains no neighboring-slice graph or CT.
    stage=dest/'graph_input';inp=stage/'outputs/p0_t1';inp.mkdir(parents=True)
    image=np.asarray(np.load(SOURCE,mmap_mode='r')[z-128]).copy()
    response,sk,chains,pieces,ids,gaps,decisions=t1.process_slice(image)
    np.savez_compressed(inp/f'slice_{z}_candidates.npz',response=response,skeleton=sk)
    dump(inp/f'slice_{z}_graph.json',{'chains':[p.tolist() for p in chains],'selected_chains':ids,'pieces':[p.tolist() for p in pieces],'gap_lengths_px':gaps,'decisions':decisions})
    Image.fromarray(image).save(dest/'ct.png')
    t2.ROOT=stage;t2.OUT=dest
    ps,orig,L,adj,edges,seeds,response,sk,tout=domain.build(z)
    np.savez_compressed(dest/'virtual_edges.npz',edges=edges)
    dump(dest/'oriented_chain_map.json',orig)
    dump(dest/'edge_columns.json',['source_oriented_chain','target_oriented_chain','distance','tangent_cos','out_chord_cos','in_chord_cos','normal_miss','spacing','turn_radians','required_bend_per_px','ridge_support','crossing_proxy','penalty','admissible'])
    top,complete,log,counts=t2.search(ps,orig,L,adj,edges,seeds,tout)
    dump(dest/'search_log.json',log)
    dump(dest/'complete_paths.json',[{'score':s[0],'path':s[1],'edge_ids':s[2]} for s in sorted(complete,key=lambda s:-s[0])])
    dump(dest/'top_paths.json',[{'rank':i+1,'score':s[0],'path':s[1],'edge_ids':s[2],'seed_piece':s[5].tolist()} for i,s in enumerate(top)])
    if not top:raise RuntimeError('No complete or partial candidate; preserved graph and search')
    segs,labs,trans,ss,xy,prov,lengths=t2.materialize(top[0],ps,edges,adj,response,z)
    src={'vertex_xy':np.concatenate(segs),'segment_offsets':np.r_[0,np.cumsum([len(p) for p in segs])],'segment_provenance':np.array(labs)}
    np.savez_compressed(dest/'precleanup_path.npz',**src)
    dump(dest/'transitions.json',trans)
    a,b,l=flatten(src);oa=a.copy();ob=b.copy();ol=l.copy();origin=np.arange(len(a));audit=[]
    before=exhaustive_topology(a,b);beforeobs=float(np.linalg.norm(b-a,axis=1)[l!=2].sum());termination='clean'
    for step in range(20):
        current=contacts.intersections_inclusive(a,b)
        if not current:break
        choices=[]
        for c in current:
            i,j=c['i'],c['j'];aa,bb,ll,oo=uncross.change(a,b,l,origin,i,j)
            after=contacts.intersections_inclusive(aa,bb)
            lost=sum(np.linalg.norm(b[k]-a[k]) for k in [i,j] if l[k]!=2)
            nb=float(np.linalg.norm(aa[i]-bb[i])+np.linalg.norm(aa[j]-bb[j]))
            rec={'crossing':c,'before_count':len(current),'after_count':len(after),'lost_observed_px':float(lost),'new_bridge_length_px':nb,
                 'old_edges_xy':[[a[i].tolist(),b[i].tolist()],[a[j].tolist(),b[j].tolist()]],'new_edges_xy':[[aa[i].tolist(),bb[i].tolist()],[aa[j].tolist(),bb[j].tolist()]],
                 'old_labels':[int(l[i]),int(l[j])],'reversed_original_edge_ids':origin[i+1:j].tolist()}
            if len(after)<len(current):choices.append((len(after),lost,nb,i,j,(aa,bb,ll,oo),rec))
        if not choices:termination='no crossing-decreasing exchange';break
        best=min(choices,key=lambda q:q[:5]);a,b,l,origin=best[5];audit.append(best[6])
    else:termination='20-exchange cap'
    # Zero-length geometric edges have no arc or signal; omit only these to avoid 0/0.
    keep=np.linalg.norm(b-a,axis=1)>1e-9;zero_removed=int((~keep).sum());a,b,l,origin=a[keep],b[keep],l[keep],origin[keep]
    base=l.copy();length=np.linalg.norm(b-a,axis=1);s=np.r_[0,np.cumsum(length)]
    for k in np.flatnonzero(origin<0):
        near=((s[:-1]>=s[k]-5)&(s[:-1]<=s[k+1]+5))|((s[1:]>=s[k]-5)&(s[1:]<=s[k+1]+5))
        l[near&np.isin(l,[1,3])]=3
    ss=np.arange(0,s[-1]+1e-8,.5);idx=np.minimum(np.searchsorted(s,ss,side='right')-1,len(a)-1)
    xy=a[idx]+((ss-s[idx])/length[idx])[:,None]*(b[idx]-a[idx]);prov=l[idx].copy()
    prov[np.isin(prov,[1,3])&~(t2.sample(response,xy)>3)]=0
    mask=~np.isin(prov,[1,3]);row=np.full(len(ss),np.nan);row[~mask]=t2.sample(image,xy[~mask])
    for name,x in [('row',row),('mask',mask),('provenance',prov)]:np.save(dest/f'{name}.npy',x)
    np.savez_compressed(dest/'path_and_samples.npz',edge_start_xy=a,edge_end_xy=b,edge_provenance=l,edge_original_provenance=base,original_T2_edge_id=origin,edge_s_px=s,s_px=ss,s_um=ss*SPACING,xy=xy,ct=row,mask=mask,provenance=prov)
    dump(dest/'corrections.json',audit)
    table(dest/'path_coordinates.csv',[{'edge':k,'start_x':a[k,0],'start_y':a[k,1],'end_x':b[k,0],'end_y':b[k,1],'start_s_um':s[k]*SPACING,'end_s_um':s[k+1]*SPACING,'provenance':int(l[k]),'original_edge_id':int(origin[k])} for k in range(len(a))])
    table(dest/'sample_provenance.csv',[{'sample':k,'s_um':ss[k]*SPACING,'x':xy[k,0],'y':xy[k,1],'provenance':int(prov[k]),'masked':int(mask[k]),'ct':'' if mask[k] else row[k]} for k in range(len(ss))])
    Image.fromarray(np.nan_to_num(row,nan=0).round().astype('uint8')[None,:]).save(dest/'unwrapped_row.png')
    Image.fromarray(mask.astype('uint8')[None,:]*255).save(dest/'unwrapped_mask.png')
    after=exhaustive_topology(a,b)
    fig,axs=plt.subplots(1,2,figsize=(13,6))
    for ax,aa,bb,ll,title in [(axs[0],oa,ob,ol,'T2 maximum-recovery'),(axs[1],a,b,l,'T3 finalized')]:
        ax.imshow(image,cmap='gray',vmin=0,vmax=255);ax.add_collection(LineCollection(np.stack([aa,bb],axis=1),colors=[t2.COLORS[int(v)] for v in ll],linewidths=.7))
        ax.scatter(*aa[0],c='lime',marker='*',s=60);ax.scatter(*bb[-1],c='red',s=20);ax.set(xlim=(25,475),ylim=(325,25),title=f'Z{z} {title}')
    fig.suptitle('Cyan observed; magenta unsampled bridge; orange low-confidence; star seed');fig.tight_layout();fig.savefig(dest/'geometry.png',dpi=140);plt.close(fig)
    record={'z':z,'total_length_mm':float(length.sum()*.12656),'observed_length_mm':float(length[l!=2].sum()*.12656),'bridge_length_mm':float(length[l==2].sum()*.12656),
            'observed_utilization_fraction':float(length[l!=2].sum()/length.sum()),'observed_retained_fraction':float(length[l!=2].sum()/beforeobs),
            'fragments':len(top[0][1]),'bridges':int((l==2).sum()),'low_confidence_fraction_samples':float((prov==3).mean()),'corrections':len(audit),
            'proper_crossings':after['proper_crossings'],'point_contacts':after['nonadjacent_point_contacts'],'collinear_overlaps':after['collinear_overlaps'],
            'final_defect_count':sum(after.values()),'exterior_reached':bool(t2.target_dist(b[-1])==0),'row_samples':len(row),'last_sample_s_um':float(ss[-1]*SPACING),
            'valid_samples':int((~mask).sum()),'bridge_samples':int((prov==2).sum()),'unresolved_samples':int((prov==0).sum()),'low_confidence_samples':int((prov==3).sum()),
            'seed_x':float(a[0,0]),'seed_y':float(a[0,1]),'end_x':float(b[-1,0]),'end_y':float(b[-1,1]),'cleanup_termination':termination,
            'zero_length_edges_removed':zero_removed,'complete_paths_evaluated':counts['complete_count']}
    dump(dest/'metrics.json',record);dump(dest/'topology_before.json',before)
    # Exact lineage and independent bilinear check, no other Z observation accessed.
    assert np.array_equal(b[:-1],a[1:]) and np.array_equal(a[0],oa[0]) and np.array_equal(b[-1],ob[-1])
    for k,oid in enumerate(origin):
        if oid<0:assert l[k]==2
        else:
            assert (np.array_equal(a[k],oa[oid]) and np.array_equal(b[k],ob[oid])) or (np.array_equal(a[k],ob[oid]) and np.array_equal(b[k],oa[oid]))
            assert base[k]==ol[oid] and (l[k]==base[k] or (base[k]==1 and l[k]==3))
    q=xy[~mask];x=np.floor(q[:,0]).astype(int);y=np.floor(q[:,1]).astype(int);dx=q[:,0]-x;dy=q[:,1]-y
    expected=image[y,x]*(1-dx)*(1-dy)+image[y,x+1]*dx*(1-dy)+image[y+1,x]*(1-dx)*dy+image[y+1,x+1]*dx*dy
    assert np.allclose(expected,row[~mask],rtol=0,atol=1e-10)
    assert np.array_equal(mask,np.isnan(row)) and np.array_equal(mask,np.isin(prov,[0,2]))
    dump(dest/'validation.json',{'geometry_continuity':'PASS','observed_lineage':'PASS','independent_bilinear':'PASS','topology':after,'topology_clean':not any(after.values())})
    dump(dest/'FINALIZED.json',{str(p.relative_to(OUT)):sha(p) for p in sorted(dest.rglob('*')) if p.is_file()})
    return record

def run_all():
    initialize();results=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(reconstruct,z):z for z in ZS}
        for f in as_completed(pending):
            z=pending[f]
            try:
                rec=f.result();results.append(rec);print(json.dumps(rec),flush=True)
            except Exception:
                dump(OUT/f'slice_{z}/failure.json',{'traceback':traceback.format_exc()});print(f'FAILED Z{z}',flush=True)
    dump(OUT/'reconstruction_summary.json',sorted(results,key=lambda x:x['z']))
    assert len(results)==50,'Failed slices preserved; no evaluation started'
    # All rows finalized before an evaluation function can load neighboring rows.
    dump(OUT/'ROWS_FINALIZED.json',{f'slice_{z}/FINALIZED.json':sha(OUT/f'slice_{z}/FINALIZED.json') for z in ZS})
    print('ALL 50 ROWS FINALIZED; evaluation not yet run',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');args=parser.parse_args()
    if args.run:run_all()
