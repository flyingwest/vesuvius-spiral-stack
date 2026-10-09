"""P0.3: strict validated-core arc-length rasters, no missing-data filling."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/p0_3'
PLAN={
 'hypothesis':'Geometrically anchored physical arc-length rows have coherent neighboring-Z CT structure.',
 'source':'outputs/orientation/ct_level4_zyx.npy', 'source_z_offset':128,
 'tracks':{'F02':'T0077','F03':'T0243'},
 'support':'Intersection of validated control u stations across all observations; no variable tails.',
 'core_u_px':{'F02':[1.,3.5],'F03':[1.,1.5]},
 'voxel_spacing_um_zyx':[126.56,126.56,126.56], 'ds_px':.25,'ds_um':31.64,
 'anchor':'s=0 at midpoint of common u interval on every observed polyline; geometry-only translation.',
 'sampling':'Piecewise-linear XY at physical arc length; float64 bilinear original uint8 CT; no CT smoothing.',
 'missing':'Inclusive track Z span; absent observations and outside-core samples NaN, mask false. Never interpolate Z.',
 'coherence_diagnostic':{'minimum_common_arc_span_px':2.,'minimum_samples':9,'minimum_row_std_CT':1.,
    'minimum_eligible_adjacent_fraction':.75,'minimum_mean_NCC':.5,'minimum_positive_fraction':.75},
 'null':'2000 fixed-seed permutations of row order; descriptive, not independent samples or a success gate.',
 'misregistration_probe':'Diagnostic only: shifts -0.5..0.5 px step0.25, same support at all shifts, never applied.',
 'misregistration_flag':'Eligible zero/shift correlations, nonzero optimum, NCC gain>=0.2 and optimum NCC>=0.5.',
 'no_tuning':'Exploratory gates fixed before P0.3 output; no retries to improve coherence.',
 'limits':'Small already-selected coarse CT cores; no new ridge detection, oracle, OCR, ink, mesh or language priors.'}


def read(path):return json.loads((ROOT/path).read_text())
def save(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def arc_geometry(xy,u,anchor):
    s=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(xy,axis=0),axis=1))]
    assert np.all(np.diff(s)>0)
    origin=float(np.interp(anchor,u,s))
    return s,origin


def sample_row(image,xy,s,origin,grid):
    local=grid+origin
    mask=(local>=-1e-12)&(local<=s[-1]+1e-12)
    pts=np.full((len(grid),2),np.nan);row=np.full(len(grid),np.nan)
    pts[mask]=np.c_[np.interp(local[mask],s,xy[:,0]),np.interp(local[mask],s,xy[:,1])]
    assert np.all((pts[mask,0]>=0)&(pts[mask,0]<=image.shape[1]-1))
    assert np.all((pts[mask,1]>=0)&(pts[mask,1]<=image.shape[0]-1))
    row[mask]=map_coordinates(image.astype(float),pts[mask,::-1].T,order=1,mode='constant',cval=np.nan,prefilter=False)
    return row,mask,pts


def correlation(a,b,grid):
    keep=np.isfinite(a)&np.isfinite(b);n=int(keep.sum())
    span=float(np.ptp(grid[keep])) if n else 0.
    aa=a[keep];bb=b[keep];std=[float(np.std(aa)),float(np.std(bb))] if n else [0.,0.]
    corr=float(np.corrcoef(aa,bb)[0,1]) if n>=2 and min(std)>1e-12 else None
    valid=n>=9 and span>=2.-1e-12 and min(std)>=1.
    return {'samples':n,'span_px':span,'std_CT':std,'descriptive_NCC':corr,
            'eligible':bool(valid),'NCC':corr if valid else None,
            'mean_absolute_difference_CT':float(np.mean(abs(aa-bb))) if n else None}


def compute(fid):
    nodes=read('outputs/p0_2/run_02/fragments.json')
    track=next(t for t in read('outputs/p0_2/run_02/tracks.json') if t['id']==PLAN['tracks'][fid])
    selected={n['z']:n for n in nodes if n['id'] in track['nodes']}
    volume=np.load(ROOT/PLAN['source'],mmap_mode='r')
    lo,hi=PLAN['core_u_px'][fid];anchor=(lo+hi)/2
    expected=read(f'outputs/p0_2/validation/{fid}_measurements.json')['common_u_px']
    assert [min(expected),max(expected)]==[lo,hi]
    geometries={};meta=[]
    zs=np.arange(min(selected),max(selected)+1)
    for z,node in selected.items():
        u=np.array(node['u']);keep=(u>=lo)&(u<=hi);xy=np.array(node['xy'])[keep];u=u[keep]
        assert u[0]==lo and u[-1]==hi
        s,origin=arc_geometry(xy,u,anchor)
        geometries[z]={'xy':xy,'u':u,'s':s,'origin':origin,'node':node['id']}
    ds=PLAN['ds_px']
    bounds=[(g['s'][0]-g['origin'],g['s'][-1]-g['origin']) for g in geometries.values()]
    grid=np.arange(int(np.floor(min(b[0] for b in bounds)/ds)),int(np.ceil(max(b[1] for b in bounds)/ds))+1)*ds
    stack=np.full((len(zs),len(grid)),np.nan);mask=np.zeros_like(stack,dtype=bool)
    coords=np.full(stack.shape+(2,),np.nan)
    for i,z in enumerate(zs):
        record={'z':int(z),'local_z':int(z-128),'observed':z in selected,'track':track['id'],'node':'',
                'arc_length_px':None,'arc_length_um':None,'alignment_offset_px':None,'alignment_offset_um':None,
                'adjacent_offset_change_px':None,'s_start_px':None,'s_end_px':None,'start_x':None,'start_y':None,
                'end_x':None,'end_y':None,'valid_samples':0}
        if z in geometries:
            g=geometries[z];stack[i],mask[i],coords[i]=sample_row(volume[z-128],g['xy'],g['s'],g['origin'],grid)
            record.update(node=g['node'],arc_length_px=float(g['s'][-1]),arc_length_um=float(g['s'][-1]*126.56),
                alignment_offset_px=-g['origin'],alignment_offset_um=-g['origin']*126.56,
                adjacent_offset_change_px=(-g['origin']+geometries[z-1]['origin']) if z-1 in geometries else None,
                s_start_px=-g['origin'],s_end_px=float(g['s'][-1]-g['origin']),
                start_x=float(g['xy'][0,0]),start_y=float(g['xy'][0,1]),end_x=float(g['xy'][-1,0]),end_y=float(g['xy'][-1,1]),
                valid_samples=int(mask[i].sum()))
        meta.append(record)
    pairs=[];scans=[];stretch=[]
    for i in range(len(zs)-1):
        z,zz=int(zs[i]),int(zs[i+1])
        r=correlation(stack[i],stack[i+1],grid);r.update(z_source=z,z_target=zz,
            observation_pair_present=z in geometries and zz in geometries)
        pairs.append(r)
        if not r['observation_pair_present']:continue
        ga,gb=geometries[z],geometries[zz]
        # Geometry residual at homologous window stations after pure s translation.
        sa=ga['s']-ga['origin'];sb=gb['s']-gb['origin'];ratios=np.diff(sb)/np.diff(sa)
        stretch.append({'z_source':z,'z_target':zz,'max_homologous_s_residual_px':float(max(abs(sb-sa))),
            'median_segment_length_ratio':float(np.median(ratios)), 'min_segment_length_ratio':float(min(ratios)),
            'max_segment_length_ratio':float(max(ratios))})
        shifts=np.arange(-.5,.501,.25);candidate=[]
        for shift in shifts:
            rr,mm,_=sample_row(volume[zz-128],gb['xy'],gb['s'],gb['origin'],grid+shift)
            candidate.append(rr)
        common=np.isfinite(stack[i])&np.all(np.isfinite(candidate),axis=0)
        scores=[]
        for shift,row in zip(shifts,candidate):
            aa=np.where(common,stack[i],np.nan);bb=np.where(common,row,np.nan)
            scores.append({'shift_px':float(shift),**correlation(aa,bb,grid)})
        eligible=[r for r in scores if r['eligible']]
        zero=scores[2];best=max(eligible,key=lambda r:r['NCC']) if eligible else None
        flag=bool(best and zero['eligible'] and best['shift_px']!=0 and best['NCC']-zero['NCC']>=.2 and best['NCC']>=.5)
        scans.append({'z_source':z,'z_target':zz,'scores':scores,'eligible_best':best,
                      'suspected_misregistration':flag,'applied':False})
    eligible=[r['NCC'] for r in pairs if r['eligible']]
    observed_pairs=sum(r['observation_pair_present'] for r in pairs)
    fraction=len(eligible)/observed_pairs if observed_pairs else 0.
    mean=float(np.mean(eligible)) if eligible else None
    positive=float(np.mean(np.array(eligible)>0)) if eligible else None
    status='NOT_EVALUABLE_SHORT_SUPPORT' if fraction<.75 else ('PASS_LOCAL_COHERENCE' if mean>=.5 and positive>=.75 else 'FAIL_COHERENCE_DIAGNOSTIC')
    # Nulls preserve measured rows/masks; Z shuffling does NOT supply more independent evidence.
    rng=np.random.default_rng(302);null=[]
    observed_indices=[i for i,z in enumerate(zs) if z in geometries]
    for _ in range(2000):
        perm=rng.permutation(observed_indices)
        vals=[correlation(stack[a],stack[b],grid)['NCC'] for a,b in zip(perm[:-1],perm[1:])]
        vals=[v for v in vals if v is not None]
        if vals:null.append(float(np.mean(vals)))
    lengths=[m['arc_length_px'] for m in meta if m['observed']]
    offset_changes=[m['adjacent_offset_change_px'] for m in meta if m['adjacent_offset_change_px'] is not None]
    summary={'track':track['id'],'z_inclusive':[int(zs[0]),int(zs[-1])],'observations':len(selected),
       'shape_z_s':list(stack.shape),'core_u_px':[lo,hi],'anchor_u_px':anchor,'ds_px':ds,'ds_um':31.64,
       'arc_length_px_min_max_mean':[float(min(lengths)),float(max(lengths)),float(np.mean(lengths))],
       'valid_row_samples_min_max':[min(m['valid_samples'] for m in meta if m['observed']),max(m['valid_samples'] for m in meta if m['observed'])],
       'missing_data_fraction':float(1-mask.mean()),'missing_observation_rows':[int(z) for z in zs if z not in selected],
       'missing_pixels':int((~mask).sum()),'valid_pixels':int(mask.sum()),
       'missing_observation_pixels':int(sum(z not in selected for z in zs)*len(grid)),
       'outside_core_padding_pixels':int(sum((~mask[i]).sum() for i,z in enumerate(zs) if z in selected)),
       'adjacent_geometric_offset_change_max_abs_px':max(abs(x) for x in offset_changes) if offset_changes else None,
       'observed_adjacent_pairs':observed_pairs,'eligible_correlation_pairs':len(eligible),
       'mean_adjacent_NCC':mean,'positive_adjacent_fraction':positive,
       'mean_adjacent_MAD_CT':float(np.mean([r['mean_absolute_difference_CT'] for r in pairs if r['observation_pair_present']])),
       'row_permutation_null_mean':float(np.mean(null)) if null else None,
       'row_permutation_null_95pct':np.percentile(null,[2.5,97.5]).tolist() if null else None,
       'null_fraction_ge_observed':float(np.mean(np.array(null)>=mean)) if null and mean is not None else None,
       'max_homologous_s_residual_px':max((r['max_homologous_s_residual_px'] for r in stretch),default=None),
       'misregistration_scan_eligible_pairs':sum(r['eligible_best'] is not None for r in scans),
       'misregistration_flags':sum(r['suspected_misregistration'] for r in scans),
       'status':status,'gap_identity_not_intensity_interpolated':True}
    return {'stack':stack,'mask':mask,'sample_xy':coords,'s_px':grid,'z':zs,'meta':meta,
            'pairs':pairs,'shift_scan':scans,'geometry_residuals':stretch,'summary':summary,
            'permutation_null':np.array(null)}


def checks(result):
    a=result['stack'];m=result['mask'];xy=result['sample_xy']
    assert np.array_equal(np.isfinite(a),m)
    assert np.array_equal(np.isfinite(xy).all(axis=-1),m)
    assert np.allclose(np.diff(result['s_px']),.25)
    assert np.all(np.diff(result['z'])==1)
    volume=np.load(ROOT/PLAN['source'],mmap_mode='r')
    # Independent direct four-pixel bilinear formula, not scipy resampling.
    for i,z in enumerate(result['z']):
        p=xy[i,m[i]];x=p[:,0];y=p[:,1];xx=np.floor(x).astype(int);yy=np.floor(y).astype(int)
        fx=x-xx;fy=y-yy;im=volume[z-128].astype(float)
        reference=(1-fx)*(1-fy)*im[yy,xx]+fx*(1-fy)*im[yy,xx+1]+(1-fx)*fy*im[yy+1,xx]+fx*fy*im[yy+1,xx+1]
        np.testing.assert_allclose(reference,a[i,m[i]],rtol=0,atol=1e-10)
    # Known bent-polyline arc length/anchor check exercises non-straight geometry.
    s,o=arc_geometry(np.array([[0.,0.],[3.,0.],[3.,4.]]),np.array([0.,1.,2.]),1.)
    np.testing.assert_array_equal(s,[0.,3.,7.]);assert o==3.


def plot(fid,r):
    grid=r['s_px']*126.56;zs=r['z'];extent=[grid[0]-15.82,grid[-1]+15.82,zs[-1]+.5,zs[0]-.5]
    fig,axes=plt.subplots(1,3,figsize=(13,6),gridspec_kw={'width_ratios':[2,1,2]})
    cmap=plt.get_cmap('gray').copy();cmap.set_bad('magenta')
    axes[0].imshow(np.ma.masked_invalid(r['stack']),cmap=cmap,vmin=0,vmax=255,extent=extent,aspect='auto',interpolation='nearest')
    for m in r['meta']:
        if m['observed']:axes[0].plot(np.array([m['s_start_px'],m['s_end_px']])*126.56,[m['z']]*2,'c|',ms=10)
    axes[0].axvline(0,color='yellow',lw=.7);axes[0].set(title='CT row stack; magenta = missing\ncyan ticks = exact retained endpoints',xlabel='Signed physical arc length (µm)',ylabel='Z level4')
    axes[1].imshow(r['mask'],cmap='gray',vmin=0,vmax=1,extent=extent,aspect='auto',interpolation='nearest');axes[1].set(title='Valid-data mask\nwhite = valid',xlabel='s (µm)')
    axes[2].plot([m['alignment_offset_um'] if m['observed'] else np.nan for m in r['meta']],zs,'o-')
    axes[2].invert_yaxis();axes[2].set(title='Geometry-only row offsets\nno scale/elastic correction',xlabel='Offset added to local arc length (µm)',ylabel='Z level4');axes[2].grid(alpha=.3)
    fig.suptitle(f'{fid} / {PLAN["tracks"][fid]} — panels enlarged for diagnosis; numeric arrays retain physical spacing')
    fig.tight_layout();fig.savefig(OUT/f'{fid}_diagnostic.png',dpi=150);plt.close(fig)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if not args.verify:
        OUT.mkdir(exist_ok=False);save(OUT/'plan.json',PLAN)
        prior=[p for folder in ['outputs/orientation','outputs/p0_1','outputs/p0_1b','outputs/p0_2','docs'] for p in (ROOT/folder).rglob('*') if p.is_file()]
        prior += [ROOT/'README.md',ROOT/'AGENTS.md',Path(__file__)]
        save(OUT/'preservation_before.json',{str(p.relative_to(ROOT)):sha(p) for p in prior})
    summaries={}
    for fid in PLAN['tracks']:
        r=compute(fid);checks(r);summaries[fid]=r['summary']
        for name in ['stack','mask','sample_xy','s_px','z','permutation_null']:
            path=OUT/f'{fid}_{name}.npy'
            if args.verify:np.testing.assert_array_equal(np.load(path),r[name])
            else:np.save(path,r[name],allow_pickle=False)
        for name in ['meta','pairs','shift_scan','geometry_residuals','summary']:
            path=OUT/f'{fid}_{name}.json'
            if args.verify:assert json.loads(path.read_text())==r[name],name
            else:save(path,r[name])
        if not args.verify:
            with (OUT/f'{fid}_alignment.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(r['meta'][0]));w.writeheader();w.writerows(r['meta'])
            preview=np.zeros(r['stack'].shape,dtype=np.uint8);preview[r['mask']]=np.rint(np.clip(r['stack'][r['mask']],0,255)).astype(np.uint8)
            Image.fromarray(preview).save(OUT/f'{fid}_stack.png')
            Image.fromarray(r['mask'].astype(np.uint8)*255).save(OUT/f'{fid}_stack_mask.png')
            plot(fid,r)
    assert summaries['F03']['missing_observation_rows']==[146]
    for name,digest in read('outputs/p0_3/preservation_before.json').items():assert sha(ROOT/name)==digest,name
    if args.verify:save(OUT/'verification.json',{'exact_replay':'PASS','independent_bilinear_reference':'PASS','bent_polyline_arc_length':'PASS','NaN_masks_and_Z_spacing':'PASS','prior_hashes':'unchanged'})
    else:save(OUT/'metrics.json',summaries)
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
