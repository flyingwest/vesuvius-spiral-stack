"""Independent 2-D polar ridge traversal with explicitly masked geometric gaps."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.signal import find_peaks
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_g1'
PLAN={'hypothesis':'Unwind each XY independently; same-slice credible continuation may cross masked gaps; stack rows without cross-Z identity requirements.',
 'source':'outputs/orientation/ct_level4_zyx.npy','z_inclusive':[128,177],'z_offset':128,'spacing_um':126.56,
 'reference':{'center_xy':[240.,140.],'ellipse_y_scale':.5,'start_angle_rad':0.,
     'direction':'increasing angle; clockwise in image coordinates','seed':'nearest eligible positive-X ray ridge to elliptical radius40 within20..80; no prior trace',
     'seed_tangent_min_cos':.65,'seed_coherence_min':.2},
 'polar':{'angle_step_deg':.5,'r_min_px':8.,'r_max_px':180.,'r_step_px':.25,'max_turns':3},
 'ridge':{'DoG_sigma_px':[.6,2.],'prominence_CT':2.5,'response_CT':3.,'radial_peak_separation_px':1.,
          'structure_tensor_gradient_sigma':1.,'structure_tensor_average_sigma':1.2},
 'continuation':{'prediction':'linear radial slope from last up-to7 consecutive points, clipped ±120 px/radian',
     'max_skipped_angle_deg':15.,'candidate_radial_window_px':6.,'max_chord_px':14.,
     'min_orientation_cos':.65,'max_normal_miss_spacing_fraction':.45,
     'score_weights':{'proximity':.4,'tangent':.35,'curvature':.15,'orientation':.1},
     'min_score':.68,'winner_margin':.06,'max_turn_angle_deg':65.,
     'repeat_stop':'after300 degrees, any return within1px of geometry at least240 degrees earlier; no forced closure'},
 'sampling':{'ds_px':.5,'ds_um':63.28,'method':'linear XY in Euclidean chord arc; bilinear original CT, float64',
     'gaps':'chord length between endpoints is a lower-bound placeholder, ALL samples inside blank; never sample connector',
     'terminal':'one masked sample sentinel after STOP; unknown remainder length not inferred; rectangular padding remains masked'},
 'no_cross_z':'Pure process_slice(image) function receives no Z, neighbor image, old trace or track.',
 'limits':'Fixed geometric reference is not a measured roll center. Angular-monotone bounded path can stop at folds; one selected candidate winding, not all sheets.',
 'success_evaluation':'50 independently constructed rows with nonempty recovered geometry and correctly masked gaps; inspect coverage and failure modes, no cross-Z persistence gate.'}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,value):p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def sample(im,p):
    p=np.asarray(p);return map_coordinates(np.asarray(im,float),[p[...,1].ravel(),p[...,0].ravel()],order=1,prefilter=False,mode='constant',cval=0).reshape(p.shape[:-1])
def norm(v):return v/max(float(np.linalg.norm(v)),1e-12)


def detect(image):
    center=np.array(PLAN['reference']['center_xy']);scale=.5
    theta=np.arange(720)*np.pi/360.;radii=np.arange(8,180.001,.25)
    xy=center+np.stack([np.cos(theta),scale*np.sin(theta)],axis=1)[:,None,:]*radii[None,:,None]
    im=np.asarray(image,float);resp=gaussian_filter(im,.6)-gaussian_filter(im,2.)
    smooth=gaussian_filter(im,1.);gy,gx=np.gradient(smooth)
    xx=gaussian_filter(gx*gx,1.2);yy=gaussian_filter(gy*gy,1.2);cross=gaussian_filter(gx*gy,1.2)
    tangent=.5*np.arctan2(2*cross,xx-yy)+np.pi/2
    tangent_sin=np.sin(2*tangent);tangent_cos=np.cos(2*tangent)
    coherence=np.sqrt((xx-yy)**2+4*cross**2)/(xx+yy+1e-9)
    values=sample(resp,xy);records=[];bytheta=[]
    for k,profile in enumerate(values):
        peaks,props=find_peaks(profile,prominence=2.5,distance=4)
        peaks=peaks[profile[peaks]>3.]
        rr=[]
        for j in peaks:
            den=profile[j-1]-2*profile[j]+profile[j+1]
            correction=np.clip(.5*(profile[j-1]-profile[j+1])/den,-.5,.5) if abs(den)>1e-12 else 0.
            rr.append(float(radii[j]+correction*.25))
        indices=[]
        for j,r in enumerate(rr):
            p=center+np.array([np.cos(theta[k]),scale*np.sin(theta[k])])*r
            if not (2<=p[0]<image.shape[1]-2 and 2<=p[1]<image.shape[0]-2):continue
            # Interpolate doubled angle to respect unoriented ridge tangent symmetry.
            aa=.5*np.arctan2(sample(tangent_sin,p[None])[0],sample(tangent_cos,p[None])[0])
            local=min(([r-rr[j-1]] if j else [])+([rr[j+1]-r] if j+1<len(rr) else []) or [6.])
            indices.append(len(records));records.append([k,r,p[0],p[1],aa,float(sample(coherence,p[None])[0]),float(profile[peaks[j]]),local])
        bytheta.append(indices)
    return np.array(records).reshape(-1,8),bytheta


def process_slice(image):
    candidates,bytheta=detect(image);center=np.array([240.,140.]);step=np.pi/360.
    seeds=[]
    for idx in bytheta[0]:
        c=candidates[idx];t=np.array([np.cos(c[4]),np.sin(c[4])])
        if 20<=c[1]<=80 and abs(t[1])>=.65 and c[5]>=.2:seeds.append(idx)
    if not seeds:return {'candidates':candidates,'points':[],'decisions':[],'gaps':[],'stop':'NO_ELIGIBLE_START_RIDGE','s':np.array([0.]),'row':np.array([np.nan]),'xy':np.full((1,2),np.nan),'summary':{'recovered_arc_px':0.,'gap_placeholder_arc_px':0.,'angular_coverage_deg':0.,'continuation_decisions':0,'accepted_gap_continuations':0,'interruptions':1,'accepted_local_steps':0,'stop_reason':'NO_ELIGIBLE_START_RIDGE','candidate_count':len(candidates),'seed_radius_px':None}}
    seed=min(seeds,key=lambda i:(abs(candidates[i,1]-40),-candidates[i,5],i))
    points=[{'station':0,'candidate':seed,'xy':candidates[seed,2:4].tolist(),'radius':float(candidates[seed,1]),'connected':False}]
    decisions=[];gaps=[];used={seed};stop='MAX_THREE_TURN_BOUND';maxstation=2160
    while points[-1]['station']<maxstation:
        last=points[-1];i=last['station'];p=np.array(last['xy']);lastc=candidates[last['candidate']]
        recent=[last]
        for prev in reversed(points[:-1]):
            if recent[-1]['station']-prev['station']!=1 or len(recent)>=7:break
            recent.append(prev)
        if len(recent)>=2:
            tt=np.array([q['station'] for q in recent])*step;rr=np.array([q['radius'] for q in recent]);slope=float(np.clip(np.polyfit(tt-tt[0],rr,1)[0],-120,120))
            tangent=norm(p-np.array(recent[min(3,len(recent)-1)]['xy']))
        else:
            tangent=np.array([np.cos(lastc[4]),np.sin(lastc[4])]);theta=i*step
            orient=np.array([-np.sin(theta),.5*np.cos(theta)])
            if tangent@orient<0:tangent=-tangent
            radial_unit=np.array([np.cos(theta),.5*np.sin(theta)])
            denom=np.linalg.det(np.c_[tangent,radial_unit])
            dr=np.linalg.lstsq(np.c_[radial_unit,-tangent],-last['radius']*orient,rcond=None)[0][0]
            slope=float(np.clip(dr,-120,120))
        accepted=None;first_failed=None;station_records=[]
        for j in range(i+1,min(i+30,maxstation)+1):
            theta=j*step;predr=last['radius']+slope*(j-i)*step
            pred=center+np.array([np.cos(theta),.5*np.sin(theta)])*predr
            options=[]
            for idx in bytheta[j%720]:
                c=candidates[idx]
                if idx in used or abs(c[1]-predr)>6:continue
                q=c[2:4];delta=q-p;dist=float(np.linalg.norm(delta))
                if dist<.03 or dist>14:continue
                direction=delta/dist;ct=np.array([np.cos(c[4]),np.sin(c[4])])
                align=abs(float(ct@direction));bend=float(np.arccos(np.clip(tangent@direction,-1,1)))
                radial_step=np.array([np.cos(theta),.5*np.sin(theta)])
                spacing=max(float(c[7]*np.linalg.norm(radial_step)),.5)
                miss=float(np.linalg.norm(q-pred));normal_miss=abs(float(np.linalg.det(np.c_[direction,q-pred])))
                # Predictive proximity and local tangent ordering, not previous-Z identity.
                comp={'proximity':float(np.exp(-miss/max(1.,.4*spacing))),
                      'tangent':max(float(tangent@direction),0.),'curvature':float(np.exp(-bend/.7)), 'orientation':align}
                score=sum(comp[k]*v for k,v in PLAN['continuation']['score_weights'].items())
                failures=[]
                if align<.65:failures.append('ridge_orientation')
                if bend>np.radians(65):failures.append('turn_angle')
                if normal_miss>.45*spacing:failures.append('layer_order_ambiguity')
                if c[5]<.15:failures.append('low_orientation_coherence')
                if score<.68:failures.append('score')
                options.append({'candidate':int(idx),'score':float(score),'normal_miss_px':normal_miss,'spacing_px':spacing,'prediction_error_px':miss,'angle_deg':float(np.degrees(bend)),'failures':failures})
            options.sort(key=lambda o:(-o['score'],o['candidate']))
            viable=[o for o in options if not o['failures']]
            chosen=viable[0] if viable else None
            if chosen:
                competitors=[o['score'] for o in options if o['candidate']!=chosen['candidate'] and o['angle_deg']<=65]
                chosen['margin']=chosen['score']-max(competitors,default=0.)
                if chosen['margin']<.06:chosen=None
            station_records.append({'station':j,'options':options,'decision':'CONTINUE' if chosen else 'MASK_CANDIDATE'})
            if chosen:
                idx=chosen['candidate'];q=candidates[idx,2:4]
                old=[v for v in points if j-v['station']>=480]
                if j>=600 and old and min(np.linalg.norm(q-np.array(v['xy'])) for v in old)<1:
                    stop='REVISIT_OR_CLOSED_LOOP';accepted='STOP';break
                accepted={'station':j,'candidate':idx,'xy':q.tolist(),'radius':float(candidates[idx,1]),'connected':j==i+1}
                break
            if first_failed is None:first_failed=j
        decisions.extend(station_records)
        if accepted=='STOP':break
        if accepted is None:
            stop='NO_CREDIBLE_CONTINUATION_WITHIN_15_DEG';break
        if not accepted['connected']:
            gaps.append({'from_point':len(points)-1,'to_point':len(points),'from_station':i,'to_station':accepted['station'],
                'placeholder_length_px':float(np.linalg.norm(np.array(accepted['xy'])-p)),
                'meaning':'straight-chord lower bound; unresolved connector never sampled'})
        points.append(accepted);used.add(accepted['candidate'])
    xy=np.array([p['xy'] for p in points]);lengths=np.linalg.norm(np.diff(xy,axis=0),axis=1);s=np.r_[0.,np.cumsum(lengths)]
    grid=np.arange(int(np.floor(s[-1]/.5))+2)*.5;row=np.full(len(grid),np.nan);coords=np.full((len(grid),2),np.nan)
    recovered=0.
    for j in range(1,len(points)):
        if not points[j]['connected']:continue
        recovered+=float(lengths[j-1]);valid=(grid>=s[j-1]-1e-12)&(grid<=s[j]+1e-12)
        w=(grid[valid]-s[j-1])/max(lengths[j-1],1e-12);q=xy[j-1]+w[:,None]*(xy[j]-xy[j-1])
        coords[valid]=q;row[valid]=sample(image,q)
    # Any sample strictly inside a declared gap must remain missing, regardless of overlap tolerances.
    for gap in gaps:
        valid=(grid>s[gap['from_point']]+1e-12)&(grid<s[gap['to_point']]-1e-12)
        assert np.isnan(row[valid]).all()
    summary={'recovered_arc_px':recovered,'gap_placeholder_arc_px':float(sum(g['placeholder_length_px'] for g in gaps)),
        'angular_coverage_deg':points[-1]['station']*.5,'continuation_decisions':len(gaps)+int(stop!='MAX_THREE_TURN_BOUND'),
        'accepted_gap_continuations':len(gaps),'interruptions':len(gaps)+int(stop!='MAX_THREE_TURN_BOUND'),
        'accepted_local_steps':sum(p['connected'] for p in points),'stop_reason':stop,
        'candidate_count':len(candidates),'seed_radius_px':float(candidates[seed,1]),'trace_points':len(points),
        's_extent_px':float(s[-1]),'native_row_width_samples':len(grid),
        'within_native_row_masked_fraction':float(np.mean(~np.isfinite(row)))}
    return {'candidates':candidates,'points':points,'decisions':decisions,'gaps':gaps,'stop':stop,'s':grid,'row':row,'xy':coords,'summary':summary}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if not args.verify:
        OUT.mkdir(exist_ok=False);(OUT/'all_50_slice_geometry_overlays').mkdir();(OUT/'slices').mkdir()
        save(OUT/'plan.json',PLAN)
        prior=[p for folder in ['outputs/orientation','outputs/p0_1','outputs/p0_1b','outputs/p0_2','outputs/p0_3','outputs/p0_4','docs','scripts'] for p in (ROOT/folder).rglob('*') if p.is_file()]+[ROOT/'README.md',ROOT/'AGENTS.md']
        save(OUT/'preservation_before.json',{str(p.relative_to(ROOT)):sha(p) for p in prior})
    volume=np.load(ROOT/PLAN['source'],mmap_mode='r');results=[]
    for z in range(128,178):
        image=np.asarray(volume[z-128]);r=process_slice(image);results.append(r)
        record={k:r[k] for k in ['points','decisions','gaps','summary']}
        if args.verify:
            assert json.loads((OUT/'slices'/f'z{z}.json').read_text())==record,z
        else:
            save(OUT/'slices'/f'z{z}.json',record)
            np.savez_compressed(OUT/'slices'/f'z{z}_numeric.npz',candidates=r['candidates'],row=r['row'],sample_xy=r['xy'],s_px=r['s'])
            fig,axes=plt.subplots(1,2,figsize=(13,6))
            for ax in axes:ax.imshow(image,cmap='gray',vmin=0,vmax=255,interpolation='nearest');ax.set(xlim=(45,435),ylim=(250,35))
            p=np.array([q['xy'] for q in r['points']]).reshape(-1,2)
            for j in range(1,len(p)):
                axes[1].plot(*p[j-1:j+1].T,color='lime' if r['points'][j]['connected'] else 'magenta',ls='-' if r['points'][j]['connected'] else '--',lw=1)
            if len(p):axes[1].plot(*p[0],'yo',ms=5);axes[1].plot(*p[-1],'rx',ms=6)
            axes[1].plot(240,140,'c+',ms=7);axes[1].plot([240,320],[140,140],'y:',lw=.7)
            axes[0].set_title(f'Z={z}: original CT; no neighboring slice inputs')
            axes[1].set_title(f"Independent geometry: {r['summary']['recovered_arc_px']:.1f}px, {r['summary']['angular_coverage_deg']:.1f}°\ngreen valid; magenta masked gap; yellow start; red STOP")
            fig.tight_layout();fig.savefig(OUT/'all_50_slice_geometry_overlays'/f'z{z}_geometry.png',dpi=130);plt.close(fig)
        if (z-127)%10==0:print(f'processed Z128..{z}',flush=True)
    width=max(len(r['row']) for r in results);stack=np.full((50,width),np.nan);coords=np.full((50,width,2),np.nan);rows=[]
    for i,r in enumerate(results):
        stack[i,:len(r['row'])]=r['row'];coords[i,:len(r['xy'])]=r['xy'];m=r['summary']
        rows.append({'z':128+i,'recovered_arc_px':m['recovered_arc_px'],'recovered_arc_um':m['recovered_arc_px']*126.56,
            'angular_coverage_deg':m['angular_coverage_deg'],'interruptions':m['interruptions'],
            'continuation_decisions':m['continuation_decisions'],'accepted_gap_continuations':m['accepted_gap_continuations'],
            'accepted_local_steps':m['accepted_local_steps'],'masked_fraction':float(np.mean(~np.isfinite(stack[i]))),
            'gap_placeholder_arc_px':m['gap_placeholder_arc_px'],'sampling_interval_um':63.28,
            'native_row_width_samples':len(r['row']),'row_width_samples':width,'valid_samples':int(np.isfinite(stack[i]).sum()),
            'reference':'fixed(240,140), ellipse yscale0.5; +X ray nearest r40; clockwise; s=0 at independent seed',
            'notes':m['stop_reason']})
    mask=np.isfinite(stack);data={'lossless_stack':stack,'lossless_mask':mask,'sample_xy':coords,'s_px':np.arange(width)*.5,'z':np.arange(128,178)}
    for name,array in data.items():
        if args.verify:np.testing.assert_array_equal(np.load(OUT/f'{name}.npy'),array)
        else:np.save(OUT/f'{name}.npy',array,allow_pickle=False)
    if not args.verify:
        with (OUT/'row_geometry.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        save(OUT/'row_geometry.json',rows)
        preview=np.zeros_like(stack,dtype=np.uint8);preview[mask]=np.rint(stack[mask]).astype(np.uint8)
        Image.fromarray(preview).save(OUT/'raw_stack.png');Image.fromarray(mask.astype(np.uint8)*255).save(OUT/'stack_mask.png')
        fig,axes=plt.subplots(2,2,figsize=(15,9));cmap=plt.get_cmap('gray').copy();cmap.set_bad('magenta')
        extent=[0,width*.5*.12656,177.5,127.5]
        axes[0,0].imshow(np.ma.masked_invalid(stack),cmap=cmap,vmin=0,vmax=255,aspect='auto',interpolation='nearest',extent=extent);axes[0,0].set(title='Independent raw CT stack; magenta = unresolved',xlabel='Geometric s (mm; gap lengths approximate)',ylabel='Z')
        axes[0,1].imshow(mask,cmap='gray',vmin=0,vmax=1,aspect='auto',interpolation='nearest',extent=extent);axes[0,1].set(title='Valid-data mask',xlabel='s (mm)',ylabel='Z')
        axes[1,0].plot([r['z'] for r in rows],[r['recovered_arc_px'] for r in rows],'o-');axes[1,0].set(title='Recovered same-slice arc length',xlabel='Z',ylabel='level4 pixels');axes[1,0].grid(alpha=.3)
        axes[1,1].plot([r['z'] for r in rows],[r['accepted_gap_continuations'] for r in rows],'o-',label='accepted masked-gap continuation');axes[1,1].plot([r['z'] for r in rows],[r['interruptions'] for r in rows],'x-',label='interruptions incl STOP');axes[1,1].legend();axes[1,1].set(xlabel='Z',ylabel='count')
        fig.tight_layout();fig.savefig(OUT/'geometry_diagnostics.png',dpi=150);plt.close(fig)
    for name,digest in json.loads((OUT/'preservation_before.json').read_text()).items():assert sha(ROOT/name)==digest,name
    lengths=np.array([r['recovered_arc_px'] for r in rows]);masked=np.array([r['masked_fraction'] for r in rows]);counts=np.array([r['accepted_gap_continuations'] for r in rows])
    summary={'slices_processed':50,'shape_z_s':list(stack.shape),'slices_with_valid_data':int(np.any(mask,axis=1).sum()),
        'recovered_arc_px_min_median_mean_max':[float(min(lengths)),float(np.median(lengths)),float(np.mean(lengths)),float(max(lengths))],
        'masked_fraction_min_median_mean_max':[float(min(masked)),float(np.median(masked)),float(np.mean(masked)),float(max(masked))],
        'accepted_gap_continuations_min_median_mean_max':[int(min(counts)),float(np.median(counts)),float(np.mean(counts)),int(max(counts))],
        'accepted_gap_continuations_total':int(sum(counts)),'stop_counts':{s:sum(r['notes']==s for r in rows) for s in sorted(set(r['notes'] for r in rows))},
        'max_angular_coverage_deg':float(max(r['angular_coverage_deg'] for r in rows)),
        'cross_z_input_or_identity_gate':False,'geometry_only':True}
    if args.verify:
        assert json.loads((OUT/'metrics.json').read_text())==summary
        save(OUT/'verification.json',{'exact_50_slice_replay':'PASS','prior_files':'unchanged','all_gap_interiors_masked':'PASS'})
    else:save(OUT/'metrics.json',summary)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
