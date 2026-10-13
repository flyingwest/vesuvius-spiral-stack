"""Two fixed estimators; gate entirely on CT geometry, never model output."""
from pathlib import Path
import json,csv,sys,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_t11_core import frame,render,quality,cube
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t11'
def serial(a):
    if isinstance(a,np.ndarray):return a.tolist()
    raise TypeError(type(a))
def main():
    replay='--replay' in sys.argv;dest=OUT/('frame_replay' if replay else 'frames');dest.mkdir(exist_ok=False)
    samples=json.loads((OUT/'samples.json').read_text());records=[];rows=[];comparisons=[];arrays={}
    for i,p in enumerate(samples):
        P,t,n=[np.array(p[k]) for k in ['P_xyz','tangent_xyz','normal_xyz']]
        try:old=render(P,t,[0,0,1],n);oq=quality(old)
        except FileNotFoundError as e:
            records.extend([{'sample':i,**p,'sigma':sigma,'accepted':False,'missing_data':str(e)} for sigma in [1.,2.]])
            continue
        for sigma in [1.,2.]:
            try:
                f=frame(P,t,n,sigma);a=render(P,f['t_s'],f['image_y'],f['n']);q=quality(a)
                record={'sample':i,**p,'sigma':sigma,**f,'old_quality':oq,'new_quality':q,'patch_sha256':hashlib.sha256(a.tobytes()).hexdigest()}
                np.savez_compressed(dest/f'patch_{i}_{int(sigma)}.npz',old=old,new=a)
                records.append(record);rows.append({'sample':i,'z':p['z'],'s_index':p['s_index'],'sigma':sigma,'accepted':f['accepted'],**{'old_'+k:v for k,v in oq.items()},**{'new_'+k:v for k,v in q.items()}})
            except FileNotFoundError as e:records.append({'sample':i,**p,'sigma':sigma,'accepted':False,'missing_data':str(e)})
        print('frame sample',i+1,'/',len(samples),flush=True)
    if any('missing_data' in r for r in records):
        (dest/'records.json').write_text(json.dumps(records,default=serial,indent=2)+'\n')
        (OUT/'DATA_BLOCKER.json').write_text(json.dumps({'status':'BLOCKED','reason':'Missing CT coverage, not orientation failure','missing':[r for r in records if 'missing_data' in r]},indent=2)+'\n')
        raise RuntimeError('Missing CT coverage; see DATA_BLOCKER.json; no scientific verdict')
    for sigma in [1.,2.]:
        rr=[r for r in records if r['sigma']==sigma];ok=[r for r in rr if r['accepted']];valid=[r for r in rr if 'new_quality' in r]
        def med(k,which):return float(np.median([r[which+'_quality'][k] for r in ok])) if ok else None
        confidence=[r['confidence'] for r in valid];pert=[r['perturbation_deg'] for r in valid]
        sanity=len(ok)/len(samples)>=.7 and bool(pert) and np.median(pert)<=15
        gates={'accepted_at_least70pct':sanity,'all72patch_pairs_available':len(valid)==len(samples)}
        if ok:
            gates.update({'heldout_longitudinal_gradient_reduction10pct':med('longitudinal_gradient_fraction','new')<=.9*med('longitudinal_gradient_fraction','old'),'depth_peak_spread_reduction10pct':med('depth_peak_spread','new')<=.9*med('depth_peak_spread','old'),'majority60pct_gradient_improvement':np.mean([r['new_quality']['longitudinal_gradient_fraction']<r['old_quality']['longitudinal_gradient_fraction'] for r in ok])>=.6,'empty_fraction_nondegradation':med('empty_fraction','new')<=med('empty_fraction','old')+.03,'center_distance_nondegradation':med('center_distance','new')<=med('center_distance','old')+.05,'band_crossings_nondegradation':med('band_crossings','new')<=med('band_crossings','old')})
        for k in gates:gates[k]=bool(gates[k])
        smooth=[]
        for z in sorted(set(r['z'] for r in ok)):
            ordered=sorted([r for r in ok if r['z']==z],key=lambda r:r['s_index'])
            smooth += [float(np.degrees(np.arccos(np.clip(abs(a['n']@b['n']),0,1)))) for a,b in zip(ordered,ordered[1:])]
        comparisons.append({'sigma':sigma,'accepted':len(ok),'rejected':len(rr)-len(ok),'missing_data':sum('missing_data'in r for r in rr),'sanity_pass':bool(sanity),'quality_pass':bool(all(gates.values())),'gates':gates,'normal_abs_z_quantiles':np.quantile([abs(r['normal_z']) for r in valid],[0,.5,.95,1]).tolist() if valid else [],'tilt_quantiles':np.quantile([r['tilt_deg'] for r in valid],[0,.5,.95,1]).tolist() if valid else [],'confidence_quantiles':np.quantile(confidence,[0,.5,.95,1]).tolist() if valid else [],'perturbation_quantiles':np.quantile(pert,[0,.5,.95,1]).tolist() if valid else [],'sparse_station_normal_angle_quantiles':np.quantile(smooth,[0,.5,.95,1]).tolist() if smooth else [],'quality_medians':{k:{which:med(k,which) for which in ['old','new']} for k in oq}})
    selected=next((r['sigma'] for r in comparisons if r['quality_pass']),None)
    result={'methods':comparisons,'selected_sigma':selected,'frame_estimation':'PASS' if any(r['sanity_pass'] for r in comparisons) else 'FAIL','decision_if_stop':'C' if any(r['sanity_pass'] for r in comparisons) else 'D','sample_count':len(samples),'right_handed_frame_vs_image':'t_s cross t_z = n; image Y=-t_z retains original T9 left-handed array convention; depth sign aligned with saved T9 normal','quality_scope':'CT occupancy/gradient proxies, not ink or a ground-truth surface'}
    (dest/'records.json').write_text(json.dumps(records,default=serial,indent=2)+'\n');(dest/'frame_metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    with (dest/'patch_quality_metrics.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    if replay:
        original=json.loads((OUT/'frames/records.json').read_text());assert json.loads((dest/'records.json').read_text())==original
        assert json.loads((dest/'frame_metrics.json').read_text())==json.loads((OUT/'frames/frame_metrics.json').read_text())
        (OUT/'frame_replay.json').write_text(json.dumps({'status':'PASS','sample_methods_replayed':len(records),'frames_metrics_patch_hashes_exact':True},indent=2)+'\n');return
    (OUT/'frame_metrics.json').write_text((dest/'frame_metrics.json').read_text());(OUT/'patch_quality_metrics.csv').write_text((dest/'patch_quality_metrics.csv').read_text())
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    for sigma in [1.,2.]:
        rr=[r for r in records if r['sigma']==sigma and 'confidence'in r]
        for ax,key in zip(axs,['normal_z','tilt_deg','confidence']):ax.hist([r[key] for r in rr],bins=15,alpha=.5,label=f'sigma{sigma}');ax.set_title(key);ax.legend()
    fig.tight_layout();fig.savefig(OUT/'frame_distributions.png');plt.close(fig)
    chosen=selected or 1.;fig,axs=plt.subplots(3,4,figsize=(12,8));overlay=OUT/'frame_overlays';overlay.mkdir()
    for row,i in enumerate([0,32,64]):
        p=samples[i];r=next(r for r in records if r['sample']==i and r['sigma']==chosen)
        with np.load(dest/f'patch_{i}_{int(chosen)}.npz') as a:
            for col,plane in enumerate([a['old'][8],a['new'][8],a['old'][:,64],a['new'][:,64]]):axs[row,col].imshow(plane,cmap='gray',vmin=0,vmax=255,aspect='auto');axs[row,col].set_title(f'sample{i}: '+['old model X-Y','new model X-Y','old depth-X','new depth-X'][col])
        v=cube(p['P_xyz']);ff,aa=plt.subplots(1,3,figsize=(12,4))
        for ax,plane,indices,title in zip(aa,[v[16],v[:,16],v[:,:,16]],[(0,1),(0,2),(1,2)],['XY','XZ','YZ']):
            ax.imshow(plane,cmap='gray',vmin=0,vmax=255);ax.set_title(title+' fixed center');ax.plot(16,16,'r+')
            for key,color in [('t_s','cyan'),('image_y','yellow'),('n','red')]:
                vv=np.asarray(r[key]);ax.arrow(16,16,10*vv[indices[0]],10*vv[indices[1]],color=color,width=.15,length_includes_head=True)
        ff.tight_layout();ff.savefig(overlay/f'sample_{i}.png');plt.close(ff)
    fig.tight_layout();fig.savefig(OUT/'t9_vs_t11_patch_comparison.png');plt.close(fig);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
