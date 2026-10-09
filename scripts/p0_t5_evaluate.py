"""Post-finalization evaluation only: never changes a slice path or row."""
import argparse,csv,json,hashlib
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_t5_consecutive import ROOT,OUT,ZS,DS,SPACING,dump,sha,table

def corr(a,b,minimum=32):
    valid=np.isfinite(a)&np.isfinite(b);n=int(valid.sum())
    if n<minimum:return None,n
    x=a[valid]-a[valid].mean();y=b[valid]-b[valid].mean();den=np.linalg.norm(x)*np.linalg.norm(y)
    return (float(x@y/den) if den else None),n

def barrier():
    hashes=json.loads((OUT/'ROWS_FINALIZED.json').read_text());assert len(hashes)==50
    for name,h in hashes.items():
        assert sha(OUT/name)==h,name
        for f,digest in json.loads((OUT/name).read_text()).items():assert sha(OUT/f)==digest,f

def build():
    records=[json.loads((OUT/f'slice_{z}/metrics.json').read_text()) for z in ZS]
    width=max(r['row_samples'] for r in records)
    a=np.full((50,width),np.nan);m=np.ones(a.shape,bool);p=np.zeros(a.shape,np.uint8);coverage=np.zeros(a.shape,bool)
    for i,z in enumerate(ZS):
        folder=OUT/f'slice_{z}';row=np.load(folder/'row.npy');mask=np.load(folder/'mask.npy');prov=np.load(folder/'provenance.npy');coords=np.load(folder/'path_and_samples.npz')
        assert np.array_equal(coords['s_um'],np.arange(len(row))*DS)
        assert np.array_equal(row,coords['ct'],equal_nan=True) and np.array_equal(np.isnan(row),mask)
        assert np.array_equal(mask,np.isin(prov,[0,2]))
        n=len(row);a[i,:n]=row;m[i,:n]=mask;p[i,:n]=prov;coverage[i,:n]=True
    return a,m,p,coverage,records

def evaluate():
    barrier();assert not (OUT/'stack.npy').exists(),'Never overwrite evaluated evidence'
    a,m,p,coverage,records=build();width=a.shape[1]
    for name,x in [('stack',a),('mask',m),('provenance',p),('coverage',coverage),('s_um',np.arange(width)*DS),('z_indices',np.array(ZS)),('z_relative_um',np.arange(50)*SPACING)]:np.save(OUT/f'{name}.npy',x)
    table(OUT/'path_metrics.csv',records)
    ext=[{'z':r['z'],'relative_z_um':i*SPACING,'path_length_mm':r['total_length_mm'],'row_samples':r['row_samples'],'last_sample_s_um':r['last_sample_s_um'],
          'offset_um':0,'valid_samples':r['valid_samples'],'valid_fraction_native':r['valid_samples']/r['row_samples'],'valid_fraction_rectangle':r['valid_samples']/width,'padding_samples':width-r['row_samples']} for i,r in enumerate(records)]
    table(OUT/'row_extents.csv',ext)
    pairs=[];windows=[];nulls=[];lags=[]
    for i in range(49):
        n=min(records[i]['row_samples'],records[i+1]['row_samples']);x=a[i,:n];y=a[i+1,:n]
        r,nvalid=corr(x,y);dr,nd=corr(np.diff(x),np.diff(y));null=[]
        for shift in np.unique(np.linspace(128,n-128,199).astype(int)):
            rr,nn=corr(x,np.roll(y,int(shift)));null.append(rr);nulls.append({'z_a':ZS[i],'z_b':ZS[i+1],'shift_samples':int(shift),'r':rr,'joint_valid':nn})
        for start in range(0,n-255,128):
            rr,nn=corr(x[start:start+256],y[start:start+256],64)
            windows.append({'z_a':ZS[i],'z_b':ZS[i+1],'s_start_um':start*DS,'s_center_um':(start+127.5)*DS,'r':rr,'joint_valid':nn})
        local_lags=[]
        for lag in range(-64,65):
            xx=x[max(0,-lag):min(n,n-lag)];yy=y[max(0,lag):min(n,n+lag)];rr,nn=corr(xx,yy,128)
            local_lags.append((rr if rr is not None else -2,lag));lags.append({'z_a':ZS[i],'z_b':ZS[i+1],'lag_samples':lag,'lag_um':lag*DS,'r':rr,'joint_valid':nn})
        best=max(local_lags,key=lambda v:(v[0],-abs(v[1])))
        pairs.append({'z_a':ZS[i],'z_b':ZS[i+1],'delta_z_um':SPACING,'shared_native_samples':n,'shared_last_s_um':(n-1)*DS,'joint_valid':nvalid,'r':r,'difference_r':dr,'difference_pairs':nd,
                      'null_median':float(np.median(null)),'null_p95':float(np.percentile(null,95)),'null_fraction_at_least_observed':float(np.mean(np.array(null)>=r)),
                      'best_lag_samples':best[1],'best_lag_um':best[1]*DS,'best_lag_r':best[0],
                      'either_topology_defective':bool(records[i]['final_defect_count'] or records[i+1]['final_defect_count']),
                      'path_length_change_mm':records[i+1]['total_length_mm']-records[i]['total_length_mm'],
                      'seed_distance_px':float(np.hypot(records[i+1]['seed_x']-records[i]['seed_x'],records[i+1]['seed_y']-records[i]['seed_y']))})
    for name,rows in [('adjacent_correlations',pairs),('local_windows',windows),('shift_controls',nulls),('diagnostic_lags',lags)]:table(OUT/f'{name}.csv',rows)
    dump(OUT/'adjacent_correlations.json',pairs)
    raw=np.nan_to_num(a,nan=0).round().clip(0,255).astype(np.uint8);lo,hi=np.percentile(a[~m],[1,99]);contrast=np.zeros_like(raw)
    contrast[~m]=np.rint(np.clip((a[~m]-lo)/(hi-lo),0,1)*255).astype(np.uint8)
    palette=np.array([[0,0,0],[0,191,255],[255,0,255],[255,176,0]],np.uint8)
    for name,v in [('raw_stack',raw),('raw_stack_contrast',contrast),('mask_stack',m.astype(np.uint8)*255),('provenance_stack',palette[p]),('physical_aspect_stack',np.repeat(raw,2,axis=0))]:Image.fromarray(v).save(OUT/f'{name}.png')
    dump(OUT/'display.json',{'raw':'round/clip0..255, missing black','contrast':'one global pooled 1st/99th percentile, visualization only','low':float(lo),'high':float(hi),
                           'physical_aspect':'each original row repeated2 vertically for126.56um/63.28um pixel-spacing ratio; no horizontal resizing, original50 rows in NPY',
                           'overview':'enlarged vertical display for inspection, NOT physical aspect; missing magenta'})
    fig,axs=plt.subplots(3,1,figsize=(18,10),sharex=True)
    cm=plt.cm.gray.copy();cm.set_bad('#bd40ba');extent=[-.5*DS/1000,(width-.5)*DS/1000,177.5,127.5]
    for ax,low,high,title in [(axs[0],0,255,'Fixed CT 0–255; missing magenta'),(axs[1],lo,hi,'Single global contrast window; no row registration')]:
        ax.imshow(np.ma.masked_array(a,m),cmap=cm,vmin=low,vmax=high,aspect='auto',interpolation='nearest',extent=extent);ax.set(title=title,ylabel='Z')
    axs[2].imshow(palette[p],aspect='auto',interpolation='nearest',extent=extent);axs[2].set(title='Black unresolved/padding; cyan observed; magenta bridge; orange low-confidence',ylabel='Z',xlabel='Seed-relative s (mm)')
    fig.tight_layout();fig.savefig(OUT/'stack_inspection.png',dpi=150);plt.close(fig)
    nwin=max(int(round(w['s_start_um']/DS/128)) for w in windows)+1;heat=np.full((49,nwin),np.nan)
    for w in windows:heat[w['z_a']-128,int(round(w['s_start_um']/DS/128))]=w['r']
    fig,axs=plt.subplots(3,1,figsize=(15,11))
    zz=[r['z_a']+.5 for r in pairs];rr=[r['r'] for r in pairs]
    axs[0].plot(zz,rr,'.-',label='raw');axs[0].plot(zz,[r['difference_r'] for r in pairs],'.-',label='first differences');axs[0].plot(zz,[r['null_p95'] for r in pairs],'--',label='shift-control95th percentile');axs[0].axhline(0,c='gray',lw=.5);axs[0].legend();axs[0].set(xlabel='Adjacent-pair midpoint Z',ylabel='Pearson r',title='Every consecutive pair; diagnostic only')
    im=axs[1].imshow(np.ma.masked_invalid(heat),vmin=-1,vmax=1,cmap='coolwarm',aspect='auto',interpolation='nearest',extent=[0,nwin*128*DS/1000,176.5,127.5]);fig.colorbar(im,ax=axs[1],label='local raw r');axs[1].set(xlabel='Window start s (mm); width16.20mm, stride8.10mm',ylabel='Pair starting Z',title='Local valid-overlap correlations; blank outside shared domains')
    axs[2].plot(zz,[r['best_lag_um']/1000 for r in pairs],'.-',label='best bounded diagnostic lag');axs[2].axhline(0,c='gray',lw=.5);axs[2].set(xlabel='Adjacent-pair midpoint Z',ylabel='lag mm, NOT applied',title='Offset sensitivity; weak or multiple peaks are not measured physical shear')
    fig.tight_layout();fig.savefig(OUT/'coherence_diagnostic.png',dpi=150);plt.close(fig)
    # Descriptive persistence: longest consecutive pair run with local r>=.2 in each fixed window.
    runs=[]
    for j in range(nwin):
        best=cur=0;bestend=None
        for i in range(49):
            cur=cur+1 if np.isfinite(heat[i,j]) and heat[i,j]>=.2 else 0
            if cur>best:best=cur;bestend=i+128
        runs.append({'window_start_s_mm':j*128*DS/1000,'longest_pair_run_r_ge_0_2':best,'run_start_pair_z':None if bestend is None else bestend-best+1,'run_end_pair_z':bestend})
    dump(OUT/'local_persistence.json',runs)
    r=np.array(rr);clean=[v['r'] for v in pairs if not v['either_topology_defective']]
    metrics={'row_count':50,'numeric_rows':50,'topology_clean_rows':sum(v['final_defect_count']==0 for v in records),'exterior_rows':sum(v['exterior_reached'] for v in records),
             'clean_exterior_rows':sum(v['exterior_reached'] and v['final_defect_count']==0 for v in records),'shape':list(a.shape),'ds_um':DS,'dz_um':SPACING,'z_center_span_mm':49*SPACING/1000,
             'common_last_s_mm':(width-1)*DS/1000,'all_native_overlap_last_s_mm':(min(v['row_samples'] for v in records)-1)*DS/1000,
             'valid_samples':int((~m).sum()),'valid_fraction':float((~m).mean()),'masked_fraction':float(m.mean()),'bridge_fraction':float((p==2).mean()),
             'padding_fraction':float((~coverage).mean()),'native_unresolved_fraction':float(((p==0)&coverage).mean()),'low_confidence_fraction':float((p==3).mean()),
             'valid_fraction_native':float((~m).sum()/coverage.sum()),'all50_joint_valid_columns':int((~m).all(axis=0).sum()),
             'mean_adjacent_r':float(r.mean()),'median_adjacent_r':float(np.median(r)),'positive_adjacent_fraction':float((r>0).mean()),'raw_r_quantiles':np.percentile(r,[0,10,25,50,75,90,100]).tolist(),
             'mean_difference_r':float(np.mean([v['difference_r'] for v in pairs])),'pairs_above_shift95':sum(v['r']>v['null_p95'] for v in pairs),
             'clean_pair_count':len(clean),'clean_pair_mean_r':float(np.mean(clean)) if clean else None,
             'local_window_count':len(windows),'local_positive_fraction':float(np.mean([w['r']>0 for w in windows if w['r'] is not None])),
             'max_local_pair_run_r_ge_0_2':max(v['longest_pair_run_r_ge_0_2'] for v in runs),
             'length_mm_quantiles':np.percentile([v['total_length_mm'] for v in records],[0,25,50,75,100]).tolist(),
             'observed_retention_range': [min(v['observed_retained_fraction'] for v in records),max(v['observed_retained_fraction'] for v in records)],
             'topology_defective_z':[v['z'] for v in records if v['final_defect_count']],
             'lag_samples_quantiles':np.percentile([v['best_lag_samples'] for v in pairs],[0,25,50,75,100]).tolist()}
    dump(OUT/'stack_metrics.json',metrics);print(json.dumps(metrics,indent=2))

def verify():
    barrier();a,m,p,c,recs=build()
    for name,x in [('stack',a),('mask',m),('provenance',p),('coverage',c)]:assert np.array_equal(x,np.load(OUT/f'{name}.npy'),equal_nan=True)
    assert np.array_equal(np.isnan(a),m) and np.array_equal(m,np.isin(p,[0,2]))
    raw=np.nan_to_num(a,nan=0).round().astype(np.uint8)
    assert np.array_equal(np.array(Image.open(OUT/'raw_stack.png')),raw)
    assert np.array_equal(np.array(Image.open(OUT/'physical_aspect_stack.png')),np.repeat(raw,2,axis=0))
    pairs=json.loads((OUT/'adjacent_correlations.json').read_text());assert len(pairs)==49
    for i,v in enumerate(pairs):
        n=min(recs[i]['row_samples'],recs[i+1]['row_samples']);assert corr(a[i,:n],a[i+1,:n])[0]==v['r']
    historical=[]
    for z in [128,144,146,160,170]:
        same={name:np.array_equal(np.load(OUT/f'slice_{z}/{name}.npy'),np.load(ROOT/f'outputs/p0_t3/run_04/slice_{z}/{name}.npy'),equal_nan=True) for name in ['row','mask','provenance']}
        historical.append({'z':z,**same})
    for name,h in json.loads((OUT/'preservation.json').read_text()).items():assert sha(ROOT/name)==h,name
    for f in OUT.rglob('*.png'):
        with Image.open(f) as im:im.verify()
    dump(OUT/'verification.json',{'status':'PASS','prior_preservation':'PASS','all50_finalization_hashes':'PASS','stack_replay':'exact; no tracing rerun','mask_provenance':'PASS','physical_aspect_png':'PASS','all49_correlations':'exact recomputation','historical_T3_replay':historical})
    print('PASS: row seals, exact stack replay, masks, all49 correlations, PNGs, preservation; historical replay:',historical)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify:verify()
    else:evaluate()
