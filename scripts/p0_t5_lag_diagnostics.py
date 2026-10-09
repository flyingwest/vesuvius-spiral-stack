"""Diagnostic local lag search on frozen rows; no transformed raster is produced."""
import json, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_t5_consecutive import OUT,DS,dump,table
from p0_t5_evaluate import barrier,corr

def lagbank(x,y,start):
    from numpy.lib.stride_tricks import sliding_window_view
    xx=x[start:start+256]
    yy=sliding_window_view(np.pad(y,(64,64),constant_values=np.nan),256)[start:start+129]
    good=np.isfinite(yy)&np.isfinite(xx)[None,:];n=good.sum(axis=1)
    a=np.where(good,xx[None,:],0);b=np.where(good,yy,0)
    sx=a.sum(1);sy=b.sum(1);denn=np.maximum(n,1)
    cov=(a*b).sum(1)-sx*sy/denn
    den=np.sqrt(np.maximum((a*a).sum(1)-sx*sx/denn,0)*np.maximum((b*b).sum(1)-sy*sy/denn,0))
    r=np.full(129,np.nan);valid=(n>=64)&(den>0);r[valid]=cov[valid]/den[valid]
    return r,n

def main():
    barrier();assert not (OUT/'local_lag_windows.csv').exists()
    a=np.load(OUT/'stack.npy');records=json.loads((OUT/'reconstruction_summary.json').read_text());pairs=json.loads((OUT/'adjacent_correlations.json').read_text())
    local=[];summ=[]
    for i,p in enumerate(pairs):
        n=p['shared_native_samples'];x=a[i,:n];y=a[i+1,:n];ps=[]
        shifts=np.unique(np.linspace(128,n-128,19).astype(int))
        for start in range(0,n-255,128):
            rr,counts=lagbank(x,y,start);best=min(range(129),key=lambda j:(-rr[j] if np.isfinite(rr[j]) else 2,abs(j-64),j));lag=best-64
            null=[]
            for shift in shifts:
                nr,_=lagbank(x,np.roll(y,int(shift)),start);null.append(float(np.nanmax(nr)))
            # Difference coherence at the same selected raw-signal lag, not a second optimized lag.
            padded=np.pad(y,(64,64),constant_values=np.nan);yy=padded[start+best:start+best+256]
            dr,dn=corr(np.diff(x[start:start+256]),np.diff(yy),32)
            rec={'z_a':p['z_a'],'z_b':p['z_b'],'s_start_um':start*DS,'s_end_um':(start+255)*DS,'zero_r':float(rr[64]),'best_r':float(rr[best]),
                 'best_lag_samples':lag,'best_lag_um':lag*DS,'joint_valid_at_best':int(counts[best]),'difference_r_at_raw_best_lag':dr,'difference_pairs':dn,
                 'peak_control_p95':float(np.percentile(null,95)),'peak_control_fraction_ge_best':float(np.mean(np.array(null)>=rr[best])),
                 'peak_above_control95':bool(rr[best]>np.percentile(null,95)),'lag_at_search_boundary':abs(lag)==64}
            local.append(rec);ps.append(rec)
        summ.append({'z_a':p['z_a'],'z_b':p['z_b'],'zero_r':p['r'],'best_global_r':p['best_lag_r'],'global_lag_samples':p['best_lag_samples'],
                     'local_windows':len(ps),'local_median_zero_r':float(np.median([q['zero_r'] for q in ps])),
                     'local_median_best_r':float(np.median([q['best_r'] for q in ps])),
                     'local_median_lag_samples':float(np.median([q['best_lag_samples'] for q in ps])),
                     'local_peak_above_control95':sum(q['peak_above_control95'] for q in ps)})
    table(OUT/'local_lag_windows.csv',local);table(OUT/'lag_correlations.csv',summ)
    raw=np.array([p['r'] for p in pairs]);peak=np.array([p['best_lag_r'] for p in pairs]);lraw=np.array([q['zero_r'] for q in local]);lpeak=np.array([q['best_r'] for q in local])
    histedges=np.arange(-1,1.0001,.1);rh,_=np.histogram(raw,bins=histedges);ph,_=np.histogram(peak,bins=histedges)
    nwin=max(int(round(q['s_start_um']/DS/128)) for q in local)+1
    hraw=np.full((49,nwin),np.nan);hpeak=hraw.copy();hlag=hraw.copy();hsig=hraw.copy()
    for q in local:
        i=q['z_a']-128;j=int(round(q['s_start_um']/DS/128));hraw[i,j]=q['zero_r'];hpeak[i,j]=q['best_r'];hlag[i,j]=q['best_lag_um']/1000;hsig[i,j]=q['peak_above_control95']
    runs=[]
    for j in range(nwin):
        run=best=0;end=None
        for i in range(49):
            run=run+1 if hsig[i,j]==1 and hpeak[i,j]>=.3 else 0
            if run>best:best=run;end=128+i
        runs.append({'s_start_mm':j*128*DS/1000,'longest_adjacent_pair_run':best,'first_pair_z':None if end is None else end-best+1,'last_pair_z':end})
    metrics={'pairs_evaluable':49,'global_raw_mean':float(raw.mean()),'global_raw_median':float(np.median(raw)),'global_raw_min':float(raw.min()),'global_raw_max':float(raw.max()),
             'global_best_mean':float(peak.mean()),'global_best_median':float(np.median(peak)),'global_best_min':float(peak.min()),'global_best_max':float(peak.max()),
             'global_best_positive_pairs':int((peak>0).sum()),'global_best_boundary_pairs':sum(abs(p['best_lag_samples'])==64 for p in pairs),
             'local_windows':len(local),'local_zero_mean':float(lraw.mean()),'local_best_mean':float(lpeak.mean()),'local_best_median':float(np.median(lpeak)),
             'local_peak_above_control95_count':sum(q['peak_above_control95'] for q in local),
             'local_peak_boundary_count':sum(q['lag_at_search_boundary'] for q in local),
             'difference_at_local_raw_peak_mean':float(np.mean([q['difference_r_at_raw_best_lag'] for q in local if q['difference_r_at_raw_best_lag'] is not None])),
             'histogram_edges':histedges.tolist(),'raw_histogram':rh.tolist(),'best_histogram':ph.tolist(),'persistent_control_exceeding_runs':runs,
             'strong_local_examples':sorted(local,key=lambda q:(not q['peak_above_control95'],-q['best_r']))[:12],
             'low_raw_pairs':[p['z_a'] for p in pairs if p['r']<=0],
             'large_length_discontinuities':sorted([{'z_a':p['z_a'],'z_b':p['z_b'],'length_change_mm':p['path_length_change_mm'],'seed_distance_px':p['seed_distance_px']} for p in pairs],key=lambda q:-abs(q['length_change_mm']))[:10],
             'warning':'Maxima are optimistically selected.19 circular-shift peak controls repeat identical ±64 optimization per window; descriptive, not independent inferential tests. No shifts applied.'}
    dump(OUT/'lag_metrics.json',metrics)
    fig,axs=plt.subplots(3,1,figsize=(15,10))
    z=np.arange(128.5,177)
    axs[0].plot(z,raw,'.-',label='authoritative zero shift');axs[0].plot(z,peak,'.-',label='diagnostic best lag ±64');axs[0].axhline(0,color='gray',lw=.5);axs[0].legend();axs[0].set(ylabel='Pearson r',xlabel='Z pair midpoint',title='No diagnostic shift is applied to the stack')
    axs[1].plot(z,[p['best_lag_um']/1000 for p in pairs],'.-');axs[1].set(ylabel='lag mm',xlabel='Z pair midpoint',title='Apparent lag, not certified physical shear')
    axs[2].hist([raw,peak],bins=np.arange(-.3,.51,.025),label=['zero shift','best lag']);axs[2].legend();axs[2].set(xlabel='Pearson r',ylabel='pair count')
    fig.tight_layout();fig.savefig(OUT/'lag_diagnostic.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(4,1,figsize=(16,13),sharex=True)
    extent=[0,nwin*128*DS/1000,176.5,127.5]
    for ax,arr,title,cm,low,high in [(axs[0],hraw,'Local zero-shift r','coolwarm',-1,1),(axs[1],hpeak,'Local best-lag r (selection-biased)','coolwarm',-1,1),(axs[2],hlag,'Local selected lag mm (not applied)','coolwarm',-64*DS/1000,64*DS/1000),(axs[3],hsig,'Exceeds matched shifted peak-control95: white=yes, black=no','gray',0,1)]:
        cmap=plt.get_cmap(cm).copy();cmap.set_bad('#bbb');im=ax.imshow(np.ma.masked_invalid(arr),vmin=low,vmax=high,cmap=cmap,aspect='auto',interpolation='nearest',extent=extent);fig.colorbar(im,ax=ax);ax.set(ylabel='pair starting Z',title=title)
    axs[-1].set_xlabel('Window start s (mm); 256 samples/window');fig.tight_layout();fig.savefig(OUT/'local_window_coherence.png',dpi=150);plt.close(fig)
    # Reproducible examples: compare values against s, do not write any registered stack.
    examples=metrics['strong_local_examples'][:3]
    fig,axs=plt.subplots(3,1,figsize=(14,9))
    for ax,q in zip(axs,examples):
        i=q['z_a']-128;start=round(q['s_start_um']/DS);lag=q['best_lag_samples'];xx=a[i,start:start+256];yy=np.full(256,np.nan);kk=np.arange(start,start+256)+lag;good=(kk>=0)&(kk<a.shape[1]);yy[good]=a[i+1,kk[good]];ss=(start+np.arange(256))*DS/1000
        ax.plot(ss,xx,lw=.8,label=f'Z{q["z_a"]}');ax.plot(ss,yy,lw=.8,label=f'Z{q["z_b"]} sampled at diagnostic lag{lag}');ax.legend();ax.set(title=f'EXPLORATORY selected example; r={q["best_r"]:.3f}; zero r={q["zero_r"]:.3f}; CT gaps retained',xlabel='s mm',ylabel='original CT')
    fig.tight_layout();fig.savefig(OUT/'local_lag_examples.png',dpi=150);plt.close(fig)
    print(json.dumps({k:v for k,v in metrics.items() if k not in ['strong_local_examples','persistent_control_exceeding_runs','histogram_edges','raw_histogram','best_histogram']},indent=2))

if __name__=='__main__':main()
