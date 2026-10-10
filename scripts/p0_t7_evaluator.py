"""T6 metric implementation reused for bounded T7 sets; only variable plot height adapted."""
import json,math
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_t6_features import OUT,T5,FEATURES,ZS,DS,sha,dump,table
from p0_t5_evaluate import corr
from p0_t5_lag_diagnostics import lagbank

def main():
 assert not (OUT/'feature_metrics.json').exists()
 for name,h in json.loads((OUT/'FEATURES_FINALIZED.json').read_text()).items():assert sha(OUT/name)==h,name
 mask=np.load(OUT/'mask.npy');shape=mask.shape;lengths=[len(np.load(T5/f'slice_{z}/row.npy')) for z in ZS]
 allmetrics=[];allpairs={};alllocal={};displays={};cm=plt.cm.gray.copy();cm.set_bad('#bf42ba')
 rng=np.random.default_rng(2606);starts=rng.integers(0,49,size=(2000,10));boot=((starts[:,:,None]+np.arange(5)[None,None,:])%49).reshape(2000,-1)[:,:49]
 np.save(OUT/'bootstrap_pair_indices.npy',boot)
 for name,(definition,display_range) in FEATURES.items():
  a=np.load(OUT/f'{name}_stack.npy');pairs=[];windows=[];lags=[];nulls=[];localcontrols=[]
  for i in range(49):
   n=min(lengths[i],lengths[i+1]);x=a[i,:n];y=a[i+1,:n];r,nvalid=corr(x,y);dr,dn=corr(np.diff(x),np.diff(y));null=[]
   for shift in np.unique(np.linspace(128,n-128,199).astype(int)):
    nr,nn=corr(x,np.roll(y,int(shift)));null.append(nr);nulls.append({'z_a':ZS[i],'shift_samples':int(shift),'r':nr})
   ll=[]
   for lag in range(-64,65):
    rr,nn=corr(x[max(0,-lag):min(n,n-lag)],y[max(0,lag):min(n,n+lag)],128);ll.append((rr if rr is not None else -2,lag));lags.append({'z_a':ZS[i],'lag_samples':lag,'r':rr,'joint_valid':nn})
   best=max(ll,key=lambda q:(q[0],-abs(q[1])))
   pw=[]
   for start in range(0,n-255,128):
    wr,wn=corr(x[start:start+256],y[start:start+256],64);bank,counts=lagbank(x,y,start);j=min(range(129),key=lambda k:(-bank[k] if np.isfinite(bank[k]) else 2,abs(k-64),k));nr=[]
    for shift in np.unique(np.linspace(128,n-128,19).astype(int)):
     nb,_=lagbank(x,np.roll(y,int(shift)),start);peak=float(np.nanmax(nb));nr.append(peak);localcontrols.append({'z_a':ZS[i],'start_sample':start,'shift_samples':int(shift),'peak_r':peak})
    rec={'z_a':ZS[i],'z_b':ZS[i+1],'start_sample':start,'s_start_um':start*DS,'r':wr,'joint_valid':wn,'best_r':float(bank[j]),'best_lag_samples':j-64,'peak_control_p95':float(np.percentile(nr,95)),'peak_above_control95':bool(bank[j]>np.percentile(nr,95))};windows.append(rec);pw.append(rec)
   pairs.append({'z_a':ZS[i],'z_b':ZS[i+1],'r':r,'joint_valid':nvalid,'difference_r':dr,'best_r':best[0],'best_lag_samples':best[1],'null_p95':float(np.percentile(null,95)),
    'local_mean_r':float(np.mean([q['r'] for q in pw])),'local_best_mean_r':float(np.mean([q['best_r'] for q in pw]))})
  rr=np.array([q['r'] for q in pairs]);loc=np.array([q['r'] for q in windows]);locpair=np.array([q['local_mean_r'] for q in pairs])
  if name=='center_intensity':baserr=rr.copy();baseloc=locpair.copy()
  delta=rr-baserr;ldelta=locpair-baseloc;ci=np.percentile(delta[boot].mean(1),[2.5,97.5]);lci=np.percentile(ldelta[boot].mean(1),[2.5,97.5])
  rec={'feature':name,'valid_fraction':float((~mask).mean()),'pairs_evaluable':len(rr),'mean_r':float(rr.mean()),'median_r':float(np.median(rr)),'min_r':float(rr.min()),'max_r':float(rr.max()),'positive_fraction':float((rr>0).mean()),
   'delta_mean_vs_center':float(delta.mean()),'pairs_better_than_center':int((delta>0).sum()),'paired_block_ci_low':float(ci[0]),'paired_block_ci_high':float(ci[1]),
   'difference_mean_r':float(np.mean([q['difference_r'] for q in pairs])),'best_lag_mean_r':float(np.mean([q['best_r'] for q in pairs])),
   'local_windows':len(loc),'local_mean_r':float(loc.mean()),'local_median_r':float(np.median(loc)),'local_positive_fraction':float((loc>0).mean()),
   'delta_local_pair_mean':float(ldelta.mean()),'local_block_ci_low':float(lci[0]),'local_block_ci_high':float(lci[1]),
   'local_best_mean_r':float(np.mean([q['best_r'] for q in windows])),'local_peak_above_control95':sum(q['peak_above_control95'] for q in windows),
   'pairs_above_global_shift95':sum(q['r']>q['null_p95'] for q in pairs)}
  allmetrics.append(rec);allpairs[name]=pairs;alllocal[name]=windows
  for suffix,rows in [('adjacent',pairs),('local_windows',windows),('lags',lags),('shift_controls',nulls),('local_peak_controls',localcontrols)]:table(OUT/f'{name}_{suffix}.csv',rows)
  dump(OUT/f'{name}_evaluation.json',{'metrics':rec,'pairs':pairs})
  def mapped(low,high):
   im=np.zeros(shape,np.uint8);im[~mask]=np.rint(np.clip((a[~mask]-low)/(high-low),0,1)*255).astype(np.uint8);return im
  lo,hi=np.percentile(a[~mask],[1,99]);assert hi>lo
  Image.fromarray(mapped(*display_range)).save(OUT/f'{name}_stack.png');Image.fromarray(mapped(lo,hi)).save(OUT/f'{name}_stack_contrast.png')
  displays[name]={'fixed_range':display_range,'global_contrast_percentiles':[1,99],'contrast_range':[float(lo),float(hi)],'definition':definition}
  nwin=max(q['start_sample']//128 for q in windows)+1;heat=np.full((49,nwin),np.nan)
  for q in windows:heat[q['z_a']-128,q['start_sample']//128]=q['r']
  fig,axs=plt.subplots(3,1,figsize=(14,9));z=np.arange(128.5,177)
  axs[0].plot(z,rr,'.-',label=name);axs[0].plot(z,baserr,color='gray',label='center baseline');axs[0].axhline(0,color='gray',lw=.5);axs[0].legend();axs[0].set(xlabel='Adjacent Z midpoint',ylabel='raw Pearson r',title='No shift applied')
  axs[1].plot(z,[q['best_r'] for q in pairs],'.-',label='diagnostic best lag');axs[1].plot(z,rr,label='zero shift');axs[1].legend();axs[1].set(xlabel='Adjacent Z midpoint',ylabel='Pearson r')
  im=axs[2].imshow(np.ma.masked_invalid(heat),cmap='coolwarm',vmin=-1,vmax=1,aspect='auto',interpolation='nearest',extent=[0,nwin*128*DS/1000,176.5,127.5]);fig.colorbar(im,ax=axs[2]);axs[2].set(xlabel='Window start s mm',ylabel='Pair starting Z',title='Local ZERO-SHIFT correlations')
  fig.tight_layout();fig.savefig(OUT/f'{name}_coherence.png',dpi=130);plt.close(fig)
  print(json.dumps(rec),flush=True)
 table(OUT/'feature_metrics.csv',allmetrics);dump(OUT/'feature_metrics.json',allmetrics);dump(OUT/'display.json',displays)
 # Criterion fixed by names/definitions, not visual letter appearance.
 best=max(allmetrics[1:],key=lambda r:r['mean_r'])['feature'];bestlocal=max(allmetrics[1:],key=lambda r:r['local_mean_r'])['feature']
 dump(OUT/'rankings.json',{'best_adjacent_nonbaseline':best,'best_local_zero_shift_nonbaseline':bestlocal,'warning':'ranking eight related features is exploratory; intervals unadjusted for selection'})
 fig,axs=plt.subplots(3,1,figsize=(13,11));names=list(FEATURES);x=np.arange(len(names))
 axs[0].bar(x,[r['mean_r'] for r in allmetrics]);axs[0].set(xticks=x,xticklabels=names,ylabel='mean adjacent r',title='Authoritative zero-shift feature comparison');axs[0].tick_params(axis='x',rotation=20)
 axs[1].bar(x,[r['local_mean_r'] for r in allmetrics],label='local zero shift');axs[1].set(xticks=x,xticklabels=names,ylabel='mean local r');axs[1].tick_params(axis='x',rotation=20)
 for name in ['center_intensity',best]:axs[2].plot(np.arange(128.5,177),[r['r'] for r in allpairs[name]],'.-',label=name)
 axs[2].legend();axs[2].set(xlabel='Adjacent Z midpoint',ylabel='r');fig.tight_layout();fig.savefig(OUT/'feature_comparison.png',dpi=150);plt.close(fig)
 fig,axs=plt.subplots(len(names),1,figsize=(17,2*len(names)),sharex=True)
 for ax,name in zip(axs,names):
  a=np.load(OUT/f'{name}_stack.npy');lo,hi=displays[name]['contrast_range'];ax.imshow(np.ma.masked_array(a,mask),vmin=lo,vmax=hi,cmap=cm,aspect='auto',interpolation='nearest',extent=[0,(shape[1]-1)*DS/1000,177.5,127.5]);ax.set(title=name,ylabel='Z')
 axs[-1].set_xlabel('Frozen seed-relative s mm; missing magenta; global per-feature display windows');fig.tight_layout();fig.savefig(OUT/'feature_stack_overview.png',dpi=140);plt.close(fig)
 print('ALL FEATURES EVALUATED; no row transformations applied',flush=True)

if __name__=='__main__':main()
