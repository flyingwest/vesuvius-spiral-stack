"""Saved-evaluation checks; no signal/geometry changes."""
import csv,json
import numpy as np
from PIL import Image
from p0_t6_features import OUT,T5,FEATURES,ZS,dump
from p0_t5_evaluate import corr

def main():
 metrics=json.loads((OUT/'feature_metrics.json').read_text());boot=np.load(OUT/'bootstrap_pair_indices.npy');baseline=None;baselocal=None;checks=[]
 for name in FEATURES:
  a=np.load(OUT/f'{name}_stack.npy');pairs=json.loads((OUT/f'{name}_evaluation.json').read_text())['pairs'];rr=np.array([v['r'] for v in pairs]);lp=np.array([v['local_mean_r'] for v in pairs])
  if baseline is None:baseline=rr;baselocal=lp
  delta=rr-baseline;ld=lp-baselocal;ci=np.percentile(delta[boot].mean(1),[2.5,97.5]);lci=np.percentile(ld[boot].mean(1),[2.5,97.5]);m=next(v for v in metrics if v['feature']==name)
  assert np.array_equal(ci,[m['paired_block_ci_low'],m['paired_block_ci_high']]);assert np.array_equal(lci,[m['local_block_ci_low'],m['local_block_ci_high']])
  rows=list(csv.DictReader((OUT/f'{name}_local_windows.csv').open()));errors=[]
  for r in rows:
   i=int(r['z_a'])-128;start=int(r['start_sample']);lag=int(r['best_lag_samples']);n=min(len(np.load(T5/f'slice_{ZS[i]}/row.npy')),len(np.load(T5/f'slice_{ZS[i+1]}/row.npy')))
   x=a[i,start:start+256];y=np.pad(a[i+1,:n],(64,64),constant_values=np.nan)
   raw,_=corr(x,y[start+64:start+64+256],64);peak,_=corr(x,y[start+64+lag:start+64+lag+256],64)
   assert raw==float(r['r']);errors.append(abs(peak-float(r['best_r'])))
  assert max(errors)<1e-10
  checks.append({'feature':name,'paired_bootstrap':'exact','local_zero_and_selected_peak_replay':'PASS','window_count':len(rows),'max_peak_error':max(errors)})
 im=np.array(Image.open(OUT/'center_intensity_stack.png'));old=np.array(Image.open(T5/'raw_stack.png'))
 dump(OUT/'metric_verification.json',{'status':'PASS','features':checks,'center_fixed_png_differing_pixels':int((im!=old).sum()),'center_numeric_stack':'bit-identical as verified independently','interpretation':'local selected peaks checked, not independent reselection of all129 lags'})
 print('PASS: paired bootstrap and all8008 local windows; center PNG differing pixels:',int((im!=old).sum()))
if __name__=='__main__':main()
