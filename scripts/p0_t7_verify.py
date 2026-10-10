"""Frozen-profile lineage and independent T7 feature/metric replay."""
import json,csv
import numpy as np
from PIL import Image
from p0_t7_family import ROOT,OUT,T5,T6,ZS,WIDTHS,OFF,FEATURES,SCALES,sha,dump
from p0_t5_evaluate import corr

def main():
 source=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r');mask=np.load(OUT/'mask.npy');checks=[]
 arrays={n:np.load(OUT/f'scalars/{n}_stack.npy') for n in FEATURES}
 for name in ['mask','provenance','coverage','s_um','z_indices','z_relative_um']:
  assert sha(OUT/f'{name}.npy')==sha(T5/f'{name}.npy')==sha(T6/f'{name}.npy')
 for i,z in enumerate(ZS):
  p=np.load(T6/f'slice_{z}_profiles.npz');g=np.load(T5/f'slice_{z}/path_and_samples.npz');ids=p['sample_indices'];q=p['profile_xy'];assert np.array_equal(ids,np.flatnonzero(~g['mask']))
  assert np.array_equal(p['xy'],g['xy'][ids]) and np.array_equal(q,p['xy'][:,None,:]+OFF[None,:,None]*p['normal'][:,None,:])
  x=np.floor(q[...,0]).astype(int);y=np.floor(q[...,1]).astype(int);dx=q[...,0]-x;dy=q[...,1]-y;im=source[z-128]
  bilinear=im[y,x]*(1-dx)*(1-dy)+im[y,x+1]*dx*(1-dy)+im[y+1,x]*(1-dx)*dy+im[y+1,x+1]*dx*dy
  assert np.allclose(bilinear,p['profiles'],rtol=0,atol=1e-10)
  # Independent formulas using stored profiles, preserving bit-exact source values.
  profile=p['profiles'];assert np.array_equal(arrays['center_intensity'][i,ids],profile[:,3])
  for width,ii in WIDTHS.items():
   pp=profile[:,ii];oo=OFF[ii];radius=max(abs(oo));avg=np.sum(pp,axis=1)/len(ii);minus=pp[:,oo<0];plus=pp[:,oo>0];outer=pp[:,abs(oo)>=radius/2]
   expected={'mean':avg,'median':np.sort(pp,axis=1)[:,len(ii)//2],'std':np.sqrt(np.sum((pp-avg[:,None])**2,axis=1)/len(ii)),
    'range':np.max(pp,axis=1)-np.min(pp,axis=1),'center_minus_mean':profile[:,3]-avg,
    'asymmetry':np.abs(np.sum(plus,axis=1)/plus.shape[1]-np.sum(minus,axis=1)/minus.shape[1]),
    'gradient_abs':np.abs(pp[:,-1]-pp[:,0])*.5/radius,'curvature':(profile[:,3]*2-pp[:,0]-pp[:,-1])/(radius*radius),
    'ridge_proxy':profile[:,3]-np.sum(outer,axis=1)/outer.shape[1]}
   for family,v in expected.items():assert np.allclose(arrays[f'{family}_{width}'][i,ids],v,rtol=0,atol=1e-10)
  checks.append({'z':z,'saved_profile_vs_frozen_XY_and_source':'PASS within1e-10','all28_features':'independent formulas within1e-10'})
 assert np.array_equal(arrays['center_intensity'],np.load(T5/'stack.npy'),equal_nan=True)
 assert np.array_equal(arrays['mean_wide'],np.load(T6/'profile_mean_stack.npy'),equal_nan=True)
 combos=json.loads((OUT/'combination_plan.json').read_text())['combinations'];assert len(combos)<=5
 for c in combos:
  expected=.5*(arrays[c['a']]/c['scale_a']+arrays[c['b']]/c['scale_b']);assert np.array_equal(expected,np.load(OUT/f'combinations/{c["name"]}_stack.npy'),equal_nan=True)
 lengths=[len(np.load(T5/f'slice_{z}/row.npy')) for z in ZS];metric_count=window_count=0;worst=0.
 for phase in ['scalars','combinations']:
  folder=OUT/phase
  for name,h in json.loads((folder/'FEATURES_FINALIZED.json').read_text()).items():assert sha(folder/name)==h,name
  metrics=json.loads((folder/'feature_metrics.json').read_text())
  for m in metrics:
   name=m['feature'];a=np.load(folder/f'{name}_stack.npy');assert np.array_equal(np.isnan(a),mask)
   pairs=json.loads((folder/f'{name}_evaluation.json').read_text())['pairs']
   for i,p in enumerate(pairs):
    n=min(lengths[i],lengths[i+1]);assert corr(a[i,:n],a[i+1,:n])[0]==p['r'];metric_count+=1
   for r in csv.DictReader((folder/f'{name}_local_windows.csv').open()):
    i=int(r['z_a'])-128;s=int(r['start_sample']);lag=int(r['best_lag_samples']);n=min(lengths[i],lengths[i+1]);x=a[i,s:s+256];y=np.pad(a[i+1,:n],(64,64),constant_values=np.nan)
    r0,_=corr(x,y[s+64:s+320],64);rp,_=corr(x,y[s+64+lag:s+320+lag],64)
    assert r0==float(r['r']);err=abs(rp-float(r['best_r']));assert err<1e-10;worst=max(worst,err);window_count+=1
 for f in OUT.rglob('*.png'):
  with Image.open(f) as im:im.verify()
 for name,h in json.loads((OUT/'preservation.json').read_text()).items():assert sha(ROOT/name)==h,name
 dump(OUT/'verification.json',{'status':'PASS','slices':checks,'scalar_feature_count':28,'actual_combinations':len(combos),'T5_center_and_T6_mean':'bit-identical',
  'masks_provenance_coordinates':'byte-identical to both T5 and T6','combinations':'exact equal weights and fixed scales','adjacent_correlations_replayed':metric_count,
  'local_windows_replayed':window_count,'local_selected_peak_max_error':worst,'phase_finalization_seals':'PASS','prior_preservation':'PASS','all_png_decode':'PASS','no_geometry_or_source_modification':True})
 print('PASS: frozen-profile/source replay;28 formulas;3 fixed combinations;',metric_count,'adjacent correlations;',window_count,'windows; preservation')
if __name__=='__main__':main()
