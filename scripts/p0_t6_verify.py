"""Independent bilinear/formula replay from immutable T5 coordinates and source CT."""
import json,csv
import numpy as np
from PIL import Image
from p0_t6_features import ROOT,OUT,T5,SOURCE,FEATURES,ZS,OFF,sha,dump,DS,SPACING
from p0_t5_evaluate import corr

def main():
 volume=np.load(SOURCE,mmap_mode='r');mask=np.load(OUT/'mask.npy');results=[];maximum=0.
 for name in ['mask','provenance','coverage','s_um','z_indices','z_relative_um']:
  assert sha(OUT/f'{name}.npy')==sha(T5/f'{name}.npy'),name
 stacks={k:np.load(OUT/f'{k}_stack.npy') for k in FEATURES}
 assert np.array_equal(stacks['center_intensity'],np.load(T5/'stack.npy'),equal_nan=True)
 for name,a in stacks.items():assert np.array_equal(np.isnan(a),mask),name
 for i,z in enumerate(ZS):
  saved=np.load(OUT/f'slice_{z}_profiles.npz');d=np.load(T5/f'slice_{z}/path_and_samples.npz');ii=saved['sample_indices'];xy=saved['xy'];t=saved['tangent'];normal=saved['normal'];q=saved['profile_xy'];profile=saved['profiles'];valid=~d['mask']
  assert np.array_equal(ii,np.flatnonzero(valid)) and np.array_equal(xy,d['xy'][valid])
  assert np.allclose(np.linalg.norm(t,axis=1),1,atol=1e-12) and np.allclose(np.sum(t*normal,axis=1),0,atol=1e-12)
  assert np.array_equal(normal,np.c_[-t[:,1],t[:,0]]) and np.array_equal(q,xy[:,None,:]+OFF[None,:,None]*normal[:,None,:])
  # Independently interpolate the saved frozen edge vertices at both clipped limits.
  es=d['edge_s_px'];vertices=np.vstack([d['edge_start_xy'],d['edge_end_xy'][-1]])
  low=saved['tangent_low_s_px'];high=saved['tangent_high_s_px'];ss=d['s_px'][valid]
  assert (low>=ss-2-1e-9).all() and (high<=ss+2+1e-9).all() and (high>low).all()
  for k in np.flatnonzero(d['edge_provenance']==2):assert not np.any((low<es[k+1]-1e-9)&(high>es[k]+1e-9))
  left=np.column_stack([np.interp(low,es,vertices[:,j]) for j in range(2)]);right=np.column_stack([np.interp(high,es,vertices[:,j]) for j in range(2)]);expected_t=right-left
  fallback=np.linalg.norm(expected_t,axis=1)<1e-8;assert np.array_equal(fallback,saved['fallback'])
  ix=np.minimum(np.searchsorted(es,ss,side='right')-1,len(es)-2)
  expected_t[fallback]=d['edge_end_xy'][ix[fallback]]-d['edge_start_xy'][ix[fallback]];expected_t/=np.linalg.norm(expected_t,axis=1)[:,None]
  assert np.array_equal(expected_t,t)
  x=np.floor(q[...,0]).astype(int);y=np.floor(q[...,1]).astype(int);dx=q[...,0]-x;dy=q[...,1]-y;im=volume[z-128]
  bilinear=im[y,x]*(1-dx)*(1-dy)+im[y,x+1]*dx*(1-dy)+im[y+1,x]*(1-dx)*dy+im[y+1,x+1]*dx*dy
  err=float(np.max(np.abs(bilinear-profile)));maximum=max(maximum,err);assert err<1e-10
  mean=bilinear.sum(1)/7
  independent={'center_intensity':bilinear[:,3],'profile_mean':mean,'profile_std':np.sqrt(((bilinear-mean[:,None])**2).sum(1)/7),
    'profile_range':np.max(bilinear,axis=1)-np.min(bilinear,axis=1),'normal_gradient_abs':np.abs(bilinear[:,4]-bilinear[:,2])*.5,
    'ridge_curvature':bilinear[:,3]*2-bilinear[:,2]-bilinear[:,4],
    'ridge_contrast':bilinear[:,3]-(bilinear[:,0]+bilinear[:,1]+bilinear[:,5]+bilinear[:,6])*.25,
    'profile_asymmetry':np.abs((bilinear[:,4]+bilinear[:,5]+bilinear[:,6]-bilinear[:,0]-bilinear[:,1]-bilinear[:,2])/3)}
  for name,v in independent.items():assert np.allclose(stacks[name][i,ii],v,rtol=0,atol=1e-10),name
  results.append({'z':z,'observed_samples':len(ii),'mask_and_coordinates':'exact','tangent_normal':'PASS','profiles_and_features_independent_bilinear':'within1e-10'})
 for name,h in json.loads((OUT/'FEATURES_FINALIZED.json').read_text()).items():assert sha(OUT/name)==h,name
 for name,h in json.loads((OUT/'preservation.json').read_text()).items():assert sha(ROOT/name)==h,name
 # Evaluation may still be running; this command is intended after its completion.
 metrics=json.loads((OUT/'feature_metrics.json').read_text());display=json.loads((OUT/'display.json').read_text());lengths=[len(np.load(T5/f'slice_{z}/row.npy')) for z in ZS]
 for name,a in stacks.items():
  pairs=json.loads((OUT/f'{name}_evaluation.json').read_text())['pairs'];assert len(pairs)==49
  for i,p in enumerate(pairs):
   n=min(lengths[i],lengths[i+1]);assert corr(a[i,:n],a[i+1,:n])[0]==p['r']
  for suffix,limits in [('stack',display[name]['fixed_range']),('stack_contrast',display[name]['contrast_range'])]:
   lo,hi=limits;im=np.zeros(a.shape,np.uint8);im[~mask]=np.rint(np.clip((a[~mask]-lo)/(hi-lo),0,1)*255).astype(np.uint8)
   assert np.array_equal(im,np.array(Image.open(OUT/f'{name}_{suffix}.png')))
 for f in OUT.glob('*.png'):
  with Image.open(f) as im:im.verify()
 base=json.loads((T5/'stack_metrics.json').read_text());assert metrics[0]['mean_r']==base['mean_adjacent_r']
 dump(OUT/'verification.json',{'status':'PASS','slices':results,'max_bilinear_error':maximum,'feature_formulas':'independent arithmetic within1e-10','T5_center_stack':'bit-identical','masks_provenance_coordinates':'byte-identical copies',
  'all392_adjacent_correlations':'exact replay','feature_seal':'PASS','prior_preservation_including_T5':'PASS','PNGs':'exact fixed/global mappings and all decode','no_search_rerun':True,'no_T5_artifact_modified':True})
 print('PASS:50x8 feature replay;392 correlations; exact T5 baseline/masks/provenance; feature and prior hashes; all PNGs')
if __name__=='__main__':main()
