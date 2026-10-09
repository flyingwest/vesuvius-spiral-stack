"""Frozen-geometry local CT profiles; no reconstruction or T5 writes."""
import argparse,csv,hashlib,json,shutil,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t6';T5=ROOT/'outputs/p0_t5';SOURCE=ROOT/'outputs/orientation/ct_level4_zyx.npy'
ZS=list(range(128,178));DS=63.28;SPACING=126.56;OFF=np.array([-4,-2,-1,0,1,2,4.])
FEATURES={
 'center_intensity':('P(0)',(0,255)),
 'profile_mean':('arithmetic mean of all7 profile samples; equal point weights, not spatial integral',(0,255)),
 'profile_std':('population standard deviation of all7 profile samples',(0,127.5)),
 'profile_range':('maximum minus minimum of all7 profile samples',(0,255)),
 'normal_gradient_abs':('abs(P(+1)-P(-1))/2, CT units/pixel',(0,127.5)),
 'ridge_curvature':('2P(0)-P(-1)-P(+1), negative normal second derivative, CT units/pixel^2',(-510,510)),
 'ridge_contrast':('P(0)-mean(P(-4),P(-2),P(+2),P(+4))',(-255,255)),
 'profile_asymmetry':('abs(mean(P(+1),P(+2),P(+4))-mean(P(-1),P(-2),P(-4)))',(0,255))}

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def table(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def values(p):
 return {'center_intensity':p[:,3],'profile_mean':p.mean(1),'profile_std':p.std(1),'profile_range':np.ptp(p,axis=1),
  'normal_gradient_abs':np.abs(p[:,4]-p[:,2])/2,'ridge_curvature':2*p[:,3]-p[:,2]-p[:,4],
  'ridge_contrast':p[:,3]-p[:,[0,1,5,6]].mean(1),'profile_asymmetry':np.abs(p[:,4:].mean(1)-p[:,:3].mean(1))}

def geometry(z):
 d=np.load(T5/f'slice_{z}/path_and_samples.npz');a=d['edge_start_xy'];b=d['edge_end_xy'];es=d['edge_s_px'];labels=d['edge_provenance'];valid=~d['mask'];ss=d['s_px'][valid];xy=d['xy'][valid]
 idx=np.minimum(np.searchsorted(es,ss,side='right')-1,len(a)-1)
 observed=np.isin(labels,[1,3]);assert observed[idx].all()
 starts=np.flatnonzero(observed&~np.r_[False,observed[:-1]]);ends=np.flatnonzero(observed&~np.r_[observed[1:],False])
 runstart=np.zeros(len(a));runend=np.zeros(len(a))
 for first,last in zip(starts,ends):runstart[first:last+1]=es[first];runend[first:last+1]=es[last+1]
 low=np.maximum(ss-2,runstart[idx]);high=np.minimum(ss+2,runend[idx]);verts=np.vstack([a,b[-1]])
 left=np.c_[np.interp(low,es,verts[:,0]),np.interp(low,es,verts[:,1])];right=np.c_[np.interp(high,es,verts[:,0]),np.interp(high,es,verts[:,1])]
 tangent=right-left;n=np.linalg.norm(tangent,axis=1);fallback=n<1e-8;tangent[fallback]=b[idx[fallback]]-a[idx[fallback]];tangent/=np.linalg.norm(tangent,axis=1)[:,None]
 normal=np.c_[-tangent[:,1],tangent[:,0]];coords=xy[:,None,:]+OFF[None,:,None]*normal[:,None,:]
 return d,valid,xy,tangent,normal,coords,fallback,low,high

def extract():
 OUT.mkdir(exist_ok=False)
 for line in (T5/'SHA256SUMS').read_text().splitlines():
  h,p=line.split('  ',1);assert sha(ROOT/p)==h,p
 prior=[p for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents]+[ROOT/'README.md',ROOT/'AGENTS.md']
 dump(OUT/'preservation.json',{str(p.relative_to(ROOT)):sha(p) for p in prior})
 dump(OUT/'plan.json',{'experiment':'P0-T6 signal extraction ONLY','slices':ZS,'source':str(SOURCE.relative_to(ROOT)),'source_sha256':sha(SOURCE),
  'frozen_paths':'outputs/p0_t5/slice_<z>/path_and_samples.npz','ds_um':DS,'voxel_spacing_um':SPACING,'normal_offsets_px':OFF.tolist(),
  'normal_offsets_um':(OFF*SPACING).tolist(),'tangent':'chord between s±2px clipped to contiguous observed edge-run; fallback to containing edge if chord norm<1e-8; normal=(-ty,tx)',
  'interpolation':'bilinear same-Z original CT, no Z neighbors; assert all offsets in bounds, no extrapolation',
  'mask':'exact T5 mask and provenance; query CT only at samples1/3; never query inferred centers or fill bridges',
  'features':{k:{'definition':v[0],'fixed_display_range':v[1]} for k,v in FEATURES.items()},
  'evaluation':'49 pairs Pearson; 256-sample windows stride128 min64 valid; diagnostic lags±64; 199 global circular controls; 19 local peak controls; no shifts applied',
  'comparison':'paired gains vs center; circular moving-block bootstrap length5 adjacent pairs,2000 replicates fixed seed2606; descriptive95% intervals;8 features exploratory multiplicity',
  'interpretation':'materiality guides: mean adjacent gain~0.03 or local mean gain~0.05 plus consistency across pairs/controls; not ink identification; no feature retuning',
  'scope':'do not edit T5, no search, cross-Z geometry, learned model, OCR, text priors, mesh; freeze definitions before feature evidence'})
 volume=np.load(SOURCE,mmap_mode='r');shape=np.load(T5/'stack.npy',mmap_mode='r').shape;stacks={k:np.full(shape,np.nan) for k in FEATURES};records=[]
 for name in ['mask','provenance','coverage','s_um','z_indices','z_relative_um']:shutil.copyfile(T5/f'{name}.npy',OUT/f'{name}.npy')
 for i,z in enumerate(ZS):
  d,valid,xy,tangent,normal,coords,fallback,low,high=geometry(z);image=volume[z-128]
  assert (coords[...,0]>=0).all() and (coords[...,0]<image.shape[1]-1).all() and (coords[...,1]>=0).all() and (coords[...,1]<image.shape[0]-1).all()
  profile=map_coordinates(image.astype(float),[coords[...,1].ravel(),coords[...,0].ravel()],order=1,mode='nearest').reshape(-1,7)
  assert np.array_equal(profile[:,3],d['ct'][valid]),'center must exactly reproduce T5'
  for name,v in values(profile).items():stacks[name][i,np.flatnonzero(valid)]=v
  np.savez_compressed(OUT/f'slice_{z}_profiles.npz',sample_indices=np.flatnonzero(valid),xy=xy,tangent=tangent,normal=normal,profile_xy=coords,profiles=profile,tangent_low_s_px=low,tangent_high_s_px=high,fallback=fallback)
  records.append({'z':z,'valid_samples':int(valid.sum()),'profile_queries':int(valid.sum()*7),'fallback_tangents':int(fallback.sum()),'clipped_tangent_windows':int(((low>d['s_px'][valid]-2)|(high<d['s_px'][valid]+2)).sum())})
 mask=np.load(OUT/'mask.npy')
 for name,a in stacks.items():
  assert np.array_equal(np.isnan(a),mask);np.save(OUT/f'{name}_stack.npy',a)
 table(OUT/'sampling_metrics.csv',records)
 dump(OUT/'FEATURES_FINALIZED.json',{str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.glob('*.np*'))})
 # Examples fixed by positions, not by coherence or feature values.
 fig,axs=plt.subplots(2,3,figsize=(15,8));fig2,bxs=plt.subplots(2,3,figsize=(15,9));examples=[]
 for ax,bx,(z,fraction) in zip(axs.flat,bxs.flat,[(128,.1),(128,.5),(152,.1),(152,.5),(177,.1),(177,.5)]):
  p=np.load(OUT/f'slice_{z}_profiles.npz');k=int(fraction*(len(p['sample_indices'])-1));q=p['profile_xy'][k];center=p['xy'][k];t=p['tangent'][k];n=p['normal'][k]
  ax.plot(OFF*SPACING,p['profiles'][k],'o-');ax.axvline(0,c='gray',lw=.5);ax.set(title=f'Z{z}, sample{p["sample_indices"][k]}',xlabel='normal offset µm',ylabel='original CT',ylim=(0,255))
  bx.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest');bx.plot(q[:,0],q[:,1],'m.-',label='normal profile');bx.plot([center[0]-2*t[0],center[0]+2*t[0]],[center[1]-2*t[1],center[1]+2*t[1]],'c-',label='tangent');bx.scatter(*center,c='lime',s=20)
  bx.set(xlim=(center[0]-10,center[0]+10),ylim=(center[1]+10,center[1]-10),title=f'Z{z}, sampled contour neighborhood');bx.legend(fontsize=7)
  examples.append({'z':z,'valid_index_fraction':fraction,'sample_index':int(p['sample_indices'][k]),'xy':center.tolist()})
 fig.tight_layout();fig.savefig(OUT/'profile_examples.png',dpi=140);plt.close(fig);fig2.tight_layout();fig2.savefig(OUT/'normal_sampling_diagnostic.png',dpi=140);plt.close(fig2);dump(OUT/'example_selection.json',examples)
 print('PASS:8 features from frozen50 paths; exact center baseline; masks/provenance unchanged; features finalized')

if __name__=='__main__':extract()
