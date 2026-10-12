"""Evaluate finalized probabilities without changing geometry or predictions."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np
from PIL import Image
from scipy.ndimage import label,find_objects
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import p0_t7_evaluator as evaluator
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t9'
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 p=np.load(OUT/'ink_probability.npy');mask=np.load(OUT/'mask.npy');prov=np.load(OUT/'provenance.npy');assert np.array_equal(np.isnan(p),mask)
 dest=OUT/'evaluation';dest.mkdir(exist_ok=False)
 shutil.copyfile(ROOT/'outputs/p0_t5/stack.npy',dest/'center_intensity_stack.npy');shutil.copyfile(OUT/'ink_probability.npy',dest/'ink_probability_stack.npy')
 for name in ['mask','provenance','s_um','z_indices','z_relative_um']:shutil.copyfile(OUT/f'{name}.npy',dest/f'{name}.npy')
 dump(dest/'FEATURES_FINALIZED.json',{a.name:hashlib.sha256(a.read_bytes()).hexdigest() for a in dest.glob('*.npy')})
 evaluator.OUT=dest;evaluator.FEATURES={'center_intensity':('frozen T5 baseline',(0,255)),'ink_probability':('fixed official pretrained model',(0,1))};evaluator.main()
 raw=np.zeros(p.shape,np.uint8);raw[~mask]=np.clip(p[~mask]*255,0,255).astype(np.uint8);Image.fromarray(raw).save(OUT/'ink_probability_stack_raw.png')
 display=np.zeros(p.shape,np.uint8);display[~mask]=np.clip((p[~mask]-.25)/.5*255,0,255).astype(np.uint8);Image.fromarray(display).save(OUT/'ink_probability_stack_display.png');Image.fromarray(display).resize((p.shape[1],p.shape[0]*2),Image.Resampling.NEAREST).save(OUT/'physical_aspect_stack.png')
 Image.fromarray((mask*255).astype(np.uint8)).save(OUT/'mask_stack.png');palette=np.array([[0,0,0],[60,190,110],[220,80,210],[240,185,50]],np.uint8);Image.fromarray(palette[prov]).save(OUT/'provenance_stack.png')
 def components(a):
  binary=np.isfinite(a)&(a>.5);labs,n=label(binary);objects=find_objects(labs);sizes=np.bincount(labs.ravel())[1:];spans=np.array([s[0].stop-s[0].start for s in objects],dtype=int)
  return {'above_half_samples':int(binary.sum()),'components':n,'components_spanning3rows':int((spans>=3).sum()),'samples_in_components_spanning3rows':int(sizes[spans>=3].sum()),'largest_vertical_span':int(spans.max()) if n else 0,'largest_component_samples':int(sizes.max()) if n else 0}
 observed=components(p);rng=np.random.default_rng(2909);null=[]
 for _ in range(199):
  shuffled=np.stack([np.roll(row,int(rng.integers(128,p.shape[1]-128))) for row in p]);null.append(components(shuffled))
 dump(OUT/'component_controls.json',{'threshold':.5,'connectivity':4,'observed':observed,'independent_row_circular_shift_controls':null,'warning':'Overlapping native CT neighborhoods can induce model persistence; components are not proof of ink or text.'})
 old=json.loads((ROOT/'outputs/p0_t6/feature_metrics.json').read_text());base=[r for r in old if r['feature']=='profile_mean'][0];new=json.loads((dest/'feature_metrics.json').read_text());ink=new[1]
 metrics={'shape':list(p.shape),'valid_inference_count':int((~mask).sum()),'valid_inference_fraction':float((~mask).mean()),'fraction_valid_T5_centers_inferred':1.,'masked_fraction':float(mask.mean()),'probability_quantiles':np.quantile(p[~mask],[0,.01,.25,.5,.75,.99,1]).tolist(),'T5_center':new[0],'T6_profile_mean':base,'T9_ink':ink,'components':observed,'component_null95':{k:float(np.percentile([a[k] for a in null],95)) for k in observed},'Z_indices':[128,177],'Z_spacing_um':126.56,'s_spacing_um':63.28,'scientific_result':'PENDING_VISUAL_REVIEW','ink_structure_visible':'PENDING_VISUAL_REVIEW'}
 dump(OUT/'metrics.json',metrics)
 fig,axs=plt.subplots(3,1,figsize=(17,8))
 for ax,name,arr in [(axs[0],'T5 raw CT',np.load(ROOT/'outputs/p0_t5/stack.npy')/255),(axs[1],'T6 profile mean',np.load(ROOT/'outputs/p0_t6/profile_mean_stack.npy')/255),(axs[2],'T9 official display (p−0.25)/0.5',(p-.25)/.5)]:
  ax.imshow(np.ma.masked_invalid(arr),cmap='gray',vmin=0,vmax=1,aspect='auto',interpolation='nearest');ax.set_ylabel('Z row');ax.set_title(name)
 fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=160);plt.close(fig)
 print(json.dumps(metrics,indent=2),flush=True)
if __name__=='__main__':main()
