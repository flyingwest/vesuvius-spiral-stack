"""Fixed pretrained readout at frozen T5 centers; no path search or cross-Z fitting."""
from pathlib import Path
from functools import lru_cache
from types import SimpleNamespace
import json,time,hashlib,concurrent.futures,sys
import numpy as np
import torch
from scipy.ndimage import map_coordinates
from numcodecs import Blosc
from vesuvius.ink_detection.inference.infer import configure_model,normalize_flat_patch,logits_to_probabilities
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t9';DATA=ROOT/'data/raw/p0_t9/native_chunks'
STEP=9.6/7.91
@lru_cache(maxsize=512)
def chunk(k,j,i):
 return np.frombuffer(Blosc().decode((DATA/f'{k}_{j}_{i}').read_bytes()),dtype=np.uint8).reshape(128,128,128)
def patch(z,xy,t,n,offset=0):
 d=(np.arange(17,dtype=float)-8+offset)[:,None,None]*STEP
 s=(np.arange(128,dtype=float)-63.5)[None,None,:]*STEP
 h=(np.arange(128,dtype=float)-63.5)[None,:,None]*STEP
 Z=np.broadcast_to(z*16+h,(17,128,128));Y=np.broadcast_to(xy[1]*16+s*t[1]+d*n[1],Z.shape);X=np.broadcast_to(xy[0]*16+s*t[0]+d*n[0],Z.shape)
 coords=[Z,Y,X];lo=np.array([int(np.floor(c.min())) for c in coords]);hi=np.array([int(np.floor(c.max()))+2 for c in coords]);volume=np.empty(hi-lo,dtype=np.uint8)
 for k in range(lo[0]//128,(hi[0]-1)//128+1):
  for j in range(lo[1]//128,(hi[1]-1)//128+1):
   for i in range(lo[2]//128,(hi[2]-1)//128+1):
    origin=np.array([k,j,i])*128;start=np.maximum(lo,origin);stop=np.minimum(hi,origin+128)
    volume[tuple(slice(a,b) for a,b in zip(start-lo,stop-lo))]=chunk(k,j,i)[tuple(slice(a,b) for a,b in zip(start-origin,stop-origin))]
 sampled=map_coordinates(volume,np.array([c-a for c,a in zip(coords,lo)]),order=1,prefilter=False,mode='constant',cval=np.nan,output=np.float32)
 assert np.isfinite(sampled).all()
 return np.clip(np.rint(sampled),0,255).astype(np.uint8)
def main():
 # Partial acquisition is allowed only for an explicitly bounded row batch;
 # each required chunk must already exist, otherwise sampling fails visibly.
 stop=int(sys.argv[1]) if len(sys.argv)>1 else 178
 if stop==178:assert (OUT/'acquisition.json').exists()
 torch.set_num_threads(4);torch.manual_seed(2909);torch.backends.cudnn.benchmark=False
 c=configure_model(SimpleNamespace(checkpoint=ROOT/'data/raw/p0_t9/step-075000.pth',amp_dtype='auto'));model=c.model.cuda().eval()
 parameters={'model_spacing_um':9.6,'native_spacing_um':7.91,'native_step':STEP,'source_transform':'scale-only native=16*level4, no fitted offset','patch_dyx':[17,128,128],'axes':'D=XY normal; Y=scan Z; X=XY tangent','tangent_normal':'exact saved T6 local frames, no cross-Z fitting','sampling':'scipy trilinear order1 float32, round half-even to uint8; no additional anti-alias filter','depth_center':0,'amp':'float16','normalization':c.preprocessing,'prediction':'official sigmoid/interpolate; mean of four central output pixels at even-grid center','batch_size':8,'row_sampling_um':63.28,'normal_sign_primary':'T6 saved normal=(-ty,tx)','sensitivity':'32 evenly spaced valid centers on each ofZ128,152,177; normal reversal and depth offsets -2,+2; no selection or replacement of primary result','display':'clip((p-0.25)/0.5,0,1)','limitations':'straight tangent-plane approximation over1.22mm; no longitudinal tilt estimation; native center registration not calibrated; no proof centers lie on ink-bearing surface'}
 param=OUT/'inference_parameters.json'
 if param.exists():assert json.loads(param.read_text())==parameters
 else:param.write_text(json.dumps(parameters,indent=2)+'\n')
 shape=np.load(ROOT/'outputs/p0_t5/mask.npy').shape;stack=np.full(shape,np.nan,dtype=np.float32)
 receipt=OUT/'orientation_depth_sensitivity.json'
 probe=json.loads(receipt.read_text()) if receipt.exists() else []
 start=time.perf_counter()
 def predict(raw):
  x=np.stack([normalize_flat_patch(p,c.preprocessing) for p in raw])[:,None]
  with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):
   q=logits_to_probabilities(model(torch.from_numpy(x).cuda()),image_hw=(128,128))
  return q[:,0,63:65,63:65].mean((1,2)).cpu().numpy()
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for z in range(128,stop):
   dest=OUT/f'slice_{z}_predictions.npz'
   # Materialize the unchanged frame once: shared NPZ decompression is not thread-safe.
   with np.load(ROOT/f'outputs/p0_t6/slice_{z}_profiles.npz') as archive:
    p={key:archive[key] for key in archive.files}
   idx=p['sample_indices'];saved=[]
   if dest.exists():
    row=np.load(dest);assert np.array_equal(row['sample_indices'],idx);stack[z-128,idx]=row['probability'];print('resume',z,flush=True);continue
   check=set(np.linspace(0,len(idx)-1,32,dtype=int)) if z in [128,152,177] else set()
   for first in range(0,len(idx),8):
    ids=list(range(first,min(first+8,len(idx))));raw=list(pool.map(lambda j:patch(z,p['xy'][j],p['tangent'][j],p['normal'][j]),ids));val=predict(raw);saved.extend(val.tolist())
    for j,r in zip(ids,raw):
     if j in check:
      alt=[r[::-1].copy(),patch(z,p['xy'][j],p['tangent'][j],p['normal'][j],-2),patch(z,p['xy'][j],p['tangent'][j],p['normal'][j],2)]
      q=predict(alt);probe.append({'z':z,'sample_index':int(idx[j]),'primary':float(val[ids.index(j)]),'reverse':float(q[0]),'minus2':float(q[1]),'plus2':float(q[2]),'patch_sha256':hashlib.sha256(r.tobytes()).hexdigest()})
      np.savez_compressed(OUT/f'probe_{z}_{idx[j]}.npz',raw_patch=r,xy=p['xy'][j],tangent=p['tangent'][j],normal=p['normal'][j],probabilities=np.r_[val[ids.index(j)],q])
   values=np.array(saved,dtype=np.float32);stack[z-128,idx]=values;np.savez_compressed(dest,sample_indices=idx,probability=values,xy=p['xy'],tangent=p['tangent'],normal=p['normal'])
   (OUT/'orientation_depth_sensitivity.json').write_text(json.dumps(probe,indent=2)+'\n')
   print('completed Z',z,'valid',len(idx),'elapsed',round(time.perf_counter()-start,1),flush=True)
 if stop<178:
  print('Bounded row batch complete; full stack not yet finalized',flush=True)
  return
 np.save(OUT/'ink_probability.npy',stack)
 mask=np.load(ROOT/'outputs/p0_t5/mask.npy');assert np.array_equal(np.isnan(stack),mask);assert np.all((stack[~mask]>=0)&(stack[~mask]<=1))
 for name in ['mask','provenance','s_um','z_indices','z_relative_um']:
  import shutil;shutil.copyfile(ROOT/f'outputs/p0_t5/{name}.npy',OUT/f'{name}.npy')
 print('ALL50 COMPLETE',flush=True)
if __name__=='__main__':main()
