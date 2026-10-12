"""Three real-CT smoke probes after complete download verification."""
from p0_t9_infer import *
def main():
 acquisition=json.loads((OUT/'acquisition.json').read_text());assert len(acquisition['chunks'])==5523
 for r in acquisition['chunks']:
  f=ROOT/r['filename'];assert f.stat().st_size==r['expected_bytes'];assert hashlib.sha256(f.read_bytes()).hexdigest()==r['sha256']
 checkpoint=ROOT/'data/raw/p0_t9/step-075000.pth';assert hashlib.sha256(checkpoint.read_bytes()).hexdigest()=='e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab'
 c=configure_model(SimpleNamespace(checkpoint=checkpoint,amp_dtype='auto'));c.model.cuda().eval();torch.set_num_threads(4)
 probes=[]
 for z,jj in [(128,0),(152,.5),(177,1)]:
  p=np.load(ROOT/f'outputs/p0_t6/slice_{z}_profiles.npz');j=int(jj*(len(p['xy'])-1));raw=patch(z,p['xy'][j],p['tangent'][j],p['normal'][j]);norm=normalize_flat_patch(raw,c.preprocessing);assert np.isfinite(norm).all()
  # Independent eight-corner trilinear checks at fixed patch grid positions.
  maxerr=0
  for dd,hh,ss in [(0,0,0),(8,63,63),(8,64,64),(16,127,127),(4,32,96)]:
   xyz=np.r_[p['xy'][j]*16,z*16];xyz[:2]+=(ss-63.5)*STEP*p['tangent'][j]+(dd-8)*STEP*p['normal'][j];xyz[2]+=(hh-63.5)*STEP
   coord=xyz[::-1];low=np.floor(coord).astype(int);fraction=coord-low;value=0.
   for dz in [0,1]:
    for dy in [0,1]:
     for dx in [0,1]:
      corner=np.array([dz,dy,dx]);q=low+corner;weight=np.prod(np.where(corner,fraction,1-fraction));value+=weight*float(chunk(*(q//128))[tuple(q%128)])
   err=abs(float(raw[dd,hh,ss])-float(np.rint(np.float32(value))));maxerr=max(maxerr,err)
  assert maxerr<=1
  with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):q=logits_to_probabilities(c.model(torch.from_numpy(norm[None,None]).cuda()),image_hw=(128,128)).cpu().numpy()
  assert np.isfinite(q).all() and q.min()>=0 and q.max()<=1
  np.savez_compressed(OUT/f'smoke_{z}.npz',raw_patch=raw,normalized=norm,probability_map=q,xy=p['xy'][j],tangent=p['tangent'][j],normal=p['normal'][j],sample_index=p['sample_indices'][j])
  probes.append({'z':z,'sample_index':int(p['sample_indices'][j]),'ct_min':int(raw.min()),'ct_max':int(raw.max()),'prob_min':float(q.min()),'prob_max':float(q.max()),'prob_center':float(q[0,0,63:65,63:65].mean()),'independent_interpolation_max_uint8_error':maxerr})
 result={'status':'PASS','download_files_verified':len(acquisition['chunks']),'checkpoint_hash_verified':True,'shape':[17,128,128],'model_spacing_um':9.6,'native_spacing_um':7.91,'probes':probes,'warning':'Finite outputs establish runtime compatibility, not ink localization or correct physical surface placement.'}
 (OUT/'smoke_test.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
