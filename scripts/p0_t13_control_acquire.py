"""Selection/acquisition only: trusted geometry may select bounds, never trace inputs."""
from pathlib import Path
import json,hashlib,urllib.request,itertools,concurrent.futures
import numpy as np,tifffile
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t13_control';D=R/'data/raw/p0_t13_control'
BASE='https://vesuvius-challenge-open-data.s3.amazonaws.com/PHercParis4/volumes/20260411134726-2.400um-0.2m-78keV-masked.zarr/'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 O.mkdir(exist_ok=False);D.mkdir(exist_ok=False)
 prior=json.loads((R/'outputs/p0_t13/preservation.json').read_text())
 for folder in ['outputs/p0_t13','docs']:
  for p in (R/folder).rglob('*'):
   if p.is_file():prior[str(p.relative_to(R))]=sha(p)
 for p,h in prior.items():assert sha(R/p)==h,p
 dump(O/'preservation.json',prior)
 p=np.stack([tifffile.imread(R/f'data/raw/p0_t10/geometry_{a}.tif') for a in 'xyz'],-1)
 # Same documented approximate UV bounds as T10, including one grid boundary.
 region=p[819:922,1049:1153];valid=region[(region>0).all(-1)];center=int(np.floor(np.median(valid[:,2])+.5));zs=list(range(center-3,center+4))
 low=(np.floor((valid.min(0)[:2]-128)/128)*128).astype(int);high=(np.ceil((valid.max(0)[:2]+128)/128)*128).astype(int)
 metadata={}
 for name in ['.zattrs','0/.zarray']:
  with urllib.request.urlopen(BASE+name,timeout=30) as r:payload=r.read()
  dest=D/('zattrs.json' if name=='.zattrs' else 'zarray.json');dest.write_bytes(payload);metadata[name]=json.loads(payload)
 assert metadata['0/.zarray']['chunks']==[128]*3 and metadata['0/.zarray']['compressor'] is None
 selection={'status':'FROZEN_BEFORE_RECONSTRUCTION','source':BASE,'source_scan_id':'20260411134726','t10_source_internal_path':'/volumes/esrf/20260311/2.4um_PHerc-Paris4_masked.zarr/','axes':'native ZYX, transverse XY, positive Z slice order','native_spacing_um':2.4,'native_z_indices':zs,'count':7,'central_evaluable_z':zs[2:-2],'native_xy_bounds_half_open':[low.tolist(),high.tolist()],'trusted_crop_grid_bounds_yx':[819,922,1049,1153],'trusted_crop_xyz_bounds':[valid.min(0).tolist(),valid.max(0).tolist()],'rule':'7 consecutive native slices centered at rounded median Z of sealed T10 crop; entire crop projected XY bbox plus128 native pixels, rounded outward to128. Geometry selects acquisition bounds ONLY.','xy_preprocessing':'exact4x4 mean with half-even uint8 rounding; pixel(i,j) maps native XY=(x0+4*j+1.5,y0+4*i+1.5); no Z pooling','trace_spacing_xy_um':9.6,'trace_slice_spacing_um':2.4,'labels_used':False,'ct_acquisition_cap_bytes':128_000_000,'registration':'same scan-ID native XYZ identity; physical error not independently calibrated'}
 dump(O/'selected_slice_range.json',selection)
 jobs=[]
 for k in range(zs[0]//128,zs[-1]//128+1):
  a=max(zs[0],k*128)-k*128;b=min(zs[-1]+1,(k+1)*128)-k*128
  for j,i in itertools.product(range(low[1]//128,high[1]//128),range(low[0]//128,high[0]//128)):
   jobs.append({'chunk_zyx':[k,j,i],'url':BASE+f'0/{k}/{j}/{i}','range_start':int(a*128**2),'range_end_inclusive':int(b*128**2-1),'local_depth':[a,b]})
 assert sum(q['range_end_inclusive']-q['range_start']+1 for q in jobs)<selection['ct_acquisition_cap_bytes'];dump(O/'acquisition_plan.json',jobs)
 def get(q):
  req=urllib.request.Request(q['url'],headers={'Range':f"bytes={q['range_start']}-{q['range_end_inclusive']}"})
  with urllib.request.urlopen(req,timeout=45) as r:
   assert r.status==206;payload=r.read(128**3+1);headers=dict(r.headers)
  assert len(payload)==q['range_end_inclusive']-q['range_start']+1
  name='chunk_'+'_'.join(map(str,q['chunk_zyx']))+'.bin';path=D/name;path.write_bytes(payload)
  return {**q,'file':str(path.relative_to(R)),'bytes':len(payload),'sha256':sha(path),'headers':headers}
 with concurrent.futures.ThreadPoolExecutor(8) as pool:receipt=list(pool.map(get,jobs))
 dump(O/'acquisition.json',receipt)
 volume=np.empty((7,int(high[1]-low[1]),int(high[0]-low[0])),np.uint8)
 for q in receipt:
  k,j,i=q['chunk_zyx'];a,b=q['local_depth'];arr=np.frombuffer((R/q['file']).read_bytes(),np.uint8).reshape(b-a,128,128)
  volume[k*128+a-zs[0]:k*128+b-zs[0],j*128-low[1]:(j+1)*128-low[1],i*128-low[0]:(i+1)*128-low[0]]=arr
 dest=O/'ct_only_inputs';dest.mkdir();files=[]
 for z,a in zip(zs,volume):
  pooled=np.rint(a.reshape(a.shape[0]//4,4,a.shape[1]//4,4).mean((1,3))).astype(np.uint8);path=dest/f'slice_{z}.npy';np.save(path,pooled);files.append({'file':path.name,'sha256':sha(path),'shape':list(pooled.shape)})
 dump(O/'ct_input_manifest.json',{'files':files,'source_payload_bytes':sum(q['bytes'] for q in receipt),'selection_sha256':sha(O/'selected_slice_range.json')});print(json.dumps(selection,indent=2));print('download bytes',sum(q['bytes'] for q in receipt))
if __name__=='__main__':main()
