"""Read-only native chunk budgeting and frozen-input preservation."""
from pathlib import Path
import numpy as np,json,urllib.request,concurrent.futures,hashlib,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t9'
BASE='https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def dump(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2)+'\n')
def main():
 prior={str(p.relative_to(ROOT)):sha(p) for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents and not p.name.startswith('p0_t9_')}
 prior.update({p:sha(ROOT/p) for p in ['README.md','AGENTS.md']});dump('preservation.json',prior)
 keys=set();counts=[]
 for z in range(128,178):
  p=np.load(ROOT/f'outputs/p0_t6/slice_{z}_profiles.npz');xy=p['xy']*16;b=abs(p['tangent'])*78+abs(p['normal'])*13+2;counts.append(len(xy))
  for (x,y),(bx,by) in zip(xy,b):
   for k in range(int((z*16-79)//128),int((z*16+79)//128)+1):
    for j in range(int((y-by)//128),int((y+by)//128)+1):
     for i in range(int((x-bx)//128),int((x+bx)//128)+1):keys.add((k,j,i))
 dump('sampling_preflight.json',{'source':BASE,'native_spacing_um':7.91,'frozen_geometry_spacing_um':126.56,'geometry_native_transform':'native index =16*level4 index, as published scale-only .zattrs; voxel-center offset is not independently calibrated','model_spacing_um':9.6,'model_native_step':9.6/7.91,'patch_dyx':[17,128,128],'patch_axes':'depth=normal in scanXY; modelY=scanZ; modelX=frozen tangent in scanXY','normal':'reuse exact T6 normal=(-ty,tx), centered depth0 primary; no cross-Z tangent estimation','center':'even128 patch geometric center63.5, center prediction bilinear mean of four central output pixels','small_diagnostic':'fixed subset only: reverse normal and depth+/-2 model voxels; no image-driven selection','native_chunk_count':len(keys),'decoded_chunk_bytes':len(keys)*128**3,'valid_centers':sum(counts),'no_CT_chunks_downloaded_yet':True})
 for suffix in ['.zattrs','0/.zarray']:
  (OUT/('native_'+suffix.replace('/','_').replace('.','')+'.json')).write_bytes(urllib.request.urlopen(BASE+suffix).read())
 def head(k):
  url=BASE+'0/'+'/'.join(map(str,k))
  try:
   with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=30) as r:return {'key':list(k),'status':r.status,'bytes':int(r.headers['Content-Length']),'etag':r.headers.get('ETag')}
  except Exception as e:return {'key':list(k),'error':str(e)}
 records=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=16) as ex:
  for rec in ex.map(head,sorted(keys)):
   records.append(rec)
   if len(records)%500==0:print('HEAD',len(records),flush=True)
 dump('native_chunk_budget.json',{'chunks':records,'compressed_bytes':sum(r.get('bytes',0) for r in records),'errors':sum('error' in r for r in records)})
 print('Budget complete',sum(r.get('bytes',0) for r in records),flush=True)
if __name__=='__main__':main()
