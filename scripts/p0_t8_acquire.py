"""Bounded HTTP-range acquisition; no label-dependent crop/depth selection."""
from pathlib import Path
import urllib.request, json, hashlib, io, concurrent.futures
import numpy as np
import tifffile
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs/p0_t8'; DATA=ROOT/'data/raw/p0_t8'
BASE='https://dl.ash2txt.org/fragments/Frag3/PHercParis1Fr34.volpkg/working/54keV_exposed_surface/'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def dump(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def fetch(url,path,span=None):
 headers={} if span is None else {'Range':f'bytes={span[0]}-{span[1]-1}'}
 with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=120) as r:
  assert r.status==(200 if span is None else 206),(url,r.status)
  b=r.read();meta={'url':url,'status':r.status,'headers':dict(r.headers),'range_half_open':span,'bytes':len(b)}
 if span:assert len(b)==span[1]-span[0]
 path.write_bytes(b);meta.update(path=str(path.relative_to(ROOT)),sha256=sha(path));return meta

def main():
 OUT.mkdir(exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
 prior={str(p.relative_to(ROOT)):sha(p) for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents and 'p0_t8' not in p.name}
 prior.update({x:sha(ROOT/x) for x in ['README.md','AGENTS.md']})
 if not (OUT/'preservation.json').exists():dump(OUT/'preservation.json',prior)
 plan={'dataset':'EduceLab-Scrolls Frag3 PHercParis1Fr34; 54keV exposed surface; uuid20230301213423','base_url':BASE,'shape_dyx':[65,7606,5249],'voxel_um':3.24,'crop_yx_half_open':[[2784,4832],[1600,3648]],'center_depth':32,'offsets':[-4,-2,-1,0,1,2,4],'max_payload_bytes':220000000,
 'features':['center_intensity','mean_wide','median_wide','std_wide','range_wide','center_minus_mean_wide','asymmetry_wide','gradient_abs_narrow','curvature_narrow','ridge_proxy_wide'],
 'definitions':'unchanged scripts/p0_t7_family.py compute; gradient/curvature narrow reproduce T6; intensity native uint16 without label-based normalization; same pixel offsets not same physical scale as T5',
 'grid':'16px pitch, centers at crop+8; exclude papyrus mask distance<=16 and either class-boundary distance<=2; nonink = label0 inside mask, potentially incomplete labels',
 'discovery':'top1024 crop rows; only select score sign per feature by AUC; no fitted formula, cutoff, depth or feature selection for confirmatory tests',
 'confirmation':'bottom1024 crop rows; evaluate every feature with discovery sign; four fixed512x512 subregions',
 'matching':'greedy deterministic positive grid order, negative without replacement, same discovery/confirmation half and same256x256 tile, distance<=128px; second match additionally normals dot>=0.98, baseline gap<=0.25 pooled SD, ridge-strength gap<=0.5 pooled SD; closest eligible',
 'nuisances':'source result.ppm surface normals; baseline Gaussian sigma32 center image; spatial gradient magnitude Gaussian sigma2 center image; controls only not new features',
 'metrics':'ROC AUC, average precision (step-integral PR AUC), standardized mean difference pooled SD, class quantiles; matched paired direction concordance and matched-pool ROC/AP; spatial tile bootstrap256px 2000 seed2808, 95%CI and99.5%CI for10-feature Bonferroni family; no iid pixel confidence intervals',
 'registration':'documented official surface registration, exact dimensions/meta, label containment, IR overlay, mask support; quantify IR darkness contrast sanity check only, no registration optimization against CT/ink metrics',
 'decision':'A: discovery-direction confirmation AUC>=0.60 with99.5% lowerCI>0.5, at least3/4 confirm regions above chance, nuisance-matched concordance95%lower>0.5; B: same significance/region/control support but effect below0.60; C: no feature meets replicated controlled evidence. Insufficient sampling/registration blocks completion rather than being biological C.',
 'limitations':'single fragment; chosen surface depth may miss ink; original label registration uncertainty; no inference about T5 coarsened geometry or actual ink chemistry; no neural training'}
 if (OUT/'plan.json').exists():assert json.loads((OUT/'plan.json').read_text())==plan
 else:
  dump(OUT/'plan.json',plan);dump(OUT/'PLAN_FROZEN.json',{'plan_sha256':sha(OUT/'plan.json'),'feature_code_sha256':sha(ROOT/'scripts/p0_t7_family.py')})
 records=[]
 for name in ['surface_volume/meta.json','inklabels.png','mask.png','ir.png']:
  records.append(fetch(BASE+name,DATA/name.replace('/','_')))
 for name,url in [('registration_tutorial.html','https://scrollprize.org/tutorial3'),('dataset_paper.html','https://arxiv.org/html/2304.02084v1')]:
  try:records.append(fetch(url,OUT/name))
  except urllib.error.HTTPError as e:
   dump(OUT/(name+'.error.json'),{'url':url,'status':e.code,'note':'Documentation accessible via web tool; data acquisition continues.'})
 y0,y1=plan['crop_yx_half_open'][0];x0,x1=plan['crop_yx_half_open'][1];width=5249
 # Every TIFF is uncompressed single strip. Read first8 bytes plus complete trailing IFD, verify layout before payload.
 tasks=[]
 for d in [28,30,31,32,33,34,36]:
  url=BASE+f'surface_volume/{d:02}.tif'
  with urllib.request.urlopen(urllib.request.Request(url,method='HEAD')) as r: length=int(r.headers['Content-Length'])
  a=fetch(url,DATA/f'{d}_header.bin',(0,8));records.append(a)
  header=(DATA/f'{d}_header.bin').read_bytes();assert header[:4]==b'II*\x00';ifd=int.from_bytes(header[4:],'little')
  records.append(fetch(url,DATA/f'{d}_ifd.bin',(ifd,length)))
  # sparse local representation only for metadata parsing
  with open(DATA/f'{d}_metadata.tif','wb') as f:f.write(header);f.seek(ifd);f.write((DATA/f'{d}_ifd.bin').read_bytes())
  with tifffile.TiffFile(DATA/f'{d}_metadata.tif') as t:
   p=t.pages[0];assert p.shape==(7606,5249) and p.dtype==np.dtype('uint16') and p.compression==1 and p.dataoffsets==(8,)
  tasks.append((url,DATA/f'{d}_rows.bin',(8+2*width*y0,8+2*width*y1)))
 ppm=BASE+'result.ppm';records.append(fetch(ppm,DATA/'ppm_header.bin',(0,200)));b=(DATA/'ppm_header.bin').read_bytes();start=b.index(b'<>\n')+3
 assert b'width: 5249' in b and b'height: 7606' in b
 for y in range(y0+8,y1,16):tasks.append((ppm,DATA/f'ppm_row_{y}.bin',(start+48*(width*y+x0),start+48*(width*y+x1))))
 assert sum(t[2][1]-t[2][0] for t in tasks)+sum(r['bytes'] for r in records)<plan['max_payload_bytes']
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  for rec in ex.map(lambda t:fetch(*t),tasks):records.append(rec)
 ct=np.stack([np.fromfile(DATA/f'{d}_rows.bin','<u2').reshape(y1-y0,width)[:,x0:x1] for d in [28,30,31,32,33,34,36]])
 np.save(DATA/'ct_profiles_dyx.npy',ct)
 geom=np.stack([np.fromfile(DATA/f'ppm_row_{y}.bin','<f8').reshape(2048,6)[8::16] for y in range(y0+8,y1,16)])
 np.save(DATA/'grid_geometry.npy',geom)
 dump(OUT/'acquisition.json',{'records':records,'network_payload_bytes':sum(r['bytes'] for r in records),'crop_shape':list(ct.shape),'grid_geometry_shape':list(geom.shape),'data_hashes':{str(p.relative_to(ROOT)):sha(p) for p in DATA.iterdir() if p.is_file()}})
 print('Acquired',ct.shape,'bytes',sum(r['bytes'] for r in records),flush=True)
if __name__=='__main__':main()
