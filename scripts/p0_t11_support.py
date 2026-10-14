"""Infrastructure-only support inventory/acquisition and isolated frame retry."""
from pathlib import Path
import hashlib,json,sys,urllib.request,itertools
import numpy as np
from numcodecs import Blosc
import p0_t11_core as c
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t11';DATA=ROOT/'data/raw/p0_t11/native_chunks'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    stage=sys.argv[1]
    if stage=='details':
        plan=json.loads((OUT/'support_plan.json').read_text());items={tuple(a['key']):{**a,'reason':'rectangular trilinear support of unchanged rotated patch','dependent_samples':[]} for a in plan['items']};bounds=[]
        samples=json.loads((OUT/'samples.json').read_text());records=json.loads((OUT/'frames/records.json').read_text())
        for i,p in enumerate(samples):
            P,t,n=[np.array(p[k]) for k in ['P_xyz','tangent_xyz','normal_xyz']]
            for sigma in [1.,2.]:
                f=next(r for r in records if r['sample']==i and r['sigma']==sigma)
                if 'n' not in f:f=c.frame(P,t,n,sigma)
                corners=np.array([P+np.asarray(f['n'])*d*c.STEP+np.asarray(f['t_s'])*s*c.STEP+np.asarray(f['image_y'])*h*c.STEP for d,s,h in itertools.product([-8,8],[-63.5,63.5],[-63.5,63.5])])
                lo=np.floor(corners.min(0)).astype(int)[::-1];hi=(np.floor(corners.max(0)).astype(int)+2)[::-1]
                keys=list(itertools.product(*[range(a//128,(b-1)//128+1) for a,b in zip(lo,hi)]));missing=[]
                for key in keys:
                    if not (ROOT/'data/raw/p0_t9/native_chunks'/'_'.join(map(str,key))).exists():
                        assert key in items;missing.append(list(key));items[key]['dependent_samples'].append({'sample':i,'z':p['z'],'s_index':p['s_index'],'sigma':sigma})
                bounds.append({'sample':i,'z':p['z'],'s_index':p['s_index'],'sigma':sigma,'native_P_xyz':p['P_xyz'],'bounds_zyx_half_open':[lo.tolist(),hi.tolist()],'required_chunks_zyx':[list(k) for k in keys],'missing_from_T9':missing})
        assert all(r['dependent_samples'] for r in items.values())
        source=json.loads((ROOT/'outputs/p0_t9/acquisition.json').read_text())['source']
        result={'source_volume':source,'native_voxel_um':7.91,'level':0,'chunk_shape_zyx':[128,128,128],'required_patches':len(bounds),'already_present_chunks':len(list((ROOT/'data/raw/p0_t9/native_chunks').iterdir())),'additional_chunks':len(items),'additional_bytes_expected':plan['total_bytes'],'missing_chunks':[{'source_volume':source,**r} for r in items.values()],'patches':bounds,'support_semantics':'Exact bounding CT reads including trilinear upper-neighbor halo; no clipping or padding'}
        (OUT/'missing_ct_inventory.json').write_text(json.dumps(result,indent=2)+'\n');print('inventory',len(bounds),'patches',len(items),'missing chunks');return
    if stage=='inventory':
        needed=set();samples=json.loads((OUT/'samples.json').read_text())
        for p in samples:
            P,t,n=[np.array(p[k]) for k in ['P_xyz','tangent_xyz','normal_xyz']]
            for sigma in [1.,2.]:
                f=c.frame(P,t,n,sigma)
                # Bounding support exactly matches renderer's rectangular CT read.
                extent=(8*abs(f['n'])+63.5*abs(f['t_s'])+63.5*abs(f['image_y']))*c.STEP
                lo=np.floor(P-extent).astype(int)[::-1];hi=(np.floor(P+extent).astype(int)+2)[::-1]
                for key in itertools.product(*[range(a//128,(b-1)//128+1) for a,b in zip(lo,hi)]):
                    if not (ROOT/'data/raw/p0_t9/native_chunks'/'_'.join(map(str,key))).exists():needed.add(key)
        base=json.loads((ROOT/'outputs/p0_t9/acquisition.json').read_text())['source'];items=[]
        for key in sorted(needed):
            url=base+'0/'+'/'.join(map(str,key))
            with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=60) as r:headers=dict(r.headers);size=int(r.headers['Content-Length'])
            items.append({'key':list(key),'url':url,'bytes':size,'headers':headers})
        total=sum(r['bytes'] for r in items);assert total<=20_000_000
        (OUT/'support_plan.json').write_text(json.dumps({'cap_bytes':20_000_000,'total_bytes':total,'items':items},indent=2)+'\n');print(json.dumps({'chunks':len(items),'bytes':total}));return
    if stage=='acquire':
        plan=json.loads((OUT/'support_plan.json').read_text());DATA.mkdir(parents=True,exist_ok=False);receipt=[]
        for item in plan['items']:
            path=DATA/'_'.join(map(str,item['key']))
            with urllib.request.urlopen(item['url'],timeout=60) as r:payload=r.read(plan['cap_bytes']+1)
            assert len(payload)==item['bytes'];path.write_bytes(payload);decoded=Blosc().decode(payload);assert len(decoded)==128**3
            receipt.append({**item,'filename':str(path.relative_to(ROOT)),'sha256':sha(path)})
        (OUT/'support_acquisition.json').write_text(json.dumps(receipt,indent=2)+'\n')
        (OUT/'supplemental_ct_manifest.json').write_text(json.dumps({'source_volume':json.loads((OUT/'missing_ct_inventory.json').read_text())['source_volume'],'chunks':receipt,'additional_chunks':len(receipt),'additional_bytes':sum(r['bytes'] for r in receipt),'verification':'Expected HEAD size, Blosc decoded shape, local SHA256. Server ETags are recorded opaque validators, not published cryptographic digests.'},indent=2)+'\n');print('acquired',len(receipt),'supplemental chunks');return
    if stage=='validate':
        inventory=json.loads((OUT/'missing_ct_inventory.json').read_text());receipt=json.loads((OUT/'supplemental_ct_manifest.json').read_text())
        for item in receipt['chunks']:
            p=ROOT/item['filename'];assert p.stat().st_size==item['bytes'] and sha(p)==item['sha256']
        uncovered=[]
        for patch in inventory['patches']:
            missing=[key for key in patch['required_chunks_zyx'] if not (ROOT/'data/raw/p0_t9/native_chunks'/'_'.join(map(str,key))).exists() and not (DATA/'_'.join(map(str,key))).exists()]
            if missing:uncovered.append({'sample':patch['sample'],'sigma':patch['sigma'],'missing':missing})
        result={'status':'PASS' if not uncovered else 'BLOCKED','REQUIRED_PATCHES':inventory['required_patches'],'FULLY_COVERED_PATCHES':inventory['required_patches']-len(uncovered),'uncovered':uncovered,'supplemental_hashes_verified':len(receipt['chunks']),'no_clipping_or_padding':True}
        (OUT/'coverage_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));assert not uncovered;return
    original=c.chunk
    def supplemented(k,j,i):
        path=ROOT/f'data/raw/p0_t9/native_chunks/{k}_{j}_{i}'
        if path.exists():return original(k,j,i)
        return np.frombuffer(Blosc().decode((DATA/f'{k}_{j}_{i}').read_bytes()),np.uint8).reshape(128,128,128)
    c.chunk=supplemented
    if stage in ['frames','replay']:
        import p0_t11_frames as runner
        runner.OUT=OUT/'complete_support_attempt'
        if stage=='frames':
            runner.OUT.mkdir(exist_ok=False)
            (runner.OUT/'samples.json').write_bytes((OUT/'samples.json').read_bytes())
        sys.argv=['frames']+(['--replay'] if stage=='replay' else [])
        runner.main();return
    if stage in ['smoke','smoke_replay','full','full_replay']:
        import p0_t11_infer as runner
        runner.OUT=OUT/'complete_support_attempt'
        sys.argv=['infer']+(['--full'] if stage.startswith('full') else [])+(['--replay'] if stage.endswith('replay') else [])
        runner.main();return
    raise ValueError(stage)
if __name__=='__main__':main()
