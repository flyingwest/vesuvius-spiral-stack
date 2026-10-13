"""Gate-controlled smoke/full inference, no geometry or model optimization."""
from pathlib import Path
from types import SimpleNamespace
import json,sys,hashlib
import numpy as np
import torch
from vesuvius.ink_detection.inference.infer import configure_model,normalize_flat_patch,logits_to_probabilities
from p0_t11_core import frame,render
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t11'
def main():
    full='--full' in sys.argv;replay='--replay' in sys.argv
    meta=json.loads((OUT/'frame_metrics.json').read_text());sigma=meta['selected_sigma'];assert sigma is not None,'geometry gate closed'
    if full:assert json.loads((OUT/'smoke_metrics.json').read_text())['status']=='PASS'
    checkpoint=ROOT/'data/raw/p0_t9/step-075000.pth'
    assert hashlib.sha256(checkpoint.read_bytes()).hexdigest()=='e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab'
    torch.set_num_threads(4);torch.manual_seed(2911);torch.backends.cudnn.benchmark=False
    free,total=torch.cuda.mem_get_info();c=configure_model(SimpleNamespace(checkpoint=checkpoint,amp_dtype='auto'));model=c.model.cuda().eval()
    def predict(raw):
        x=np.stack([normalize_flat_patch(a,c.preprocessing) for a in raw])[:,None]
        with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):p=logits_to_probabilities(model(torch.from_numpy(x).cuda()),image_hw=(128,128))[:,0].cpu().numpy()
        assert np.isfinite(p).all() and p.min()>=0 and p.max()<=1;return p
    if not full:
        samples=json.loads((OUT/'samples.json').read_text());raw=[]
        records=json.loads((OUT/'frames/records.json').read_text());accepted={r['sample']:r['accepted'] for r in records if r['sigma']==sigma}
        for i in range(len(samples)):
            with np.load(OUT/f'frames/patch_{i}_{int(sigma)}.npz') as a:raw.extend([a['old'],a['new'] if accepted[i] else a['old']])
        maps=np.concatenate([predict(raw[i:i+8]) for i in range(0,len(raw),8)]);centers=maps[:,63:65,63:65].mean((1,2)).reshape(-1,2)
        if replay:
            original=np.load(OUT/'smoke_probability_maps.npy');error=float(abs(maps-original).max());assert error<=.002
            (OUT/'inference_replay.json').write_text(json.dumps({'status':'PASS','patches':len(raw),'maximum_probability_error':error,'scope':'all smoke patches; full-stack predictions, if present, are separately mapped to saved coordinates'},indent=2)+'\n');return
        np.save(OUT/'smoke_probability_maps.npy',maps);np.save(OUT/'smoke_centers.npy',centers)
        adjacent={}
        for k,name in enumerate(['T9','T11']):
            corrs=[]
            for group in range(3):
                a=centers[group*24:(group+1)*24,k].reshape(3,8)
                for j in range(2):corrs.append(float(np.corrcoef(a[j],a[j+1])[0,1]))
            adjacent[name]={'mean':float(np.mean(corrs)),'median':float(np.median(corrs)),'pairs':corrs}
        result={'status':'PASS','finite_fraction':1.,'cuda_free_before_load_bytes':free,'probability_range':[float(maps.min()),float(maps.max())],'center_quantiles':{name:np.quantile(centers[:,k],[0,.25,.5,.75,1]).tolist() for k,name in enumerate(['T9','T11'])},'adjacent_Z':adjacent,'selection':'geometry-only; model output does not choose estimator','patches':len(raw)}
        assert np.isfinite(centers).all() and np.std(centers[:,1])>0
        (OUT/'smoke_metrics.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result),flush=True);return
    if replay:
        maximum=0.;replayed=0;stack=np.load(OUT/'ink_probability.npy');mask=np.load(ROOT/'outputs/p0_t5/mask.npy');assert np.array_equal(np.isnan(stack),mask)
        for z in range(128,178):
            with np.load(ROOT/f'outputs/p0_t6/slice_{z}_profiles.npz') as a:r={k:a[k] for k in a.files}
            with np.load(OUT/f'full/slice_{z}.npz') as a:saved={k:a[k] for k in a.files}
            assert np.array_equal(saved['xy'],r['xy']) and np.array_equal(saved['sample_indices'],r['sample_indices'])
            assert np.array_equal(saved['probability'].astype(np.float32),stack[z-128,r['sample_indices']])
            choices=np.flatnonzero(saved['accepted']);j=int(choices[0]) if len(choices) else 0
            P=np.r_[r['xy'][j]*16,z*16];t=np.r_[r['tangent'][j],0.];n=np.r_[r['normal'][j],0.];f=frame(P,t,n,sigma)
            y=f['image_y'] if f['accepted'] else np.array([0,0,1]);nn=f['n'] if f['accepted'] else n
            assert np.array_equal(np.stack([t,y,nn]),saved['frames'][j])
            prob=float(predict([render(P,t,y,nn)])[:,63:65,63:65].mean());error=abs(prob-float(saved['probability'][j]));maximum=max(maximum,error);replayed+=1
        assert maximum<=.002
        (OUT/'full_replay.json').write_text(json.dumps({'status':'PASS','coordinate_rows_exact':50,'centers_replayed':replayed,'maximum_probability_error':maximum},indent=2)+'\n');return
    dest=OUT/'full';dest.mkdir(exist_ok=False);mask=np.load(ROOT/'outputs/p0_t5/mask.npy');stack=np.full(mask.shape,np.nan,np.float32)
    for z in range(128,178):
        with np.load(ROOT/f'outputs/p0_t6/slice_{z}_profiles.npz') as a:r={k:a[k] for k in a.files}
        frames=[];probs=[];accepted=[];raw=[]
        for j in range(len(r['sample_indices'])):
            P=np.r_[r['xy'][j]*16,z*16];t=np.r_[r['tangent'][j],0.];n=np.r_[r['normal'][j],0.];f=frame(P,t,n,sigma)
            y=f['image_y'] if f['accepted'] else np.array([0,0,1]);nn=f['n'] if f['accepted'] else n
            raw.append(render(P,t,y,nn));frames.append(np.stack([t,y,nn]));accepted.append(f['accepted'])
            if len(raw)==8 or j==len(r['sample_indices'])-1:
                p=predict(raw);probs.extend(p[:,63:65,63:65].mean((1,2)).tolist());raw=[]
        stack[z-128,r['sample_indices']]=probs
        np.savez_compressed(dest/f'slice_{z}.npz',sample_indices=r['sample_indices'],xy=r['xy'],frames=frames,accepted=accepted,probability=probs)
        print('completed Z',z,'centers',len(probs),flush=True)
    assert np.array_equal(np.isnan(stack),mask);np.save(OUT/'ink_probability.npy',stack)
if __name__=='__main__':main()
