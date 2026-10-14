"""Compare every old-frame control to the unchanged sealed T9 renderer."""
from pathlib import Path
import json
import numpy as np
import p0_t9_infer as reference
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t11'
def main():
    samples=json.loads((OUT/'samples.json').read_text());maximum=0
    for i,p in enumerate(samples):
        with np.load(ROOT/f"outputs/p0_t6/slice_{p['z']}_profiles.npz") as a:
            j=p['j'];expected=reference.patch(p['z'],a['xy'][j],a['tangent'][j],a['normal'][j])
        for sigma in [1,2]:
            with np.load(OUT/f'complete_support_attempt/frames/patch_{i}_{sigma}.npz') as a:error=int(abs(expected.astype(int)-a['old'].astype(int)).max())
            maximum=max(maximum,error)
    assert maximum==0
    (OUT/'old_frame_verification.json').write_text(json.dumps({'status':'PASS','sample_points':len(samples),'paired_patch_arrays_compared':len(samples)*2,'sealed_T9_renderer_max_uint8_error':maximum,'model_inference':False},indent=2)+'\n');print('T9 renderer exact',len(samples),'points',flush=True)
if __name__=='__main__':main()
