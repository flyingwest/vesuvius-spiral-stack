"""Preserve predecessor evidence and select geometry-only sample identities."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs/p0_t11'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    assert not OUT.exists()
    prior=json.loads((ROOT/'outputs/p0_t10/preservation.json').read_text())
    for line in (ROOT/'outputs/p0_t10/SHA256SUMS').read_text().splitlines():
        h,p=line.split('  ',1);prior[p]=h
    prior['outputs/p0_t10/SHA256SUMS']=sha(ROOT/'outputs/p0_t10/SHA256SUMS')
    for p,h in prior.items():assert sha(ROOT/p)==h,p
    samples=[]
    for group in [[128,129,130],[151,152,153],[175,176,177]]:
        rows={}
        for z in group:
            with np.load(ROOT/f'outputs/p0_t6/slice_{z}_profiles.npz') as f:rows[z]={k:f[k] for k in ['sample_indices','xy','tangent','normal']}
        common=sorted(set.intersection(*[set(r['sample_indices'].tolist()) for r in rows.values()]))
        chosen=np.array(common)[np.linspace(0,len(common)-1,8,dtype=int)]
        for z in group:
            r=rows[z]
            for s in chosen:
                j=int(np.flatnonzero(r['sample_indices']==s)[0]);samples.append({'z':z,'s_index':int(s),'j':j,'P_xyz':[*list(r['xy'][j]*16),z*16],'tangent_xyz':[*list(r['tangent'][j]),0.],'normal_xyz':[*list(r['normal'][j]),0.]})
    OUT.mkdir();(OUT/'preservation.json').write_text(json.dumps(prior,indent=2)+'\n');(OUT/'samples.json').write_text(json.dumps(samples,indent=2)+'\n')
    print('predecessors',len(prior),'fixed samples',len(samples))
if __name__=='__main__':main()
