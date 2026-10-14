"""Preserve all predecessor evidence before T12 fitting."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t12'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    assert not OUT.exists()
    prior=json.loads((ROOT/'outputs/p0_t11/preservation.json').read_text())
    for line in (ROOT/'outputs/p0_t11/SHA256SUMS').read_text().splitlines():
        h,p=line.split('  ',1);prior[p]=h
    prior['outputs/p0_t11/SHA256SUMS']=sha(ROOT/'outputs/p0_t11/SHA256SUMS')
    for p,h in prior.items():assert sha(ROOT/p)==h,p
    OUT.mkdir();(OUT/'preservation.json').write_text(json.dumps(prior,indent=2)+'\n')
    (OUT/'samples.json').write_bytes((ROOT/'outputs/p0_t11/samples.json').read_bytes())
    print('preserved',len(prior),'files; same72sample identities')
if __name__=='__main__':main()
