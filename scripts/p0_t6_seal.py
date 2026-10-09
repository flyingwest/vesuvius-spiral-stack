"""Preserve prior evidence, seal T6, immediately verify every manifest entry."""
import json,subprocess,platform,shutil,sys
import numpy,scipy,matplotlib,PIL
from p0_t6_features import ROOT,OUT,T5,FEATURES,sha,dump

def main():
 manifest=OUT/'SHA256SUMS';assert not manifest.exists(),'Already sealed; use sha256sum -c'
 for name in ['verification.json','metric_verification.json']:assert json.loads((OUT/name).read_text())['status']=='PASS'
 for name,h in json.loads((OUT/'FEATURES_FINALIZED.json').read_text()).items():assert sha(OUT/name)==h,name
 for name,h in json.loads((OUT/'preservation.json').read_text()).items():assert sha(ROOT/name)==h,name
 for line in (T5/'SHA256SUMS').read_text().splitlines():
  h,name=line.split('  ',1);assert sha(ROOT/name)==h,name
 for feature in FEATURES:
  for suffix in ['_stack.png','_stack_contrast.png','_coherence.png','_stack.npy']:assert (OUT/(feature+suffix)).is_file()
 for name in ['feature_comparison.png','feature_metrics.csv','feature_metrics.json','profile_examples.png','normal_sampling_diagnostic.png','verification.json','mask.npy','provenance.npy']:assert (OUT/name).is_file()
 dump(OUT/'environment.json',{'python':sys.version,'platform':platform.platform(),'numpy':numpy.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'PIL':PIL.__version__})
 shutil.copyfile('/home/wpotech/.hermes/cache/scratch/p0_t6_extraction.log',OUT/'extraction.log')
 files=sorted(p for p in OUT.rglob('*') if p.is_file());names=[str(p.relative_to(ROOT)) for p in files]
 r=subprocess.run(['git','check-ignore','--stdin'],input='\n'.join(names)+'\n',capture_output=True,text=True,cwd=ROOT);assert r.returncode in (0,1);ignored=set(r.stdout.splitlines())
 inventory=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'git':'excluded; preserved locally' if str(p.relative_to(ROOT)) in ignored else 'included in milestone'} for p in files]
 dump(OUT/'git_artifact_inventory.json',{'files':inventory,'note':'Inventory, seal receipt and manifest are also committed; excluded arrays/profiles/logs remain local and fully checksummed'})
 dump(OUT/'seal_receipt.json',{'P0_T6_COMPLETE':'YES','P0_T6_RESULT':'PASS','BEST_FEATURE':'profile_mean','READY_FOR_NEXT_STEP':'YES',
  'scope':'narrow CT repeatability PASS; ready for bounded deterministic feature study, not ink deployment','prior_T5_and_earlier_preservation':'PASS','feature_array_seal':'PASS','independent_replay':'PASS','all_metrics_audit':'PASS',
  'checksum':'script requires immediate full manifest validation; no tracing or geometry changes'})
 files=sorted([p for p in OUT.rglob('*') if p.is_file() and p!=manifest]+list((ROOT/'scripts').glob('p0_t6_*.py'))+[ROOT/'docs/p0_t6_local_ct_features.md'])
 manifest.write_text(''.join(sha(p)+'  '+str(p.relative_to(ROOT))+'\n' for p in files))
 check=subprocess.run(['sha256sum','-c',str(manifest)],capture_output=True,text=True,cwd=ROOT);assert check.returncode==0,check.stdout[-3000:]+check.stderr
 print('PASS: all',len(files),'manifest entries; frozen features; prior T5 seal and preservation; all required outputs')
 print('Git included artifact bytes',sum(v['bytes'] for v in inventory if v['git']=='included in milestone'),'local excluded bytes',sum(v['bytes'] for v in inventory if v['git']!='included in milestone'))
if __name__=='__main__':main()
