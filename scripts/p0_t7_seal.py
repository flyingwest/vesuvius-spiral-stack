"""Seal the bounded T7 search and verify all included evidence."""
import json,subprocess,shutil,sys,platform
import numpy,scipy,matplotlib,PIL
from PIL import Image
from p0_t7_family import ROOT,OUT,T6,sha,dump
from p0_t7_select import enrich

def main():
 manifest=OUT/'SHA256SUMS';assert not manifest.exists()
 assert json.loads((OUT/'verification.json').read_text())['status']=='PASS'
 # Reproduce joint/regional summaries without changing source arrays.
 for phase,target in [('scalars','feature_metrics.json'),('combinations','combination_metrics.json')]:
  rows,regions,ranked=enrich(OUT/phase)
  if phase=='combinations':rows=[r for r in rows if r['feature'].startswith('combo_')]
  assert rows==json.loads((OUT/target).read_text())
  for name,h in json.loads((OUT/phase/'FEATURES_FINALIZED.json').read_text()).items():assert sha(OUT/phase/name)==h,name
 for name,h in json.loads((OUT/'preservation.json').read_text()).items():assert sha(ROOT/name)==h,name
 for line in (T6/'SHA256SUMS').read_text().splitlines():
  h,name=line.split('  ',1);assert sha(ROOT/name)==h,name
 required=['baseline_center_stack.png','baseline_profile_mean_stack.png','best_individual_feature_stack.png','best_combined_feature_stack.png','feature_comparison.png','coherence_comparison.png','per_z_comparison.png','feature_metrics.csv','feature_metrics.json','combination_metrics.csv','verification.json','mask.npy','provenance.npy']
 assert all((OUT/name).is_file() for name in required)
 for f in OUT.rglob('*.png'):
  with Image.open(f) as im:im.verify()
 decision=json.loads((OUT/'decision.json').read_text());chosen=decision['BEST_INDIVIDUAL_FEATURE'];combo=decision['best_tested_combination']
 for root,phase,name in [('baseline_center_stack.png','scalars','center_intensity'),('baseline_profile_mean_stack.png','scalars','mean_wide'),('best_individual_feature_stack.png','scalars',chosen),('best_combined_feature_stack.png','combinations',combo)]:
  assert sha(OUT/root)==sha(OUT/phase/f'{name}_stack_contrast.png')
 shutil.copyfile('/home/wpotech/.hermes/cache/scratch/p0_t7_extraction.log',OUT/'extraction.log')
 dump(OUT/'environment.json',{'python':sys.version,'platform':platform.platform(),'numpy':numpy.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'PIL':PIL.__version__})
 files=sorted(p for p in OUT.rglob('*') if p.is_file());names=[str(p.relative_to(ROOT)) for p in files]
 r=subprocess.run(['git','check-ignore','--stdin'],input='\n'.join(names)+'\n',capture_output=True,text=True,cwd=ROOT);assert r.returncode in (0,1);ignored=set(r.stdout.splitlines())
 inventory=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'git':'excluded, preserved locally' if str(p.relative_to(ROOT)) in ignored else 'included'} for p in files]
 dump(OUT/'git_artifact_inventory.json',{'files':inventory,'note':'Inventory, receipt and manifest are additionally included; all excluded numeric arrays/logs remain locally sealed'})
 dump(OUT/'seal_receipt.json',{'P0_T7_COMPLETE':'YES','P0_T7_RESULT':'C','BEST_INDIVIDUAL_FEATURE':chosen,'BEST_COMBINATION':'NONE','best_tested_combination':combo,'READY_FOR_LEARNED_MODEL':'NO',
  'scalar_count':28,'combination_count':3,'prior_preservation':'PASS','T6_seal':'PASS','phase_array_seals':'PASS','independent_replay':'PASS','joint_ranking_and_regions_replay':'PASS','required_visuals':'PASS'})
 files=sorted([p for p in OUT.rglob('*') if p.is_file() and p!=manifest]+list((ROOT/'scripts').glob('p0_t7_*.py'))+[ROOT/'docs/p0_t7_deterministic_feature_search.md'])
 manifest.write_text(''.join(sha(p)+'  '+str(p.relative_to(ROOT))+'\n' for p in files))
 r=subprocess.run(['sha256sum','-c',str(manifest)],capture_output=True,text=True,cwd=ROOT);assert r.returncode==0,r.stdout[-2000:]+r.stderr
 print('PASS: all',len(files),'checksums;28 scalar+3 combination seals; preserved prior evidence; joint rankings/regions and visuals verified')
 print('Git payload',sum(x['bytes'] for x in inventory if x['git']=='included'),'local excluded bytes',sum(x['bytes'] for x in inventory if x['git']!='included'))
if __name__=='__main__':main()
