"""Predeclared joint ranking and at most three fixed equal-weight combinations."""
import json,shutil,itertools
import numpy as np
from scipy.stats import rankdata
from p0_t7_family import ROOT,OUT,FEATURES,SCALES,sha,dump,table

REGIONS=[(128,137),(138,147),(148,157),(158,167),(168,176)]
def enrich(folder,reference=None):
 metrics=json.loads((folder/'feature_metrics.json').read_text())
 basefolder=OUT/'scalars'
 cp=json.loads((basefolder/'center_intensity_evaluation.json').read_text())['pairs'];mp=json.loads((basefolder/'mean_wide_evaluation.json').read_text())['pairs']
 center=np.array([p['r'] for p in cp]);mean=np.array([p['r'] for p in mp]);ml=np.array([p['local_mean_r'] for p in mp]);boot=np.load(folder/'bootstrap_pair_indices.npy')
 results=[];regionrows=[]
 for m in metrics:
  p=json.loads((folder/f'{m["feature"]}_evaluation.json').read_text())['pairs'];r=np.array([x['r'] for x in p]);l=np.array([x['local_mean_r'] for x in p]);regions=0
  for lo,hi in REGIONS:
   sel=np.array([lo<=x['z_a']<=hi for x in p]);dr=float((r-mean)[sel].mean());dl=float((l-ml)[sel].mean());improve=dr>1e-12 and dl>1e-12;regions+=int(improve)
   regionrows.append({'feature':m['feature'],'z_start':lo,'z_end_pair_start':hi,'pairs':int(sel.sum()),'mean_r':float(r[sel].mean()),'local_pair_mean_r':float(l[sel].mean()),'delta_r_vs_profile_mean':dr,'delta_local_vs_profile_mean':dl,'both_improved':improve})
  ci=np.percentile((r-mean)[boot].mean(1),[2.5,97.5]);lci=np.percentile((l-ml)[boot].mean(1),[2.5,97.5])
  rec={**m,'delta_vs_T6_mean':float((r-mean).mean()),'pairs_better_than_T6_mean':int((r-mean>1e-12).sum()),'delta_local_pair_vs_T6_mean':float((l-ml).mean()),
       'regions_both_improved':regions,'T6_mean_gain_ci_low':float(ci[0]),'T6_mean_gain_ci_high':float(ci[1]),'T6_local_gain_ci_low':float(lci[0]),'T6_local_gain_ci_high':float(lci[1])}
  rec['broad_gain_guide']=bool(rec['delta_vs_T6_mean']>=.02 and rec['local_mean_r']-next(x['local_mean_r'] for x in json.loads((basefolder/'feature_metrics.json').read_text()) if x['feature']=='mean_wide')>=.01 and rec['pairs_better_than_T6_mean']>=30 and regions>=4 and rec['local_peak_above_control95']>=next(x['local_peak_above_control95'] for x in json.loads((basefolder/'feature_metrics.json').read_text()) if x['feature']=='mean_wide'))
  results.append(rec)
 candidates=[m for m in results if m['feature']!='center_intensity']
 for m in results:m['joint_rank_score']=None
 keys=['mean_r','local_mean_r','local_peak_above_control95','regions_both_improved']
 ranks=np.array([rankdata([m[k] for m in candidates],method='average')/len(candidates) for k in keys]).mean(0)
 for m,score in zip(candidates,ranks):m['joint_rank_score']=float(score)
 ranked=sorted(candidates,key=lambda m:(-m['joint_rank_score'],m['feature']))
 return results,regionrows,ranked

def main():
 folder=OUT/'scalars';results,regions,ranked=enrich(folder)
 dump(OUT/'feature_metrics.json',results);table(OUT/'feature_metrics.csv',results);table(OUT/'scalar_region_metrics.csv',regions)
 selected=[];dedup=[];mask=np.load(OUT/'mask.npy')
 for m in ranked:
  name=m['feature'];a=np.load(folder/f'{name}_stack.npy')[~mask]/SCALES[name]
  if any(abs(float(np.corrcoef(a,b)[0,1]))>1-1e-12 for _,b in selected):dedup.append(name);continue
  selected.append((name,a))
  if len(selected)==3:break
 assert len(selected)==3
 dest=OUT/'combinations';dest.mkdir(exist_ok=False)
 specs={'center_intensity':{'definition':FEATURES['center_intensity'][0],'display_range':[0,255]},'mean_wide':{'definition':FEATURES['mean_wide'][0],'display_range':[0,255]}}
 combos=[]
 for i,j in itertools.combinations(range(3),2):
  na=selected[i][0];nb=selected[j][0];name=f'combo_{i+1}_{j+1}';aa=np.load(folder/f'{na}_stack.npy');bb=np.load(folder/f'{nb}_stack.npy');a=.5*(aa/SCALES[na]+bb/SCALES[nb])
  specs[name]={'definition':f'0.5*({na}/{SCALES[na]} + {nb}/{SCALES[nb]})','display_range':[-1,1]};combos.append({'name':name,'a':na,'b':nb,'scale_a':SCALES[na],'scale_b':SCALES[nb],'weights':[.5,.5]})
  np.save(dest/f'{name}_stack.npy',a)
 for name in ['center_intensity','mean_wide']:shutil.copyfile(folder/f'{name}_stack.npy',dest/f'{name}_stack.npy')
 for name in ['mask','provenance','coverage','s_um','z_indices','z_relative_um']:shutil.copyfile(OUT/f'{name}.npy',dest/f'{name}.npy')
 dump(OUT/'combination_plan.json',{'policy':'three pairwise fixed equal-weight means of joint-ranked top3 non-affine-equivalent signals; no further search','ranked_scalars':[m['feature'] for m in ranked],
  'selected':[name for name,_ in selected],'affine_duplicate_skips':dedup,'best_individual':ranked[0]['feature'],'combinations':combos,'specifications':specs})
 dump(dest/'FEATURES_FINALIZED.json',{str(p.relative_to(dest)):sha(p) for p in sorted(dest.glob('*.npy'))})
 print(json.dumps({'selected':[name for name,_ in selected],'best':ranked[0],'combination_count':3},indent=2))
if __name__=='__main__':main()
