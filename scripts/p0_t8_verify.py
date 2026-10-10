"""Independent saved-byte/profile/metric checks; no path reconstruction."""
from p0_t8_acquire import *
from scipy.stats import rankdata
from PIL import Image
import csv

def auc(y,v):
 n=int(y.sum());m=len(y)-n
 return float((rankdata(v)[y].sum()-n*(n+1)/2)/(n*m))
def ap(y,v):
 # Independently accumulate tied thresholds from high score to low.
 vals=np.unique(v)[::-1];order=np.argsort(-v);ys=y[order];ss=v[order];ends=np.r_[np.flatnonzero(np.diff(ss)),len(y)-1];tp=np.cumsum(ys)[ends];counts=ends+1
 return float(np.sum(np.diff(np.r_[0,tp])*tp/counts)/y.sum())
def main():
 prior=json.loads((OUT/'preservation.json').read_text());assert all(sha(ROOT/p)==h for p,h in prior.items())
 acq=json.loads((OUT/'acquisition.json').read_text());assert all(sha(ROOT/p)==h for p,h in acq['data_hashes'].items())
 a=np.load(DATA/'ct_profiles_dyx.npy');source=np.stack([np.fromfile(DATA/f'{d}_rows.bin','<u2').reshape(2048,5249)[:,1600:3648] for d in [28,30,31,32,33,34,36]])
 assert np.array_equal(a,source)
 s=np.load(OUT/'sampled_features.npz');p=s['profiles'];xy=s['xy_local'];assert np.array_equal(p,a[:,xy[:,1],xy[:,0]].T)
 expected={'center_intensity':p[:,3],'mean_wide':np.sum(p,axis=1)/7,'median_wide':np.sort(p,axis=1)[:,3],'std_wide':np.sqrt(((p-p.mean(1)[:,None])**2).sum(1)/7),'range_wide':p.max(1)-p.min(1),'center_minus_mean_wide':p[:,3]-p.mean(1),'asymmetry_wide':abs(p[:,4:].mean(1)-p[:,:3].mean(1)),'gradient_abs_narrow':abs(p[:,4]-p[:,2])/2,'curvature_narrow':2*p[:,3]-p[:,2]-p[:,4],'ridge_proxy_wide':p[:,3]-p[:,[0,1,5,6]].mean(1)}
 assert all(np.array_equal(s[k],v) for k,v in expected.items())
 frozen=json.loads((OUT/'SAMPLES_FROZEN.json').read_text());assert sha(OUT/'sampled_features.npz')==frozen['sample_hash'];assert sha(OUT/'matched_pairs.npz')==frozen['pair_hash']
 plan=json.loads((OUT/'PLAN_FROZEN.json').read_text());assert sha(OUT/'plan.json')==plan['plan_sha256'];assert sha(ROOT/'scripts/p0_t7_family.py')==plan['feature_code_sha256']
 af=json.loads((OUT/'ANALYSIS_PLAN_FROZEN.json').read_text());assert sha(OUT/'analysis_plan_amendment.json')==af['amendment_sha256'];assert sha(ROOT/'scripts/p0_t8_evaluate.py')==af['evaluation_code_sha256']
 labels=s['labels'];orig=np.array(Image.open(DATA/'inklabels.png'))>0;assert np.array_equal(labels,orig[s['xy_global'][:,1],s['xy_global'][:,0]])
 pairs=np.load(OUT/'matched_pairs.npz');test=s['confirmation'];assert not set(test)&set(s['discovery'])
 for name,pair in pairs.items():
  pos,neg=pair.T;assert np.all(labels[pos]) and np.all(~labels[neg]);assert len(set(neg))==len(neg);assert set(pair.ravel()).issubset(set(test));assert np.all(s['blocks'][pos]==s['blocks'][neg]);assert np.all(np.linalg.norm(xy[pos]-xy[neg],axis=1)<=128)
  if name=='nuisance':
   assert np.all(np.sum(s['normals'][pos]*s['normals'][neg],axis=1)>=.98)
   assert np.all(abs(s['baseline'][pos]-s['baseline'][neg])<=.25*s['baseline'][test].std())
   assert np.all(abs(s['strength'][pos]-s['strength'][neg])<=.5*s['strength'][test].std())
 metrics=json.loads((OUT/'feature_metrics.json').read_text());directions=json.loads((OUT/'DIRECTIONS_FROZEN.json').read_text())['signs'];boot=np.load(OUT/'bootstrap_metrics.npz');count=0
 for k in expected:
  discovery=s['discovery'];assert directions[k]==(1 if auc(labels[discovery],s[k][discovery])>=.5 else -1)
  for cohort,ids in [('discovery',discovery),('confirmation',test)]+[(n,q.ravel()) for n,q in pairs.items()]:
   y=labels[ids];v=s[k][ids];r=metrics[k][cohort];sign=directions[k]
   assert abs(auc(y,sign*v)-r['auc'])<1e-12;assert abs(ap(y,sign*v)-r['average_precision'])<1e-12
   d=(v[y].mean()-v[~y].mean())/np.sqrt((v[y].var()+v[~y].var())/2);assert abs(d-r['effect_d'])<1e-12
   for j,n in enumerate(['auc','average_precision','effect_d']):
    assert np.allclose(np.quantile(boot[k+'_'+cohort][:,j],[.025,.975]),r[n+'_ci95'],rtol=0,atol=1e-12)
    assert np.allclose(np.quantile(boot[k+'_'+cohort][:,j],[.0025,.9975]),r[n+'_ci995'],rtol=0,atol=1e-12)
   # Independently replay ten deterministic bootstrap draws by literal block replication.
   blocks=s['blocks'][ids];uniq,inv=np.unique(blocks,return_inverse=True);rng=np.random.default_rng(2808)
   for rep in range(10):
    weights=np.bincount(rng.integers(len(uniq),size=len(uniq)),minlength=len(uniq))[inv];ix=np.repeat(np.arange(len(ids)),weights)
    assert abs(auc(y[ix],sign*v[ix])-boot[k+'_'+cohort][rep,0])<1e-12
    assert abs(ap(y[ix],sign*v[ix])-boot[k+'_'+cohort][rep,1])<1e-12
   count+=1
 for pth in OUT.glob('*.png'):
  with Image.open(pth) as im:im.verify()
 result={'status':'PASS','prior_files_unchanged':len(prior),'source_files_hash_verified':len(acq['data_hashes']),'ct_crop_replayed_exact':True,'all_sample_profiles_replayed_exact':True,'ten_feature_formulas_replayed_exact':True,'all_labels_xy_replayed_exact':True,'metric_cohorts_independently_replayed':count,'bootstrap_replicates_independently_replayed':count*10,'all_saved_bootstrap_confidence_intervals_verified':True,'bootstrap_scope':'10/2000 draws per feature/cohort independently regenerated; all saved draw quantiles checked','matched_constraints_verified':True,'plan_and_source_formula_hashes_verified':True,'no_T5_geometry_or_rows_modified':True,'limitations':'Registration subpixel error not independently measured. C is lack of established signal, not proof of equivalence to chance.'}
 dump(OUT/'verification.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
