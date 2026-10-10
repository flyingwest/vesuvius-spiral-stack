"""Frozen ten-feature evaluation with spatial confirmation and matched controls."""
from p0_t8_acquire import *
from p0_t7_family import compute
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
import csv
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def table(p,rows):
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def metric_setup(y,score,values=None):
 order=np.argsort(score,kind='stable');s=score[order];yy=y[order].astype(float);v=(score if values is None else values)[order]
 starts=np.r_[0,np.flatnonzero(np.diff(s))+1]
 def calc(w=None):
  w=np.ones(len(y)) if w is None else w[order]
  p=np.add.reduceat(w*yy,starts);n=np.add.reduceat(w*(1-yy),starts);P=p.sum();N=n.sum()
  if P==0 or N==0:return np.full(3,np.nan)
  auc=np.sum(p*(np.cumsum(n)-.5*n))/(P*N)
  prec=np.cumsum(p[::-1])/np.maximum(np.cumsum((p+n)[::-1]),1e-100)
  ap=np.sum(p[::-1]*prec)/P
  mp=np.sum(w*yy*v)/P;mn=np.sum(w*(1-yy)*v)/N
  vp=np.sum(w*yy*(v-mp)**2)/P;vn=np.sum(w*(1-yy)*(v-mn)**2)/N
  effect=(mp-mn)/max(np.sqrt((vp+vn)/2),1e-100)
  return np.array([auc,ap,effect])
 return calc

def summarize(y,v,blocks,sign,seed=2808):
 f=metric_setup(y,sign*v,v);point=f();ids,inv=np.unique(blocks,return_inverse=True);rng=np.random.default_rng(seed)
 boot=np.array([f(np.bincount(rng.integers(len(ids),size=len(ids)),minlength=len(ids))[inv]) for _ in range(2000)])
 assert np.isfinite(boot).all()
 r={'n':len(y),'positive':int(y.sum()),'negative':int((~y).sum()),'blocks':len(ids),'prevalence':float(y.mean()),'sign_from_discovery':sign}
 for j,name in enumerate(['auc','average_precision','effect_d']):
  r[name]=float(point[j]);r[name+'_ci95']=np.quantile(boot[:,j],[.025,.975]).tolist();r[name+'_ci995']=np.quantile(boot[:,j],[.0025,.9975]).tolist()
 r['positive_quantiles']=np.quantile(v[y],[.05,.25,.5,.75,.95]).tolist();r['negative_quantiles']=np.quantile(v[~y],[.05,.25,.5,.75,.95]).tolist()
 return r,boot

def match(ids,y,xy,blocks,normal,baseline,strength,controlled):
 pos=ids[y[ids]];neg=ids[~y[ids]];tree=cKDTree(xy[neg]);used=set();pairs=[]
 sb=baseline[ids].std();sg=strength[ids].std()
 for p in pos:
  c=neg[tree.query_ball_point(xy[p],128)];c=c[blocks[c]==blocks[p]]
  if controlled:
   c=c[(normal[c]@normal[p]>=.98)&(abs(baseline[c]-baseline[p])<=.25*sb)&(abs(strength[c]-strength[p])<=.5*sg)]
  if not len(c):continue
  dist=((xy[c]-xy[p])**2).sum(1);c=c[np.lexsort((c,dist))]
  for n in c:
   if int(n) not in used:pairs.append((p,n));used.add(int(n));break
 return np.array(pairs,dtype=int).reshape(-1,2)

def main():
 plan=json.loads((OUT/'plan.json').read_text());assert sha(OUT/'plan.json')==json.loads((OUT/'PLAN_FROZEN.json').read_text())['plan_sha256']
 a=np.load(DATA/'ct_profiles_dyx.npy');reg=np.load(OUT/'registration_arrays.npz');geom=np.load(DATA/'grid_geometry.npy')
 yy,xx=np.mgrid[8:2048:16,8:2048:16];eligible=reg['usable'][yy,xx]
 xy=np.column_stack([xx[eligible],yy[eligible]]);profiles=a[:,yy[eligible],xx[eligible]].T.astype(float)
 labels=reg['labels'][yy[eligible],xx[eligible]];blocks=(xy[:,1]//256)*8+xy[:,0]//256
 normals=geom[:,:,3:][eligible];normals=normals/np.linalg.norm(normals,axis=1)[:,None]
 baseline=gaussian_filter(a[3].astype(float),32)[yy[eligible],xx[eligible]]
 sm=gaussian_filter(a[3].astype(float),2);gy,gx=np.gradient(sm);strength=np.hypot(gx,gy)[yy[eligible],xx[eligible]]
 features={k:v for k,v in compute(profiles).items() if k in plan['features']}
 amendment=json.loads((OUT/'analysis_plan_amendment.json').read_text())
 assert amendment['before_ct_discrimination_metrics']
 discovery=np.flatnonzero(xy[:,0]<1024);test=np.flatnonzero(xy[:,0]>=1024)
 assert labels[discovery].sum()>100 and labels[test].sum()>100
 signs={k:1 if metric_setup(labels[discovery],v[discovery])()[0]>=.5 else -1 for k,v in features.items()}
 dump(OUT/'DIRECTIONS_FROZEN.json',{'signs':signs,'discovery_only':True,'plan_sha256':sha(OUT/'plan.json')})
 np.savez_compressed(OUT/'sampled_features.npz',xy_local=xy,xy_global=xy+[1600,2784],profiles=profiles,labels=labels,blocks=blocks,normals=normals,baseline=baseline,strength=strength,discovery=discovery,confirmation=test,**features)
 matches={name:match(test,labels,xy,blocks,normals,baseline,strength,ctrl) for name,ctrl in [('nearby',False),('nuisance',True)]}
 np.savez_compressed(OUT/'matched_pairs.npz',**matches)
 dump(OUT/'SAMPLES_FROZEN.json',{'sample_hash':sha(OUT/'sampled_features.npz'),'pair_hash':sha(OUT/'matched_pairs.npz'),'counts':{'discovery':len(discovery),'confirmation':len(test),**{k:len(p) for k,p in matches.items()}}})
 rows=[];details={};boots={};regional=[];balance=[]
 for m,pairs in matches.items():
  assert len(pairs)>100,(m,len(pairs))
  p,n=pairs.T
  balance.append({'matching':m,'pairs':len(p),'median_distance_px':float(np.median(np.linalg.norm(xy[p]-xy[n],axis=1))),'min_normal_dot':float(np.min(np.sum(normals[p]*normals[n],axis=1))),
   'baseline_signed_difference_sd':float((baseline[p]-baseline[n]).mean()/baseline[test].std()),'baseline_absolute_difference_sd':float(abs(baseline[p]-baseline[n]).mean()/baseline[test].std()),'strength_signed_difference_sd':float((strength[p]-strength[n]).mean()/strength[test].std()),'strength_absolute_difference_sd':float(abs(strength[p]-strength[n]).mean()/strength[test].std())})
 for name,v in features.items():
  sign=signs[name];result={}
  for cohort,ids in [('discovery',discovery),('confirmation',test)]+[(k,p.ravel()) for k,p in matches.items()]:
   r,boot=summarize(labels[ids],v[ids],blocks[ids],sign);result[cohort]=r;boots[name+'_'+cohort]=boot
   flat={'feature':name,'cohort':cohort,**{k:val for k,val in r.items() if not isinstance(val,list)}}
   for k in ['auc','average_precision','effect_d']:
    for ci in ['ci95','ci995']:
     flat[k+'_'+ci+'_low'],flat[k+'_'+ci+'_high']=r[k+'_'+ci]
   if cohort in matches:
    p,n=matches[cohort].T;win=(sign*v[p]>sign*v[n]).astype(float)+.5*(v[p]==v[n]);b=blocks[p];idsb,inv=np.unique(b,return_inverse=True)
    sums=np.bincount(inv,weights=win);nums=np.bincount(inv);rng=np.random.default_rng(2808)
    cb=[]
    for _ in range(2000):
     ix=rng.integers(len(idsb),size=len(idsb));cb.append(sums[ix].sum()/nums[ix].sum())
    r['paired_concordance']=float(win.mean());r['paired_concordance_ci95']=np.quantile(cb,[.025,.975]).tolist();r['paired_concordance_ci995']=np.quantile(cb,[.0025,.9975]).tolist();boots[name+'_'+cohort+'_paired']=np.array(cb)
    flat['paired_concordance']=r['paired_concordance'];flat['paired_ci95_low'],flat['paired_ci95_high']=r['paired_concordance_ci95']
   else:flat.update(paired_concordance='',paired_ci95_low='',paired_ci95_high='')
   rows.append(flat)
  regional_pass=0
  # Four fixed512x512 corner blocks in the right-half confirmation area.
  for region in range(4):
   rx=1024+(region%2)*512;ry=(region//2)*1536
   ids=test[(xy[test,0]>=rx)&(xy[test,0]<rx+512)&(xy[test,1]>=ry)&(xy[test,1]<ry+512)]
   vals=metric_setup(labels[ids],sign*v[ids])();regional_pass+=int(vals[0]>.5)
   regional.append({'feature':name,'region':region,'x0':rx,'y0':ry,'n':len(ids),'positive':int(labels[ids].sum()),'auc':float(vals[0])})
  result['regions_above_chance']=regional_pass;details[name]=result
  print(name,'confirm AUC',round(result['confirmation']['auc'],4),'controlled paired',round(result['nuisance']['paired_concordance'],4),flush=True)
 np.savez_compressed(OUT/'bootstrap_metrics.npz',**boots)
 dump(OUT/'feature_metrics.json',details);table(OUT/'feature_metrics.csv',rows);table(OUT/'regional_metrics.csv',regional);table(OUT/'matching_balance.csv',balance)
 qualified=[k for k,d in details.items() if d['confirmation']['auc_ci995'][0]>.5 and d['regions_above_chance']>=3 and d['nuisance']['paired_concordance_ci95'][0]>.5]
 best=max(qualified,key=lambda k:details[k]['confirmation']['auc']) if qualified else 'NONE'
 decision='C' if not qualified else 'A' if details[best]['confirmation']['auc']>=.60 else 'B'
 dump(OUT/'decision.json',{'result':decision,'best_feature':best,'qualifying_features':qualified,'learned_model_justified':'NO' if decision=='C' else 'MAYBE','meaning':'single-fragment frozen-feature association; not an ink detector, not proof of chemistry or transfer to downsampled T5'})
 fig,axs=plt.subplots(5,2,figsize=(13,17))
 for ax,(name,v) in zip(axs.ravel(),features.items()):
  lo,hi=np.quantile(v[test],[.005,.995]);bins=np.linspace(lo,hi,50)
  for value,color,label in [(False,'steelblue','nonink'),(True,'tomato','ink')]:ax.hist(v[test][labels[test]==value],bins=bins,density=True,histtype='step',color=color,label=label)
  ax.set_title(name+f" | confirmed AUC {details[name]['confirmation']['auc']:.3f}");ax.legend();ax.set_xlabel('Native CT-derived units; tails outside0.5–99.5% not displayed')
 fig.tight_layout();fig.savefig(OUT/'feature_distributions.png',dpi=150);plt.close(fig)
 fig,ax=plt.subplots(figsize=(12,6));names=list(features);x=np.arange(len(names))
 for j,(cohort,color) in enumerate([('confirmation','black'),('nearby','tab:blue'),('nuisance','tab:orange')]):
  val=np.array([details[k][cohort]['auc'] for k in names]);ci=np.array([details[k][cohort]['auc_ci95'] for k in names]);ax.errorbar(x+(j-1)*.2,val,yerr=np.abs(ci.T-val),fmt='o',color=color,label=cohort)
 ax.axhline(.5,color='gray');ax.set_xticks(x,names,rotation=40,ha='right');ax.set_ylabel('ROC AUC; discovery-selected sign; block95%CI');ax.legend();fig.tight_layout();fig.savefig(OUT/'auc_comparison.png',dpi=160);plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(12,5))
 for name,v in features.items():
  score=signs[name]*v[test];order=np.argsort(-score,kind='stable');y=labels[test][order];ends=np.r_[np.flatnonzero(np.diff(score[order])),len(y)-1];tp=np.cumsum(y)[ends];fp=np.cumsum(~y)[ends]
  axs[0].plot(np.r_[0,fp/fp[-1]],np.r_[0,tp/tp[-1]],label=name);axs[1].plot(tp/tp[-1],tp/(tp+fp))
 axs[0].plot([0,1],[0,1],'k--');axs[0].set(xlabel='False positive rate',ylabel='True positive rate',title='Confirmation ROC');axs[0].legend(fontsize=7)
 axs[1].axhline(labels[test].mean(),color='gray',ls='--');axs[1].set(xlabel='Recall',ylabel='Precision',title='Confirmation PR; dashed=prevalence')
 fig.tight_layout();fig.savefig(OUT/'roc_pr_curves.png',dpi=160);plt.close(fig)
 print('Decision',decision,best,flush=True)
if __name__=='__main__':main()
