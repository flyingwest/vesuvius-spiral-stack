"""Post-search summaries and comparison visuals; no further candidate generation."""
import json,csv,shutil
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_t7_family import OUT,FEATURES,DS,dump,table
from p0_t7_select import enrich,REGIONS

def main():
 combo,regions,ranked=enrich(OUT/'combinations');actual=[r for r in combo if r['feature'].startswith('combo_')];table(OUT/'combination_metrics.csv',actual);dump(OUT/'combination_metrics.json',actual);table(OUT/'combination_region_metrics.csv',regions)
 selection=json.loads((OUT/'combination_plan.json').read_text());individual=selection['best_individual'];bestcombo=next(r['feature'] for r in ranked if r['feature'].startswith('combo_'))
 scalars=json.loads((OUT/'feature_metrics.json').read_text());by={r['feature']:r for r in scalars};bc=next(r for r in actual if r['feature']==bestcombo)
 picked=[('scalars','center_intensity','center baseline'),('scalars','mean_wide','T6 profile mean'),('scalars',individual,'joint-ranked individual'),('combinations',bestcombo,'best tested combination')]
 # Native size, one global contrast per feature; comparison figures disclose mapping.
 for phase,name,label,outname in [(a,b,c,d) for (a,b,c),d in zip(picked,['baseline_center_stack.png','baseline_profile_mean_stack.png','best_individual_feature_stack.png','best_combined_feature_stack.png'])]:
  shutil.copyfile(OUT/phase/f'{name}_stack_contrast.png',OUT/outname)
 measures=[by['center_intensity'],by['mean_wide'],by[individual],bc];labels=[r['feature'] for r in measures]
 fig,axs=plt.subplots(2,2,figsize=(14,9));x=np.arange(4)
 for ax,key,title in [(axs[0,0],'mean_r','Mean adjacent ZERO-SHIFT r'),(axs[0,1],'local_mean_r','Mean local ZERO-SHIFT r'),(axs[1,0],'local_peak_above_control95','Local peaks exceeding matched controls /1001'),(axs[1,1],'regions_both_improved','Regions improving BOTH metrics vs T6 mean /5')]:
  ax.bar(x,[r[key] for r in measures]);ax.set(xticks=x,xticklabels=labels,title=title);ax.tick_params(axis='x',rotation=15)
 fig.suptitle('Fixed joint selection is not a claim of superiority; no feature beats T6 across criteria');fig.tight_layout();fig.savefig(OUT/'feature_comparison.png',dpi=145);plt.close(fig)
 fig,axs=plt.subplots(2,1,figsize=(15,9),sharex=True)
 for phase,name,label in picked:
  p=json.loads((OUT/phase/f'{name}_evaluation.json').read_text())['pairs'];axs[0].plot([r['z_a']+.5 for r in p],[r['r'] for r in p],label=name,lw=1);axs[1].plot([r['z_a']+.5 for r in p],[r['local_mean_r'] for r in p],label=name,lw=1)
 for ax in axs:ax.axhline(0,color='gray',lw=.5);ax.legend();ax.set_ylabel('Pearson r')
 axs[0].set_title('Every adjacent pair: authoritative zero-shift');axs[1].set(title='Within-pair mean local zero-shift correlation',xlabel='Adjacent-pair midpoint Z');fig.tight_layout();fig.savefig(OUT/'per_z_comparison.png',dpi=145);plt.close(fig)
 fig,axs=plt.subplots(4,2,figsize=(17,12));mask=np.load(OUT/'mask.npy');cm=plt.cm.gray.copy();cm.set_bad('#be40b5')
 for (phase,name,label),axes in zip(picked,axs):
  a=np.load(OUT/phase/f'{name}_stack.npy');lo,hi=np.percentile(a[~mask],[1,99]);axes[0].imshow(np.ma.masked_array(a,mask),cmap=cm,vmin=lo,vmax=hi,aspect='auto',interpolation='nearest',extent=[0,(a.shape[1]-1)*DS/1000,177.5,127.5]);axes[0].set(title=f'{label}: {name}; global contrast',ylabel='Z')
  rows=list(csv.DictReader((OUT/phase/f'{name}_local_windows.csv').open()));nwin=max(int(r['start_sample'])//128 for r in rows)+1;heat=np.full((49,nwin),np.nan)
  for r in rows:heat[int(r['z_a'])-128,int(r['start_sample'])//128]=float(r['r'])
  axes[1].imshow(np.ma.masked_invalid(heat),cmap='coolwarm',vmin=-1,vmax=1,aspect='auto',interpolation='nearest',extent=[0,nwin*128*DS/1000,176.5,127.5]);axes[1].set(title='Local ZERO-SHIFT r; red positive / blue negative',ylabel='Pair start Z')
 for ax in axs[-1]:ax.set_xlabel('s mm')
 fig.tight_layout();fig.savefig(OUT/'coherence_comparison.png',dpi=145);plt.close(fig)
 # Compare combination to selected individual using exactly the same paired blocks.
 ip=json.loads((OUT/f'scalars/{individual}_evaluation.json').read_text())['pairs'];cp=json.loads((OUT/f'combinations/{bestcombo}_evaluation.json').read_text())['pairs'];boot=np.load(OUT/'combinations/bootstrap_pair_indices.npy')
 dr=np.array([a['r']-b['r'] for a,b in zip(cp,ip)]);dl=np.array([a['local_mean_r']-b['local_mean_r'] for a,b in zip(cp,ip)])
 region_summary=[]
 sr=list(csv.DictReader((OUT/'scalar_region_metrics.csv').open()));cr=list(csv.DictReader((OUT/'combination_region_metrics.csv').open()))
 for phase,name,label in picked:
  for r in (sr if phase=='scalars' else cr):
   if r['feature']==name:region_summary.append(r)
 table(OUT/'selected_region_comparison.csv',region_summary)
 result={'P0_T7_COMPLETE':'YES','P0_T7_RESULT':'C','result_name':'SIGNAL_INSUFFICIENT','BEST_INDIVIDUAL_FEATURE':individual,'BEST_COMBINATION':'NONE','best_tested_combination':bestcombo,'READY_FOR_LEARNED_MODEL':'NO',
  'reason':'No scalar or tested combination clearly improves both adjacent and local repeatability across Z beyond T6 mean; local-control gains trade off against raw correlations.',
  'selected_individual':by[individual],'best_tested_combination_metrics':bc,
  'combination_vs_individual':{'mean_adjacent_gain':float(dr.mean()),'paired95':np.percentile(dr[boot].mean(1),[2.5,97.5]).tolist(),'mean_pair_local_gain':float(dl.mean()),'paired_local95':np.percentile(dl[boot].mean(1),[2.5,97.5]).tolist()},
  'broad_guide_satisfying_scalars':[r['feature'] for r in scalars if r['broad_gain_guide']], 'broad_guide_satisfying_combinations':[r['feature'] for r in actual if r['broad_gain_guide']],
  'warning':'C refers only to gain beyond T6 in this bounded family, not rejection of all CT signals or geometry. Learned ink readiness requires target-specific data and validation, not merely exhausted summaries.'}
 assert not result['broad_guide_satisfying_scalars'] and not result['broad_guide_satisfying_combinations']
 dump(OUT/'decision.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
