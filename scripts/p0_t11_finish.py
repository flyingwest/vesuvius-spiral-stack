"""Seal-stage diagnostics for the frozen early-stop branch; no model inference."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np
from PIL import Image,ImageDraw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_t11_core import cube
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t11';RUN=OUT/'complete_support_attempt'
def read(p):return json.loads(p.read_text())
def write(p,a):p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    m=read(RUN/'frame_metrics.json');assert m['selected_sigma'] is None
    assert not any(a['quality_pass'] for a in m['methods'])
    result=m['decision_if_stop'];records=read(RUN/'frames/records.json');samples=read(OUT/'samples.json')
    preserve=read(OUT/'preservation.json');blocked=read(OUT/'blocked_attempt_preservation.json');frozen=read(OUT/'PRE_EVIDENCE_FREEZE.json')['source_sha256'];support=read(OUT/'SUPPORT_FREEZE.json')['sha256']
    for manifest in [preserve,blocked,frozen,support]:
        for name,h in manifest.items():assert sha(ROOT/name)==h,name
    original_ct=read(ROOT/'outputs/p0_t9/acquisition.json')['chunks']
    for item in original_ct:
        p=ROOT/item['filename'];assert p.stat().st_size==item['bytes'] and sha(p)==item['sha256']
    receipt=read(OUT/'supplemental_ct_manifest.json')
    for item in receipt['chunks']:assert sha(ROOT/item['filename'])==item['sha256'] and (ROOT/item['filename']).stat().st_size==item['bytes']
    replay=read(RUN/'frame_replay.json');assert replay['status']=='PASS' and replay['sample_methods_replayed']==144
    assert read(RUN/'frames/records.json')==read(RUN/'frame_replay/records.json')
    oldcheck=read(OUT/'old_frame_verification.json');assert oldcheck['status']=='PASS'
    for r in records:
        i=r['sample'];s=int(r['sigma']);p=samples[i];assert r['P_xyz']==p['P_xyz'] and r['s_index']==p['s_index'] and r['z']==p['z']
        assert np.allclose(np.cross(r['t_s'],r['t_z']),r['n'],atol=1e-6)
        assert np.dot(r['n'],r['normal_xyz'])>=0
        with np.load(RUN/f'frames/patch_{i}_{s}.npz') as a,np.load(RUN/f'frame_replay/patch_{i}_{s}.npz') as b:
            assert np.array_equal(a['old'],b['old']) and np.array_equal(a['new'],b['new'])
            assert hashlib.sha256(a['new'].tobytes()).hexdigest()==r['patch_sha256']
            prior=OUT/f'frames/patch_{i}_{s}.npz'
            if prior.exists():
                with np.load(prior) as c:assert np.array_equal(a['old'],c['old']) and np.array_equal(a['new'],c['new'])
    for name in ['frame_metrics.json','patch_quality_metrics.csv']:shutil.copyfile(RUN/name,OUT/name)
    shutil.copyfile(RUN/'frame_replay.json',OUT/'frame_replay.json')
    diagnostics=[]
    for method in m['methods']:
        rr=[r for r in records if r['sigma']==method['sigma']]
        diagnostics.append({'sigma':method['sigma'],'accepted_fraction':method['accepted']/72,'rejection_reason_counts_overlapping':{'confidence':sum(r['confidence']<.25 for r in rr),'perturbation':sum(r['perturbation_deg']>20 for r in rr),'tangent_incompatibility':sum(r['tangent_incompatibility_deg']>30 for r in rr),'tilt':sum(r['tilt_deg']>75 for r in rr)},'paired_median_relative_reduction':{k:(v['old']-v['new'])/v['old'] if v['old'] else None for k,v in method['quality_medians'].items()}})
    t10=read(ROOT/'outputs/p0_t10/control_vs_t9_geometry.json')['measured_official_frame'];base=read(ROOT/'outputs/p0_t9/metrics.json')
    decision={'P0_T11_COMPLETE':'YES','P0_T11_RESULT':result,'MISSING_CT_COVERAGE':'RESOLVED','ADDITIONAL_CHUNKS':receipt['additional_chunks'],'ADDITIONAL_BYTES':receipt['additional_bytes'],'PATCH_COVERAGE':'144/144','FRAME_ESTIMATION':m['frame_estimation'],'PATCH_GEOMETRY_IMPROVED':'WEAK','SMOKE_INFERENCE':'NOT_RUN','FULL_INFERENCE':'NOT_RUN','INK_STRUCTURE_VISIBLE':'NO','READY_FOR_NEXT_STEP':'NO','ink_visibility_scope':'Not assessed: no T11 probabilities exist; NO means not established, not a negative ink observation.','local_llm_dependency':'NONE','stop_reason':'No estimator passes unchanged geometry gate; no smoke or full inference authorized by original progression.'}
    metrics={'decision':decision,'frame_methods':m['methods'],'diagnostics':diagnostics,'trusted_control_frame':t10,'baseline_metrics':{k:base[k] for k in ['T5_center','T6_profile_mean','T9_ink']},'T11_probability_metrics':None,'T11_valid_inference_fraction_of_frozen_centers':0.,'expected_full_stack_shape':base['shape'],'actual_probability_stack_shape':None,'full_target_coverage_verified':False,'qualitative':'Repeated bands and empty regions remain at fixed samples; some orientation changes improve CT proxies but not the frozen geometry gate.'}
    write(OUT/'metrics.json',metrics);write(OUT/'decision.json',decision)
    write(OUT/'smoke_metrics.json',{'status':'NOT_RUN','reason':'Original geometry gate closed','probability_metrics':None,'not_a_CUDA_or_model_failure':True})
    # These visibly labeled status placeholders contain no probability raster.
    for name in ['ink_probability_stack_raw.png','ink_probability_stack_display.png']:
        image=Image.new('RGB',(1000,160),'white');draw=ImageDraw.Draw(image)
        draw.text((30,30),'NOT RUN: P0-T11 geometry gate closed. No model inference or probability stack.',fill='black')
        draw.text((30,65),'This is a status placeholder, NOT scientific probability data.',fill='black');image.save(OUT/name)
    fig,axs=plt.subplots(1,3,figsize=(14,4))
    for sigma in [1,2]:
        rr=[r for r in records if r['sigma']==sigma]
        for ax,key,values in [(axs[0],'|normal Z|',[abs(r['normal_z']) for r in rr]),(axs[1],'longitudinal tilt degrees',[r['tilt_deg'] for r in rr]),(axs[2],'orientation confidence',[r['confidence'] for r in rr])]:ax.hist(values,bins=15,alpha=.5,label=f'sigma{sigma}');ax.set_title(key);ax.legend()
    axs[0].axvline(t10['surface_normal_abs_z_quantiles'][1],color='black',ls='--',label='T10 median');axs[1].axvline(t10['V_deviation_from_scanZ_degrees_quantiles'][1],color='black',ls='--',label='T10 median')
    for ax in axs[:2]:ax.legend()
    fig.tight_layout();fig.savefig(OUT/'frame_distributions.png');plt.close(fig)
    fig,axs=plt.subplots(3,6,figsize=(18,9));overlay=OUT/'frame_overlays';overlay.mkdir(exist_ok=True)
    for row,i in enumerate([0,32,64]):
        with np.load(RUN/f'frames/patch_{i}_1.npz') as a,np.load(RUN/f'frames/patch_{i}_2.npz') as b:
            for col,(title,plane) in enumerate([('old model X-Y',a['old'][8]),('sigma1 model X-Y',a['new'][8]),('sigma2 model X-Y',b['new'][8]),('old depth-X',a['old'][:,64]),('sigma1 depth-X',a['new'][:,64]),('sigma2 depth-X',b['new'][:,64])]):axs[row,col].imshow(plane,cmap='gray',vmin=0,vmax=255,aspect='auto');axs[row,col].set_title(f'sample{i} {title}')
        v=cube(samples[i]['P_xyz'])
        for sigma in [1,2]:
            r=next(r for r in records if r['sample']==i and r['sigma']==sigma);ff,aa=plt.subplots(1,3,figsize=(12,4))
            for ax,plane,indices,title in zip(aa,[v[16],v[:,16],v[:,:,16]],[(0,1),(0,2),(1,2)],['XY','XZ','YZ']):
                ax.imshow(plane,cmap='gray',vmin=0,vmax=255);ax.set_title(f'{title} sample{i} sigma{sigma}');ax.plot(16,16,'r+')
                for key,color in [('t_s','cyan'),('image_y','yellow'),('n','red')]:
                    vec=np.array(r[key]);ax.arrow(16,16,10*vec[indices[0]],10*vec[indices[1]],color=color,width=.15,length_includes_head=True)
            ff.suptitle('Cyan: frozen tangent; yellow: image longitudinal axis; red: depth normal');ff.tight_layout();ff.savefig(overlay/f'sample_{i}_sigma{sigma}.png');plt.close(ff)
    fig.tight_layout();fig.savefig(OUT/'t9_vs_t11_patch_comparison.png');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4));names=['T5','T6','T9'];keys=['T5_center','T6_profile_mean','T9_ink']
    for ax,key,title in zip(axs,['mean_r','local_mean_r'],['Sealed adjacent-Z correlation','Sealed local-window correlation']):ax.bar(names,[base[k][key] for k in keys]);ax.set_title(title+'; T11 NOT RUN');ax.axhline(0,color='gray')
    fig.tight_layout();fig.savefig(OUT/'coherence_comparison.png');plt.close(fig)
    verification={'status':'PASS','predecessor_files_unchanged':len(preserve),'blocked_attempt_and_freeze_files_unchanged':len(blocked),'original_pre_evidence_hashes_unchanged':len(frozen),'support_freeze_hashes_unchanged':len(support),'supplemental_hashes_verified':len(receipt['chunks']),'coverage':read(OUT/'coverage_validation.json'),'frame_replay':replay,'old_frame_T9_renderer_check':oldcheck,'original_completed_patch_arrays_match_retry':141,'centers_and_sample_indices_exact':True,'right_handed_geometric_frames_and_original_depth_sign':True,'smoke_inference_replay':'NOT_APPLICABLE_NOT_RUN','full_inference_replay':'NOT_APPLICABLE_NOT_RUN','probability_pngs':'Explicit NOT RUN status placeholders; no probability data generated','original_scientific_gates_unchanged':True}
    # Count actual preserved original arrays rather than infer a declared total.
    verification['original_completed_patch_arrays_match_retry']=len(list((OUT/'frames').glob('patch_*.npz')))
    verification['original_T9_CT_chunks_hash_verified_again']=len(original_ct)
    write(OUT/'verification.json',verification);print(json.dumps(decision),flush=True)
if __name__=='__main__':main()
