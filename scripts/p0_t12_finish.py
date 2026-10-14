"""Report and verify the stopped T12 geometry branch without model inference."""
from pathlib import Path
import json,hashlib,shutil,collections,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t12'
def read(p):return json.loads(p.read_text())
def dump(p,a):p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    gate=read(OUT/'geometry/geometry_gate.json');assert gate['status']=='FAIL' and gate['stop_result']=='D'
    records=read(OUT/'fit/records.json');controls=read(OUT/'fit/control.json');samples=read(OUT/'samples.json');assert len(records)==72
    assert not any(r['plane']['accepted'] or r['quadratic']['accepted'] for r in records)
    preserved=read(OUT/'preservation.json');frozen=read(OUT/'PRE_EVIDENCE_FREEZE.json')['sha256']
    for manifest in [preserved,frozen]:
        for p,h in manifest.items():assert sha(ROOT/p)==h,p
    assert records==read(OUT/'fit_replay/records.json') and controls==read(OUT/'fit_replay/control.json')
    assert read(OUT/'geometry/metrics.json')==read(OUT/'geometry_replay/metrics.json')
    for i,p in enumerate(samples):
        r=records[i]
        for k in ['P_xyz','z','s_index','tangent_xyz','normal_xyz']:assert r[k]==p[k]
        with np.load(OUT/f'geometry/patch_{i}.npz') as a,np.load(OUT/f'geometry_replay/patch_{i}.npz') as b:
            for name in a.files:assert np.array_equal(a[name],b[name])
            assert np.array_equal(a['plane'],a['T11']) and np.array_equal(a['plane_then_quadratic'],a['T11'])
    for name in ['local_fit_metrics.csv','centering_offsets.csv','geometry_gate.json']:shutil.copyfile(OUT/'geometry'/name,OUT/name)
    with (OUT/'centering_offsets.csv').open() as f:offsets=list(csv.DictReader(f))
    assert all(float(r['displacement_um'])==0 and r['accepted']=='False' for r in offsets)
    med=gate['methods'][0]['medians'];deltas={k:{'T9':med['T9'][k],'T11':med['T11'][k],'T12':med['plane'][k],'T12_minus_T9':med['plane'][k]-med['T9'][k],'relative_reduction_percent':100*(med['T9'][k]-med['plane'][k])/med['T9'][k] if med['T9'][k] else None,'T12_minus_T11':med['plane'][k]-med['T11'][k]} for k in med['T9']}
    control_summary={k:{'values':[c['quality'][k] for c in controls],'median':float(np.median([c['quality'][k] for c in controls]))} for k in med['T9']}
    result={'P0_T12_COMPLETE':'YES','P0_T12_RESULT':'D','LOCAL_SHEET_FIT':'FAIL','PATCH_GEOMETRY_GATE':'FAIL','SMOKE_INFERENCE':'NOT_RUN','FULL_INFERENCE':'NOT_RUN','INK_STRUCTURE_VISIBLE':'NOT_ASSESSED','READY_FOR_NEXT_STEP':'NO','scope':'Operational D: this bounded detector/placement method cannot establish identity; it also rejects all3trusted controls. Not proof of intrinsic CT ambiguity.','selected_method':None,'accepted_target_fits':0,'accepted_control_fits':0,'target_quadratic_attempts':0,'median_applied_displacement_um':0.,'maximum_applied_displacement_um':0.,'fitted_displacement_distribution':None,'fit_residual_distribution':None,'rejection_counts':dict(collections.Counter(r['plane'].get('reason','checks') for r in records)),'deltas':deltas,'control_metrics':control_summary,'control_rejections':[c['fit']['plane'].get('reason','checks') for c in controls],'expected_fitted_patch_support':0,'new_CT_chunks':0,'inherited_supply':'All T12 patch metrics are exactly T11 fallback performance, not T12 local placement improvement.'}
    dump(OUT/'patch_quality_metrics.json',{'summary':result,'all_samples':read(OUT/'geometry/metrics.json'),'control_diagnostics':controls})
    dump(OUT/'decision.json',result);dump(OUT/'sheet_candidate_audit.json',{'target':records,'controls':controls,'interpretation':'All candidates/rejections retained; accepted fit residuals are undefined, not zero.'})
    fig,axs=plt.subplots(2,3,figsize=(14,8));keys=['depth_peak_spread','empty_fraction','band_crossings','center_distance','longitudinal_gradient_fraction','center_slab_occupancy']
    for ax,k in zip(axs.ravel(),keys):
        ax.scatter(np.ones(3),[c['quality'][k] for c in controls],label='trusted controls (n=3)')
        for x,name in enumerate(['T9','T11','plane'],2):ax.scatter(x,med[name][k],label=name if name!='plane' else 'T12 fallback')
        ax.set_xticks([1,2,3,4],['control','T9','T11','T12']);ax.set_title(k)
    fig.suptitle('Same CT-quality metrics; target medians over72; all local fits rejected');fig.tight_layout();fig.savefig(OUT/'control_metric_distributions.png');plt.close(fig)
    fig,axs=plt.subplots(3,6,figsize=(18,9))
    for row,i in enumerate([0,32,64]):
        with np.load(OUT/f'geometry/patch_{i}.npz') as a:
            for col,(name,key,cut) in enumerate([('T9','T9',0),('T11','T11',0),('T12 fallback','plane',0),('T9','T9',1),('T11','T11',1),('T12 fallback','plane',1)]):
                plane=a[key][8] if cut==0 else a[key][:,64];axs[row,col].imshow(plane,cmap='gray',vmin=0,vmax=255,aspect='auto');axs[row,col].set_title(f'sample{i} {name} '+('model XY' if cut==0 else 'depth-X'))
    fig.suptitle('No accepted T12 placement; unchanged T11 fallback shown explicitly');fig.tight_layout();fig.savefig(OUT/'t9_t11_t12_patch_comparison.png');plt.close(fig)
    verification={'status':'PASS','predecessor_files_byte_identical':len(preserved),'pre_evidence_hashes_unchanged':len(frozen),'same_frozen_targets':72,'same_control_locations':3,'fit_replay':read(OUT/'fit_replay.json'),'patch_replay':read(OUT/'patch_replay.json'),'frozen_coordinates_exact':True,'all_T12_patches_equal_T11_fallback':True,'no_applied_center_displacement':True,'coverage':read(OUT/'coverage.json'),'no_model_loaded_or_inference_run':True,'model_replay':'NOT_APPLICABLE_NOT_RUN','no_new_CT_acquisition':True,'geometry_gate_unchanged':True}
    dump(OUT/'verification.json',verification);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
