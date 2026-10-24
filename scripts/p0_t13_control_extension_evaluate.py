"""Evaluate fixed15 slices only after new threads freeze; no reconstruction changes."""
import numpy as np,tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_t13_control_extension import R,O,OLD,NEW,ALL,read,dump,sha,table,thread_dir,freeze_check
from p0_t13_control_evaluate import section,closest,unit

def main():
 freeze_check()
 for manifest in [O/'new_independent_threads_FROZEN.json',OLD/'independent_threads_FROZEN.json']:
  for p,h in read(manifest).items():assert sha(R/p)==h,p
 assert len(read(O/'new_independent_threads/runs.json'))==8 and all(r['exit_code']==0 for r in read(O/'new_independent_threads/runs.json'))
 # Trusted coordinates first opened after all8 new threads frozen and checked.
 p=np.stack([tifffile.imread(R/f'data/raw/p0_t10/geometry_{a}.tif') for a in 'xyz'],-1).astype(float);region=p[818:923,1048:1154];du=(region[1:-1,2:]-region[1:-1,:-2])/2;dv=(region[2:,1:-1]-region[:-2,1:-1])/2;P=region[1:-1,1:-1];N=unit(np.cross(du,dv));V=unit(dv)
 origin=np.array(read(O/'extension_plan.json')['native_xy_bounds_half_open'][0]);dest=O/'registered_evaluation';dest.mkdir(exist_ok=False);rows=[];metrics=[];roundtrip=0.
 for z in ALL:
  td=thread_dir(z)
  with np.load(td/'thread.npz') as a:d={k:a[k].copy() for k in a.files}
  m=read(td/'metrics.json');xy=d['xy']*4+origin+1.5;res=float(abs((xy-origin-1.5)/4-d['xy']).max());roundtrip=max(roundtrip,res);segments,normals,vs=section(z,P,N,V);j,f,dist=closest(xy,segments);points=segments[j,0]+f[:,None]*(segments[j,1]-segments[j,0]);normal=unit(normals[j,0]+f[:,None]*(normals[j,1]-normals[j,0]));vdir=unit(vs[j,0]+f[:,None]*(vs[j,1]-vs[j,0]));long=np.array([0.,0.,1.])-normal*normal[:,2:3];longnorm=np.linalg.norm(long,axis=1);long=unit(long);spread=np.degrees(np.arccos(np.clip(abs((unit(normals[j,0])*unit(normals[j,1])).sum(1)),0,1)))
  valid=~d['mask'];near=valid&(dist*2.4<=19.2)&(longnorm>=.25);idx=np.flatnonzero(near);span=float(np.ptp(d['s_px'][idx])*9.6) if len(idx)>1 else 0.
  arrays={'frozen_xyz':np.column_stack([xy,np.full(len(xy),z)]),'source_s_px':d['s_px'],'source_provenance':d['provenance'],'source_mask':d['mask'],'trusted_xyz':points,'trusted_normal':normal,'trusted_longitudinal':long,'trusted_v':vdir,'distance_um':dist*2.4,'near':near,'normal_endpoint_variation_deg':spread,'trusted_section_segments':segments}
  if z not in NEW:
   with np.load(OLD/f'registered_evaluation/slice_{z}.npz') as old:
    for key,a in arrays.items():assert np.array_equal(a,old[key],equal_nan=True),(z,key)
  np.savez_compressed(dest/f'slice_{z}.npz',**arrays)
  rows.append({'z':z,'source':'new' if z in NEW else 'sealed_reuse','observed_sample_count':int(valid.sum()),'qualifying_count':int(near.sum()),'qualifying_fraction':float(near.sum()/valid.sum()),'distance_only_count':int((valid&(dist*2.4<=19.2)).sum()),'near_arc_span_um':span,'observed_distance_median_um':float(np.median(dist[valid]*2.4)),'qualifying_distance_median_um':float(np.median(dist[near]*2.4)) if near.any() else None,'qualifying_distance_max_um':float(np.max(dist[near]*2.4)) if near.any() else None,'coordinate_roundtrip_native_voxel':res})
  metrics.append({'z':z,'source':'new' if z in NEW else 'sealed_reuse','status':m['status'],'observed_length_um':m['observed_length_px']*9.6,'bridge_length_um':m['bridge_length_px']*9.6,'observed_samples':m['observed_samples'],'low_confidence_samples':m['low_confidence_samples'],'bridges':m['bridges'],'score':m['score'],**m['final_topology']})
 windows=[]
 for start in range(44449,44460):
  sample=rows[start-44449:start-44449+5];counts=[r['qualifying_count'] for r in sample];counts_pass=min(counts)>=25;span_pass=sample[2]['near_arc_span_um']>=100
  assert not counts_pass or span_pass
  windows.append({'start_z':start,'end_z':start+4,**{f'count_{j}':n for j,n in enumerate(counts)},'minimum_count':min(counts),'total_count':sum(counts),'total_deficit':sum(max(25-n,0) for n in counts),'failing_slices':sum(n<25 for n in counts),'central_near_span_um':sample[2]['near_arc_span_um'],'per_slice_count_gate':counts_pass,'inherited_span_gate':span_pass,'pass':counts_pass and span_pass})
 passing=[w for w in windows if w['pass']];best=passing[0] if passing else min(windows,key=lambda w:(-w['minimum_count'],w['total_deficit'],-w['total_count'],w['start_z']));bounds={}
 for side,ordered in [('lower',rows[:4]),('upper',rows[-4:][::-1])]:
  run=0
  for row in ordered:
   if row['qualifying_count']<25:break
   run+=1
  bounds[side]={'boundary_inward_qualifying_run':run,'boundary_counts_outside_to_inside':[r['qualifying_count'] for r in ordered],'additional_slices_needed_if_all_qualify':5-run if run else None,'prospective_promising_rule_at_least3':run>=3}
 result={'CONTROL_EXTENSION_COMPLETE':'YES','CONTROL_EXTENSION_RESULT':'PASS' if passing else 'FAIL','TOTAL_CONTROL_SLICES':15,'NEW_SLICES_RECONSTRUCTED':8,'BEST_5_SLICE_WINDOW':[best['start_z'],best['end_z']],'BEST_WINDOW_COUNTS':[best[f'count_{j}'] for j in range(5)],'TRUSTED_OVERLAP':'PASS' if passing else 'PARTIAL' if any(r['qualifying_count'] for r in rows) else 'FAIL','READY_TO_RESUME_P0_T13':'YES' if passing else 'NO','windows':windows,'best_window':best,'per_slice':rows,'passing_windows':len(passing),'boundary_continuation':bounds,'further_extension_executed':False,'T13_executed':False,'distance_definition':'unchanged nearest trusted SAME-Z triangle-section point; upper bound on unrestricted3Dsurface distance','frame_eligibility':'unchanged scanZ-projection norm>=.25; observed only; no bridge samples','original7_registered_arrays_exact':True,'coordinate_roundtrip_max_native_voxel':roundtrip,'physical_registration_error_um':None,'physical_registration_uncertainty':'Uncalibrated; no independent calibration supplied in sealed source. Arithmetic residual is NOT physical accuracy.','uncertainty':read(OLD/'trusted_overlap_audit.json')['uncertainty'],'trusted_geometry_read_after_new_threads_frozen':True}
 table(O/'thread_metrics.csv',metrics);table(O/'per_slice_overlap.csv',rows);table(O/'five_slice_windows.csv',windows);dump(O/'trusted_overlap_audit.json',result)
 if passing:
  center=(best['start_z']+best['end_z'])//2
  with np.load(dest/f'slice_{center}.npz') as a:indices=np.flatnonzero(a['near']).tolist()
  dump(O/'control_input_manifest.json',{'selected_window':[best['start_z'],best['end_z']],'center_z':center,'target_sample_indices':indices,'thread_files':{str(z):str((thread_dir(z)/'thread.npz').relative_to(R)) for z in range(best['start_z'],best['end_z']+1)},'evaluation_files':{str(z):str((dest/f'slice_{z}.npz').relative_to(R)) for z in range(best['start_z'],best['end_z']+1)},'trusted_frames_evaluation_only':True,'physical_registration_uncertainty_um':None,'correspondence_not_selected':True})
 fig,axs=plt.subplots(2,1,figsize=(12,7));axs[0].bar(ALL,[r['qualifying_count'] for r in rows],color=['tab:blue' if z in NEW else 'gray' for z in ALL]);axs[0].axhline(25,color='red',ls='--');axs[0].set(title='Qualifying observed points: blue new / gray unchanged; red minimum25',ylabel='Count');axs[1].bar([w['start_z'] for w in windows],[w['minimum_count'] for w in windows]);axs[1].axhline(25,color='red',ls='--');axs[1].set(title='EVERY consecutive five-slice window: minimum per-slice count',xlabel='Window start native Z',ylabel='Minimum count');fig.tight_layout();fig.savefig(O/'overlap_windows.png',dpi=140);plt.close(fig)
 print({k:result[k] for k in ['CONTROL_EXTENSION_RESULT','BEST_5_SLICE_WINDOW','BEST_WINDOW_COUNTS','boundary_continuation']},flush=True)
if __name__=='__main__':main()
