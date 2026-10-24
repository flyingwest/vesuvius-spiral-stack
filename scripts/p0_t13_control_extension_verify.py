"""Preservation, exact new-thread replay, CT reassembly, isolation and gate checks."""
import json,re,subprocess,hashlib
from pathlib import Path
import numpy as np
from p0_t13_control_extension import R,O,OLD,NEW,ALL,IMAGE,WORKER,read,dump,sha,freeze_check,thread_dir
from p0_t3_verify import exhaustive_topology

def main():
 freeze_check();prior=read(O/'preservation.json')
 for p,h in prior.items():assert sha(R/p)==h,p
 for stage in ['new_independent_threads','replay','audit_replay']:
  for p,h in read(O/(stage+'_FROZEN.json')).items():assert sha(R/p)==h,p
 plan=read(O/'extension_plan.json');origin=np.array(plan['native_xy_bounds_half_open'][0]);end=np.array(plan['native_xy_bounds_half_open'][1]);receipt=read(O/'acquisition.json')
 for q in receipt:
  assert sha(R/q['file'])==q['sha256'];assert (R/q['file']).stat().st_size==q['bytes'];assert q['headers']['Content-Range']==f"bytes {q['range_start']}-{q['range_end_inclusive']}/2097152"
 checks=[];leak=[];files={};runtime=set()
 for z in NEW:
  native=np.empty((end[1]-origin[1],end[0]-origin[0]),np.uint8)
  for q in receipt:
   if not q['z_first']<=z<=q['z_last']:continue
   _,j,i=q['chunk_zyx'];a=np.frombuffer((R/q['file']).read_bytes(),np.uint8).reshape(4,128,128)[z-q['z_first']];native[j*128-origin[1]:(j+1)*128-origin[1],i*128-origin[0]:(i+1)*128-origin[0]]=a
  calc=np.rint(native.reshape(native.shape[0]//4,4,native.shape[1]//4,4).astype(float).sum((1,3))/16).astype(np.uint8);image=np.load(O/f'ct_only_inputs/slice_{z}.npy');assert np.array_equal(calc,image)
  for name in ['graph.npz','thread.npz']:
   with np.load(thread_dir(z)/name) as a:
    for stage in ['replay','audit_replay']:
     with np.load(O/stage/f'slice_{z}'/name) as b:
      assert a.files==b.files
      for key in a.files:assert np.array_equal(a[key],b[key],equal_nan=True),(z,name,key)
  for name in ['links.json','search.json','decisions.json','metrics.json']:
   for stage in ['replay','audit_replay']:assert read(thread_dir(z)/name)==read(O/stage/f'slice_{z}'/name)
  with np.load(thread_dir(z)/'thread.npz') as a:d={k:a[k].copy() for k in a.files}
  vertices=d['vertices_xy'];topo=exhaustive_topology(vertices[:-1],vertices[1:]);assert not any(topo.values());s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(vertices,axis=0),axis=1))];assert np.array_equal(s,d['edge_s_px']);xy=np.stack([np.interp(d['s_px'],s,vertices[:,k]) for k in [0,1]],1);assert np.array_equal(xy,d['xy']);mask=d['mask'];assert np.array_equal(mask,d['provenance']==2);assert np.isnan(d['ct'][mask]).all();q=xy[~mask];x=q[:,0].astype(int);y=q[:,1].astype(int);dx=q[:,0]-x;dy=q[:,1]-y;bilinear=image[y,x]*(1-dx)*(1-dy)+image[y,x+1]*dx*(1-dy)+image[y+1,x]*(1-dx)*dy+image[y+1,x+1]*dx*dy;error=float(abs(bilinear-d['ct'][~mask]).max());assert error<1e-10
  checks.append({'z':z,'replay':'EXACT_GRAPH_SEARCH_DECISIONS_ARRAYS_BOTH_REPLAYS','CT_pooling_replay':'EXACT','topology':topo,'bilinear_max_error':error})
  for stage in ['new_independent_threads','replay','audit_replay']:
   conf=read(O/stage/f'container_{z}.json');assert conf['network']=='none' and conf['readonly'] and conf['cap_drop']==['ALL'];assert 'no-new-privileges' in conf['security_opt'];mounts={m['Destination']:m for m in conf['mounts']};assert set(mounts)==({'/input.npy','/worker.py','/packages','/output','/strace'} if stage=='audit_replay' else {'/input.npy','/worker.py','/packages','/output'});assert all(not m['RW'] for k,m in mounts.items() if k!='/output');assert mounts['/worker.py']['Source']==str(WORKER) and mounts['/input.npy']['Source']==str(O/f'ct_only_inputs/slice_{z}.npy')
   events=read(O/stage/f'slice_{z}/access_log.json');outside=[e for e in events if not e['path'].startswith(('/packages/','/usr/','/output/')) and e['path']!='/input.npy'];assert not outside,outside;assert sum(e['path']=='/input.npy' for e in events)==1
  opened=[]
  for line in (O/f'audit_replay/syscalls_{z}.txt').read_text().splitlines():
   if ('openat(' in line or 'open(' in line) and 'O_RDONLY' in line and 'O_DIRECTORY' not in line:
    match=re.search(r'= \d+<([^>]+)>',line)
    if match:opened.append(match.group(1))
  assert '/input.npy' in opened and '/worker.py' in opened
  unexpected=[p for p in opened if p not in ['/input.npy','/worker.py'] and not p.startswith(('/packages/','/usr/','/lib/','/etc/','/proc/','/sys/'))];assert not unexpected,unexpected
  for p in set(opened):
   host=R/'.venv/lib/python3.12/site-packages'/p[len('/packages/'):] if p.startswith('/packages/') else WORKER if p=='/worker.py' else O/f'ct_only_inputs/slice_{z}.npy' if p=='/input.npy' else None
   if host:files[f'{z}:{p}' if p=='/input.npy' else p]={'file':str(host.relative_to(R)),'sha256':sha(host)}
   else:runtime.add(p)
  leak.append({'z':z,'isolation':'PASS_ALL3RUNS','Python_file_audit':'PASS_ALL3RUNS','syscall_log':f'audit_replay/syscalls_{z}.txt','successful_readonly_open_paths':sorted(set(opened)),'trusted_or_label_reads':0,'input_sha256':sha(O/f'ct_only_inputs/slice_{z}.npy')})
 cmd=['docker','run','--rm','-i','--entrypoint','python','--network','none','--read-only',IMAGE,'-c',"import sys,json,hashlib,os; paths=json.load(sys.stdin); print(json.dumps({p:hashlib.sha256(open(p,'rb').read()).hexdigest() for p in paths if os.path.isfile(p) and not p.startswith(('/proc/','/sys/'))}))"]
 runtimehash=json.loads(subprocess.check_output(cmd,input=json.dumps(sorted(runtime)),text=True));dump(O/'reconstruction_accessed_file_hashes.json',{'host_files':files,'runtime_image_files':runtimehash,'dynamic_runtime_files_unhashed':[p for p in runtime if p.startswith(('/proc/','/sys/'))],'image':IMAGE,'worker_sha256':sha(WORKER),'strace_sha256':sha(Path('/usr/bin/strace'))})
 dump(O/'leakage_audit.json',{'status':'PASS','new_slices':8,'exact_replays':16,'worker_identical_to576b562':True,'trusted_geometry_label_probability_overlap_reads_during_reconstruction':0,'audits':leak,'scope':'New worker runs see exactly one CT slice; no trusted data or neighboring slices mounted. All24 successful runs have Python-open logs;8instrumented replays have syscall logs. Host evaluation begins only after new freeze. Original7 are reused without reconstruction.','input_and_code_hashes':'reconstruction_accessed_file_hashes.json'})
 overlap=read(O/'trusted_overlap_audit.json');counts={}
 for z in ALL:
  with np.load(O/f'registered_evaluation/slice_{z}.npz') as a,np.load(thread_dir(z)/'thread.npz') as d:
   assert np.array_equal(a['frozen_xyz'][:,:2],d['xy']*4+origin+1.5);assert (a['frozen_xyz'][:,2]==z).all();dist=np.linalg.norm(a['frozen_xyz']-a['trusted_xyz'],axis=1)*2.4;assert np.allclose(dist,a['distance_um'],atol=1e-9);assert np.allclose(np.linalg.norm(a['trusted_normal'],axis=1),1);assert np.allclose((a['trusted_normal']*a['trusted_longitudinal']).sum(1),0,atol=1e-9);counts[z]=int(a['near'].sum())
   eligible=(~d['mask'])&(a['distance_um']<=19.2)&(np.sqrt(np.maximum(0,1-a['trusted_normal'][:,2]**2))>=.25);assert np.array_equal(eligible,a['near'])
   if z not in NEW:
    with np.load(OLD/f'registered_evaluation/slice_{z}.npz') as b:
     for key in a.files:assert np.array_equal(a[key],b[key],equal_nan=True)
 assert len(overlap['windows'])==11
 for w in overlap['windows']:
  expected=[counts[z] for z in range(w['start_z'],w['end_z']+1)];assert expected==[w[f'count_{j}'] for j in range(5)];assert w['pass']==all(n>=25 for n in expected);assert w['total_deficit']==sum(max(25-n,0) for n in expected)
 assert overlap['CONTROL_EXTENSION_RESULT']==('PASS' if any(w['pass'] for w in overlap['windows']) else 'FAIL')
 # Explicitly verify every prior seal, including manifests not transitively covered.
 seals=sorted((R/'outputs').glob('p0_t*/SHA256SUMS'));sealcounts={}
 for seal in seals:
  if seal.parent==O:continue
  n=0
  for line in seal.read_text().splitlines():
   h,p=line.split('  ',1);assert sha(R/p)==h,(str(seal),p);n+=1
  sealcounts[str(seal.relative_to(R))]=n
 result={k:overlap[k] for k in ['CONTROL_EXTENSION_COMPLETE','CONTROL_EXTENSION_RESULT','TOTAL_CONTROL_SLICES','NEW_SLICES_RECONSTRUCTED','BEST_5_SLICE_WINDOW','BEST_WINDOW_COUNTS','TRUSTED_OVERLAP','READY_TO_RESUME_P0_T13']};result.update({'status':'PASS','meaning':'Reproducibility/preservation verification distinct from scientific extension result','LEAKAGE_AUDIT':'PASS','predecessor_entries_unchanged':len(prior),'prior_seals_verified':sealcounts,'old7_threads_and_original_control_evidence_byte_identical':True,'old7_registered_arrays_exact':True,'new_slice_checks':checks,'acquisition_ranges_verified':len(receipt),'all11windows_checked':True,'gate_unchanged':True,'T13_executed':False,'further_expansion':False,'physical_registration_error_um':None})
 dump(O/'verification.json',result);print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
