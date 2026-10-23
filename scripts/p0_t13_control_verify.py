"""Independent saved-data verification and leakage audit; never retraces targets."""
from pathlib import Path
import json,hashlib,re,subprocess,csv
import numpy as np
from p0_t3_verify import exhaustive_topology
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t13_control'
IMAGE='nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a'
def read(p):return json.loads(p.read_text())
def dump(p,a):p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 prior=read(O/'preservation.json')
 for p,h in prior.items():assert sha(R/p)==h,p
 freeze=read(O/'PRE_RECONSTRUCTION_FREEZE.json')
 for p,h in freeze.items():
  source=O/'attempt_01_permission_denied/launcher_original.py' if p=='scripts/p0_t13_control_run.py' else R/p
  assert sha(source)==h,p
 assert sha(R/'scripts/p0_t13_control_run.py')==read(O/'infrastructure_amendment.json')['launcher_sha256']
 for stage in ['independent_threads','replay','audit_replay']:
  for p,h in read(O/(stage+'_FROZEN.json')).items():assert sha(R/p)==h,p
 select=read(O/'selected_slice_range.json');zs=select['native_z_indices'];origin=np.array(select['native_xy_bounds_half_open'][0]);end=np.array(select['native_xy_bounds_half_open'][1]);volume=np.empty((len(zs),int(end[1]-origin[1]),int(end[0]-origin[0])),np.uint8)
 receipts=read(O/'acquisition.json')
 for q in receipts:
  p=R/q['file'];assert sha(p)==q['sha256'] and p.stat().st_size==q['bytes'];assert q['headers']['Content-Range']==f"bytes {q['range_start']}-{q['range_end_inclusive']}/2097152"
  k,j,i=q['chunk_zyx'];a,b=q['local_depth'];arr=np.frombuffer(p.read_bytes(),np.uint8).reshape(b-a,128,128);volume[k*128+a-zs[0]:k*128+b-zs[0],j*128-origin[1]:(j+1)*128-origin[1],i*128-origin[0]:(i+1)*128-origin[0]]=arr
 checks=[];leak=[];files={};runtime=set();totalnear=0;observed=0
 for k,z in enumerate(zs):
  image=np.load(O/f'ct_only_inputs/slice_{z}.npy');a=volume[k];expected=np.rint(a.reshape(a.shape[0]//4,4,a.shape[1]//4,4).astype(float).sum((1,3))/16).astype(np.uint8);assert np.array_equal(image,expected)
  primary=O/f'independent_threads/slice_{z}'
  for name in ['graph.npz','thread.npz']:
   with np.load(primary/name) as a:
    for stage in ['replay','audit_replay']:
     with np.load(O/stage/f'slice_{z}'/name) as b:
      assert a.files==b.files
      for key in a.files:assert np.array_equal(a[key],b[key],equal_nan=True),(z,name,key)
  for name in ['links.json','search.json','decisions.json','metrics.json']:
   for stage in ['replay','audit_replay']:assert read(primary/name)==read(O/stage/f'slice_{z}'/name),(z,name)
  with np.load(primary/'thread.npz') as a:d={key:a[key].copy() for key in a.files}
  p=d['vertices_xy'];topo=exhaustive_topology(p[:-1],p[1:]);assert not any(topo.values());ss=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))];assert np.array_equal(ss,d['edge_s_px']);xy=np.stack([np.interp(d['s_px'],ss,p[:,j]) for j in [0,1]],1);assert np.array_equal(xy,d['xy']);mask=d['mask'];assert np.array_equal(mask,d['provenance']==2) and np.isnan(d['ct'][mask]).all()
  q=xy[~mask];x=q[:,0].astype(int);y=q[:,1].astype(int);dx=q[:,0]-x;dy=q[:,1]-y;calc=image[y,x]*(1-dx)*(1-dy)+image[y,x+1]*dx*(1-dy)+image[y+1,x]*(1-dx)*dy+image[y+1,x+1]*dx*dy;error=float(abs(calc-d['ct'][~mask]).max());assert error<1e-10
  with np.load(O/f'registered_evaluation/slice_{z}.npz') as ev:
   assert np.array_equal(ev['frozen_xyz'][:,:2],xy*4+origin+1.5);assert (ev['frozen_xyz'][:,2]==z).all();assert np.allclose(ev['trusted_xyz'][:,2],z,atol=1e-9);dist=np.linalg.norm(ev['trusted_xyz']-ev['frozen_xyz'],axis=1)*2.4;assert np.allclose(dist,ev['distance_um'],atol=1e-9);assert np.allclose(np.linalg.norm(ev['trusted_normal'],axis=1),1);assert np.allclose((ev['trusted_normal']*ev['trusted_longitudinal']).sum(1),0,atol=1e-9);totalnear+=int(ev['near'].sum());observed+=int((~mask).sum())
  checks.append({'z':z,'full_graph_search_thread_replay':'EXACT_BOTH_REPLAYS','independent_topology':topo,'bilinear_max_error':error,'native_CT_pooling_replay':'EXACT','registration_arrays':'PASS'})
  for stage in ['independent_threads','replay','audit_replay']:
   conf=read(O/stage/f'container_{z}.json');assert conf['network']=='none' and conf['readonly'] and conf['cap_drop']==['ALL'];mounts={m['Destination']:m for m in conf['mounts']};assert set(mounts)==({'/packages','/worker.py','/input.npy','/output','/strace'} if stage=='audit_replay' else {'/packages','/worker.py','/input.npy','/output'});assert all(not m['RW'] for k,m in mounts.items() if k!='/output');assert mounts['/input.npy']['Source']==str(O/f'ct_only_inputs/slice_{z}.npy')
   events=read(O/stage/f'slice_{z}/access_log.json');outside=[e for e in events if not e['path'].startswith(('/packages/','/usr/','/output/')) and e['path']!='/input.npy'];assert not outside,outside
   inputreads=[e for e in events if e['path']=='/input.npy'];assert len(inputreads)==1
  trace=(O/f'audit_replay/syscalls_{z}.txt').read_text();opened=[]
  for line in trace.splitlines():
   if ('openat(' in line or 'open(' in line) and 'O_RDONLY' in line and 'O_DIRECTORY' not in line:
    match=re.search(r'= \d+<([^>]+)>',line)
    if match:opened.append(match.group(1))
  assert '/input.npy' in opened and '/worker.py' in opened
  allowed=('/packages/','/usr/','/lib/','/etc/','/proc/','/sys/')
  unexpected=[p for p in opened if p not in ['/input.npy','/worker.py'] and not p.startswith(allowed)];assert not unexpected,unexpected
  for p in set(opened):
   host=R/'.venv/lib/python3.12/site-packages'/p[len('/packages/'):] if p.startswith('/packages/') else (R/'scripts/p0_t13_control_reconstruct.py' if p=='/worker.py' else O/f'ct_only_inputs/slice_{z}.npy' if p=='/input.npy' else None)
   if host:files[f'{z}:{p}' if p=='/input.npy' else p]={'host_file':str(host.relative_to(R)),'sha256':sha(host)}
   else:runtime.add(p)
  leak.append({'z':z,'container_isolation':'PASS_3RUNS','python_open_audit':'PASS_3RUNS','syscall_log':f'audit_replay/syscalls_{z}.txt','successful_readonly_open_paths':sorted(set(opened)),'trusted_or_label_reads':0,'exact_CT_input_sha256':sha(O/f'ct_only_inputs/slice_{z}.npy')})
 # Hash actual runtime image files in the same pinned image; no scientific data.
 cmd=['docker','run','--rm','-i','--entrypoint','python','--network','none','--read-only',IMAGE,'-c',"import sys,json,hashlib,os; paths=json.load(sys.stdin); print(json.dumps({p:hashlib.sha256(open(p,'rb').read()).hexdigest() for p in paths if os.path.isfile(p) and not p.startswith(('/proc/','/sys/'))}))"]
 runtimehash=json.loads(subprocess.check_output(cmd,input=json.dumps(sorted(runtime)),text=True))
 dump(O/'reconstruction_accessed_file_hashes.json',{'scientific_and_library_files':files,'pinned_image_runtime_files':runtimehash,'dynamic_runtime_files_unhashed':[p for p in runtime if p.startswith(('/proc/','/sys/'))],'image':IMAGE,'strace_binary_sha256':sha(Path('/usr/bin/strace'))})
 leakage={'status':'PASS','reconstruction_runs':7,'exact_replays':14,'scientific_input_reads':'one CT-only slice per worker, no neighboring slice read','trusted_geometry_or_label_reads':0,'audits':leak,'hash_inventory':'reconstruction_accessed_file_hashes.json','source_sha256':sha(R/'scripts/p0_t13_control_reconstruct.py'),'scope':'Python open logs for primary and both replays; syscall file-access logs for second exact replay include loader/native library accesses. Read-only mount/network isolation applies to all runs.','selection_disclosure':'Host read trusted geometry solely to select slice/crop bounds before tracing, as Task1. CT intensity files contain no trusted coordinates. Evaluation read geometry only after all threads hash-frozen. Inherited T10 crop originally selected using labels; no labels read in this experiment.','failed_attempt':'Initial permission failure before CT read preserved; no contaminated reconstruction attempt.'}
 dump(O/'leakage_audit.json',leakage)
 overlap=read(O/'trusted_overlap_audit.json');gates={'at_least5_consecutive_independent_threads':len(checks)>=5,'central_slice_has_Zpm2_neighbors':len(zs)>=5,'meaningful_trusted_overlap':overlap['status']=='PASS','trusted_frames_associated_postfreeze':totalnear>0,'leakage_audit_pass':True}
 verification={'status':'PASS','meaning':'execution/replay/preservation verification, not scientific prerequisite acceptance','CONTROL_THREAD_PREREQUISITE_COMPLETE':'YES','INDEPENDENT_CONTROL_THREADS':'PASS','CONTROL_SLICE_COUNT':len(zs),'TRUSTED_OVERLAP':overlap['status'],'LEAKAGE_AUDIT':'PASS','READY_TO_RESUME_P0_T13':'YES' if all(gates.values()) else 'NO','prerequisite_gate':'PASS' if all(gates.values()) else 'FAIL','gate_checks':gates,'near_samples':totalnear,'observed_samples':observed,'near_fraction':totalnear/observed,'predecessor_files_unchanged':len(prior),'coordinate_transform_roundtrip_native_voxel':overlap['coordinate_transform_roundtrip_max_native_voxel'],'physical_registration_error_um':None,'checks':checks,'T13_executed':False,'model_loaded':False,'acquisition_payloads_verified':len(receipts),'all_CT_inputs_replayed_from_range_payloads':True,'pre_reconstruction_freeze':'PASS; launcher-only infrastructure amendment checked against preserved original'}
 dump(O/'verification.json',verification);print(json.dumps(verification,indent=2))
if __name__=='__main__':main()
