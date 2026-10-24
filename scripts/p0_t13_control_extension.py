"""Separate bounded extension; frozen control worker imported only through Docker."""
from pathlib import Path
import hashlib,json,sys,os,subprocess,itertools,urllib.request,concurrent.futures
import numpy as np
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t13_control_extension';OLD=R/'outputs/p0_t13_control';DATA=R/'data/raw/p0_t13_control_extension'
IMAGE='nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a'
NEW=list(range(44449,44453))+list(range(44460,44464));ALL=list(range(44449,44464))
WORKER=R/'scripts/p0_t13_control_reconstruct.py'
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def table(p,rows):
 import csv
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def thread_dir(z):return (O/'new_independent_threads' if z in NEW else OLD/'independent_threads')/f'slice_{z}'
def setup():
 assert not O.exists();prior=read(OLD/'preservation.json')
 for line in (OLD/'SHA256SUMS').read_text().splitlines():
  h,p=line.split('  ',1);prior[p]=h
 prior['outputs/p0_t13_control/SHA256SUMS']=sha(OLD/'SHA256SUMS')
 for p,h in prior.items():assert sha(R/p)==h,p
 frozen=subprocess.check_output(['git','show','576b562:scripts/p0_t13_control_reconstruct.py'],cwd=R)
 assert hashlib.sha256(frozen).hexdigest()==sha(WORKER)
 O.mkdir();DATA.mkdir(exist_ok=False);dump(O/'preservation.json',prior)
 old=read(OLD/'selected_slice_range.json');plan={'experiment':'P0-T13 control prerequisite bounded extension ONLY','hypothesis':'Nearby independent CT-only threads may extend trusted overlap support; the specified15-slice test cannot itself pass given sealed zero endpoints, but measures boundary continuation prospectively.','base_commit':'576b56209607260c0d4377b2ebd5e95d9a1d140c','total_native_z':ALL,'new_native_z':NEW,'reuse_native_z':list(range(44453,44460)),'source':old['source'],'native_xy_bounds_half_open':old['native_xy_bounds_half_open'],'native_spacing_um':2.4,'xy_spacing_um':9.6,'xy_preprocessing':old['xy_preprocessing'],'source_worker':str(WORKER.relative_to(R)),'worker_sha256':sha(WORKER),'frozen_method':read(OLD/'plan.json')['method'],'gate':{'distance_um':19.2,'minimum_qualifying_points_per_slice':25,'all5_slices_required':True,'frame_projection_min_norm':.25,'central_near_arc_span_um':100,'central_span_note':'Inherited safeguard; any25 distinct0.5px samples span at least115.2um, so redundant with count gate.','no_averaging':True},'window_selection':'earliest passing window; if none, maximize minimum per-slice count, then minimize total point deficit, then maximize total qualifying count, then earliest start','prospective_bound':'Sealed endpoints44453,44459 each have0 qualifying points. Each new window in44449..44463 contains one endpoint. All3old windows failed; no PASS possible without changing sealed evidence/gate.','leakage':'identical isolated Docker and Python open auditing; syscall-instrumented reverse-order replay; all new threads frozen before trusted evaluation','acquisition_cap_bytes':128000000,'prospective_promising_boundary_rule':'At least3 consecutive boundary slices each qualifying>=25 supports a separately authorized test of1-2additional slices; not automatic expansion or a prediction.', 'registration_uncertainty':'No independent calibration in sealed evidence; physical error unknown; numerical roundtrip not physical accuracy','stop':'No expansion beyond44449..44463; no T13 execution regardless of result'}
 dump(O/'extension_plan.json',plan)
 deps=[R/'scripts/p0_t13_control_extension.py',R/'scripts/p0_t13_control_extension_evaluate.py',R/'scripts/p0_t13_control_extension_verify.py',WORKER,R/'scripts/p0_t13_control_evaluate.py',O/'extension_plan.json']
 dump(O/'PRE_EVIDENCE_FREEZE.json',{str(p.relative_to(R)):sha(p) for p in deps});print('preserved',len(prior),'worker equals576b562',flush=True)
def freeze_check():
 for p,h in read(O/'PRE_EVIDENCE_FREEZE.json').items():assert sha(R/p)==h,p

def acquire():
 freeze_check();plan=read(O/'extension_plan.json');low=np.array(plan['native_xy_bounds_half_open'][0]);high=np.array(plan['native_xy_bounds_half_open'][1]);jobs=[]
 for first,last in [(44449,44452),(44460,44463)]:
  assert first//128==last//128;k=first//128;a=first-k*128;b=last-k*128+1
  for j,i in itertools.product(range(low[1]//128,high[1]//128),range(low[0]//128,high[0]//128)):
   jobs.append({'chunk_zyx':[k,int(j),int(i)],'z_first':first,'z_last':last,'url':plan['source']+f'0/{k}/{j}/{i}','range_start':a*128**2,'range_end_inclusive':b*128**2-1})
 assert sum(q['range_end_inclusive']-q['range_start']+1 for q in jobs)<=plan['acquisition_cap_bytes'];dump(O/'acquisition_plan.json',jobs)
 def fetch(q):
  req=urllib.request.Request(q['url'],headers={'Range':f"bytes={q['range_start']}-{q['range_end_inclusive']}"})
  with urllib.request.urlopen(req,timeout=60) as r:
   assert r.status==206;headers=dict(r.headers);payload=r.read(128**3+1)
  assert headers['Content-Range']==f"bytes {q['range_start']}-{q['range_end_inclusive']}/2097152";assert len(payload)==q['range_end_inclusive']-q['range_start']+1
  path=DATA/(f"z{q['z_first']}_"+'_'.join(map(str,q['chunk_zyx']))+'.bin');assert not path.exists();path.write_bytes(payload)
  return {**q,'file':str(path.relative_to(R)),'sha256':sha(path),'bytes':len(payload),'headers':headers}
 with concurrent.futures.ThreadPoolExecutor(8) as pool:receipt=list(pool.map(fetch,jobs))
 dump(O/'acquisition.json',receipt);dest=O/'ct_only_inputs';dest.mkdir();files=[]
 for z in NEW:
  a=np.empty((high[1]-low[1],high[0]-low[0]),np.uint8)
  for q in receipt:
   if not q['z_first']<=z<=q['z_last']:continue
   _,j,i=q['chunk_zyx'];raw=np.frombuffer((R/q['file']).read_bytes(),np.uint8).reshape(4,128,128)[z-q['z_first']];a[j*128-low[1]:(j+1)*128-low[1],i*128-low[0]:(i+1)*128-low[0]]=raw
  pooled=np.rint(a.reshape(a.shape[0]//4,4,a.shape[1]//4,4).mean((1,3))).astype(np.uint8);path=dest/f'slice_{z}.npy';np.save(path,pooled);files.append({'z':z,'file':str(path.relative_to(R)),'sha256':sha(path),'shape':list(pooled.shape)})
 dump(O/'CT_INPUTS_FROZEN.json',{'files':files,'payload_bytes':sum(q['bytes'] for q in receipt)});print('new CT inputs frozen',len(files),'bytes',sum(q['bytes'] for q in receipt),flush=True)
def run(stage):
 freeze_check();assert stage in ['new_independent_threads','replay','audit_replay'];dest=O/stage;dest.mkdir(exist_ok=False);runs=[]
 for item in read(O/'CT_INPUTS_FROZEN.json')['files']:assert sha(R/item['file'])==item['sha256']
 audit=stage=='audit_replay'
 for z in NEW if stage=='new_independent_threads' else NEW[::-1]:
  name=f't13ext-{stage.replace("_","-")}-{z}';cmd=['docker','create','--name',name,'--user',f'{os.getuid()}:{os.getgid()}','--entrypoint','/strace' if audit else 'python','--network','none','--read-only','--cap-drop','ALL','--security-opt','no-new-privileges','--tmpfs','/tmp:rw,noexec,nosuid,size=128m','-e','PYTHONPATH=/packages','-e','OPENBLAS_NUM_THREADS=1','-e','OMP_NUM_THREADS=1','-e','MKL_NUM_THREADS=1','-v',str(R/'.venv/lib/python3.12/site-packages')+':/packages:ro','-v',str(WORKER)+':/worker.py:ro','-v',str(O/f'ct_only_inputs/slice_{z}.npy')+':/input.npy:ro','-v',str(dest)+':/output:rw']
  if audit:cmd+=['-v','/usr/bin/strace:/strace:ro']
  cmd+=[IMAGE]
  if audit:cmd+=['-f','-yy','-s','4096','-e','trace=%file','-o',f'/output/syscalls_{z}.txt','python']
  cmd+=['-B','/worker.py','/input.npy',f'/output/slice_{z}'];cid=subprocess.check_output(cmd,text=True).strip();ins=json.loads(subprocess.check_output(['docker','inspect',cid],text=True))[0]
  dump(dest/f'container_{z}.json',{'command':cmd,'id':cid,'image':ins['Image'],'mounts':ins['Mounts'],'network':ins['HostConfig']['NetworkMode'],'readonly':ins['HostConfig']['ReadonlyRootfs'],'cap_drop':ins['HostConfig']['CapDrop'],'security_opt':ins['HostConfig']['SecurityOpt']})
  proc=subprocess.run(['docker','start','-a',cid],capture_output=True,text=True);(dest/f'log_{z}.txt').write_text(proc.stdout+proc.stderr);state=json.loads(subprocess.check_output(['docker','inspect',cid],text=True))[0]['State'];subprocess.run(['docker','rm',cid],check=True,stdout=subprocess.DEVNULL);runs.append({'z':z,'exit_code':state['ExitCode']});dump(dest/'runs.json',runs);print(stage,z,state['ExitCode'],flush=True)
  if state['ExitCode']!=0:raise RuntimeError(proc.stdout+proc.stderr)
 dump(O/(stage+'_FROZEN.json'),{str(p.relative_to(R)):sha(p) for p in sorted(dest.rglob('*')) if p.is_file()})
if __name__=='__main__':
 stage=sys.argv[1]
 if stage=='setup':setup()
 elif stage=='acquire':acquire()
 else:run(stage)
