"""Acquire and launch only the preregistered alternate region, with frozen worker."""
from pathlib import Path
import sys,os,json,hashlib,subprocess,itertools,urllib.request,concurrent.futures
import numpy as np
from p0_t13_control_extension import IMAGE,WORKER,sha,read,dump,table
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t13_control_alt';OLD=R/'outputs/p0_t13_control';DATA=R/'data/raw/p0_t13_control_alt'
PLAN=read(O/'selected_region.json') if (O/'selected_region.json').exists() else {}
NEW=ALL=PLAN.get('native_z_indices',[])
def thread_dir(z):return O/'independent_threads'/f'slice_{z}'
def freeze_check():
 for p,h in read(O/'PRE_RECONSTRUCTION_FREEZE.json').items():assert sha(R/p)==h,p

def acquire():
 freeze_check();DATA.mkdir(exist_ok=False);low=np.array(PLAN['native_xy_bounds_half_open'][0]);high=np.array(PLAN['native_xy_bounds_half_open'][1]);jobs=[]
 for k in range(ALL[0]//128,ALL[-1]//128+1):
  first=max(ALL[0],k*128);last=min(ALL[-1],k*128+127)
  for j,i in itertools.product(range(low[1]//128,high[1]//128),range(low[0]//128,high[0]//128)):
   jobs.append({'chunk_zyx':[k,int(j),int(i)],'z_first':first,'z_last':last,'url':PLAN['source']+f'0/{k}/{j}/{i}','range_start':(first-k*128)*128**2,'range_end_inclusive':(last-k*128+1)*128**2-1})
 assert sum(q['range_end_inclusive']-q['range_start']+1 for q in jobs)<=PLAN['acquisition_cap_bytes'];dump(O/'acquisition_plan.json',jobs)
 def fetch(q):
  req=urllib.request.Request(q['url'],headers={'Range':f"bytes={q['range_start']}-{q['range_end_inclusive']}"})
  with urllib.request.urlopen(req,timeout=60) as r:
   assert r.status==206;headers=dict(r.headers);payload=r.read(128**3+1)
  assert headers['Content-Range']==f"bytes {q['range_start']}-{q['range_end_inclusive']}/2097152";assert len(payload)==q['range_end_inclusive']-q['range_start']+1
  path=DATA/(f"z{q['z_first']}_"+'_'.join(map(str,q['chunk_zyx']))+'.bin');assert not path.exists();path.write_bytes(payload)
  return {**q,'file':str(path.relative_to(R)),'sha256':sha(path),'bytes':len(payload),'headers':headers}
 with concurrent.futures.ThreadPoolExecutor(8) as pool:receipt=list(pool.map(fetch,jobs))
 dump(O/'acquisition.json',receipt);dest=O/'ct_only_inputs';dest.mkdir();files=[]
 for z in ALL:
  a=np.empty((high[1]-low[1],high[0]-low[0]),np.uint8)
  for q in receipt:
   if not q['z_first']<=z<=q['z_last']:continue
   _,j,i=q['chunk_zyx'];raw=np.frombuffer((R/q['file']).read_bytes(),np.uint8).reshape(q['z_last']-q['z_first']+1,128,128)[z-q['z_first']];a[j*128-low[1]:(j+1)*128-low[1],i*128-low[0]:(i+1)*128-low[0]]=raw
  pooled=np.rint(a.reshape(a.shape[0]//4,4,a.shape[1]//4,4).mean((1,3))).astype(np.uint8);path=dest/f'slice_{z}.npy';np.save(path,pooled);files.append({'z':z,'file':str(path.relative_to(R)),'sha256':sha(path),'shape':list(pooled.shape)})
 dump(O/'CT_INPUTS_FROZEN.json',{'files':files,'payload_bytes':sum(q['bytes'] for q in receipt)});print('CT inputs frozen',len(files),'bytes',sum(q['bytes'] for q in receipt),flush=True)
def run(stage):
 freeze_check();assert stage in ['independent_threads','replay','audit_replay'];dest=O/stage;dest.mkdir(exist_ok=False);runs=[]
 for item in read(O/'CT_INPUTS_FROZEN.json')['files']:assert sha(R/item['file'])==item['sha256']
 audit=stage=='audit_replay'
 for z in NEW if stage=='independent_threads' else NEW[::-1]:
  name=f't13alt-{stage.replace("_","-")}-{z}';cmd=['docker','create','--name',name,'--user',f'{os.getuid()}:{os.getgid()}','--entrypoint','/strace' if audit else 'python','--network','none','--read-only','--cap-drop','ALL','--security-opt','no-new-privileges','--tmpfs','/tmp:rw,noexec,nosuid,size=128m','-e','PYTHONPATH=/packages','-e','OPENBLAS_NUM_THREADS=1','-e','OMP_NUM_THREADS=1','-e','MKL_NUM_THREADS=1','-v',str(R/'.venv/lib/python3.12/site-packages')+':/packages:ro','-v',str(WORKER)+':/worker.py:ro','-v',str(O/f'ct_only_inputs/slice_{z}.npy')+':/input.npy:ro','-v',str(dest)+':/output:rw']
  if audit:cmd+=['-v','/usr/bin/strace:/strace:ro']
  cmd+=[IMAGE]
  if audit:cmd+=['-f','-yy','-s','4096','-e','trace=%file','-o',f'/output/syscalls_{z}.txt','python']
  cmd+=['-B','/worker.py','/input.npy',f'/output/slice_{z}'];cid=subprocess.check_output(cmd,text=True).strip();ins=json.loads(subprocess.check_output(['docker','inspect',cid],text=True))[0]
  dump(dest/f'container_{z}.json',{'command':cmd,'id':cid,'image':ins['Image'],'mounts':ins['Mounts'],'network':ins['HostConfig']['NetworkMode'],'readonly':ins['HostConfig']['ReadonlyRootfs'],'cap_drop':ins['HostConfig']['CapDrop'],'security_opt':ins['HostConfig']['SecurityOpt']})
  proc=subprocess.run(['docker','start','-a',cid],capture_output=True,text=True);(dest/f'log_{z}.txt').write_text(proc.stdout+proc.stderr);state=json.loads(subprocess.check_output(['docker','inspect',cid],text=True))[0]['State'];subprocess.run(['docker','rm',cid],check=True,stdout=subprocess.DEVNULL);runs.append({'z':z,'exit_code':state['ExitCode']});dump(dest/'runs.json',runs);print(stage,z,state['ExitCode'],flush=True)
  if state['ExitCode']!=0:raise RuntimeError(proc.stdout+proc.stderr)
 dump(O/(stage+'_FROZEN.json'),{str(p.relative_to(R)):sha(p) for p in sorted(dest.rglob('*')) if p.is_file()})

if __name__=='__main__':
 if sys.argv[1]=='acquire':acquire()
 else:run(sys.argv[1])
