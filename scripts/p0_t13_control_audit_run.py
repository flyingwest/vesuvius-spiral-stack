"""Launch isolated, single-slice workers; only CT, code, libraries, output mounted."""
from pathlib import Path
import subprocess,json,hashlib,sys,os
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t13_control'
IMAGE='nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a'
def main():
 stage=sys.argv[1];assert stage == 'audit_replay';dest=O/stage;dest.mkdir(exist_ok=False)
 zs=json.loads((O/'selected_slice_range.json').read_text())['native_z_indices'];runs=[]
 for z in (zs if stage=='independent_threads' else zs[::-1]):
  name=f't13-control-{stage.replace("_","-")}-{z}'
  cmd=['docker','create','--name',name,'--user',f'{os.getuid()}:{os.getgid()}','--entrypoint','/strace','--network','none','--read-only','--cap-drop','ALL','--security-opt','no-new-privileges','--tmpfs','/tmp:rw,noexec,nosuid,size=128m','-e','PYTHONPATH=/packages','-e','OPENBLAS_NUM_THREADS=1','-e','OMP_NUM_THREADS=1','-e','MKL_NUM_THREADS=1','-v','/usr/bin/strace:/strace:ro','-v',str(R/'.venv/lib/python3.12/site-packages')+':/packages:ro','-v',str(R/'scripts/p0_t13_control_reconstruct.py')+':/worker.py:ro','-v',str(O/f'ct_only_inputs/slice_{z}.npy')+':/input.npy:ro','-v',str(dest)+':/output:rw',IMAGE,'-f','-yy','-s','4096','-e','trace=%file','-o',f'/output/syscalls_{z}.txt','python','-B','/worker.py','/input.npy',f'/output/slice_{z}']
  cid=subprocess.check_output(cmd,text=True).strip();inspect=json.loads(subprocess.check_output(['docker','inspect',cid],text=True))[0]
  (dest/f'container_{z}.json').write_text(json.dumps({'command':cmd,'id':cid,'image':inspect['Image'],'mounts':inspect['Mounts'],'network':inspect['HostConfig']['NetworkMode'],'readonly':inspect['HostConfig']['ReadonlyRootfs'],'cap_drop':inspect['HostConfig']['CapDrop'],'security_opt':inspect['HostConfig']['SecurityOpt']},indent=2)+'\n')
  proc=subprocess.run(['docker','start','-a',cid],capture_output=True,text=True);(dest/f'log_{z}.txt').write_text(proc.stdout+proc.stderr)
  state=json.loads(subprocess.check_output(['docker','inspect',cid],text=True))[0]['State'];subprocess.run(['docker','rm',cid],check=True,stdout=subprocess.DEVNULL)
  runs.append({'z':z,'exit_code':state['ExitCode']});(dest/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(stage,z,'exit',state['ExitCode'],flush=True)
  if state['ExitCode']!=0:raise RuntimeError(proc.stdout+proc.stderr)
 files=[p for p in dest.rglob('*') if p.is_file()];(O/(stage+'_FROZEN.json')).write_text(json.dumps({str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
if __name__=='__main__':main()
