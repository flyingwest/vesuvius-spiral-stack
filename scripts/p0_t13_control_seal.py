"""Final visualization, inventory and checksum sealing; refuses an existing seal."""
from pathlib import Path
import json,hashlib,sys,subprocess
import numpy as np,scipy,skimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t13_control'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 assert not (O/'SHA256SUMS').exists();v=json.loads((O/'verification.json').read_text());assert v['status']=='PASS'
 for p,h in json.loads((O/'preservation.json').read_text()).items():assert sha(R/p)==h,p
 select=json.loads((O/'selected_slice_range.json').read_text());zs=select['native_z_indices'];fig,axs=plt.subplots(4,2,figsize=(16,13))
 for ax,z in zip(axs.flat,zs):
  image=np.load(O/f'ct_only_inputs/slice_{z}.npy')
  with np.load(O/f'independent_threads/slice_{z}/thread.npz') as a:p=a['vertices_xy'];labels=a['edge_provenance']
  ax.imshow(image,cmap='gray',vmin=0,vmax=255);ax.add_collection(LineCollection(np.stack([p[:-1],p[1:]],1),colors=[{1:'cyan',2:'magenta',3:'orange'}[int(k)] for k in labels],linewidths=.7));ax.scatter(*p[0],s=16,color='lime');ax.set_title(f'Z{z}; cyan observed / magenta inferred / orange low confidence')
 axs.flat[-1].axis('off');fig.tight_layout();fig.savefig(O/'independent_threads_overview.png',dpi=140);plt.close(fig)
 m=json.loads((O/'control_input_manifest.json').read_text());m.update({'prerequisite_gate':v['prerequisite_gate'],'READY_TO_RESUME_P0_T13':v['READY_TO_RESUME_P0_T13'],'partial_evaluation_subset_available':True});dump(O/'control_input_manifest.json',m)
 dump(O/'environment.json',{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'skimage':skimage.__version__,'image':'nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a','library_hash_inventory':'reconstruction_accessed_file_hashes.json','scientific_thread_env':{'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}})
 files=[p for p in O.rglob('*') if p.is_file()]+list((R/'scripts').glob('p0_t13_control_*.py'))+[R/'docs/p0_t13_control_thread_prerequisite.md']+list((R/'data/raw/p0_t13_control').glob('*'))
 # Published small raw CT-only inputs and primary/registered arrays are explicitly
 # staged despite generic ignore rules. Redundant replay arrays remain local.
 def published(p):
  s=str(p.relative_to(R))
  return not s.startswith('data/') and not (p.suffix=='.npz' and any(q in p.parts for q in ['replay','audit_replay']))
 inventory=[{'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'planned_git':published(p)} for p in sorted(files)]
 dump(O/'git_artifact_inventory.json',inventory);files.append(O/'git_artifact_inventory.json')
 (O/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(R)}\n' for p in sorted(files)))
 proc=subprocess.run(['sha256sum','-c',str(O/'SHA256SUMS')],cwd=R,capture_output=True,text=True);assert proc.returncode==0,proc.stdout+proc.stderr
 print('sealed and verified',len(files),'files; report,source,all CT payloads and replay evidence')
if __name__=='__main__':main()
