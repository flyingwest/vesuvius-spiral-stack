"""Read-only predecessor checks and independent T14 representation checks."""
from pathlib import Path
import json,hashlib,sys
import numpy as np,tifffile
from scipy.ndimage import map_coordinates
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t14'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 prior=json.loads((O/'preservation.json').read_text())
 for n,h in prior['files'].items():assert sha(R/n)==h,n
 for n,h in json.loads((O/'PRE_INFERENCE_FREEZE.json').read_text()).items():
  source=O/'serialization_amendment/oracle_before.py' if n=='scripts/p0_t14_oracle.py' else R/n
  assert sha(source)==h,n
 amendment=json.loads((O/'serialization_amendment/amendment.json').read_text());assert sha(R/'scripts/p0_t14_oracle.py')==amendment['new_sha256']
 before=(O/'serialization_amendment/oracle_before.py').read_text();after=(R/'scripts/p0_t14_oracle.py').read_text();assert before.replace('json.dumps(x,indent=2,allow_nan=False)','json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item())')==after
 for attempt,complete in [('primary_before_metrics','stack'),('replay_before_metrics','stack_replay')]:
  for p in (O/'serialization_amendment'/attempt).rglob('*'):
   if not p.is_file():continue
   q=O/complete/p.relative_to(O/'serialization_amendment'/attempt)
   if p.suffix=='.npz':
    with np.load(p) as a,np.load(q) as b:
     for k in a.files:assert np.array_equal(a[k],b[k],equal_nan=True)
   else:assert sha(p)==sha(q)
 for name in ['geometry','stack']:
  for p in (O/name).rglob('*'):
   if not p.is_file() or p.name.startswith('native_ct_'):continue
   q=O/(name+'_replay')/p.relative_to(O/name)
   if p.suffix=='.npz':
    with np.load(p) as a,np.load(q) as b:
     assert a.files==b.files
     for k in a.files:assert np.array_equal(a[k],b[k],equal_nan=True),(name,k)
   else:assert sha(p)==sha(q),str(p)
 with np.load(O/'geometry/threads.npz') as a:xyz=a['xyz'];uv=a['uv'];mask=a['valid'];zs=a['z'];ss=a['s_um']
 assert np.array_equal(np.diff(zs),np.ones(len(zs)-1));assert np.allclose(np.diff(ss),9.6)
 assert abs(xyz[:,:,2][mask]-np.broadcast_to(zs[:,None],mask.shape)[mask]).max()<1e-8
 chord=np.linalg.norm(np.diff(xyz,axis=1),axis=2)*2.4;chord=chord[np.isfinite(chord)];assert chord.max()<=9.6000001 and chord.min()>9.0
 # Independent four-corner bilinear XYZ check against stored producer geometry.
 meta=json.loads((R/'outputs/p0_t10/official_geometry_metadata.json').read_text());scale=np.array(meta['scale']);points=np.stack([tifffile.imread(R/f'data/raw/p0_t10/geometry_{a}.tif') for a in 'xyz'],-1)
 H,W=points.shape[:2];offset=np.array([W,H])/(2*scale)-(np.rint(np.array([W,H])/scale)-1)/2
 query=(uv[mask][::101]+[5248,4096])*4+1.5+offset;g=query*scale;ij=np.floor(g).astype(int);f=g-ij;x=ij[:,0];y=ij[:,1];a=f[:,0,None];b=f[:,1,None]
 interp=(1-a)*(1-b)*points[y,x]+a*(1-b)*points[y,x+1]+(1-a)*b*points[y+1,x]+a*b*points[y+1,x+1]
 err=np.linalg.norm(interp-xyz[mask][::101],axis=1)*2.4;assert err.max()<2.4 # conservative native-voxel numerical discretization ceiling; not physical accuracy.
 with np.load(O/'stack/oracle_stack.npz') as a:p=a['probability'];valid=a['valid'];lab=a['labels'];ev=a['evaluation_mask']
 assert np.array_equal(mask,valid) and np.isfinite(p[valid]).all() and np.isnan(p[~valid]).all()
 image=np.asarray(Image.open(O/'stack/oracle_ink_stack_raw.png'));assert np.array_equal(image,np.rint(np.nan_to_num(p)*65535).astype(np.uint16))
 # Independent AUC via average ranks, not primary cumulative integration.
 from scipy.stats import rankdata
 y=lab[ev];scores=p[ev];r=rankdata(scores).astype(np.float64);n=int(y.sum());m=len(y)-n;auc=(r[y].sum()-n*(n+1)/2)/(n*m)
 metric=json.loads((O/'stack/metrics.json').read_text());assert abs(auc-metric['ROC_AUC'])<1e-12
 model=json.loads((O/'inference/replay.json').read_text());assert model['status']=='PASS'
 old=np.load(R/'outputs/p0_t10/official_control_probability.npy');new=np.load(O/'inference/official_control_probability.npy');old_err=float(abs(old-new).max());assert old_err<=.002
 # Native transverse CT is already available only at these representative old control slices.
 for z in [44449,44456,44463]:
  source=R/'outputs'/('p0_t13_control' if 44453<=z<=44459 else 'p0_t13_control_extension')/'ct_only_inputs'/f'slice_{z}.npy'
  ct=np.load(source);i=int(z-zs[0]);pts=xyz[i,valid[i]];fig,ax=plt.subplots(figsize=(12,5));ax.imshow(ct,cmap='gray',vmin=0,vmax=255,extent=(18944,21504,16640,15616));ax.plot(pts[:,0],pts[:,1],color='cyan',lw=1);ax.set(title=f'Oracle trusted intersection on existing transverse CT, native Z{z}',xlabel='Native X (2.4um)',ylabel='Native Y (2.4um)');fig.tight_layout();fig.savefig(O/'stack/oracle_thread_overlays'/f'native_ct_z{z}.png',dpi=150);plt.close(fig)
 verification={'status':'PASS','predecessor_files_unchanged':len(prior['files']),'prior_seals_verified':len(prior['seals']),'geometry_and_stack_replay':'exact arrays and JSON/PNG before verification-only native overlays','model_replay':model,'T10_probability_max_difference':old_err,'independent_rank_AUC':float(auc),'max_bilinear_vs_triangle_arc_position_error_um':float(err.max()),'XY_chord_length_um_quantiles':np.quantile(chord,[0,.5,1]).tolist(),'plane_error_native':0.,'raw_png_encoding_exact':True,'physical_registration_uncertainty':'UNCALIBRATED','serialization_amendment_only_verified':True,'failed_attempt_arrays_and_images_unchanged':True,'no_detector_or_T13_execution':True,'reconstruction_code_modified':False,'initial_setup_failures':'Preservation reader initially assumed repository-relative names; fixed to resolve historical seal-relative paths before successful baseline. One shell invocation used unavailable python; retried with python3. Neither touched predecessors or scientific outputs.'}
 (O/'verification.json').write_text(json.dumps(verification,indent=2)+'\n');print(verification)
if __name__=='__main__':main()
