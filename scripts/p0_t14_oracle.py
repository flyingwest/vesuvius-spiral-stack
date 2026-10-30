"""Oracle-only transverse intersections and physical arc-length stack. Never traces CT."""
from pathlib import Path
import json,hashlib,csv,sys
import numpy as np
from scipy.ndimage import map_coordinates
import tifffile
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t14';T=R/'outputs/p0_t10'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item())+'\n')
def preserve():
 entries={}; seals=[]
 for seal in sorted((R/'outputs').rglob('SHA256SUMS')):
  if O in seal.parents:continue
  # Every existing seal is checked, including nested original audits.
  for line in seal.read_text().splitlines():
   if not line.strip():continue
   h,n=line.split(maxsplit=1);n=n.lstrip('* ');p=R/n
   if not p.is_file():p=seal.parent/n
   assert p.is_file(),n
   n=str(p.relative_to(R))
   assert sha(p)==h,n
   if n in entries:assert entries[n]==h
   entries[n]=h
  entries[str(seal.relative_to(R))]=sha(seal);seals.append(str(seal.relative_to(R)))
 dump(O/'preservation.json',{'files':entries,'seals':seals})
def geometry(dest):
 dest.mkdir(exist_ok=False)
 meta=json.loads((T/'official_geometry_metadata.json').read_text());scale=np.array(meta['scale']);axes=[tifffile.imread(R/f'data/raw/p0_t10/geometry_{a}.tif') for a in 'xyz'];H,W=axes[0].shape
 # Current producer canvas convention: center=count/(2*scale), origin=-(canvas-1)/2.
 canvas=np.rint(np.array([W,H])/scale).astype(int);offset=np.array([W,H])/(2*scale)-(canvas-1)/2
 uv=np.stack(np.meshgrid(np.arange(64,448),np.arange(64,448)),axis=-1).astype(float)
 nativeuv=(uv+np.array([5248,4096]))*4+1.5+offset
 grid=nativeuv*scale
 P=np.stack([map_coordinates(a.astype(float),[grid[:,:,1],grid[:,:,0]],order=1,prefilter=False) for a in axes],-1)
 assert (P>0).all() and (np.diff(P[:,:,2],axis=0)<0).all()
 zlo=int(np.ceil(P[-1,:,2].max()));zhi=int(np.floor(P[0,:,2].min()));assert zhi-zlo+1>=300
 # Consistent two-triangle tessellation; UV is carried only to query producer CT/probability.
 def tris(a):return np.concatenate([np.stack([a[:-1,:-1],a[:-1,1:],a[1:,1:]],-2).reshape(-1,3,a.shape[-1]),np.stack([a[:-1,:-1],a[1:,1:],a[1:,:-1]],-2).reshape(-1,3,a.shape[-1])])
 tri=tris(P);tuv=tris(uv);low=tri[:,:,2].min(1);high=tri[:,:,2].max(1);rows=[];coords=[];uvs=[];lengths=[];maxerr=0.; maxgap=0.
 for z in range(zlo,zhi+1):
  chosen=(low<=z)&(high>z);p=tri[chosen];u=tuv[chosen];pp=[];uu=[]
  for a,b in [(0,1),(1,2),(2,0)]:
   hit=((p[:,a,2]<=z)&(p[:,b,2]>z))|((p[:,b,2]<=z)&(p[:,a,2]>z));pa=p[hit,a];pb=p[hit,b];f=(z-pa[:,2])/(pb[:,2]-pa[:,2]);pp.append(pa+f[:,None]*(pb-pa));uu.append(u[hit,a]+f[:,None]*(u[hit,b]-u[hit,a]))
  pts=np.concatenate(pp);us=np.concatenate(uu);order=np.argsort(us[:,0],kind='stable');pts=pts[order];us=us[order];keep=np.r_[True,np.diff(us[:,0])>1e-8];pts=pts[keep];us=us[keep]
  assert abs(us[0,0]-64)<1e-7 and abs(us[-1,0]-447)<1e-7
  # Strict monotone Z(V) makes the intersection a single U-monotone component.
  maxgap=max(maxgap,float(np.diff(us[:,0]).max()));assert np.diff(us[:,0]).max()<=1.000001
  s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(pts,axis=0),axis=1)*2.4)];ss=np.arange(0,s[-1]+1e-9,9.6)
  xyz=np.column_stack([np.interp(ss,s,pts[:,k]) for k in range(3)]);outuv=np.column_stack([np.interp(ss,s,us[:,k]) for k in range(2)])
  maxerr=max(maxerr,float(abs(xyz[:,2]-z).max()));coords.append(xyz);uvs.append(outuv);lengths.append(s[-1]);rows.append({'z':z,'length_um':float(s[-1]),'samples':len(ss),'components':1,'gaps':0,'start_u':64,'start_v':float(us[0,1]),'end_u':447,'end_v':float(us[-1,1])})
 width=max(map(len,coords));xyz=np.full((len(rows),width,3),np.nan);uvout=np.full((len(rows),width,2),np.nan);mask=np.zeros((len(rows),width),bool)
 for i,(p,u) in enumerate(zip(coords,uvs)):xyz[i,:len(p)]=p;uvout[i,:len(p)]=u;mask[i,:len(p)]=True
 np.savez_compressed(dest/'threads.npz',xyz=xyz,uv=uvout,valid=mask,z=np.arange(zlo,zhi+1),s_um=np.arange(width)*9.6)
 with (dest/'thread_metrics.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 dump(dest/'geometry.json',{'z_range':[zlo,zhi],'slices':len(rows),'height_um':(zhi-zlo)*2.4,'arc_step_um':9.6,'z_step_um':2.4,'length_um_quantiles':np.quantile(lengths,[0,.5,1]).tolist(),'shape':list(mask.shape),'valid_samples':int(mask.sum()),'missing_fraction':float(1-mask.mean()),'max_plane_error_native':maxerr,'maximum_intersection_U_gap_pixels':maxgap,'components_each_row':1,'genuine_gaps':0,'padding_masked':True,'producer_canvas_offset_native_uv':offset.tolist(),'grid_scale':scale.tolist(),'historical_renderer_flags_unknown':True,'normal_sign':'delivered producer depth order preserved','origin':'s=0 at fixed left UV boundary U=5312; no image/label-based row alignment','interpolation':'bilinear stored trusted grid evaluated at level2 pooling centers; piecewise-linear triangle intersections; arc of resulting XYZ polyline','crop':'T10 UV level2 Y4096:4608 X5248:5760, interior64:448','physical_registration_uncertainty':'UNCALIBRATED'})
 print(json.loads((dest/'geometry.json').read_text()),flush=True)
def stack(dest):
 from p0_t10_evaluate import ranking
 dest.mkdir(exist_ok=False)
 with np.load(O/'geometry/threads.npz') as a:uv=a['uv'];valid=a['valid'];zs=a['z'];xyz=a['xyz']
 prob=np.load(O/'inference/official_control_probability.npy');q=uv.copy();q[~valid]=0
 vals=map_coordinates(prob,[q[:,:,1],q[:,:,0]],order=1,prefilter=False);vals[~valid]=np.nan
 with np.load(T/'selected_labels.npz') as a:lab=a['inklabels']>0;sup=a['supervision_mask']>0
 label=map_coordinates(lab.astype(np.uint8),[q[:,:,1],q[:,:,0]],order=0,prefilter=False)>0
 ev=map_coordinates(np.load(T/'evaluation_mask.npy').astype(np.uint8),[q[:,:,1],q[:,:,0]],order=0,prefilter=False)>0;ev &=valid
 auc,ap=ranking(label[ev].astype(int),vals[ev]);pred=vals[ev]>.5;truth=label[ev];tp=int((pred&truth).sum());fp=int((pred&~truth).sum());fn=int((~pred&truth).sum())
 # No claims that oversampled adjacent native planes are independent.
 pair=valid[:-1]&valid[1:];corr=float(np.corrcoef(vals[:-1][pair],vals[1:][pair])[0,1]);stride=4;pair4=valid[:-stride]&valid[stride:];corr4=float(np.corrcoef(vals[:-stride][pair4],vals[stride:][pair4])[0,1])
 np.savez_compressed(dest/'oracle_stack.npz',probability=vals,valid=valid,labels=label,evaluation_mask=ev,z=zs)
 raw=np.rint(np.nan_to_num(vals)*65535).astype(np.uint16);display=np.rint(np.clip((np.nan_to_num(vals)-.25)/.5,0,1)*255).astype(np.uint8)
 Image.fromarray(raw).save(dest/'oracle_ink_stack_raw.png');Image.fromarray(display).save(dest/'oracle_ink_stack_display.png');Image.fromarray((valid*255).astype(np.uint8)).save(dest/'mask_stack.png');Image.fromarray(valid.astype(np.uint8)).save(dest/'provenance_stack.png');Image.fromarray(display).resize((display.shape[1],round((len(zs)-1)*2.4/9.6)+1),Image.Resampling.BOX).save(dest/'physical_aspect_stack.png')
 Image.fromarray((label&ev).astype(np.uint8)*255).save(dest/'ground_truth_stack.png')
 fig,axs=plt.subplots(1,3,figsize=(15,5));extent=[0,display.shape[1]*9.6/1000,zs[-1]*2.4/1000,zs[0]*2.4/1000]
 for ax,a,title in zip(axs,[np.ma.masked_where(~valid,vals),np.ma.masked_where(~ev,label),np.ma.masked_where(~ev,(vals>.5)!=label)],['Oracle probabilities [.25,.75]','Registered ink labels','Threshold disagreement']):
  ax.imshow(a,cmap='gray',vmin=.25 if 'probabilities' in title else 0,vmax=.75 if 'probabilities' in title else 1,extent=extent,aspect='equal');ax.set(title=title,xlabel='Thread arc s (mm)',ylabel='Native Z (mm)')
 fig.tight_layout();fig.savefig(dest/'registered_comparison.png',dpi=160);plt.close(fig)
 overlays=dest/'oracle_thread_overlays';overlays.mkdir();ct=np.load(T/'official_input.npy')[8]
 for i in [0,len(zs)//4,len(zs)//2,3*len(zs)//4,len(zs)-1]:
  fig,axs=plt.subplots(1,2,figsize=(11,4));axs[0].imshow(ct,cmap='gray',vmin=0,vmax=255);axs[0].plot(uv[i,valid[i],0],uv[i,valid[i],1],color='cyan',lw=1);axs[0].set(title=f'Trusted rendered CT with XY-plane intersection Z{zs[i]}',xlabel='Producer U (crop pixel)',ylabel='Producer V (crop pixel)');p=xyz[i,valid[i]];axs[1].plot(p[:,0]*2.4/1000,p[:,1]*2.4/1000);axs[1].set_aspect('equal');axs[1].set(title='Oracle thread in native transverse XY',xlabel='X (mm)',ylabel='Y (mm)');fig.tight_layout();fig.savefig(overlays/f'z{zs[i]}.png',dpi=140);plt.close(fig)
 metrics={'evaluated_samples':int(ev.sum()),'ink_samples':int(label[ev].sum()),'prevalence':float(label[ev].mean()),'ROC_AUC':auc,'average_precision':ap,'IoU_at_05':tp/(tp+fp+fn),'Dice_at_05':2*tp/(2*tp+fp+fn),'mean_ink_probability':float(vals[ev&label].mean()),'mean_nonink_probability':float(vals[ev&~label].mean()),'adjacent_Z_correlation':corr,'four_Z_9p6um_correlation':corr4,'zero_or_nonfinite_valid_probabilities':int((~np.isfinite(vals[valid])).sum()),'unique_rounded_label_pixels':int(len(np.unique(np.rint(q[ev]).astype(int),axis=0))),'no_independent_sample_count_claim':True,'metric_gate_frozen_before_run':{'AUC_ge_075':auc>=.75,'AP_minus_prevalence_ge_015':ap-label[ev].mean()>=.15,'Dice_ge_05':2*tp/(2*tp+fp+fn)>=.5},'scope':'training-corpus oracle upper bound; no detector validation; labels and probabilities share trusted UV association; no physical registration accuracy claim','T13_executed':False}
 dump(dest/'metrics.json',metrics);print(metrics,flush=True)
if __name__=='__main__':
 if sys.argv[1]=='preserve':preserve()
 elif sys.argv[1]=='geometry':geometry(O/sys.argv[2])
 elif sys.argv[1]=='stack':stack(O/sys.argv[2])
