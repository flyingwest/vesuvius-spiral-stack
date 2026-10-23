"""Post-freeze evaluation only. No edits to independent CT-derived threads."""
from pathlib import Path
import json,hashlib,csv
import numpy as np,tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t13_control'
def read(p):return json.loads(p.read_text())
def dump(p,a):p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def table(p,rows):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def unit(a):return a/np.maximum(np.linalg.norm(a,axis=-1,keepdims=True),1e-12)
def section(z,P,N,V):
 # Triangles are evaluation-only tessellation of already trusted local grid.
 # No reconstruction triangulation, no reconstructed point influences triangles.
 faces=[];norm=[];vv=[]
 for iy in range(P.shape[0]-1):
  for ix in range(P.shape[1]-1):
   for off in [[(0,0),(0,1),(1,1)],[(0,0),(1,1),(1,0)]]:
    coords=np.array([P[iy+a,ix+b] for a,b in off]);ns=np.array([N[iy+a,ix+b] for a,b in off]);vs=np.array([V[iy+a,ix+b] for a,b in off])
    if not (coords>0).all() or not (coords[:,2].min()<=z<coords[:,2].max()):continue
    points=[];nsout=[];vsout=[]
    for a,b in [(0,1),(1,2),(2,0)]:
     if (coords[a,2]<=z<coords[b,2]) or (coords[b,2]<=z<coords[a,2]):
      f=(z-coords[a,2])/(coords[b,2]-coords[a,2]);points.append(coords[a]+f*(coords[b]-coords[a]));nsout.append(ns[a]+f*(ns[b]-ns[a]));vsout.append(vs[a]+f*(vs[b]-vs[a]))
    if len(points)==2 and np.linalg.norm(points[0]-points[1])>1e-8:faces.append(points);norm.append(nsout);vv.append(vsout)
 return np.array(faces),np.array(norm),np.array(vv)
def closest(xy,segments):
 a=segments[:,0,:2];b=segments[:,1,:2];d=b-a;den=(d*d).sum(1);near=[];fra=[];dist=[]
 for start in range(0,len(xy),256):
  q=xy[start:start+256];f=np.clip(((q[:,None]-a)*d).sum(2)/den,0,1);point=a+d*f[:,:,None];err=((point-q[:,None])**2).sum(2);j=err.argmin(1);near.extend(j.tolist());fra.extend(f[np.arange(len(q)),j].tolist());dist.extend(np.sqrt(err[np.arange(len(q)),j]).tolist())
 return np.array(near),np.array(fra),np.array(dist)
def main():
 frozen=read(O/'independent_threads_FROZEN.json')
 for p,h in frozen.items():assert sha(R/p)==h,p
 assert all(r['exit_code']==0 for r in read(O/'independent_threads/runs.json'))
 # First trusted-geometry reads in this stage occur AFTER hash freeze check.
 p=np.stack([tifffile.imread(R/f'data/raw/p0_t10/geometry_{a}.tif') for a in 'xyz'],-1).astype(float)
 region=p[818:923,1048:1154];du=(region[1:-1,2:]-region[1:-1,:-2])/2;dv=(region[2:,1:-1]-region[:-2,1:-1])/2
 P=region[1:-1,1:-1];N=unit(np.cross(du,dv));V=unit(dv)
 select=read(O/'selected_slice_range.json');zs=select['native_z_indices'];origin=np.array(select['native_xy_bounds_half_open'][0]);rows=[];threads=[];inputs=[];dest=O/'registered_evaluation';dest.mkdir(exist_ok=False);examples=[];roundtrip=0.
 for z in zs:
  with np.load(O/f'independent_threads/slice_{z}/thread.npz') as a:d={k:a[k].copy() for k in a.files}
  m=read(O/f'independent_threads/slice_{z}/metrics.json');xy=d['xy']*4+origin+1.5;roundtrip=max(roundtrip,float(abs((xy-origin-1.5)/4-d['xy']).max()));segments,normals,vs=section(z,P,N,V);assert len(segments)
  j,f,dist=closest(xy,segments);points=segments[j,0]+f[:,None]*(segments[j,1]-segments[j,0]);normal=unit(normals[j,0]+f[:,None]*(normals[j,1]-normals[j,0]));vdir=unit(vs[j,0]+f[:,None]*(vs[j,1]-vs[j,0]));long=np.array([0.,0.,1.])-normal*normal[:,2:3];longnorm=np.linalg.norm(long,axis=1);long=unit(long)
  normal_endspread=np.degrees(np.arccos(np.clip(abs((unit(normals[j,0])*unit(normals[j,1])).sum(1)),0,1)))
  valid=~d['mask'];near=valid&(dist*2.4<=19.2)&(longnorm>=.25)
  idx=np.flatnonzero(near);span=float(np.ptp(d['s_px'][idx])*9.6) if len(idx)>1 else 0.
  # Full coordinate evidence and evaluation-only associations, no matched cross-Z track.
  xyz=np.column_stack([xy,np.full(len(xy),z)])
  np.savez_compressed(dest/f'slice_{z}.npz',frozen_xyz=xyz,source_s_px=d['s_px'],source_provenance=d['provenance'],source_mask=d['mask'],trusted_xyz=points,trusted_normal=normal,trusted_longitudinal=long,trusted_v=vdir,distance_um=dist*2.4,near=near,normal_endpoint_variation_deg=normal_endspread,trusted_section_segments=segments)
  rows.append({'z':z,'observed_points':int(valid.sum()),'near_points':int(near.sum()),'near_fraction':float(near.sum()/max(valid.sum(),1)),'near_arc_span_um':span,'distance_median_um':float(np.median(dist[valid]*2.4)),'distance_p95_um':float(np.quantile(dist[valid]*2.4,.95)),'near_distance_median_um':float(np.median(dist[near]*2.4)) if near.any() else None,'near_distance_max_um':float(np.max(dist[near]*2.4)) if near.any() else None,'normal_endpoint_variation_median_deg':float(np.median(normal_endspread[near])) if near.any() else None,'coordinate_roundtrip_max_native_voxel':roundtrip,'surface_section_segments':len(segments),'proper_crossings':m['final_topology']['proper'],'nonadjacent_contacts':m['final_topology']['contacts'],'collinear_overlaps':m['final_topology']['overlaps']})
  threads.append({'z':z,'status':m['status'],'observed_length_um':m['observed_length_px']*9.6,'bridge_length_um':m['bridge_length_px']*9.6,'observed_samples':m['observed_samples'],'low_confidence_samples':m['low_confidence_samples'],'bridges':m['bridges'],'score':m['score'],**m['final_topology']})
  if z in select['central_evaluable_z']:
   inputs.append({'z':z,'target_sample_indices':idx.tolist(),'registered_file':str((dest/f'slice_{z}.npz').relative_to(R)),'frozen_thread':f'outputs/p0_t13_control/independent_threads/slice_{z}/thread.npz','neighbor_independent_threads':{str(zz):f'outputs/p0_t13_control/independent_threads/slice_{zz}/thread.npz' for zz in range(z-2,z+3)},'neighbors_have_no_selected_correspondence':True,'trusted_fields_evaluation_only':True})
  examples.append((z,xy,points,near,segments))
 table(O/'registration_metrics.csv',rows);table(O/'thread_metrics.csv',threads)
 central_pass=[]
 for z in select['central_evaluable_z']:
  support=[r for r in rows if abs(r['z']-z)<=2];center=next(r for r in rows if r['z']==z);central_pass.append({'z':z,'pass':all(r['near_points']>=25 for r in support) and center['near_arc_span_um']>=100,'near_counts':[r['near_points'] for r in support]})
 overlap={'status':'PASS' if any(r['pass'] for r in central_pass) else 'PARTIAL' if any(r['near_points'] for r in rows) else 'FAIL','distance_definition':'exact nearest point on piecewise-linear trusted surface intersection with the SAME native Z plane; upper bound on unconstrained3D nearest-surface distance','threshold_um':19.2,'central_windows':central_pass,'rows':rows,'transform':'nativeXYZ=(4*x_trace+x0+1.5,4*y_trace+y0+1.5,nativeZ); same native scan asset ID. No fitted transform or corrected thread.','coordinate_transform_roundtrip_max_native_voxel':roundtrip,'physical_registration_error_um':None,'physical_registration_error_reason':'No publisher registration-error bound or independent fiducials. Arithmetic residual is NOT physical registration accuracy.','uncertainty':{'xy_pooling_cell_half_diagonal_um':float(np.sqrt(2)*4.8),'native_voxel_um':2.4,'trusted_mesh_grid_pitch_parameter_native':20,'normal_variation':'perpoint endpoint-normal variation saved; not a calibrated uncertainty bound','producer_registration_uncertainty_um':None},'gross_topology':'Independent paths have exhaustive zero proper crossings, contacts, overlaps. Spatial overlap with one trusted sheet does not certify entire path identity.','evaluation_reads':[f'data/raw/p0_t10/geometry_{a}.tif' for a in 'xyz'],'thread_freeze_verified_before_geometry_reads':True}
 dump(O/'trusted_overlap_audit.json',overlap)
 dump(O/'control_input_manifest.json',{'status':overlap['status'],'source_selection':'selected_slice_range.json','units':'native2.4umXYZ; source local XY9.6um','frames_evaluation_only':True,'central_inputs':inputs,'uncertainty':overlap['uncertainty'],'registration_error_um':None,'explicit_limit':'No T13 matching, ribbons, frame-error experiment or model inference performed. Full neighboring independent paths supplied without forcing correspondence.','gate_before_leakage_and_replay':central_pass})
 fig,axs=plt.subplots(3,1,figsize=(13,10))
 for ax,z in zip(axs,select['central_evaluable_z']):
  _,xy,points,near,segments=next(q for q in examples if q[0]==z);im=np.load(O/f'ct_only_inputs/slice_{z}.npy');ax.imshow(im,cmap='gray',vmin=0,vmax=255,extent=(origin[0],origin[0]+im.shape[1]*4,origin[1]+im.shape[0]*4,origin[1]));ax.plot(xy[:,0],xy[:,1],color='cyan',lw=.7,label='frozen independent thread');
  for k,seg in enumerate(segments):ax.plot(seg[:,0],seg[:,1],color='yellow',lw=1,label='trusted section (post-freeze only)' if k==0 else None)
  ax.scatter(xy[near,0],xy[near,1],s=3,color='lime',label='within19.2um');ax.set_title(f'Native Z{z}, independent path and evaluation-only surface');ax.legend(fontsize=7)
 fig.tight_layout();fig.savefig(O/'trusted_overlap_overview.png',dpi=140);plt.close(fig)
 print(json.dumps({'overlap':overlap['status'],'central':central_pass,'counts':[r['near_points'] for r in rows]},indent=2))
if __name__=='__main__':main()
