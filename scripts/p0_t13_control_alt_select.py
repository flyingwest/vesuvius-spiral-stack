"""Neutral geometry/availability selection, before ANY alternate reconstruction."""
from pathlib import Path
import json,hashlib,urllib.request,itertools,concurrent.futures,subprocess
import numpy as np,tifffile
from p0_t13_control_evaluate import section,unit
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t13_control_alt';BASE='https://vesuvius-challenge-open-data.s3.amazonaws.com/PHercParis4/volumes/20260411134726-2.400um-0.2m-78keV-masked.zarr/'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,a):p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
def main():
 O.mkdir(exist_ok=False)
 prior=json.loads((R/'outputs/p0_t13_control_extension/preservation.json').read_text())
 for line in (R/'outputs/p0_t13_control_extension/SHA256SUMS').read_text().splitlines():
  h,p=line.split('  ',1);prior[p]=h
 prior['outputs/p0_t13_control_extension/SHA256SUMS']=sha(R/'outputs/p0_t13_control_extension/SHA256SUMS')
 for p,h in prior.items():assert sha(R/p)==h,p
 dump(O/'preservation.json',prior)
 policy={'stage':'neutral selection BEFORE reconstruction','candidate_order':'centers at (H//4,W//2), (3H//4,W//2), (H//2,W//4), (H//2,3W//4); first eligible. No selection after thread evidence.','grid_half_size_yx':[51,52],'slice_count':7,'z_center':'rounded median native Z of candidate grid','excluded_z':'candidate7-slice range must be separated from44449..44463 by at least1024 native slices','geometry_gates':{'all_grid_coordinates_positive':True,'minimum_3D_edge_length_native':.01,'maximum_3D_edge_length_native':200,'maximum_neighbor_normal_angle_deg':45,'minimum_each_section_arc_native':500,'single_connected_section':True,'minimum_section_to_CT_crop_boundary_native':64},'crop':'full candidate projectedXYbbox plus128native pixels, round outward128; same4x4XY pooling as frozen worker','coverage':'HEAD every native chunk touched by selected7planes; require size2097152, uncompressed128cubed; failures preserved; no intensity/label/model inspection','selection_only_geometry_access':True,'thread_worker':'exact576b562 worker','gate':'each of5consecutive slices has>=25 observed points within19.2um; inherited frame-validity unchanged','maximum_regions_reconstructed':1,'registration_uncertainty':'physical error remains uncalibrated','stop':'No T13; no second reconstructed region even if FAIL'}
 dump(O/'selection_policy.json',policy)
 dump(O/'SELECTION_POLICY_FROZEN.json',{str(p.relative_to(R)):sha(p) for p in [O/'selection_policy.json',Path(__file__).resolve()]})
 points=np.stack([tifffile.imread(R/f'data/raw/p0_t10/geometry_{a}.tif') for a in 'xyz'],-1).astype(float);H,W=points.shape[:2];candidates=[(H//4,W//2),(3*H//4,W//2),(H//2,W//4),(H//2,3*W//4)];audit=[];selected=None
 for cy,cx in candidates:
  y0,y1=cy-51,cy+52;x0,x1=cx-52,cx+53;reg=points[y0-1:y1+1,x0-1:x1+1];P=reg[1:-1,1:-1];du=(reg[1:-1,2:]-reg[1:-1,:-2])/2;dv=(reg[2:,1:-1]-reg[:-2,1:-1])/2;N=unit(np.cross(du,dv));V=unit(dv);edges=np.r_[np.linalg.norm(np.diff(reg,axis=0),axis=-1).ravel(),np.linalg.norm(np.diff(reg,axis=1),axis=-1).ravel()];dots=np.r_[(N[1:]*N[:-1]).sum(-1).ravel(),(N[:,1:]*N[:,:-1]).sum(-1).ravel()];normaljump=float(np.degrees(np.arccos(np.clip(dots,-1,1))).max());center=int(np.floor(np.median(P[:,:,2])+.5));zs=list(range(center-3,center+4));low=(np.floor((P.min((0,1))[:2]-128)/128)*128).astype(int);high=(np.ceil((P.max((0,1))[:2]+128)/128)*128).astype(int)
  checks={'positive_coordinates':bool((reg>0).all()),'edge_lengths_bounded':bool(edges.min()>=.01 and edges.max()<=200),'normal_continuity':normaljump<=45,'separate_from_old_range':zs[-1]<44449-1024 or zs[0]>44463+1024};sections=[]
  if all(checks.values()):
   for z in zs:
    seg,_,_=section(z,P,N,V)
    if len(seg):
     # Shared triangle intersections form one connected polyline if graph has
     # exactly two degree1 endpoints, all others degree2, and one component.
     keys=[tuple(np.round(p,6)) for p in seg.reshape(-1,3)];adj={}
     for a,b in zip(keys[::2],keys[1::2]):adj.setdefault(a,set()).add(b);adj.setdefault(b,set()).add(a)
     seen=set();todo=[next(iter(adj))]
     while todo:
      a=todo.pop()
      if a in seen:continue
      seen.add(a);todo.extend(adj[a]-seen)
     degrees=[len(v) for v in adj.values()];continuous=len(seen)==len(adj) and degrees.count(1)==2 and all(d in [1,2] for d in degrees);arc=float(np.linalg.norm(seg[:,1]-seg[:,0],axis=1).sum());xy=seg[:,:,:2];margin=float(min((xy-low).min(),(high-xy).min()))
    else:continuous=False;arc=0.;margin=0.
    sections.append({'z':z,'segments':len(seg),'single_connected_section':continuous,'arc_native':arc,'crop_margin_native':margin})
   checks['all7_continuous_usable_interior_sections']=all(s['single_connected_section'] and s['arc_native']>=500 and s['crop_margin_native']>=64 for s in sections)
  record={'grid_bounds_yx':[y0,y1,x0,x1],'native_z':zs,'native_xy_bounds_half_open':[low.tolist(),high.tolist()],'edge_length_minmax_native':[float(edges.min()),float(edges.max())],'maximum_neighbor_normal_angle_deg':normaljump,'sections':sections,'checks':checks};audit.append(record);dump(O/'selection_audit.json',audit)
  if not all(checks.values()):continue
  jobs=[]
  for k,j,i in itertools.product(range(zs[0]//128,zs[-1]//128+1),range(low[1]//128,high[1]//128),range(low[0]//128,high[0]//128)):
   jobs.append({'chunk_zyx':[int(k),int(j),int(i)],'url':BASE+f'0/{k}/{j}/{i}'})
  def head(q):
   try:
    with urllib.request.urlopen(urllib.request.Request(q['url'],method='HEAD'),timeout=40) as r:h=dict(r.headers);size=int(r.headers['Content-Length'])
    return {**q,'bytes':size,'headers':h,'pass':size==2097152}
   except Exception as e:return {**q,'pass':False,'error':str(e)}
  with concurrent.futures.ThreadPoolExecutor(8) as pool:coverage=list(pool.map(head,jobs))
  record['coverage_head']=coverage;record['checks']['full_CT_coverage']=all(q['pass'] for q in coverage);dump(O/'selection_audit.json',audit)
  if not record['checks']['full_CT_coverage']:continue
  selected={'source':BASE,'native_z_indices':zs,'native_xy_bounds_half_open':[low.tolist(),high.tolist()],'trusted_grid_bounds_yx':[y0,y1,x0,x1],'trusted_grid_derivative_halo_bounds_yx':[y0-1,y1+1,x0-1,x1+1],'native_spacing_um':2.4,'xy_sampling_um':9.6,'arc_spacing_um':4.8,'coordinate_transform':'nativeXY=lowXY+4*localXY+1.5; nativeZ unchanged','selection_candidate_rank':len(audit),'geometry_checks':record['checks'],'selection_used_thread_overlap_labels_models':False,'source_worker_sha256':sha(R/'scripts/p0_t13_control_reconstruct.py'),'gate':{'distance_um':19.2,'minimum_each_slice_count':25,'frame_projection_min_norm':.25},'window_selection':'earliest passing; otherwise largest minimum count, smallest summed deficit, greatest total, earliest start','acquisition_cap_bytes':128000000,'registration_error_um':None,'maximum_reconstructed_regions':1};break
 assert selected is not None,'No neutral candidate qualified; do not trace or silently relax criteria'
 assert sha(R/'scripts/p0_t13_control_reconstruct.py')==hashlib.sha256(subprocess.check_output(['git','show','576b562:scripts/p0_t13_control_reconstruct.py'],cwd=R)).hexdigest()
 dump(O/'selected_region.json',selected)
 sources=list((R/'scripts').glob('p0_t13_control_alt*.py'))+[R/'scripts/p0_t13_control_reconstruct.py',R/'scripts/p0_t13_control_evaluate.py',O/'selected_region.json',O/'selection_policy.json']
 dump(O/'PRE_RECONSTRUCTION_FREEZE.json',{str(p.relative_to(R)):sha(p) for p in sources});print(json.dumps(selected,indent=2),flush=True)
if __name__=='__main__':main()
