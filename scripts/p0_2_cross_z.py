"""Independent 2-D ridge fragments and bounded cross-Z associations; no surface mesh."""
import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d, map_coordinates
from scipy.signal import find_peaks
from scipy.spatial.distance import cdist
from skimage.morphology import skeletonize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
ZS=list(range(137,152))
ROI=[160,296,100,174]  # x0,x1,y0,y1, half-open; includes all previous endpoints
PLAN={'hypothesis':'Neighboring Z supports fragment identity and may disambiguate single-slice continuations.',
 'scope':'15 existing real CT XY slices, independent extraction, no in-plane joins, no mesh or final stacked raster',
 'z_inclusive':[137,151],'z_center':144,'source':'outputs/orientation/ct_level4_zyx.npy',
 'local_slab_z_origin':128,'roi_xy_half_open':ROI,'spacing_um':126.56,
 'extraction':{'DoG_sigma':[.6,2.],'threshold_CT_units':3.,'min_vertices':7,'min_length_px':6.,
               'branch_rule':'Remove degree>=3 skeleton nodes; no closing or gap filling; split into remaining paths'},
 'registration':{'max_shift_px':2.,'grid_step_px':.5,'objective':'NCC of sigma-1 image at fixed ROI interior'},
 'association':{'candidate_distance_px':4.,'near_distance_px':1.5,'score_weights':
  {'proximity':.2,'overlap':.25,'texture':.3,'orientation':.1,'curvature':.05,'spacing':.05,'radial_order':.05},
  'accept_score':.72,'accept_margin':.06,'accept_median_distance_px':1.25,'accept_overlap':.5,
  'accept_angle_deg':25.,'accept_NCC':.65,'normal_spacing_fraction':.45,'radial_spacing_fraction':.5,
  'hard_reject':{'distance_px':2.5,'overlap':.25,'angle_deg':45.,'NCC':.15}},
 'gap_reappearance':{'max_missing_slices':1,'min_score':.78,'min_NCC':.75,'min_overlap':.65,'margin':.10},
 'controls':{'ids':['F02','F03','F04','F05','F06'],'anchor':'middle CSV point of each original fragment',
             'success':'>=3 of 5 controls attached uniquely within 2 px and tracked across >=7 observed slices'},
 'endpoint_support':{'source_attachment_px':2.,'candidate_distance_px':1.5,'source_distance_px':2.,
                     'min_consecutive_neighbor_votes':3,'min_vote_fraction':.6,'min_runner_up_vote_margin':2},
 'status':'Prospective parameters for an exploratory pilot; not an independently reviewed preregistration',
 'failure_rule':'Preserve failures; if controls fail, diagnose gates/extraction before claiming any success.'}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(out,name,data):(out/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def sample(image,xy):
    p=np.asarray(xy);return map_coordinates(image,[p[...,1].ravel(),p[...,0].ravel()],order=1,mode='nearest').reshape(p.shape[:-1])
def ncc(a,b):
    a=np.asarray(a).ravel();b=np.asarray(b).ravel();a=a-a.mean();b=b-b.mean()
    return float(np.clip(a@b/max(np.linalg.norm(a)*np.linalg.norm(b),1e-9),-1,1))
def unit(a):return a/np.maximum(np.linalg.norm(a,axis=-1,keepdims=True),1e-9)


def features(points,image,response):
    p=np.asarray(points,float)
    smooth=gaussian_filter1d(p,1.,axis=0,mode='nearest')
    tangent=unit(np.gradient(smooth,axis=0))
    angles=np.unwrap(np.arctan2(tangent[:,1],tangent[:,0]))
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
    kappa=np.gradient(angles)/np.maximum(np.gradient(arc),.1)
    normal=np.c_[-tangent[:,1],tangent[:,0]]
    mid=len(p)//2;off=np.arange(-8,8.01,.25)
    profile=sample(response,p[mid]+off[:,None]*normal[mid])
    peaks,_=find_peaks(profile,prominence=2,distance=4);peaks=peaks[profile[peaks]>1]
    spacing=float(np.median(np.diff(off[peaks]))) if len(peaks)>1 else None
    texture=sample(image,p[:,None,:]+np.arange(-2,3)[None,:,None]*normal[:,None,:])
    return {'xy':p.tolist(),'tangent':tangent.tolist(),'curvature':kappa.tolist(),
            'spacing_px':spacing,'texture':texture.tolist(),'centroid':p.mean(axis=0).tolist(),
            'arc_length_px':float(arc[-1]),'mean_ridge_response':float(sample(response,p).mean())}


def extract(image,z):
    x0,x1,y0,y1=ROI
    response=gaussian_filter(image,.6)-gaussian_filter(image,2.)
    sk=skeletonize(response[y0:y1,x0:x1]>3.)
    pixels={tuple(p) for p in np.argwhere(sk)}
    neighbors={}
    for p in sorted(pixels):
        ns=[];y,x=p
        for dy in [-1,0,1]:
            for dx in [-1,0,1]:
                q=(y+dy,x+dx)
                if (dy==0 and dx==0) or q not in pixels:continue
                # Avoid diagonal shortcuts across existing orthogonal skeleton edges.
                if dy and dx and ((y+dy,x) in pixels or (y,x+dx) in pixels):continue
                ns.append(q)
        neighbors[p]=sorted(ns)
    branch={p for p in pixels if len(neighbors[p])>=3}
    remain=pixels-branch
    graph={p:[q for q in neighbors[p] if q in remain] for p in remain}
    seen=set();fragments=[]
    for seed in sorted(remain):
        if seed in seen:continue
        component=set();todo=[seed]
        while todo:
            p=todo.pop()
            if p in component:continue
            component.add(p);todo.extend(q for q in graph[p] if q not in component)
        seen.update(component)
        ends=sorted(p for p in component if len(graph[p])<=1)
        start=ends[0] if ends else min(component)
        ordered=[];visited=set();current=start
        while current not in visited:
            ordered.append(current);visited.add(current)
            nxt=[q for q in graph[current] if q not in visited]
            if not nxt:break
            current=nxt[0]
        assert len(ordered)==len(component)
        if len(ordered)<7:continue
        xy=np.array([[x+x0,y+y0] for y,x in ordered],float)
        if np.linalg.norm(np.diff(xy,axis=0),axis=1).sum()<6:continue
        descriptor=features(xy,image,response)
        fragments.append({'id':f'z{z}_f{len(fragments):03}','z':z,
                          'cycle_cut':not bool(ends),**descriptor})
    return fragments,{'z':z,'fragments':len(fragments),'skeleton_pixels':int(sk.sum()),
                      'removed_branch_nodes':len(branch),'ridge_mask_pixels':int(np.sum(response[y0:y1,x0:x1]>3))}


def registration(a,b):
    x0,x1,y0,y1=ROI
    yy,xx=np.mgrid[y0+4:y1-4,x0+4:x1-4]
    points=np.stack([xx,yy],axis=-1)
    aa=gaussian_filter(a,1.);bb=gaussian_filter(b,1.)
    source=sample(aa,points)
    choices=[]
    for dy in np.arange(-2,2.01,.5):
        for dx in np.arange(-2,2.01,.5):
            choices.append((ncc(source,sample(bb,points+[dx,dy])),float(dx),float(dy)))
    choices.sort(key=lambda r:(-r[0],r[1]**2+r[2]**2,r[1],r[2]))
    best=choices[0]
    return {'shift_xy':[best[1],best[2]],'NCC':best[0],'at_bound':abs(best[1])==2 or abs(best[2])==2}


def pair_evidence(a,b,shift):
    p=np.asarray(a['xy']);q=np.asarray(b['xy']);ps=p+shift
    if np.any(ps.min(axis=0)>q.max(axis=0)+4) or np.any(q.min(axis=0)>ps.max(axis=0)+4):return None
    distances=cdist(ps,q)
    near=distances.argmin(axis=1);dp=distances.min(axis=1);dq=distances.min(axis=0)
    if min(dp.min(),dq.min())>4:return None
    med=float((np.median(dp)+np.median(dq))/2)
    overlap=float(min(np.mean(dp<=1.5),np.mean(dq<=1.5)))
    tp=np.array(a['tangent']);tq=np.array(b['tangent'])[near]
    dot=np.sum(tp*tq,axis=1)
    angle=float(np.degrees(np.arccos(np.clip(np.abs(dot),0,1))).mean())
    texb=np.asarray(b['texture'])[near].copy();texb[dot<0]=texb[dot<0,::-1]
    corr=ncc(np.asarray(a['texture']),texb)
    delta=q[near]-ps
    normal=np.c_[-tp[:,1],tp[:,0]]
    normal_error=float(np.median(np.abs(np.sum(delta*normal,axis=1))))
    radial_error=float(np.median(np.abs(np.linalg.norm(q[near]-[240,140],axis=1)-np.linalg.norm(ps-[240,140],axis=1))))
    curvature=float(np.median(np.abs(np.abs(a['curvature'])-np.abs(np.array(b['curvature'])[near]))))
    la,lb=a['spacing_px'],b['spacing_px'];spacing=min(la,lb) if la and lb else None
    components={'proximity':float(np.exp(-med/1.5)),'overlap':overlap,'texture':max(0,corr),
       'orientation':float(np.mean(np.abs(dot))),'curvature':float(np.exp(-curvature/.15)),
       'spacing':float(np.exp(-abs(np.log(la/lb)))) if la and lb else .5,
       'radial_order':float(np.exp(-radial_error/max(spacing or 3.,1)))}
    score=sum(components[k]*v for k,v in PLAN['association']['score_weights'].items())
    hard=[]
    for condition,label in [(med>2.5,'distance'),(overlap<.25,'overlap'),(angle>45,'orientation'),(corr<.15,'texture')]:
        if condition:hard.append(label)
    gates=[]
    for condition,label in [(score<.72,'score'),(med>1.25,'distance'),(overlap<.5,'overlap'),
       (angle>25,'orientation'),(corr<.65,'texture'),
       (spacing is None,'spacing_unknown'),(spacing is not None and normal_error>.45*spacing,'normal_layer_risk'),
       (spacing is not None and radial_error>.5*spacing,'radial_layer_risk')]:
        if condition:gates.append(label)
    return {'source':a['id'],'target':b['id'],'z_source':a['z'],'z_target':b['z'],
      'score':float(score),'components':components,'median_distance_px':med,'overlap':overlap,
      'angle_deg':angle,'NCC':corr,'normal_error_px':normal_error,'radial_error_px':radial_error,
      'spacing_px':spacing,'curvature_difference':curvature,'hard_reject':hard,
      'gate_failures':gates,'status':'REJECTED' if hard else 'UNKNOWN','shift_xy':list(shift)}


def associations(left,right,shift,gap=False):
    records=[]
    for a in left:
        for b in right:
            record=pair_evidence(a,b,shift)
            if record is not None:records.append(record)
    by_source=defaultdict(list);by_target=defaultdict(list)
    for r in records:
        if r['hard_reject']:continue
        by_source[r['source']].append(r);by_target[r['target']].append(r)
    for group in [by_source,by_target]:
        for key in group:group[key].sort(key=lambda r:(-r['score'],r['source'],r['target']))
    for r in records:
        if r['hard_reject']:continue
        row=by_source[r['source']];col=by_target[r['target']]
        r['source_margin']=r['score']-max((t['score'] for t in row if t is not r),default=0.)
        r['target_margin']=r['score']-max((t['score'] for t in col if t is not r),default=0.)
        r['mutual_best']=row[0] is r and col[0] is r
        extra=[]
        if not r['mutual_best']:extra.append('not_mutual_best')
        if min(r['source_margin'],r['target_margin'])<(.10 if gap else .06):extra.append('competitor_margin')
        if gap:
            if r['score']<.78 or r['NCC']<.75 or r['overlap']<.65:extra.append('gap_strict_evidence')
        r['decision_failures']=r['gate_failures']+extra
        if not r['decision_failures']:r['status']='ACCEPTED'
    return records


def tracks_from(nodes,edges):
    graph={n['id']:set() for n in nodes};lookup={n['id']:n for n in nodes}
    for e in edges:
        if e['status']=='ACCEPTED':graph[e['source']].add(e['target']);graph[e['target']].add(e['source'])
    seen=set();tracks=[];membership={}
    for seed in sorted(graph):
        if seed in seen:continue
        todo=[seed];component=set()
        while todo:
            p=todo.pop()
            if p in component:continue
            component.add(p);todo.extend(graph[p]-component)
        seen.update(component);ordered=sorted(component,key=lambda x:(lookup[x]['z'],x))
        zz=[lookup[x]['z'] for x in ordered];assert len(set(zz))==len(zz)
        track={'id':f'T{len(tracks):04}','nodes':ordered,'z_observed':zz,'length_slices':len(zz),
               'span_slices':max(zz)-min(zz)+1,'missing_z':[z for z in range(min(zz),max(zz)+1) if z not in zz]}
        tracks.append(track)
        for node in ordered:membership[node]=track['id']
    return tracks,membership


def nearest_fragment(point,fragments):
    choices=sorted((float(np.linalg.norm(np.array(f['xy'])-point,axis=1).min()),f['id']) for f in fragments)
    return choices


def evaluate(nodes,by_z,edges,tracks,membership,regs):
    lookup={n['id']:n for n in nodes};trackmap={t['id']:t for t in tracks}
    controls=[]
    for fid in PLAN['controls']['ids']:
        a=np.genfromtxt(ROOT/f'outputs/p0_1b/{fid}_coordinates.csv',delimiter=',',names=True)
        point=np.array([a['x'][len(a)//2],a['y'][len(a)//2]])
        choices=nearest_fragment(point,by_z[144]);best=choices[0]
        unique=best[0]<=2 and (len(choices)==1 or choices[1][0]-best[0]>=.5)
        track=trackmap[membership[best[1]]] if unique else None
        controls.append({'id':fid,'anchor_xy':point.tolist(),'attachment_distance_px':best[0],
          'attachment_unique':unique,'fragment':best[1] if unique else None,
          'track':track['id'] if track else None,'length_slices':track['length_slices'] if track else 0,
          'z_observed':track['z_observed'] if track else [],'coverage_gate_pass':bool(track and track['length_slices']>=7)})
    offsets={144:np.zeros(2)}
    for z in range(144,151):offsets[z+1]=offsets[z]+regs[z]['shift_xy']
    for z in range(143,136,-1):offsets[z]=offsets[z+1]-regs[z]['shift_xy']
    endpoint_records=[]
    original=json.loads((ROOT/'outputs/p0_1b/endpoints.json').read_text())
    for e in original:
        point=np.array(e['xy']);choices=nearest_fragment(point,by_z[144]);best=choices[0]
        attached=best[0]<=2 and (len(choices)==1 or choices[1][0]-best[0]>=.5)
        track=trackmap[membership[best[1]]] if attached else None
        votes=[]
        for c in e['candidates']:
            yes=[];observed=[]
            for node_id in track['nodes'] if track else []:
                node=lookup[node_id];z=node['z'];q=np.array(node['xy']);shift=offsets[z]
                ds=np.linalg.norm(q-point-shift,axis=1);dc=np.linalg.norm(q-np.array(c['xy'])-shift,axis=1)
                if ds.min()>2:continue
                observed.append(z)
                ia,ib=int(ds.argmin()),int(dc.argmin())
                lo,hi=sorted([ia,ib]);arc=np.linalg.norm(np.diff(q[lo:hi+1],axis=0),axis=1).sum()
                # A same-fragment path must be local, not a circuit around a winding.
                if dc.min()<=1.5 and arc<=1.8*np.linalg.norm(np.array(c['xy'])-point)+3:yes.append(z)
            neighbor_yes=[z for z in yes if z!=144]
            longest=run=0;last=-999
            for z in neighbor_yes:
                run=run+1 if z==last+1 else 1;longest=max(longest,run);last=z
            votes.append({'candidate':c['id'],'xy':c['xy'],'z_support':yes,'observed_z':observed,
               'neighbor_support':len(neighbor_yes),'consecutive_neighbor_support':longest,
               'vote_fraction':len(yes)/len(observed) if observed else 0.,'single_slice_supported':144 in yes})
        votes.sort(key=lambda v:(-v['neighbor_support'],-v['vote_fraction'],v['candidate']))
        win=votes[0];runner=votes[1]['neighbor_support'] if len(votes)>1 else 0
        qualifies=win['consecutive_neighbor_support']>=3 and win['vote_fraction']>=.6 and win['neighbor_support']-runner>=2
        decision='SUPPORTED_CROSS_Z' if qualifies and not win['single_slice_supported'] else ('ALREADY_SUPPORTED_SINGLE_SLICE' if qualifies else 'UNKNOWN')
        endpoint_records.append({'endpoint':e['id'],'source_attachment_unique':attached,
         'source_distance_px':best[0],'track':track['id'] if track else None,'candidates':votes,
         'decision':decision,'winner':win['candidate'] if qualifies else None,'physical_oracle_correct':None})
    return controls,endpoint_records


def compute():
    volume=np.load(ROOT/PLAN['source'],mmap_mode='r')
    images={z:np.asarray(volume[z-128],dtype=float) for z in ZS}
    assert np.array_equal(images[144],np.asarray(__import__('PIL.Image',fromlist=['Image']).open(ROOT/'outputs/p0_1/source_slice.png')))
    by_z={};counts=[]
    for z in ZS:by_z[z],count=extract(images[z],z);counts.append(count)
    nodes=[n for z in ZS for n in by_z[z]]
    regs={z:registration(images[z],images[z+1]) for z in ZS[:-1]}
    adjacent=[]
    for z in ZS[:-1]:adjacent.extend(associations(by_z[z],by_z[z+1],np.array(regs[z]['shift_xy'])))
    outgoing={e['source'] for e in adjacent if e['status']=='ACCEPTED'}
    incoming={e['target'] for e in adjacent if e['status']=='ACCEPTED'}
    gap=[]
    for z in ZS[:-2]:
        left=[n for n in by_z[z] if n['id'] not in outgoing]
        right=[n for n in by_z[z+2] if n['id'] not in incoming]
        shift=np.array(regs[z]['shift_xy'])+regs[z+1]['shift_xy']
        gap.extend(associations(left,right,shift,gap=True))
    edges=adjacent+gap
    tracks,membership=tracks_from(nodes,edges)
    controls,endpoints=evaluate(nodes,by_z,edges,tracks,membership,regs)
    summary={'z_inclusive':[137,151],'slices':15,'fragment_counts':counts,'total_fragments':len(nodes),
      'adjacent_associations':dict(Counter(e['status'] for e in adjacent)),
      'gap_associations':dict(Counter(e['status'] for e in gap)),
      'association_total':len(edges),'tracks':len(tracks),
      'track_length_histogram':dict(sorted(Counter(t['length_slices'] for t in tracks).items())),
      'tracks_with_missing_slice':sum(bool(t['missing_z']) for t in tracks),
      'control_pass_count':sum(c['coverage_gate_pass'] for c in controls),'control_count':len(controls),
      'control_gate_met':sum(c['coverage_gate_pass'] for c in controls)>=3,
      'endpoint_decisions':dict(Counter(e['decision'] for e in endpoints)),
      'true_false_association_rate':None,'true_identity_oracle':'Unavailable; indicators are not ground truth'}
    return images,nodes,regs,edges,tracks,controls,endpoints,summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='outputs/p0_2/run_01');parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=ROOT/args.out
    if not args.verify:
        out.mkdir(parents=True,exist_ok=False)
        prior=[p for r in ['outputs/orientation','outputs/p0_1','outputs/p0_1b','docs','scripts','data/raw/p0_orientation'] for p in (ROOT/r).rglob('*') if p.is_file()]+[ROOT/'README.md',ROOT/'AGENTS.md']
        save(out,'preservation_before.json',{str(p.relative_to(ROOT)):sha(p) for p in prior})
        save(out,'plan.json',PLAN)
        save(out,'environment.json',{p:importlib.metadata.version(p) for p in ['numpy','scipy','scikit-image','matplotlib']})
        save(out,'source_identity.json',{'script_sha256':sha(Path(__file__)),'source_sha256':sha(ROOT/PLAN['source'])})
    images,nodes,regs,edges,tracks,controls,endpoints,summary=compute()
    artifacts={'fragments.json':nodes,'registration.json':{str(k):v for k,v in regs.items()},'associations.json':edges,
       'tracks.json':tracks,'controls.json':controls,'endpoint_reconsideration.json':endpoints,'metrics.json':summary}
    if args.verify:
        for name,data in artifacts.items():assert json.loads((out/name).read_text())==json.loads(json.dumps(data)),name
        for name,digest in json.loads((out/'preservation_before.json').read_text()).items():assert sha(ROOT/name)==digest,name
        save(out,'verification.json',{'numerical_replay':'exact','prior_hashes':'unchanged','source_scope':'15 XY slices only'})
        print('PASS: exact replay and preservation');return
    for name,data in artifacts.items():save(out,name,data)
    look={n['id']:n for n in nodes};membership={n:t for t in tracks for n in t['nodes']}
    palette=plt.get_cmap('tab20')
    for z in [137,140,144,148,151]:
        fig,ax=plt.subplots(figsize=(14,8));ax.imshow(images[z],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
        for node in [n for n in nodes if n['z']==z]:
            track=membership[node['id']];xy=np.array(node['xy'])
            if track['length_slices']>=3:
                ax.plot(*xy.T,lw=1,color=palette(int(track['id'][1:])%20))
            else:ax.plot(*xy.T,lw=.5,color='red',alpha=.35)
        for c in controls:
            if c['track']:
                track=next(t for t in tracks if t['id']==c['track'])
                for n in track['nodes']:
                    if look[n]['z']==z:
                        xy=np.array(look[n]['xy']);ax.plot(*xy.T,'cyan',lw=2)
                        ax.text(*xy[len(xy)//2],c['id'],color='yellow',fontsize=9)
        ax.set(xlim=(ROI[0],ROI[1]),ylim=(ROI[3],ROI[2]),title=f'Z={z}: independent fragments; color = provisional cross-Z track, cyan = controls')
        fig.tight_layout();fig.savefig(out/f'overlay_z{z}.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,4))
    for i,c in enumerate(controls):ax.plot(c['z_observed'],[i]*len(c['z_observed']),'o-',label=c['id'])
    ax.set(xlabel='Level-4 Z index',yticks=range(len(controls)),yticklabels=[c['id'] for c in controls],title='Control track coverage (not an unwrapped image)');ax.set_xticks(ZS);ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(out/'control_tracks.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2));print(json.dumps(controls,indent=2))


if __name__=='__main__':main()
