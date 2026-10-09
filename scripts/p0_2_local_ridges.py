"""P0.2 engineering revision: independent local ridge fragments, bounded cross-Z matching."""
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_2_cross_z import ROOT,ZS,ROI,sha,save,sample,ncc,unit,tracks_from

PLAN={'revision':'run_02; designed after preserved run_01 control failure, not independent confirmation',
 'z_inclusive':[137,151],'center_z':144,'spacing_um':126.56,
 'source':'outputs/orientation/ct_level4_zyx.npy','extraction_independence':
 'Each XY independently supplies all local peaks. Fixed inspection windows come from P0.1b; no preceding-Z initialization.',
 'windows':{'controls':5,'endpoints':14,'u_control':[-4,4,.5],'u_endpoint':[-4,10,.5],'v':[-8,8,.25]},
 'ridge':{'sigma':[.6,2.],'peak_prominence':1.,'min_response':1.,'peak_separation_px':1.,
          'neighbor_u_max_v_change':.75,'neighbor_u_competitor_margin':.125,'min_points':7,'min_arc_px':3.},
 'association':{'min_common_stations':5,'accept_normal_residual_px':1.,'accept_NCC':.55,
   'accept_overlap':.6,'accept_angle_deg':25.,'accept_score':.70,'margin':.06,
   'layer_spacing_fraction':.45,'weights':{'position':.2,'overlap':.15,'texture':.35,'orientation':.1,
                                        'curvature':.05,'spacing':.05,'radial_order':.1}},
 'gap':{'missing_slices_max':1,'score':.76,'NCC':.7,'margin':.10},
 'controls':'Same F02,F03,F04,F05,F06 midpoints, >=3/5 attached within 2 px and observed in >=7 slices',
 'endpoint_resolution':'Unique candidate with >=3 consecutive neighboring support slices, >=0.6 supported fraction, >=2 vote lead; compare central-slice baseline',
 'limits':'Window-local proposals may overlap; counts are not unique anatomical sheets. No in-plane gap filling or mesh.'}


def windows():
    result=[]
    for fid in ['F02','F03','F04','F05','F06']:
        a=np.genfromtxt(ROOT/f'outputs/p0_1b/{fid}_coordinates.csv',delimiter=',',names=True)
        xy=np.c_[a['x'],a['y']];mid=len(xy)//2
        t=unit(xy[min(mid+5,len(xy)-1)]-xy[max(mid-5,0)])
        result.append({'id':'control_'+fid,'kind':'control','origin':xy[mid].tolist(),
                       'tangent':t.tolist(),'u':np.arange(-4,4.01,.5).tolist(),'reference':fid})
    for e in json.loads((ROOT/'outputs/p0_1b/endpoints.json').read_text()):
        result.append({'id':e['id'],'kind':'endpoint','origin':e['xy'],
                       'tangent':e['outward_tangent'],'u':np.arange(-4,10.01,.5).tolist(),'reference':e['id']})
    return result


def extract_window(image,z,w):
    response=gaussian_filter(image,.6)-gaussian_filter(image,2.)
    origin=np.array(w['origin']);t=np.array(w['tangent']);normal=np.array([-t[1],t[0]])
    us=np.array(w['u']);vs=np.arange(-8,8.01,.25)
    grid=origin+us[:,None,None]*t+vs[None,:,None]*normal
    values=sample(response,grid)
    peaks=[];all_spacing=[]
    for i,profile in enumerate(values):
        idx,_=find_peaks(profile,prominence=1.,distance=4)
        idx=idx[profile[idx]>1.]
        nodes=[]
        for k in idx:
            denom=profile[k-1]-2*profile[k]+profile[k+1]
            correction=np.clip(.5*(profile[k-1]-profile[k+1])/denom,-.5,.5) if abs(denom)>1e-9 else 0.
            nodes.append((i,float(vs[k]+correction*.25)))
        peaks.append(nodes)
        all_spacing.extend(np.diff([p[1] for p in nodes]).tolist())
    local_spacing=float(np.median(all_spacing)) if all_spacing else None
    graph={p:[] for row in peaks for p in row}
    for a,b in zip(peaks[:-1],peaks[1:]):
        if not a or not b:continue
        d=np.abs(np.array([p[1] for p in a])[:,None]-np.array([p[1] for p in b]))
        for i,p in enumerate(a):
            j=int(d[i].argmin())
            if d[:,j].argmin()!=i or d[i,j]>.75:continue
            row=np.sort(d[i]);col=np.sort(d[:,j])
            if (len(row)>1 and row[1]-row[0]<.125) or (len(col)>1 and col[1]-col[0]<.125):continue
            graph[p].append(b[j]);graph[b[j]].append(p)
    seen=set();fragments=[]
    for seed in sorted(graph):
        if seed in seen:continue
        comp=set();todo=[seed]
        while todo:
            p=todo.pop()
            if p in comp:continue
            comp.add(p);todo.extend(q for q in graph[p] if q not in comp)
        seen.update(comp)
        if len(comp)<7:continue
        ordered=sorted(comp);stations=np.array([p[0] for p in ordered]);v=np.array([p[1] for p in ordered]);u=us[stations]
        assert len(set(stations))==len(stations)
        xy=origin+u[:,None]*t+v[:,None]*normal
        arc=np.linalg.norm(np.diff(xy,axis=0),axis=1).sum()
        if arc<3:continue
        slope=np.gradient(v,u);angle=np.arctan(slope);curvature=np.gradient(angle,u)/np.sqrt(1+slope**2)
        # Common window-normal profiles avoid reversing/rotating textured patches inconsistently.
        texture=sample(image,xy[:,None,:]+np.arange(-2,2.01,.5)[None,:,None]*normal)
        fragments.append({'id':f"z{z}_{w['id']}_r{len(fragments):02}",'z':z,'window':w['id'],
          'stations':stations.tolist(),'u':u.tolist(),'v':v.tolist(),'xy':xy.tolist(),
          'texture':texture.tolist(),'angle':angle.tolist(),'curvature':curvature.tolist(),
          'spacing_px':local_spacing,'arc_length_px':float(arc),'normal':normal.tolist()})
    return fragments


def evidence(a,b,shift):
    common,ia,ib=np.intersect1d(a['stations'],b['stations'],return_indices=True)
    if len(common)<3:return None
    normal=np.array(a['normal']);p=np.array(a['xy'])[ia];q=np.array(b['xy'])[ib]
    dv=np.array(b['v'])[ib]-np.array(a['v'])[ia]-np.array(shift)@normal
    distance=float(np.median(np.abs(dv)))
    if distance>4:return None
    overlap=float(len(common)/min(len(a['stations']),len(b['stations'])))
    corr=ncc(np.array(a['texture'])[ia],np.array(b['texture'])[ib])
    angle=float(np.degrees(np.abs(np.array(a['angle'])[ia]-np.array(b['angle'])[ib])).mean())
    dk=float(np.median(np.abs(np.array(a['curvature'])[ia]-np.array(b['curvature'])[ib])))
    la,lb=a['spacing_px'],b['spacing_px'];spacing=min(la,lb) if la and lb else None
    radial=float(np.median(np.abs(np.linalg.norm(q-[240,140],axis=1)-np.linalg.norm(p+shift-[240,140],axis=1))))
    components={'position':float(np.exp(-distance/1.5)),'overlap':overlap,'texture':max(corr,0),
      'orientation':float(np.cos(np.radians(min(angle,90)))),'curvature':float(np.exp(-dk/.3)),
      'spacing':float(np.exp(-abs(np.log(la/lb)))) if la and lb else .5,
      'radial_order':float(np.exp(-radial/max(spacing or 3,1)))}
    score=float(sum(components[k]*v for k,v in PLAN['association']['weights'].items()))
    hard=[];fail=[]
    for cond,name in [(distance>2.5,'distance'),(corr<.1,'texture'),(angle>50,'orientation')]:
        if cond:hard.append(name)
    for cond,name in [(len(common)<5,'short_common_support'),(distance>1,'distance'),(corr<.55,'texture'),
                      (angle>25,'orientation'),(overlap<.6,'overlap'),(score<.70,'score'),
                      (spacing is not None and distance>.45*spacing,'layer_order_risk'),
                      (spacing is None and distance>.5,'spacing_unknown_and_displaced')]:
        if cond:fail.append(name)
    return {'source':a['id'],'target':b['id'],'z_source':a['z'],'z_target':b['z'],'window':a['window'],
      'score':score,'components':components,'common_stations':len(common),'normal_error_px':distance,
      'overlap':overlap,'NCC':corr,'angle_deg':angle,'curvature_difference':dk,'radial_error_px':radial,
      'spacing_px':spacing,'hard_reject':hard,'gate_failures':fail,'status':'REJECTED' if hard else 'UNKNOWN'}


def associate(left,right,shift,gap=False):
    records=[]
    for a in left:
        for b in right:
            if a['window']!=b['window']:continue
            r=evidence(a,b,shift)
            if r is not None:records.append(r)
    source=defaultdict(list);target=defaultdict(list)
    for r in records:
        if not r['hard_reject']:source[r['source']].append(r);target[r['target']].append(r)
    for group in [source,target]:
        for k in group:group[k].sort(key=lambda r:(-r['score'],r['source'],r['target']))
    for r in records:
        if r['hard_reject']:continue
        row,col=source[r['source']],target[r['target']]
        mutual=row[0] is r and col[0] is r
        margin=min(r['score']-max((q['score'] for q in row if q is not r),default=0),
                   r['score']-max((q['score'] for q in col if q is not r),default=0))
        r['mutual_best']=mutual;r['minimum_margin']=margin
        fail=list(r['gate_failures'])
        if not mutual:fail.append('not_mutual_best')
        if margin<(.10 if gap else .06):fail.append('competitor_margin')
        if gap and (r['score']<.76 or r['NCC']<.7):fail.append('gap_strict_evidence')
        r['decision_failures']=fail
        if not fail:r['status']='ACCEPTED'
    return records


def attach(point,nodes):
    options=sorted((float(np.linalg.norm(np.array(n['xy'])-point,axis=1).min()),n['id']) for n in nodes)
    if not options:return None,None
    best=options[0];unique=best[0]<=2 and (len(options)==1 or options[1][0]-best[0]>=.5)
    return best[1] if unique else None,best[0]


def evaluate(nodes,by_z,tracks,membership,ws,regs):
    look={n['id']:n for n in nodes};tm={t['id']:t for t in tracks}
    controls=[];endpoints=[]
    offsets={144:np.zeros(2)}
    for z in range(144,151):offsets[z+1]=offsets[z]+regs[str(z)]['shift_xy']
    for z in range(143,136,-1):offsets[z]=offsets[z+1]-regs[str(z)]['shift_xy']
    old={e['id']:e for e in json.loads((ROOT/'outputs/p0_1b/endpoints.json').read_text())}
    for w in ws:
        point=np.array(w['origin']);nid,distance=attach(point,[n for n in by_z[144] if n['window']==w['id']])
        track=tm[membership[nid]] if nid else None
        if w['kind']=='control':
            controls.append({'id':w['reference'],'anchor_xy':w['origin'],'attachment':nid,'distance_px':distance,
               'track':track['id'] if track else None,'length_slices':track['length_slices'] if track else 0,
               'z_observed':track['z_observed'] if track else [],
               'coverage_gate_pass':bool(track and track['length_slices']>=7)})
            continue
        records=[]
        for c in old[w['id']]['candidates']:
            yes=[];observed=[]
            for node_id in track['nodes'] if track else []:
                node=look[node_id];z=node['z'];xy=np.array(node['xy']);shift=offsets[z]
                d0=np.linalg.norm(xy-point-shift,axis=1);d1=np.linalg.norm(xy-np.array(c['xy'])-shift,axis=1)
                if d0.min()>2:continue
                observed.append(z)
                ia,ib=int(d0.argmin()),int(d1.argmin());lo,hi=sorted([ia,ib])
                length=np.linalg.norm(np.diff(xy[lo:hi+1],axis=0),axis=1).sum()
                if d1.min()<=1.5 and length<=1.8*np.linalg.norm(np.array(c['xy'])-point)+3:yes.append(z)
            neighbor=[z for z in yes if z!=144];last=-999;run=longest=0
            for z in neighbor:run=run+1 if z==last+1 else 1;longest=max(longest,run);last=z
            records.append({'candidate':c['id'],'xy':c['xy'],'observed_z':observed,'z_support':yes,
              'neighbor_support':len(neighbor),'longest_consecutive_neighbors':longest,
              'support_fraction':len(yes)/len(observed) if observed else 0,'single_slice_supported':144 in yes})
        records.sort(key=lambda r:(-r['neighbor_support'],-r['support_fraction'],r['candidate']))
        win=records[0];runner=records[1]['neighbor_support'] if len(records)>1 else 0
        qualifies=win['longest_consecutive_neighbors']>=3 and win['support_fraction']>=.6 and win['neighbor_support']-runner>=2
        endpoints.append({'endpoint':w['id'],'attachment':nid,'distance_px':distance,'track':track['id'] if track else None,
          'candidates':records,'winner':win['candidate'] if qualifies else None,
          'decision':('ALREADY_SUPPORTED_SINGLE_SLICE' if win['single_slice_supported'] else 'SUPPORTED_CROSS_Z') if qualifies else 'UNKNOWN',
          'physical_oracle_correct':None})
    return controls,endpoints


def compute():
    ws=windows();volume=np.load(ROOT/PLAN['source'],mmap_mode='r');images={z:np.asarray(volume[z-128],float) for z in ZS}
    by_z={z:[n for w in ws for n in extract_window(images[z],z,w)] for z in ZS}
    nodes=[n for z in ZS for n in by_z[z]]
    regs=json.loads((ROOT/'outputs/p0_2/run_01/registration.json').read_text())
    adjacent=[]
    for z in ZS[:-1]:adjacent.extend(associate(by_z[z],by_z[z+1],np.array(regs[str(z)]['shift_xy'])))
    outgoing={e['source'] for e in adjacent if e['status']=='ACCEPTED'};incoming={e['target'] for e in adjacent if e['status']=='ACCEPTED'}
    gap=[]
    for z in ZS[:-2]:
        shift=np.array(regs[str(z)]['shift_xy'])+regs[str(z+1)]['shift_xy']
        gap.extend(associate([n for n in by_z[z] if n['id'] not in outgoing],
          [n for n in by_z[z+2] if n['id'] not in incoming],shift,gap=True))
    edges=adjacent+gap;tracks,membership=tracks_from(nodes,edges)
    controls,endpoints=evaluate(nodes,by_z,tracks,membership,ws,regs)
    summary={'z_inclusive':[137,151],'slices':15,'total_window_local_fragments':len(nodes),
      'fragment_counts':[{'z':z,'all':len(by_z[z]),'controls':sum(n['window'].startswith('control') for n in by_z[z]),
                          'endpoint_windows':sum(not n['window'].startswith('control') for n in by_z[z])} for z in ZS],
      'adjacent_associations':dict(Counter(e['status'] for e in adjacent)),
      'gap_associations':dict(Counter(e['status'] for e in gap)),'association_total':len(edges),
      'track_length_histogram':dict(sorted(Counter(t['length_slices'] for t in tracks).items())),
      'tracks':len(tracks),'tracks_with_missing_slice':sum(bool(t['missing_z']) for t in tracks),
      'control_pass_count':sum(c['coverage_gate_pass'] for c in controls),'control_count':5,
      'control_gate_met':sum(c['coverage_gate_pass'] for c in controls)>=3,
      'endpoint_decisions':dict(Counter(e['decision'] for e in endpoints)),
      'true_false_association_rate':None,'oracle':'No registered physical layer identity in this region'}
    return images,ws,nodes,edges,tracks,controls,endpoints,summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=ROOT/'outputs/p0_2/run_02'
    if not args.verify:
        out.mkdir(parents=True,exist_ok=False);save(out,'plan.json',PLAN)
        prior=[p for r in ['outputs/orientation','outputs/p0_1','outputs/p0_1b','outputs/p0_2/run_01','docs','scripts','data/raw/p0_orientation'] for p in (ROOT/r).rglob('*') if p.is_file()]+[ROOT/'README.md',ROOT/'AGENTS.md']
        save(out,'preservation_before.json',{str(p.relative_to(ROOT)):sha(p) for p in prior})
        save(out,'script_identity.json',{'script_sha256':sha(Path(__file__))})
    images,ws,nodes,edges,tracks,controls,endpoints,summary=compute()
    artifacts={'windows.json':ws,'fragments.json':nodes,'associations.json':edges,'tracks.json':tracks,
               'controls.json':controls,'endpoint_reconsideration.json':endpoints,'metrics.json':summary}
    if args.verify:
        for name,data in artifacts.items():assert json.loads((out/name).read_text())==json.loads(json.dumps(data)),name
        for name,digest in json.loads((out/'preservation_before.json').read_text()).items():assert sha(ROOT/name)==digest,name
        save(out,'verification.json',{'replay':'exact','prior_files':'unchanged','extraction':'no neighboring-Z inputs'})
        print('PASS: exact replay and preservation');return
    for name,data in artifacts.items():save(out,name,data)
    look={n['id']:n for n in nodes}
    selected=[t for t in tracks if t['id'] in {c['track'] for c in controls}]
    palette=plt.get_cmap('tab10')
    for z in [137,140,144,148,151]:
        fig,ax=plt.subplots(figsize=(14,8));ax.imshow(images[z],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
        for node in [n for n in nodes if n['z']==z]:ax.plot(*np.array(node['xy']).T,color='cyan',alpha=.3,lw=.6)
        for i,c in enumerate(controls):
            if c['track']:
                track=next(t for t in selected if t['id']==c['track'])
                for nid in track['nodes']:
                    if look[nid]['z']==z:
                        q=np.array(look[nid]['xy']);ax.plot(*q.T,color=palette(i),lw=2)
                        ax.text(*q[len(q)//2],c['id'],color='yellow',fontsize=10)
        ax.set(xlim=(160,296),ylim=(174,100),title=f'Z={z}: independently extracted window-local ridges; labeled cross-Z controls')
        fig.tight_layout();fig.savefig(out/f'overlay_z{z}.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,4))
    for i,c in enumerate(controls):ax.plot(c['z_observed'],[i]*len(c['z_observed']),'o-',color=palette(i))
    ax.set(xlabel='Z level 4',yticks=range(5),yticklabels=[c['id'] for c in controls],title='Control track observations (no stacked CT raster)');ax.set_xticks(ZS);ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(out/'control_tracks.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2));print(json.dumps(controls,indent=2))


if __name__=='__main__':main()
