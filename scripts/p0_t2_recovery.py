"""Bounded prize-collecting beam search on independent 2-D CT fragment graphs."""
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree
from scipy.signal import find_peaks
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t2'
ZS=[128,144,146,160,170];COLORS=['black','#00bfff','#ff00ff','#ffb000']
PLAN={'hypothesis':'Maximum-recovery continuous inferred XY thread with local uncertainty and honest CT provenance',
 'slices':ZS,'source':'outputs/orientation/ct_level4_zyx.npy','z_origin':128,'spacing_um':126.56,'ds_px':.5,
 'input_graph':'unchanged P0-T1 maximal skeleton chains; retain arc >=3px plus seed chain',
 'seed':'exact P0-T1 seed, explore both directions on its parent chain',
 'exterior_target_xy':[430,460,150,210],
 'target_meaning':'right peripheral papyrus band, not outermost material boundary; fixed using existing CT context',
 'max_bridge_px':12.,'min_tangent_cos':-.5,'min_connector_cos':-.35,
 'beam_width':256,'max_chain_depth':180,'outgoing_cap':12,
 'objective':'observed arc - bridge penalty + 1.5*(initial target distance-current target distance)',
 'bridge_penalty':'2*d + 3*turn + 2*required_bend + 4*normal_miss/spacing + 3*crossing_proxy + 2*(1-gap_ridge_support)',
 'low_confidence':'turn cosine<.7, normal miss>.45spacing, crossing proxy>0, local penalty margin<1, or within 5px of T1 suspect zone',
 'constraints':'no repeated chain; no repeated observed pixel except immediate shared endpoint; no bridge crossing prior sampled geometry',
 'search':'diverse beam by endpoint 12px cell/orientation quadrant, then top score fill; top3 distinct complete paths if available',
 'provenance':{'0':'unresolved/no usable data','1':'observed CT-supported','2':'inferred bridge (always unsampled)','3':'low-confidence observed, CT sampled only if ridge-supported'},
 'scope':'exploratory bounded search, no claim of global optimum, no observed-percentage gate'}

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def arc(p):return np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
def unit(v):return v/np.maximum(np.linalg.norm(v,axis=-1,keepdims=True),1e-9)
def sample(a,p):return map_coordinates(a.astype(float),[p[...,1].ravel(),p[...,0].ravel()],order=1,mode='nearest').reshape(p.shape[:-1])
def target_dist(p):return np.linalg.norm(np.maximum([430,150]-p,0)+np.maximum(p-[460,210],0),axis=-1)
def codes(p):return frozenset((np.rint(p[:,1]).astype(int)*506+np.rint(p[:,0]).astype(int)).tolist())

def build(z):
    old=json.loads((ROOT/f'outputs/p0_t1/slice_{z}_graph.json').read_text())
    arrays=np.load(ROOT/f'outputs/p0_t1/slice_{z}_candidates.npz');response=arrays['response'];sk=arrays['skeleton']
    chains=[np.array(p,float) for p in old['chains']];sid=old['selected_chains'][0];seed=np.array(old['pieces'][0][0])
    ids=[i for i,p in enumerate(chains) if arc(p)[-1]>=3 or i==sid]
    ps=[];original=[]
    for i in ids:
        ps.extend([chains[i],chains[i][::-1]]);original.extend([i,i])
    starts=np.array([p[0] for p in ps]);ends=np.array([p[-1] for p in ps])
    tin=unit(np.array([p[min(5,len(p)-1)]-p[0] for p in ps]));tout=unit(np.array([p[-1]-p[max(0,len(p)-6)] for p in ps]))
    lengths=np.array([arc(p)[-1] for p in ps]);tree=cKDTree(starts)
    offsets=np.arange(-10,10.01,.5);normal=np.c_[-tout[:,1],tout[:,0]]
    profiles=sample(response,ends[:,None,:]+offsets[None,:,None]*normal[:,None,:])
    spacing=[]
    for profile in profiles:
        peaks,_=find_peaks(profile,prominence=2,distance=2)
        spacing.append(float(np.clip(np.median(np.diff(offsets[peaks])),2,10)) if len(peaks)>1 else 4.)
    adjacency=[];edge_rows=[]
    for i,p in enumerate(ps):
        js=np.array([j for j in sorted(tree.query_ball_point(ends[i],12)) if original[j]!=original[i]],int)
        if not len(js):adjacency.append([]);continue
        delta=starts[js]-ends[i];d=np.linalg.norm(delta,axis=1);v=unit(delta)
        tc=tin[js]@tout[i];ca=np.sum(v*tout[i],axis=1);cb=np.sum(v*tin[js],axis=1)
        ca[d==0]=1;cb[d==0]=1
        miss=np.abs(delta@normal[i]);turn=np.arccos(np.clip(tc,-1,1));bend=np.arccos(np.clip(ca,-1,1))+np.arccos(np.clip(cb,-1,1))
        mid=ends[i]+delta[:,None,:]*np.linspace(.2,.8,5)[None,:,None]
        evidence=np.mean(np.clip(sample(response,mid)/12,0,1),axis=1)
        hits=np.mean(sk[np.rint(mid[...,1]).astype(int),np.rint(mid[...,0]).astype(int)],axis=1)
        cross=hits*np.minimum(miss/spacing[i],3)*(d>2)
        penalty=2*d+3*turn+2*bend+4*miss/spacing[i]+3*cross+2*(1-evidence)*(d>0)
        admissible=(tc>=-.5)&(ca>=-.35)&(cb>=-.35)
        edges=[]
        for k,j in enumerate(js):
            e=[int(i),int(j),float(d[k]),float(tc[k]),float(ca[k]),float(cb[k]),float(miss[k]),spacing[i],
               float(turn[k]),float(bend[k]/max(d[k],1)),float(evidence[k]),float(cross[k]),float(penalty[k]),int(admissible[k])]
            edge_rows.append(e)
            if admissible[k]:edges.append((float(penalty[k]),int(j),len(edge_rows)-1))
        adjacency.append(sorted(edges))
    edge_array=np.array(edge_rows,float)
    seed_states=[]
    for orientation in [0,1]:
        idx=2*ids.index(sid)+orientation;p=ps[idx];j=int(np.argmin(np.linalg.norm(p-seed,axis=1)));part=p[j:]
        if len(part)>1:seed_states.append((idx,part))
    return ps,original,lengths,adjacency,edge_array,seed_states,response,sk,tout

def search(ps,original,lengths,adj,edges,seeds,tout):
    pixels=[codes(p) for p in ps];dtarget=target_dist(np.array([p[-1] for p in ps]));initial=float(target_dist(seeds[0][1][0]))
    # state: score, oriented chain sequence, edge IDs, visited original chains, occupied pixels, first piece
    beam=[(float(arc(p)[-1])+1.5*(initial-dtarget[i]),(i,),(),frozenset([original[i]]),codes(p),p) for i,p in seeds]
    complete=[];bestpartial=[];log=[];expanded=0;history_rejected=0
    for depth in range(180):
        new=[]
        for state in beam:
            score,path,eids,used,occupied,first=state;i=path[-1];endcode=next(iter(codes(ps[i][-1:])))
            for penalty,j,eid in adj[i][:12]:
                if original[j] in used:history_rejected+=1;continue
                overlap=pixels[j]&occupied
                if overlap-{endcode}:history_rejected+=1;continue
                e=edges[eid];d=e[2]
                bridge=ps[i][-1]+np.linspace(0,1,max(2,int(d*2)+1))[:,None]*(ps[j][0]-ps[i][-1])
                interior=bridge[(np.linalg.norm(bridge-bridge[0],axis=1)>1.5)&(np.linalg.norm(bridge-bridge[-1],axis=1)>1.5)]
                if codes(interior)&occupied:history_rejected+=1;continue
                nscore=score+lengths[j]-penalty+1.5*(dtarget[i]-dtarget[j])
                s=(float(nscore),path+(j,),eids+(eid,),used|{original[j]},occupied|pixels[j]|codes(bridge),first)
                expanded+=1
                if dtarget[j]==0:complete.append(s)
                else:new.append(s)
        if not new:break
        new.sort(key=lambda s:(-s[0],s[1]))
        # Coverage diversity keeps a long locally rewarding detour from monopolizing the beam.
        chosen=[];keys=set();chosen_paths=set()
        for s in new:
            i=s[1][-1];x,y=ps[i][-1];v=tout[i];key=(int(x//12),int(y//12),int((math.atan2(v[1],v[0])+math.pi)/(math.pi/2))%4)
            if key not in keys:
                chosen.append(s);keys.add(key);chosen_paths.add(s[1])
            if len(chosen)>=192:break
        for s in new:
            if len(chosen)>=256:break
            if s[1] not in chosen_paths:chosen.append(s);chosen_paths.add(s[1])
        beam=chosen;bestpartial=sorted(bestpartial+beam,key=lambda s:(-s[0],s[1]))[:50]
        log.append({'depth':depth+1,'expanded_total':expanded,'beam':len(beam),'complete_total':len(complete),
                    'closest_target_px':float(min(dtarget[s[1][-1]] for s in beam)),'best_score':beam[0][0]})
    pool=sorted(complete or bestpartial,key=lambda s:(-s[0],s[1]));top=[]
    for s in pool:
        # Report alternatives with different observed chain sets; exact paths also retained in complete archive.
        if all(len(s[3]&t[3])/len(s[3]|t[3])<.9 for t in top):top.append(s)
        if len(top)==3:break
    if len(top)<3:
        for s in pool:
            if all(s[1]!=t[1] for t in top):top.append(s)
            if len(top)==3:break
    return top,complete,log,{'expanded':expanded,'history_rejected':history_rejected,'complete_count':len(complete),'exterior_reached':bool(complete)}

def materialize(state,ps,edges,adj,response,z):
    score,path,eids,used,occupied,first=state;segments=[];labels=[];transition=[]
    suspect={146:[269,109],160:[300,112]}.get(z)
    for k,i in enumerate(path):
        p=first if k==0 else ps[i];low=False
        if k:
            e=edges[eids[k-1]];prev=path[k-1];gap=np.array([ps[prev][-1],p[0]])
            ranked=adj[prev];margin=(ranked[1][0]-ranked[0][0]) if len(ranked)>1 else 999.
            low=bool(e[3]<.7 or e[6]>.45*e[7] or e[11]>0 or margin<1)
            if suspect is not None and np.linalg.norm(gap[0]-suspect)<5:low=True
            transition.append({'edge_id':eids[k-1],'low_confidence':low,'local_margin':margin})
            if e[2]>0:segments.append(gap);labels.append(2)
        # Local uncertainty only: first 5px of this observed chain, not its downstream suffix.
        if low and arc(p)[-1]>5:
            a=arc(p);j=np.searchsorted(a,5);cut=np.array([np.interp(5,a,p[:,0]),np.interp(5,a,p[:,1])])
            segments.extend([np.vstack([p[:j],cut]),np.vstack([cut,p[j:]])]);labels.extend([3,1])
        else:segments.append(p);labels.append(3 if low else 1)
    lengths=np.array([arc(p)[-1] for p in segments]);starts=np.r_[0,np.cumsum(lengths)]
    ss=np.arange(0,starts[-1]+1e-8,.5);xy=np.zeros((len(ss),2));provenance=np.zeros(len(ss),np.uint8);ct=np.full(len(ss),np.nan)
    for k,(p,label) in enumerate(zip(segments,labels)):
        sel=(ss>=starts[k]-1e-9)&(ss<=starts[k+1]+1e-9);a=arc(p)
        xy[sel]=np.c_[np.interp(ss[sel]-starts[k],a,p[:,0]),np.interp(ss[sel]-starts[k],a,p[:,1])];provenance[sel]=label
    # Label0 means no usable ridge evidence even if graph construction supplied a polyline.
    support=sample(response,xy)>3
    provenance[np.isin(provenance,[1,3])&~support]=0
    return segments,labels,transition,ss,xy,provenance,lengths

def render(z,image,paths,ps,edges,adj,response):
    fig,axs=plt.subplots(1,len(paths),figsize=(6*len(paths),6),squeeze=False)
    for rank,(ax,state) in enumerate(zip(axs[0],paths)):
        segs,labs,*_=materialize(state,ps,edges,adj,response,z)
        ax.imshow(image,cmap='gray',vmin=0,vmax=255);ax.add_collection(LineCollection(segs,colors=[COLORS[l] for l in labs],linewidths=.9))
        ax.add_patch(Rectangle((430,150),30,60,fill=False,edgecolor='lime'));ax.scatter(*segs[0][0],c='lime',marker='*',s=60)
        ax.set(xlim=(150,475),ylim=(320,25),title=f'Z{z} rank{rank+1}, score={state[0]:.1f}')
    fig.tight_layout();fig.savefig(OUT/f'slice_{z}_candidate_paths.png',dpi=160);plt.close(fig)
    state=paths[0];segs,labs,trans,ss,xy,prov,lengths=materialize(state,ps,edges,adj,response,z)
    fig,ax=plt.subplots(figsize=(12,9));ax.imshow(image,cmap='gray',vmin=0,vmax=255)
    for i in state[1]:
        for _,j,eid in adj[i][:3]:
            if eid not in state[2]:ax.plot(*np.array([ps[i][-1],ps[j][0]]).T,color='gray',ls=':',lw=.5,alpha=.4)
    ax.add_collection(LineCollection(segs,colors=[COLORS[l] for l in labs],linewidths=1.2))
    ax.add_patch(Rectangle((430,150),30,60,fill=False,edgecolor='lime'));ax.scatter(*segs[0][0],c='lime',marker='*',s=90)
    ax.set(xlim=(150,475),ylim=(320,25),title=f'Z{z}: cyan observed; magenta inferred; orange low confidence; dotted gray alternatives')
    fig.tight_layout();fig.savefig(OUT/f'slice_{z}_best_continuous_thread.png',dpi=180);plt.close(fig)
    palette=np.array([[0,0,0],[0,191,255],[255,0,255],[255,176,0]],np.uint8)
    Image.fromarray(palette[prov][None,:,:]).save(OUT/f'slice_{z}_provenance.png')
    ct=np.full(len(ss),np.nan);valid=np.isin(prov,[1,3]);ct[valid]=sample(image,xy[valid])
    Image.fromarray(np.nan_to_num(ct,nan=0).round().astype('uint8')[None,:]).save(OUT/f'slice_{z}_unwrapped_signal.png')
    fig,axs=plt.subplots(3,1,figsize=(14,5),sharex=True)
    axs[0].plot(ss*.12656,xy[:,0],lw=.7,label='X');axs[0].plot(ss*.12656,xy[:,1],lw=.7,label='Y');axs[0].legend();axs[0].set_ylabel('coordinate px')
    axs[1].imshow(palette[prov][None,:,:],aspect='auto',extent=[0,ss[-1]*.12656,0,1],interpolation='nearest');axs[1].set_ylabel('provenance')
    axs[2].plot(ss*.12656,ct,lw=.6);axs[2].set(xlabel='physical arc coordinate s (mm)',ylabel='CT uint8',ylim=(0,255))
    fig.tight_layout();fig.savefig(OUT/f'slice_{z}_unwrapped_geometry.png',dpi=150);plt.close(fig)
    np.savez_compressed(OUT/f'slice_{z}_samples.npz',s_px=ss,s_um=ss*126.56,xy=xy,provenance=prov,ct=ct,
                        vertex_xy=np.concatenate(segs),vertex_provenance=np.concatenate([np.full(len(p),l,np.uint8) for p,l in zip(segs,labs)]),
                        segment_offsets=np.r_[0,np.cumsum([len(p) for p in segs])],segment_provenance=labs)
    return segs,labs,trans,ss,xy,prov,lengths

def main():
    OUT.mkdir(exist_ok=False)
    prior=[p for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents]
    prior.extend(ROOT/p for p in ['README.md','AGENTS.md','spiral.png'])
    dump(OUT/'preservation.json',{str(p.relative_to(ROOT)):sha(p) for p in prior});dump(OUT/'plan.json',PLAN)
    volume=np.load(ROOT/PLAN['source'],mmap_mode='r');metrics=[]
    baseline={r['z']:r for r in json.loads((ROOT/'outputs/p0_t1/visual_assessment.json').read_text())}
    for z in ZS:
        ps,orig,L,adj,edges,seeds,response,sk,tout=build(z)
        np.savez_compressed(OUT/f'slice_{z}_virtual_edges.npz',edges=edges)
        dump(OUT/f'slice_{z}_edge_columns.json',['source_oriented_chain','target_oriented_chain','distance','tangent_cos','out_chord_cos','in_chord_cos','normal_miss','spacing','turn_radians','required_bend_per_px','ridge_support','crossing_proxy','penalty','admissible'])
        dump(OUT/f'slice_{z}_oriented_chain_map.json',orig)
        top,complete,log,counts=search(ps,orig,L,adj,edges,seeds,tout)
        dump(OUT/f'slice_{z}_search_log.json',log)
        dump(OUT/f'slice_{z}_complete_paths.json',[{'score':s[0],'path':s[1],'edge_ids':s[2]} for s in sorted(complete,key=lambda s:-s[0])])
        dump(OUT/f'slice_{z}_top_paths.json',[{'rank':i+1,'score':s[0],'path':s[1],'edge_ids':s[2],'seed_piece':s[5].tolist()} for i,s in enumerate(top)])
        if not top:raise RuntimeError(f'No candidate path for {z}; evidence retained')
        segs,labs,trans,ss,xy,prov,lengths=render(z,volume[z-128],top,ps,edges,adj,response)
        labs=np.array(labs);bridges=lengths[labs==2];obs=float(lengths[labs!=2].sum());total=float(lengths.sum())
        m={'z':z,**counts,'ranked_alternatives':len(top),'total_arc_px':total,'total_arc_mm':total*.12656,
           'observed_graph_arc_px':obs,'inferred_arc_px':float(bridges.sum()),'observed_graph_percent':100*obs/total,
           'usable_ct_length_px':float(np.isin(prov,[1,3]).sum()*.5),'usable_ct_sample_count':int(np.isin(prov,[1,3]).sum()),
           'fragments':len(top[0][1]),'bridges':len(bridges),'longest_bridge_px':float(bridges.max()) if len(bridges) else 0.,
           'median_bridge_px':float(np.median(bridges)) if len(bridges) else 0.,'low_confidence_transitions':sum(t['low_confidence'] for t in trans),
           'seed_xy':segs[0][0].tolist(),'end_xy':segs[-1][-1].tolist(),'target_distance_reduction_px':float(target_dist(segs[0][0])-target_dist(segs[-1][-1])),
           'radial_progress_px':float(np.linalg.norm(segs[-1][-1]-[240,140])-np.linalg.norm(segs[0][0]-[240,140])),
           'T1_raw_observed_arc_px':baseline[z]['visible_arc_px'],'T1_audited_observed_arc_px':baseline[z]['pre_flag_visible_arc_px'],
           'T1_radial_progress_px':baseline[z]['radial_progress_px'],'observed_gain_vs_T1_raw':obs/baseline[z]['visible_arc_px'],
           'sample_counts_by_label':np.bincount(prov,minlength=4).tolist(),'transitions':trans}
        metrics.append(m);dump(OUT/'metrics.json',metrics);print(json.dumps({k:v for k,v in m.items() if k!='transitions'}),flush=True)
    simple=[{k:v for k,v in m.items() if k not in ['transitions']} for m in metrics]
    with (OUT/'measurements.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(simple[0]));w.writeheader();w.writerows(simple)

if __name__=='__main__':main()
