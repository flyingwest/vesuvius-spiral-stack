"""Independent XY ridge graph traversal. No polar fit or neighboring slice input."""
import json, hashlib, csv, argparse
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.signal import find_peaks
from skimage.morphology import skeletonize
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/p0_t1'
ZS=[128,144,146,160,170]
SPACING=126.56
PLAN={'slices':ZS,'source':'outputs/orientation/ct_level4_zyx.npy','z_origin':128,
      'spacing_um':SPACING,'ds_px':0.5,'roi_xy':[25,475,25,325],
      'seed_reference_xy':[246,112],'interior_reference_xy':[240,140],
      'detection':'DoG sigma 0.6 minus 2, threshold >3; skeleton; maximal graph chains; no branch removal',
      'seed':'nearest point to reference on chain >=20 px; direction toward farther endpoint from interior reference',
      'continuation':'endpoints within 8 px; score .3 tangent + .25 chord + .2 proximity + .15 normal consistency + .1 ridge strength',
      'gates':'score >=.72; margin >=.08 among eligible; tangent cosine >=.5; chord cosines >=.65 for nonzero gaps; normal miss <=.45 local layer spacing; no intervening skeleton crossing',
      'scope':'exploratory fixed-parameter 2D geometry; no physical identity oracle; no cross-Z input',
      'gap':'all nonzero connectors masked; chord length is lower bound, not true missing arc',
      'limit':'200 chains; revisits forbidden; terminal unknown one masked sample, not estimated remaining sheet length'}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x): p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def unit(v): return v/max(np.linalg.norm(v),1e-9)
def arc(p): return np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
def sample(a,p): return map_coordinates(a.astype(float),[p[:,1],p[:,0]],order=1,mode='nearest')

def graph(image):
    response=gaussian_filter(image.astype(float),.6)-gaussian_filter(image.astype(float),2.)
    mask=response>3
    roi=np.zeros_like(mask);roi[25:325,25:475]=True
    sk=skeletonize(mask & roi)
    pixels=set(map(tuple,np.argwhere(sk)))
    adj={}
    for y,x in sorted(pixels):
        adj[y,x]=[(y+dy,x+dx) for dy in [-1,0,1] for dx in [-1,0,1]
                  if (dy or dx) and (y+dy,x+dx) in pixels
                  and not (dy and dx and ((y+dy,x) in pixels or (y,x+dx) in pixels))]
    used=set(); chains=[]
    # Start at junctions/endpoints, then consume remaining cycles deterministically.
    starts=sorted(adj,key=lambda p:(len(adj[p])==2,p))
    for p in starts:
        for q in adj[p]:
            edge=tuple(sorted((p,q)))
            if edge in used:continue
            line=[p,q];used.add(edge)
            while len(adj[q])==2 and q!=p:
                nxt=next(v for v in adj[q] if v!=line[-2])
                edge=tuple(sorted((q,nxt)))
                if edge in used:break
                used.add(edge);line.append(nxt);q=nxt
            chains.append(np.array([(x,y) for y,x in line],float))
    return response,sk,chains

def process_slice(image):
    response,sk,chains=graph(image)
    center=np.array([240,140]);reference=np.array([246,112])
    eligible=[(float(np.linalg.norm(p-reference,axis=1).min()),i) for i,p in enumerate(chains) if arc(p)[-1]>=20]
    if not eligible:raise RuntimeError('No seed chain; retain failure rather than fabricate a seed')
    _,sid=min(eligible);p=chains[sid];j=int(np.argmin(np.linalg.norm(p-reference,axis=1)))
    # Choose direction locally once; no monotonic radial or angular objective thereafter.
    if np.linalg.norm(p[0]-center)>np.linalg.norm(p[-1]-center):p=p[::-1];j=len(p)-1-j
    first=p[j:]
    if len(first)<2:first=p[max(0,j-1):]
    pieces=[first];ids=[sid];gaps=[];decisions=[];visited={sid}
    for step in range(200):
        old=pieces[-1];end=old[-1];t=unit(old[-1]-old[max(0,len(old)-6)])
        normal=np.array([-t[1],t[0]])
        offsets=np.arange(-10,10.01,.25);profile=sample(response,end+offsets[:,None]*normal)
        peaks,_=find_peaks(profile,prominence=2,distance=4)
        spacing=float(np.median(np.diff(offsets[peaks]))) if len(peaks)>1 else 4.
        spacing=float(np.clip(spacing,2,10))
        candidates=[]
        for cid,q0 in enumerate(chains):
            if cid in visited:continue
            for reverse in [False,True]:
                q=q0[::-1] if reverse else q0
                delta=q[0]-end;d=float(np.linalg.norm(delta))
                if d>8 or arc(q)[-1]<3:continue
                u=unit(q[min(5,len(q)-1)]-q[0]);alignment=float(t@u)
                chord=min(float(t@unit(delta)),float(u@unit(delta))) if d else 1.
                miss=abs(float(delta@normal));strength=float(np.clip(sample(response,q[:6]).mean()/12,0,1))
                connector=end+np.linspace(0,1,max(3,int(d*4)+1))[:,None]*delta
                interior=connector[(np.linalg.norm(connector-end,axis=1)>1.5)&(np.linalg.norm(connector-q[0],axis=1)>1.5)]
                crossing=bool(len(interior) and np.any(sk[np.rint(interior[:,1]).astype(int),np.rint(interior[:,0]).astype(int)]))
                score=.3*max(0,alignment)+.25*max(0,chord)+.2*np.exp(-d/8)+.15*np.exp(-miss/spacing)+.1*strength
                reasons=[]
                if alignment<.5:reasons.append('tangent')
                if d and chord<.65:reasons.append('chord')
                if miss>.45*spacing:reasons.append('normal/layer-spacing')
                if crossing:reasons.append('intervening-ridge')
                if score<.72:reasons.append('score')
                candidates.append({'chain':cid,'reverse':reverse,'xy':q[0].tolist(),'distance':d,'tangent_cos':alignment,
                                   'chord_cos':chord,'normal_miss':miss,'spacing':spacing,'ridge_strength':strength,
                                   'score':float(score),'rejected_by':reasons})
        candidates.sort(key=lambda c:(-c['score'],c['chain'],c['reverse']))
        viable=[c for c in candidates if not c['rejected_by']]
        winner=viable[0] if viable else None
        margin=winner['score']-viable[1]['score'] if len(viable)>1 else 1.
        accept=bool(winner and margin>=.08)
        decisions.append({'endpoint':end.tolist(),'candidates':candidates,'margin':margin,
                          'decision':'CONTINUE' if accept else 'UNKNOWN','selected':winner['chain'] if accept else None})
        if not accept:break
        q=chains[winner['chain']];q=q[::-1] if winner['reverse'] else q
        gaps.append(float(np.linalg.norm(q[0]-end)));pieces.append(q);ids.append(winner['chain']);visited.add(winner['chain'])
    return response,sk,chains,pieces,ids,gaps,decisions

def raster(image,pieces,gaps):
    lengths=[float(arc(p)[-1]) for p in pieces];extent=sum(lengths)+sum(gaps)
    ss=np.arange(0,extent+.000001,.5);ss=np.r_[ss,ss[-1]+.5]
    xy=np.full((len(ss),2),np.nan);row=np.full(len(ss),np.nan);start=0.
    for i,p in enumerate(pieces):
        a=arc(p);sel=(ss>=start-1e-9)&(ss<=start+a[-1]+1e-9)
        xy[sel]=np.c_[np.interp(ss[sel]-start,a,p[:,0]),np.interp(ss[sel]-start,a,p[:,1])]
        row[sel]=sample(image,xy[sel]);start+=a[-1]+(gaps[i] if i<len(gaps) else 0)
    return ss,xy,row

def render(z,image,response,sk,chains,pieces,gaps,decisions,row):
    prefix=OUT/f'slice_{z}'
    Image.fromarray(image).save(str(prefix)+'_ct.png')
    fig,ax=plt.subplots(figsize=(12,8));ax.imshow(image,cmap='gray',vmin=0,vmax=255)
    for p in chains:ax.plot(p[:,0],p[:,1],color='cyan',lw=.35,alpha=.65)
    ax.set(xlim=(25,475),ylim=(325,25),title=f'Z={z}: all detected ridge graph chains (not all are papyrus)')
    fig.savefig(str(prefix)+'_ridge_candidates.png',dpi=160);plt.close(fig)
    for with_gaps in [False,True]:
        fig,axs=plt.subplots(1,2,figsize=(15,7))
        allp=np.concatenate(pieces)
        for ax in axs:
            ax.imshow(image,cmap='gray',vmin=0,vmax=255)
            for i,p in enumerate(pieces):ax.plot(p[:,0],p[:,1],lw=1.5,label=f'visible {i+1}',color=plt.cm.tab10(i%10))
            ax.scatter(*pieces[0][0],c='lime',marker='*',s=120,label='seed',zorder=5)
            ax.scatter(*pieces[-1][-1],c='red',marker='x',s=70,label='UNKNOWN endpoint',zorder=5)
            if with_gaps:
                for i,d in enumerate(gaps):
                    if d:ax.plot(*np.array([pieces[i][-1],pieces[i+1][0]]).T,'m--',lw=1.2,label='masked connector' if i==0 else None)
                for dec in decisions:
                    for c in dec['candidates'][:5]:
                        if c['chain']!=dec['selected']:ax.plot(*np.array([dec['endpoint'],c['xy']]).T,color='yellow',ls=':',lw=.65,alpha=.6)
        axs[0].set(xlim=(25,475),ylim=(325,25),title=f'Z={z}: independent thread / provisional identity')
        lo=allp.min(0)-15;hi=allp.max(0)+15
        axs[1].set(xlim=(lo[0],hi[0]),ylim=(hi[1],lo[1]),title='Detail: dashed magenta = unsampled gap; yellow = alternatives')
        axs[0].legend(loc='lower left',fontsize=7)
        fig.tight_layout();fig.savefig(str(prefix)+('_thread_trace_with_gaps.png' if with_gaps else '_thread_trace.png'),dpi=160);plt.close(fig)
    Image.fromarray(np.nan_to_num(row,nan=0).round().clip(0,255).astype('uint8')[None,:]).save(str(prefix)+'_unwrapped_row.png')
    Image.fromarray((~np.isfinite(row)).astype('uint8')[None,:]*255).save(str(prefix)+'_row_mask.png')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    volume=np.load(ROOT/PLAN['source'],mmap_mode='r')
    if args.verify:
        results=[]
        for z in reversed(ZS):
            image=volume[z-128];r=process_slice(image);saved=json.loads((OUT/f'slice_{z}_graph.json').read_text())
            assert r[6]==saved['decisions'] and r[4]==saved['selected_chains']
            ss,xy,row=raster(image,r[3],r[5]);a=np.load(OUT/f'slice_{z}_samples.npz')
            assert np.array_equal(row,a['ct'],equal_nan=True) and np.array_equal(xy,a['xy'],equal_nan=True)
            valid=np.isfinite(row);p=xy[valid];x=np.floor(p[:,0]).astype(int);y=np.floor(p[:,1]).astype(int);dx=p[:,0]-x;dy=p[:,1]-y
            expected=image[y,x]*(1-dx)*(1-dy)+image[y,x+1]*dx*(1-dy)+image[y+1,x]*(1-dx)*dy+image[y+1,x+1]*dx*dy
            assert np.allclose(expected,row[valid],atol=1e-10,rtol=0)
            for f in OUT.glob(f'slice_{z}_*.png'):
                with Image.open(f) as im:im.verify()
            results.append({'z':z,'replay_and_bilinear':'PASS'})
        for p,h in json.loads((OUT/'preservation.json').read_text()).items():assert sha(ROOT/p)==h,p
        dump(OUT/'verification.json',{'checks':results,'prior_artifact_preservation':'PASS'})
        print('PASS: five independent reverse-order replays, bilinear sampling, PNGs, preservation');return
    OUT.mkdir(exist_ok=False)
    previous=[p for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents]
    previous += [ROOT/'README.md',ROOT/'AGENTS.md',ROOT/'spiral.png']
    dump(OUT/'preservation.json',{str(p.relative_to(ROOT)):sha(p) for p in previous})
    dump(OUT/'plan.json',PLAN);dump(OUT/'source_identity.json',{'script_sha256':sha(Path(__file__)),'volume_sha256':sha(ROOT/PLAN['source'])})
    metrics=[]
    for z in ZS:
        image=np.asarray(volume[z-128]);response,sk,chains,pieces,ids,gaps,decisions=process_slice(image)
        ss,xy,row=raster(image,pieces,gaps)
        np.savez_compressed(OUT/f'slice_{z}_candidates.npz',response=response,skeleton=sk)
        dump(OUT/f'slice_{z}_graph.json',{'chains':[p.tolist() for p in chains],'selected_chains':ids,'pieces':[p.tolist() for p in pieces],'gap_lengths_px':gaps,'decisions':decisions})
        np.savez_compressed(OUT/f'slice_{z}_samples.npz',s_px=ss,s_um=ss*SPACING,xy=xy,ct=row,missing=~np.isfinite(row))
        points=np.concatenate(pieces);radius=np.linalg.norm(points-[240,140],axis=1)
        lengths=[float(arc(p)[-1]) for p in pieces]
        m={'z':z,'seed_x':float(points[0,0]),'seed_y':float(points[0,1]),'visible_arc_px':sum(lengths),'gap_arc_lower_bound_px':sum(gaps),
           'extent_px':sum(lengths)+sum(gaps),'visible_arc_mm':sum(lengths)*SPACING/1000,'visible_fragments':1+sum(d>0 for d in gaps),
           'graph_chains_traversed':len(pieces),'interruptions':1+sum(d>0 for d in gaps),'accepted_reconnections':sum(d>0 for d in gaps),
           'unresolved_endpoints':1,'suspected_wrong_layer_jumps':None,'start_radius_px':float(radius[0]),'end_radius_px':float(radius[-1]),
           'radial_progress_px':float(radius[-1]-radius[0]),'maximum_radius_px':float(radius.max()),'sampling_um':SPACING*.5,
           'row_width':len(row),'masked_fraction':float(np.mean(~np.isfinite(row))),'candidate_chains':len(chains),
           'confidence':'provisional; requires visual audit, no oracle'}
        metrics.append(m);render(z,image,response,sk,chains,pieces,gaps,decisions,row)
        print(json.dumps(m),flush=True)
    dump(OUT/'metrics.json',metrics)
    with (OUT/'measurements.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)

if __name__=='__main__':main()
