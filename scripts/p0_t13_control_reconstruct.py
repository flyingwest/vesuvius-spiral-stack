"""CT-only single-slice reconstruction worker. No trusted paths or label imports."""
import sys,os,json,hashlib
access=[]
def audit(event,args):
 if event=='open':
  path=args[0]
  if isinstance(path,(str,bytes)):
   path=os.fsdecode(path);access.append({'path':path,'mode':str(args[1]),'flags':args[2]})
   if any(x in path.lower() for x in ['geometry','inklabel','supervision','trusted','probability']):raise RuntimeError('Forbidden data access: '+path)
sys.addaudithook(audit)
import numpy as np
from scipy.ndimage import gaussian_filter,map_coordinates
from scipy.spatial import cKDTree
from skimage.morphology import skeletonize

def unit(a):return a/max(np.linalg.norm(a),1e-9)
def arc(p):return np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
def graph(image):
 response=gaussian_filter(image.astype(float),.6)-gaussian_filter(image.astype(float),2)
 mask=response>3;mask[:4]=False;mask[-4:]=False;mask[:,:4]=False;mask[:,-4:]=False
 sk=skeletonize(mask);pixels=set(map(tuple,np.argwhere(sk)));adj={}
 for y,x in sorted(pixels):
  adj[y,x]=[(y+dy,x+dx) for dy in [-1,0,1] for dx in [-1,0,1] if (dx or dy) and (y+dy,x+dx) in pixels and not (dx and dy and ((y+dy,x) in pixels or (y,x+dx) in pixels))]
 used=set();chains=[]
 for p in sorted(adj,key=lambda p:(len(adj[p])==2,p)):
  for q in adj[p]:
   edge=tuple(sorted((p,q)))
   if edge in used:continue
   line=[p,q];used.add(edge)
   while len(adj[q])==2 and q!=p:
    nxt=next(v for v in adj[q] if v!=line[-2]);edge=tuple(sorted((q,nxt)))
    if edge in used:break
    used.add(edge);line.append(nxt);q=nxt
   a=np.array([(x,y) for y,x in line],float)
   # Open a closed cycle by dropping only its duplicate endpoint edge.
   if np.array_equal(a[0],a[-1]):a=a[:-1]
   if len(a)>1 and arc(a)[-1]>=3:chains.append(a)
 return response,sk,chains

def intersections(a,b,c,d):
 # Vectorized exact-enough pixel/polyline orientation test; includes contacts/overlap.
 def cross(v,w):return v[...,0]*w[...,1]-v[...,1]*w[...,0]
 o1=cross(b-a,c-a);o2=cross(b-a,d-a);o3=cross(d-c,a-c);o4=cross(d-c,b-c)
 box=(np.maximum(np.minimum(a,b),np.minimum(c,d))<=np.minimum(np.maximum(a,b),np.maximum(c,d))+1e-9).all(-1)
 return box&(o1*o2<=1e-9)&(o3*o4<=1e-9)
def topology(p):
 a=p[:-1];b=p[1:];counts={'proper':0,'contacts':0,'overlaps':0}
 for i in range(len(a)-2):
  js=np.flatnonzero(intersections(a[i],b[i],a[i+2:],b[i+2:]))+i+2
  for j in js:
   cross=lambda v,w:v[0]*w[1]-v[1]*w[0]
   u=cross(b[i]-a[i],a[j]-a[i]);v=cross(b[i]-a[i],b[j]-a[i]);q=cross(b[j]-a[j],a[i]-a[j]);r=cross(b[j]-a[j],b[i]-a[j])
   if u*v<0 and q*r<0:counts['proper']+=1
   elif abs(u)+abs(v)+abs(q)+abs(r)<1e-9:
    axis=int(np.argmax(abs(b[i]-a[i])));extent=min(max(a[i,axis],b[i,axis]),max(a[j,axis],b[j,axis]))-max(min(a[i,axis],b[i,axis]),min(a[j,axis],b[j,axis]))
    counts['overlaps' if extent>1e-9 else 'contacts']+=1
   else:counts['contacts']+=1
 return counts

def reconstruct(image):
 response,sk,chains=graph(image)
 # Preserve all detector chains and all eligible links. Filter internally non-simple
 # chains before search (they remain in graph archive) rather than erase final data.
 eligible=[i for i,p in enumerate(chains) if not any(topology(p).values())]
 if not eligible:return None,chains,[],{'reason':'no_simple_CT_chain'}
 oriented=[chains[i][::sign] for i in eligible for sign in [1,-1]];starts=np.array([p[0] for p in oriented]);tree=cKDTree(starts);length=np.array([arc(p)[-1] for p in oriented]);links=[]
 for i,p in enumerate(oriented):
  t=unit(p[-1]-p[max(0,len(p)-6)]);n=np.array([-t[1],t[0]]);rows=[]
  for j in sorted(tree.query_ball_point(p[-1],12)):
   if i//2==j//2:continue
   q=oriented[j];delta=q[0]-p[-1];d=float(np.linalg.norm(delta));u=unit(q[min(5,len(q)-1)]-q[0]);co=float(t@u);c1=float(t@unit(delta)) if d else 1.;c2=float(u@unit(delta)) if d else 1.
   if co<-.5 or min(c1,c2)<-.35:continue
   miss=abs(float(delta@n));line=p[-1]+np.linspace(.2,.8,5)[:,None]*delta;support=float(np.clip(map_coordinates(response,[line[:,1],line[:,0]],order=1)/12,0,1).mean());cost=2*d+3*np.arccos(np.clip(co,-1,1))+2*(np.arccos(np.clip(c1,-1,1))+np.arccos(np.clip(c2,-1,1)))+4*miss/4+2*(1-support)*(d>0)
   rows.append({'to':j,'distance':d,'cost':float(cost),'tangent_cos':co,'normal_miss':miss,'support':support,'low_confidence':bool(co<.7 or miss>1.8)})
  links.append(sorted(rows,key=lambda a:(a['cost'],a['to'])))
 seeds=sorted(range(len(oriented)),key=lambda i:(-length[i],i))[:16]
 # Both orientations of the eight longest eligible source chains.
 seedids=sorted(range(len(eligible)),key=lambda i:(-length[2*i],i))[:8];seeds=[2*i+s for i in seedids for s in [0,1]]
 states=[{'score':float(length[i]),'path':[i],'poly':oriented[i].copy(),'labels':np.ones(len(oriented[i])-1,dtype=np.uint8),'used':frozenset([i//2]),'decisions':[]} for i in seeds];best=max(states,key=lambda s:s['score']);trace=[];rejected_topology=0
 for depth in range(64):
  expanded=[]
  for state in states:
   i=state['path'][-1];p=state['poly'];a=p[:-1];b=p[1:]
   for link in links[i][:12]:
    j=link['to']
    if j//2 in state['used']:continue
    q=oriented[j];d=link['distance'];new=np.vstack([p[-1],q]) if d else q
    ca=new[:-1];cb=new[1:];valid=True
    if d and len(q)>2 and intersections(p[-1],q[0],q[1:-1],q[2:]).any():rejected_topology+=1;continue
    # Only the immediate joining endpoint is exempt; all other contacts rejected.
    for k,(aa,bb) in enumerate(zip(ca,cb)):
     olda=a[:-1] if k==0 else a;oldb=b[:-1] if k==0 else b
     if len(olda) and intersections(aa,bb,olda,oldb).any():valid=False;break
     if k==0 and len(a) and intersections(aa,bb,a[-1],b[-1]):
      # Collinear reversal at joining edge is forbidden; straight continuation fine.
      if (bb-aa)@(b[-1]-a[-1])<0:valid=False;break
    if not valid:rejected_topology+=1;continue
    labels=np.ones(len(q)-1,np.uint8)
    if link['low_confidence']:labels[np.r_[0,np.cumsum(np.linalg.norm(np.diff(q,axis=0),axis=1))][:-1]<=5]=3
    poly=np.vstack([p,q]) if d else np.vstack([p,q[1:]])
    lab=np.r_[state['labels'],[2] if d else [],labels].astype(np.uint8)
    score=state['score']+length[j]-link['cost']
    newstate={'score':float(score),'path':state['path']+[j],'poly':poly,'labels':lab,'used':state['used']|{j//2},'decisions':state['decisions']+[link]};expanded.append(newstate)
    if score>best['score']:best=newstate
  if not expanded:break
  expanded.sort(key=lambda s:(-s['score'],s['path']));states=expanded[:128];trace.append({'depth':depth+1,'states':len(expanded),'best_score':best['score']})
 final=topology(best['poly']);assert not any(final.values()),final
 return best,chains,links,{'search_trace':trace,'eligible_source_chains':eligible,'seed_source_chains':[eligible[i] for i in seedids],'all_chains':len(chains),'topology_rejected_extensions':rejected_topology,'final_topology':final,'beam':128,'depth_cap':64,'outgoing_cap':12}

def main():
 inp,out=sys.argv[1:3];os.mkdir(out);image=np.load(inp);result,chains,links,stats=reconstruct(image)
 def write(name,obj):
  with open(out+'/'+name,'w') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
 flat=np.concatenate(chains) if chains else np.empty((0,2));offsets=np.r_[0,np.cumsum([len(p) for p in chains])];np.savez_compressed(out+'/graph.npz',xy=flat,offsets=offsets);write('links.json',links);write('search.json',stats)
 if result:
  p=result['poly'];ss=arc(p);grid=np.arange(0,ss[-1]+1e-8,.5);ids=np.minimum(np.searchsorted(ss,grid,side='right')-1,len(p)-2);xy=np.stack([np.interp(grid,ss,p[:,a]) for a in [0,1]],1);labels=result['labels'][ids];mask=labels==2;signal=map_coordinates(image.astype(float),[xy[:,1],xy[:,0]],order=1);signal[mask]=np.nan
  np.savez_compressed(out+'/thread.npz',vertices_xy=p,edge_provenance=result['labels'],edge_s_px=ss,s_px=grid,xy=xy,provenance=labels,mask=mask,ct=signal)
  metrics={'status':'RECONSTRUCTED','path':result['path'],'score':result['score'],'observed_length_px':float(np.diff(ss)[result['labels']!=2].sum()),'bridge_length_px':float(np.diff(ss)[result['labels']==2].sum()),'samples':len(grid),'observed_samples':int((~mask).sum()),'low_confidence_samples':int((labels==3).sum()),'bridges':int((result['labels']==2).sum()),'final_topology':stats['final_topology']};write('decisions.json',result['decisions']);write('metrics.json',metrics)
 else:write('metrics.json',{'status':'FAIL',**stats})
 # This hook records Python open events including imports. Container mounts enforce
 # absence of evaluation data even for native-extension reads not visible to hook.
 with open(out+'/access_log.json','w') as f:json.dump(access,f,indent=2)
 print(inp,'complete',flush=True)
if __name__=='__main__':main()
