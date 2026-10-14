"""Bounded local ridge assignment, plane-first fit and conditional quadratic."""
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from numcodecs import Blosc
import p0_t11_core as ct
ROOT=Path(__file__).resolve().parents[1]
GRID=np.arange(15,113,12)  # 15,27,...111: fixed9x9 local footprint
@lru_cache(maxsize=384)
def chunk(k,j,i):
    for phase in ['p0_t9','p0_t11','p0_t12']:
        p=ROOT/f'data/raw/{phase}/native_chunks/{k}_{j}_{i}'
        if p.exists():return np.frombuffer(Blosc().decode(p.read_bytes()),np.uint8).reshape(128,128,128)
    raise FileNotFoundError(f'Missing CT support {k}_{j}_{i}; never pad/clip')
ct.chunk=chunk
def candidates(profile):
    p=gaussian_filter1d(np.asarray(profile,float),.75,mode='reflect');idx,props=find_peaks(p,prominence=8)
    answer=[]
    for i,prom in zip(idx,props['prominences']):
        den=p[i-1]-2*p[i]+p[i+1];delta=float(np.clip(.5*(p[i-1]-p[i+1])/den,-.5,.5)) if den else 0.
        answer.append({'depth':float(i-8+delta),'prominence':float(prom),'height':float(p[i])})
    return answer
def choose(peaks,expected,limit=2.):
    ordered=sorted(peaks,key=lambda r:(abs(r['depth']-expected),r['depth']))
    if not ordered:return None,'no_peak'
    a=ordered[0];distance=abs(a['depth']-expected)
    if distance>limit:return None,'outside_displacement_bound'
    if len(ordered)>1 and abs(ordered[1]['depth']-expected)-distance<1.:return None,'competing_ridges'
    return a,'unique_nearest'
def design(u,v,degree):
    return np.stack([np.ones_like(u),u,v]+([u*u,u*v,v*v] if degree==2 else []),axis=-1)
def surface(coeff,degree):
    v,u=np.meshgrid((np.arange(128)-63.5)/63.5,(np.arange(128)-63.5)/63.5,indexing='ij')
    h=design(u,v,degree)@coeff
    du=np.full_like(h,coeff[1]);dv=np.full_like(h,coeff[2])
    if degree==2:du+=2*coeff[3]*u+coeff[4]*v;dv+=coeff[4]*u+2*coeff[5]*v
    return h,du/63.5,dv/63.5
def fit(volume):
    central=candidates(volume[:,61:67,61:67].mean((1,2)));anchor,reason=choose(central,0)
    evidence={'central_candidates':central,'anchor':anchor,'anchor_reason':reason,'sites':[]}
    if anchor is None:return {'plane':{'accepted':False,'reason':reason},'quadratic':{'attempted':False,'accepted':False,'reason':'no_unique_anchor'},'evidence':evidence}
    for iy,y in enumerate(GRID):
        for ix,x in enumerate(GRID):
            peaks=candidates(volume[:,y-1:y+2,x-1:x+2].mean((1,2)));chosen,why=choose(peaks,anchor['depth'])
            evidence['sites'].append({'iy':iy,'ix':ix,'y':int(y),'x':int(x),'candidates':peaks,'chosen':chosen,'reason':why})
    good=[s for s in evidence['sites'] if s['chosen'] is not None];fraction=len(good)/81
    if len(good)<57:return {'plane':{'accepted':False,'reason':'less_than70pct_unique_local_ridges','ridge_fraction':fraction},'quadratic':{'attempted':False,'accepted':False,'reason':'insufficient_identity_evidence'},'evidence':evidence}
    u=np.array([(s['x']-63.5)/63.5 for s in good]);v=np.array([(s['y']-63.5)/63.5 for s in good]);depth=np.array([s['chosen']['depth'] for s in good]);train=np.array([(s['ix']+s['iy'])%2==0 for s in good])
    if train.sum()<20 or (~train).sum()<20:return {'plane':{'accepted':False,'reason':'insufficient_split_support'},'quadratic':{'attempted':False,'accepted':False,'reason':'insufficient_split_support'},'evidence':evidence}
    def fit_degree(degree):
        X=design(u,v,degree);weights=np.ones(train.sum());coeff=np.zeros(X.shape[1]);base=np.sqrt(np.array([s['chosen']['prominence'] for s in good])[train])
        for _ in range(3):
            w=weights*base;coeff=np.linalg.lstsq(X[train]*w[:,None],depth[train]*w,rcond=None)[0];res=abs(X[train]@coeff-depth[train]);weights=np.minimum(1,1/np.maximum(res,1e-9))
        predicted=X@coeff;residual=predicted-depth;h,du,dv=surface(coeff,degree)
        conflict=[]
        for s,estimate in zip(good,predicted):
            d=abs(s['chosen']['depth']-estimate);other=[abs(p['depth']-estimate) for p in s['candidates'] if p is not s['chosen']]
            conflict.append(bool(other and min(other)-d<.5))
        neighbors=[];lookup={(s['iy'],s['ix']):s['chosen']['depth'] for s in good}
        for (iy,ix),d in lookup.items():
            for key in [(iy+1,ix),(iy,ix+1)]:
                if key in lookup:neighbors.append(abs(d-lookup[key]))
        checks={'center_displacement_le2':abs(float(coeff[0]))<=2,'all_patch_displacement_le2':float(abs(h).max())<=2,'heldout_residual_le1':float(np.median(abs(residual[~train])))<=1,'heldout_residual95_le2':float(np.quantile(abs(residual[~train]),.95))<=2,'competitor_conflict_le5pct':float(np.mean(conflict))<=.05,'no_large_neighbor_ridge_steps':bool(neighbors) and max(neighbors)<=2,'max_slope_le01':float(np.hypot(du,dv).max())<=.1}
        checks['center_stays_with_anchor']=abs(float(coeff[0])-anchor['depth'])<=1 and all(abs(p['depth']-float(coeff[0]))-abs(anchor['depth']-float(coeff[0]))>=.5 for p in central if p is not anchor)
        curvature=float(max(abs(2*coeff[3]),abs(coeff[4]),abs(2*coeff[5]))) if degree==2 else 0.
        checks['curvature_le1']=curvature<=1
        return {'attempted':True,'accepted':bool(all(checks.values())),'degree':degree,'coefficients':coeff.tolist(),'checks':checks,'ridge_fraction':fraction,'sheet_center_depth':float(coeff[0]),'anchor_depth':anchor['depth'],'sheet_center_residual':abs(float(coeff[0])-anchor['depth']),'fit_residual_train':float(np.median(abs(residual[train]))),'fit_residual_holdout':float(np.median(abs(residual[~train]))),'competitor_conflict_fraction':float(np.mean(conflict)),'maximum_neighbor_ridge_step':float(max(neighbors)),'maximum_surface_displacement':float(abs(h).max()),'maximum_slope':float(np.hypot(du,dv).max()),'curvature':curvature}
    plane=fit_degree(1)
    # A quadratic is allowed only when stable ridge assignments exist and the
    # plane's held-out residual fails; never rescue ambiguous sheet identity.
    identity=all(plane['checks'][k] for k in ['center_displacement_le2','center_stays_with_anchor','competitor_conflict_le5pct','no_large_neighbor_ridge_steps'])
    # A stricter adequacy trigger is distinct from admissibility. Fixed using
    # only synthetic curvature probes before any target evidence.
    plane_residual=np.abs(design(u,v,1)@np.asarray(plane['coefficients'])-depth)[~train]
    if identity and (np.median(plane_residual)>.15 or np.quantile(plane_residual,.95)>.35):
        quad=fit_degree(2)
        quad['checks']['heldout_improvement20pct']=quad['fit_residual_holdout']<=.8*plane['fit_residual_holdout']
        quad['accepted']=bool(all(quad['checks'].values()))
    else:quad={'attempted':False,'accepted':False,'reason':'plane_residual_not_demonstrably_insufficient_or_identity_failed'}
    return {'plane':plane,'quadratic':quad,'identity_reliable':bool(identity),'evidence':evidence}
def coordinates(P,t,y,n,coeff,degree):
    h,du,dv=surface(np.asarray(coeff),degree);t,y,n=map(np.asarray,[t,y,n]);nn=n[:,None,None]-t[:,None,None]*du-y[:,None,None]*dv;nn/=np.linalg.norm(nn,axis=0)
    s=(np.arange(128)-63.5)[None,:]*ct.STEP;v=(np.arange(128)-63.5)[:,None]*ct.STEP
    sheet=np.asarray(P)[:,None,None]+t[:,None,None]*s+y[:,None,None]*v+n[:,None,None]*h*ct.STEP
    d=(np.arange(17)-8)*ct.STEP
    return sheet[:,None]+nn[:,None]*d[None,:,None,None]
def render(P,t,y,n,fit):
    coords=coordinates(P,t,y,n,fit['coefficients'],fit['degree']);return np.clip(np.rint(ct.sample_xyz(coords)),0,255).astype(np.uint8)
