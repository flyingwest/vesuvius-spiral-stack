"""Frozen-point CT orientation and patch rendering, native XYZ conventions."""
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.ndimage import gaussian_filter,map_coordinates
from scipy.signal import find_peaks
from numcodecs import Blosc
ROOT=Path(__file__).resolve().parents[1]; STEP=9.6/7.91
@lru_cache(maxsize=384)
def chunk(k,j,i):
    return np.frombuffer(Blosc().decode((ROOT/f'data/raw/p0_t9/native_chunks/{k}_{j}_{i}').read_bytes()),np.uint8).reshape(128,128,128)
def sample_xyz(coords):
    c=np.asarray(coords)[::-1];lo=np.floor(c.reshape(3,-1).min(1)).astype(int);hi=np.floor(c.reshape(3,-1).max(1)).astype(int)+2
    volume=np.empty(hi-lo,np.uint8)
    for k in range(lo[0]//128,(hi[0]-1)//128+1):
        for j in range(lo[1]//128,(hi[1]-1)//128+1):
            for i in range(lo[2]//128,(hi[2]-1)//128+1):
                origin=np.array([k,j,i])*128;start=np.maximum(lo,origin);stop=np.minimum(hi,origin+128)
                volume[tuple(slice(a,b) for a,b in zip(start-lo,stop-lo))]=chunk(k,j,i)[tuple(slice(a,b) for a,b in zip(start-origin,stop-origin))]
    local=c-lo.reshape((3,)+(1,)*(c.ndim-1));v=map_coordinates(volume,local,order=1,prefilter=False,mode='constant',cval=np.nan,output=np.float32)
    assert np.isfinite(v).all();return v
def unit(v):return v/np.linalg.norm(v)
def cube(P):
    z,y,x=np.meshgrid(*[np.arange(-16,17)]*3,indexing='ij');return sample_xyz(np.stack([x,y,z])+np.asarray(P)[:,None,None,None])
def estimate(P,t,oldn,sigma):
    a=cube(P);g=np.stack([gaussian_filter(a,sigma,order=o,mode='reflect') for o in [(0,0,1),(0,1,0),(1,0,0)]],-1)
    z,y,x=np.meshgrid(*[np.arange(-16,17)]*3,indexing='ij');w=np.exp(-(x*x+y*y+z*z)/72);flat=g.reshape(-1,3).astype(float);J=flat.T@(flat*w.reshape(-1,1))/w.sum()
    t=unit(t);oldn=unit(oldn);basis=np.stack([oldn,[0,0,1]],axis=1);ev,vec=np.linalg.eigh(basis.T@J@basis);n=unit(basis@vec[:,-1]);n*=1 if n@oldn>=0 else -1
    fullval,fullvec=np.linalg.eigh(J);tz=np.cross(n,t);confidence=float((ev[-1]-ev[0])/max(ev[-1],1e-12));disagreement=float(np.degrees(np.arcsin(np.clip(abs(fullvec[:,-1]@t),0,1))))
    return {'n':n,'t_s':t,'t_z':tz,'image_y':-tz,'confidence':confidence,'tangent_incompatibility_deg':disagreement,'normal_z':float(n[2]),'tilt_deg':float(np.degrees(np.arccos(np.clip(abs(tz[2]),0,1)))),'energy':float(np.trace(J))}
def frame(P,t,n,sigma):
    f=estimate(P,t,n,sigma);pert=[]
    for shift in [-2,2]:
        q=estimate(np.asarray(P)+[0,0,shift],t,n,sigma);pert.append(float(np.degrees(np.arccos(np.clip(abs(f['n']@q['n']),0,1)))))
    f['perturbation_deg']=max(pert);f['accepted']=bool(f['confidence']>=.25 and f['perturbation_deg']<=20 and f['tangent_incompatibility_deg']<=30 and f['tilt_deg']<=75 and f['energy']>1e-6)
    assert np.allclose(np.cross(f['t_s'],f['t_z']),f['n'],atol=1e-6)
    return f
def render(P,t,y,n):
    d=(np.arange(17)-8)[:,None,None]*STEP;s=(np.arange(128)-63.5)[None,None,:]*STEP;h=(np.arange(128)-63.5)[None,:,None]*STEP
    c=np.asarray(P)[:,None,None,None]+np.asarray(n)[:,None,None,None]*d+np.asarray(t)[:,None,None,None]*s+np.asarray(y)[:,None,None,None]*h
    return np.clip(np.rint(sample_xyz(c)),0,255).astype(np.uint8)
def quality(a):
    # Hold out distant Y wings from the central native +/-16 orientation cube.
    b=gaussian_filter(a.astype(float),1);wings=np.r_[0:40,88:128];q=b[:,wings,:]
    gy=np.gradient(b,axis=1)[:,wings,:];gn=np.gradient(b,axis=0)[:,wings,:];ratio=float(np.mean(gy*gy)/(np.mean(gn*gn)+np.mean(gy*gy)+1e-9))
    peaks=np.argmax(q,axis=0);center=float(np.mean(abs(peaks-8)/8));spread=float(np.std(peaks));crossings=[]
    for y in [16,64,111]:crossings.append(len(find_peaks(b[8,y],prominence=20)[0]))
    return {'longitudinal_gradient_fraction':ratio,'depth_peak_spread':spread,'center_distance':center,'empty_fraction':float((a==0).mean()),'band_crossings':float(np.mean(crossings)),'center_slab_occupancy':float((a[6:11]>32).mean())}
