"""Exact crossing audit and bounded local loop excision, preserving edge provenance."""
import numpy as np
from scipy.spatial import cKDTree
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
def flatten(a):
    verts=a['vertex_xy'];offset=a['segment_offsets'];labels=a['segment_provenance']
    aa=[];bb=[];ll=[]
    for i,label in enumerate(labels):
        p=verts[offset[i]:offset[i+1]]
        for start,end in zip(p,p[1:]):
            if np.linalg.norm(end-start)>1e-9:aa.append(start);bb.append(end);ll.append(label)
    return np.array(aa),np.array(bb),np.array(ll)

def crossings(a,b):
    v=b-a;length=np.linalg.norm(v,axis=1);s=np.r_[0,np.cumsum(length)]
    pairs=cKDTree((a+b)/2).query_pairs(float(length.max())+1e-8,output_type='ndarray')
    if not len(pairs):return []
    i,j=pairs.T;keep=j>i+1;i=i[keep];j=j[keep]
    cross=lambda x,y:x[...,0]*y[...,1]-x[...,1]*y[...,0]
    den=cross(v[i],v[j]);good=np.abs(den)>1e-9;i=i[good];j=j[good];den=den[good]
    t=cross(a[j]-a[i],v[j])/den;u=cross(a[j]-a[i],v[i])/den
    keep=(t>1e-8)&(t<1-1e-8)&(u>1e-8)&(u<1-1e-8)
    return [{'i':int(ii),'j':int(jj),'t':float(tt),'u':float(uu),'xy':(a[ii]+tt*v[ii]).tolist(),
             's_first':float(s[ii]+tt*length[ii]),'s_second':float(s[jj]+uu*length[jj]),
             'loop_arc_px':float(s[jj]+uu*length[jj]-s[ii]-tt*length[ii])} for ii,jj,tt,uu in zip(i[keep],j[keep],t[keep],u[keep])]

if __name__=='__main__':
    report=[]
    for z in [128,144,146,160,170]:
        a,b,l=flatten(np.load(ROOT/f'outputs/p0_t2/slice_{z}_samples.npz'))
        c=crossings(a,b);report.append({'z':z,'crossings':c});print(z,json.dumps(c),flush=True)
    (ROOT/'outputs/p0_t3/local_crossing_diagnosis.json').write_text(json.dumps(report,indent=2)+'\n')
