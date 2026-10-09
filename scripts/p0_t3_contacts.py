"""Broaden T3 uncrossing to vertex-on-edge contacts without removing observed runs."""
import numpy as np
from scipy.spatial import cKDTree
import p0_t3_uncross as c

def intersections_inclusive(a,b):
    v=b-a;length=np.linalg.norm(v,axis=1);s=np.r_[0,np.cumsum(length)]
    pairs=cKDTree((a+b)/2).query_pairs(float(length.max())+1e-8,output_type='ndarray')
    if not len(pairs):return []
    i,j=pairs.T;good=j>i+1;i=i[good];j=j[good]
    cross=lambda x,y:x[...,0]*y[...,1]-x[...,1]*y[...,0]
    den=cross(v[i],v[j]);good=np.abs(den)>1e-9;i=i[good];j=j[good];den=den[good]
    t=cross(a[j]-a[i],v[j])/den;u=cross(a[j]-a[i],v[i])/den
    good=(t>=-1e-8)&(t<=1+1e-8)&(u>=-1e-8)&(u<=1+1e-8)
    return [{'i':int(ii),'j':int(jj),'t':float(tt),'u':float(uu),'xy':(a[ii]+tt*v[ii]).tolist(),
             's_first':float(s[ii]+tt*length[ii]),'s_second':float(s[jj]+uu*length[jj]),
             'loop_arc_px':float(s[jj]+uu*length[jj]-s[ii]-tt*length[ii]),
             'type':'proper' if 1e-8<tt<1-1e-8 and 1e-8<uu<1-1e-8 else 'point_contact'}
            for ii,jj,tt,uu in zip(i[good],j[good],t[good],u[good])]

if __name__=='__main__':
    c.OUT=c.ROOT/'outputs/p0_t3/run_04'
    c.PLAN={**c.PLAN,'contact_amendment':'run03 exhaustive check found two incident segment pairs at one Z146 point contact; admit inclusive vertex-on-edge intersections to the same two-edge correction; exact collinear audit remains separate'}
    c.crossings=intersections_inclusive
    c.main()
