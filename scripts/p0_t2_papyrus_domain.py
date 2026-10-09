"""Disclosed T2 domain correction: exclude holder, leave objective/gates unchanged."""
import p0_t2_recovery as m
import numpy as np
from matplotlib.path import Path

POLYGON=[[50,80],[95,52],[135,72],[180,65],[245,25],[280,33],[350,65],[415,115],
         [465,150],[470,200],[435,235],[375,280],[300,310],[230,320],[165,300],[120,250],[90,220],[50,165]]
original_build=m.build

def build(z):
    ps,orig,L,adj,edges,seeds,response,sk,tout=original_build(z)
    boundary=Path(np.array(POLYGON))
    admitted=np.array([boundary.contains_points(p,radius=1e-7).all() for p in ps])
    filtered=[]
    for i,options in enumerate(adj):
        keep=[]
        for penalty,j,eid in options:
            connector=np.linspace(ps[i][-1],ps[j][0],25)
            if admitted[i] and admitted[j] and boundary.contains_points(connector,radius=1e-7).all():keep.append((penalty,j,eid))
            else:edges[eid,13]=0
        filtered.append(keep)
    assert all(admitted[i] for i,p in seeds)
    yy,xx=np.mgrid[:493,:506]
    np.savez_compressed(m.OUT/f'slice_{z}_domain.npz',mask=boundary.contains_points(np.c_[xx.ravel(),yy.ravel()]).reshape(493,506),admitted_orientations=admitted)
    return ps,orig,L,filtered,edges,seeds,response,sk,tout

if __name__=='__main__':
    m.OUT=m.ROOT/'outputs/p0_t2/run_02'
    m.PLAN={**m.PLAN,'domain_polygon_xy':POLYGON,'amendment':'After run01 holder-edge excursion observed: exclude external support with one coarse same-XY polygon; no objective/bridge/beam changes; exploratory, not preregistered'}
    m.build=build
    m.main()
