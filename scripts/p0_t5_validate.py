"""Saved-geometry replay and source checks; no search rerun or path modification."""
import json
import numpy as np
from PIL import Image
from p0_t5_consecutive import ROOT,OUT,ZS,SOURCE,SPACING,DS,t1,t2,uncross,flatten,dump,exhaustive_topology
from p0_t5_evaluate import barrier

def main():
    barrier();volume=np.load(SOURCE,mmap_mode='r');checks=[]
    for z in ZS:
        dest=OUT/f'slice_{z}';image=volume[z-128];p=np.load(dest/'path_and_samples.npz')
        assert np.array_equal(image,np.asarray(Image.open(dest/'ct.png')))
        response,sk,chains=t1.graph(image)
        saved=np.load(dest/f'graph_input/outputs/p0_t1/slice_{z}_candidates.npz')
        graph=json.loads((dest/f'graph_input/outputs/p0_t1/slice_{z}_graph.json').read_text())
        assert np.array_equal(response,saved['response']) and np.array_equal(sk,saved['skeleton'])
        assert len(chains)==len(graph['chains']) and all(np.array_equal(a,b) for a,b in zip(chains,graph['chains']))
        aa,bb,ll=flatten(np.load(dest/'precleanup_path.npz'));origin=np.arange(len(aa))
        for correction in json.loads((dest/'corrections.json').read_text()):
            cr=correction['crossing'];aa,bb,ll,origin=uncross.change(aa,bb,ll,origin,cr['i'],cr['j'])
        keep=np.linalg.norm(bb-aa,axis=1)>1e-9;aa,bb,ll,origin=aa[keep],bb[keep],ll[keep],origin[keep]
        assert np.array_equal(aa,p['edge_start_xy']) and np.array_equal(bb,p['edge_end_xy'])
        assert np.array_equal(ll,p['edge_original_provenance']) and np.array_equal(origin,p['original_T2_edge_id'])
        assert np.array_equal(bb[:-1],aa[1:]);length=np.linalg.norm(bb-aa,axis=1);s=np.r_[0,np.cumsum(length)]
        ss=np.arange(0,s[-1]+1e-8,.5);xy=np.column_stack([np.interp(ss,s,np.r_[aa[:,j],bb[-1,j]]) for j in range(2)])
        assert np.array_equal(s,p['edge_s_px']) and np.array_equal(ss*SPACING,p['s_um'])
        assert np.allclose(xy,p['xy'],rtol=0,atol=1e-10)
        for k in np.flatnonzero(origin<0):
            near=((s[:-1]>=s[k]-5)&(s[:-1]<=s[k+1]+5))|((s[1:]>=s[k]-5)&(s[1:]<=s[k+1]+5));ll[near&np.isin(ll,[1,3])]=3
        assert np.array_equal(ll,p['edge_provenance'])
        idx=np.minimum(np.searchsorted(s,ss,side='right')-1,len(aa)-1);labels=ll[idx].copy()
        labels[np.isin(labels,[1,3])&~(t2.sample(response,p['xy'])>3)]=0
        assert np.array_equal(labels,p['provenance'])
        mask=np.isin(labels,[0,2]);assert np.array_equal(mask,p['mask']) and np.isnan(p['ct'][mask]).all()
        q=p['xy'][~mask];x=np.floor(q[:,0]).astype(int);y=np.floor(q[:,1]).astype(int);dx=q[:,0]-x;dy=q[:,1]-y
        bilinear=image[y,x]*(1-dx)*(1-dy)+image[y,x+1]*dx*(1-dy)+image[y+1,x]*(1-dx)*dy+image[y+1,x+1]*dx*dy
        assert np.allclose(bilinear,p['ct'][~mask],rtol=0,atol=1e-10)
        topo=exhaustive_topology(aa,bb);assert topo==json.loads((dest/'validation.json').read_text())['topology']
        checks.append({'z':z,'source_plane':'exact','ridge_graph_replay':'exact','correction_replay':'exact','XY_interpolation':'within1e-10','bilinear':'within1e-10','masks_and_provenance':'exact','topology':topo})
    dump(OUT/'geometry_replay_validation.json',{'status':'PASS','meaning':'all recorded values reproduced, not a claim every topology is clean','slices':checks})
    print('PASS: 50 source planes, independent ridge graphs, saved correction replay, cumulative arc, XY, bilinear, provenance and exhaustive topology')
if __name__=='__main__':main()
