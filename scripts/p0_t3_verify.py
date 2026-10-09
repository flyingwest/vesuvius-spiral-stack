"""Independent full segment-pair audit, provenance lineage and sampling validation."""
import csv,json
from pathlib import Path
import numpy as np
from PIL import Image
import p0_t2_recovery as m
from p0_t3_local_topology import flatten

ROOT=m.ROOT;OUT=ROOT/'outputs/p0_t3';RUN=OUT/'run_03'

def exhaustive_topology(a,b):
    proper=touch=overlap=0
    cross=lambda x,y:x[...,0]*y[...,1]-x[...,1]*y[...,0]
    for i in range(len(a)-2):
        c=a[i+2:];d=b[i+2:];v=b[i]-a[i];w=d-c
        box=(np.maximum(np.minimum(a[i],b[i]),np.minimum(c,d))<=np.minimum(np.maximum(a[i],b[i]),np.maximum(c,d))+1e-10).all(axis=1)
        c=c[box];d=d[box];w=w[box]
        o1=cross(v,c-a[i]);o2=cross(v,d-a[i]);o3=cross(w,a[i]-c);o4=cross(w,b[i]-c)
        p=(o1*o2<-1e-12)&(o3*o4<-1e-12);proper+=int(p.sum())
        inclusive=(o1*o2<=1e-12)&(o3*o4<=1e-12)
        col=(np.abs(o1)<1e-10)&(np.abs(o2)<1e-10)&(np.abs(o3)<1e-10)&(np.abs(o4)<1e-10)
        dim=int(np.argmax(np.abs(v)))
        ol=np.minimum(max(a[i,dim],b[i,dim]),np.maximum(c[:,dim],d[:,dim]))-np.maximum(min(a[i,dim],b[i,dim]),np.minimum(c[:,dim],d[:,dim]))
        ov=col&(ol>1e-9);overlap+=int(ov.sum());touch+=int((inclusive&~p&~ov).sum())
    return {'proper_crossings':proper,'nonadjacent_point_contacts':touch,'collinear_overlaps':overlap}

def main():
    volume=np.load(ROOT/m.PLAN['source'],mmap_mode='r');results=[]
    metrics={r['z']:r for r in json.loads((RUN/'metrics.json').read_text())}
    for z in m.ZS:
        dest=RUN/f'slice_{z}';p=np.load(dest/'path_and_samples.npz');a=p['edge_start_xy'];b=p['edge_end_xy'];l=p['edge_provenance']
        olda,oldb,oldl=flatten(np.load(ROOT/f'outputs/p0_t2/slice_{z}_samples.npz'));origin=p['original_T2_edge_id'];base=p['edge_original_provenance']
        topology=exhaustive_topology(a,b);assert topology['proper_crossings']==0
        assert np.array_equal(b[:-1],a[1:]) and np.array_equal(a[0],olda[0]) and np.array_equal(b[-1],oldb[-1])
        # Every retained CT edge is the exact original edge, possibly reversed; no synthetic observed geometry.
        for i,oid in enumerate(origin):
            if oid<0:assert l[i]==2 and base[i]==2
            else:
                assert (np.array_equal(a[i],olda[oid]) and np.array_equal(b[i],oldb[oid])) or (np.array_equal(a[i],oldb[oid]) and np.array_equal(b[i],olda[oid]))
                assert base[i]==oldl[oid]
                assert l[i]==base[i] or (base[i]==1 and l[i]==3)
        idx=origin[origin>=0];assert len(np.unique(idx))==len(idx)
        lengths=np.linalg.norm(b-a,axis=1);arc=np.r_[0,np.cumsum(lengths)]
        assert np.array_equal(arc,p['edge_s_px']) and np.allclose(np.diff(p['s_um']),63.28,atol=1e-8)
        ss=p['s_px'];xy=p['xy'];prov=p['provenance'];ct=p['ct'];mask=p['mask']
        assert np.array_equal(np.isnan(ct),mask) and np.array_equal(mask,np.isin(prov,[0,2]))
        assert np.array_equal(np.load(dest/'row.npy'),ct,equal_nan=True) and np.array_equal(np.load(dest/'mask.npy'),mask)
        # Check inferred intervals against s, not just trusting stored provenance.
        gap_samples=0
        for i in np.flatnonzero(l==2):
            inside=(ss>arc[i]+1e-8)&(ss<arc[i+1]-1e-8)
            assert np.all(prov[inside]==2) and np.isnan(ct[inside]).all();gap_samples+=int(inside.sum())
        q=xy[~mask];x=np.floor(q[:,0]).astype(int);y=np.floor(q[:,1]).astype(int);dx=q[:,0]-x;dy=q[:,1]-y;im=volume[z-128].astype(float)
        expected=im[y,x]*(1-dx)*(1-dy)+im[y,x+1]*dx*(1-dy)+im[y+1,x]*(1-dx)*dy+im[y+1,x+1]*dx*dy
        assert np.allclose(expected,ct[~mask],rtol=0,atol=1e-10)
        assert len(list(csv.DictReader((dest/'sample_provenance.csv').open())))==len(ss)
        png=np.array(Image.open(dest/'unwrapped_row.png'));mpng=np.array(Image.open(dest/'unwrapped_mask.png'))
        assert np.array_equal(png[0],np.nan_to_num(ct,nan=0).round().astype('uint8')) and np.array_equal(mpng[0],mask.astype('uint8')*255)
        for f in dest.glob('*.png'):
            with Image.open(f) as image:image.verify()
        rec={'z':z,**topology,'lineage':'all retained observed geometry exact, only local1→3 upgrades','gap_interior_samples_checked':gap_samples,
             'total_samples':len(ss),'independent_bilinear':'PASS','continuous_s_and_XY':'PASS','source_T2_edges_retained':len(idx)}
        results.append(rec);print(json.dumps(rec),flush=True)
    for name,h in json.loads((RUN/'preservation.json').read_text()).items():assert m.sha(ROOT/name)==h,name
    m.dump(OUT/'verification.json',{'slices':results,'preservation':'all run03 pre-existing hashes PASS, including T2 evidence/report and README/AGENTS',
                                  'method':'exhaustive orientation/bounding-box pair tests independent of KD-tree crossing detector'})

if __name__=='__main__':main()
