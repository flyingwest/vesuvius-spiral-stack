"""Exploratory, manually guided 2-D ridge tracing on one existing CT slice."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter, map_coordinates

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_1'
GUIDE = np.array([[246,111],[258,120],[271,131],[281,141],[281,147],
                 [270,155],[245,159],[210,158],[183,154],[175,150],
                 [187,141],[211,132],[229,122],[241,114]], dtype=float)
OFFSETS = np.arange(-3,3.01,.25)


def sample(image, xy):
    return map_coordinates(image.astype(float), [xy[:,1],xy[:,0]],
                           order=1, mode='nearest', prefilter=False)


def geometry():
    s = np.r_[0,np.cumsum(np.linalg.norm(np.diff(GUIDE,axis=0),axis=1))]
    u = np.linspace(0,s[-1],int(np.ceil(s[-1]/.5))+1)
    spline = CubicSpline(s,GUIDE,axis=0,bc_type='natural')
    p = spline(u)
    tangent = spline(u,1)
    tangent /= np.linalg.norm(tangent,axis=1)[:,None]
    normal = np.c_[-tangent[:,1],tangent[:,0]]
    return p,normal


def trace(image, guide, normal, sigma=.6, shift=0., penalty=.01):
    response = gaussian_filter(image.astype(float),sigma)-gaussian_filter(image.astype(float),2.)
    candidates = guide[:,None,:]+normal[:,None,:]*(OFFSETS[None,:,None]+shift)
    values = sample(response,candidates.reshape(-1,2)).reshape(candidates.shape[:2])
    unary = -values/255 + .001*OFFSETS[None,:]**2
    score=unary[0].copy()
    back=np.zeros(values.shape,dtype=np.int32)
    for k in range(1,len(guide)):
        delta=OFFSETS[:,None]-OFFSETS[None,:]
        distance=np.linalg.norm(candidates[k-1,:,None,:]-candidates[k,None,:,:],axis=2)
        cost=score[:,None]+penalty*delta**2
        cost[(np.abs(delta)>.5)|(distance>1.25)]=np.inf
        back[k]=np.argmin(cost,axis=0)
        score=cost[back[k],np.arange(len(OFFSETS))]+unary[k]
    assert np.isfinite(score.min())
    states=np.zeros(len(guide),dtype=int);states[-1]=np.argmin(score)
    for k in range(len(guide)-1,0,-1):states[k-1]=back[k,states[k]]
    return candidates[np.arange(len(guide)),states], response, states


def resample(path):
    cumulative=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
    assert np.all(np.diff(cumulative)>0)
    s=np.linspace(0,cumulative[-1],int(np.ceil(cumulative[-1]/.5))+1)
    xy=np.column_stack([np.interp(s,cumulative,path[:,d]) for d in range(2)])
    return s,xy


def main():
    assert not (OUT/'metrics.json').exists(), 'Preserve prior runs, do not overwrite'
    image=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')[16].copy()
    assert np.array_equal(image,np.asarray(Image.open(OUT/'source_slice.png')))
    plan={'slice_local_z':16,'slice_global_level4_z':144,'native_nominal_z':2304,
          'guide_xy':GUIDE.tolist(),'guide_step_max_px':.5,'offsets_px':OFFSETS.tolist(),
          'sigma_px':.6,'background_sigma_px':2.,'guide_offset_penalty':.001,
          'transition_penalty':.01,'max_offset_change_px':.5,'max_step_px':1.25,
          'row_step_max_px':.5,'scope':'Exploratory, open near-turn candidate; no forced closure or gap repair.',
          'variants':['sigma_0.4','sigma_0.8','guide_shift_-0.5','guide_shift_+0.5'],
          'interpretation':'Manual guide and mathematical continuity do not establish one physical sheet.'}
    (OUT/'trace_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    guide,normal=geometry()
    path,response,states=trace(image,guide,normal)
    s,xy=resample(path);intensity=sample(image,xy)
    row=np.clip(np.rint(intensity),0,255).astype(np.uint8)[None,:]
    Image.fromarray(row).save(OUT/'unwrapped_row.png')
    Image.fromarray(np.repeat(row,48,axis=0)).save(OUT/'unwrapped_row_preview.png')
    np.save(OUT/'sampling_values.npy',intensity)
    np.save(OUT/'optimizer_path_xy.npy',path)
    np.save(OUT/'ridge_response.npy',response)
    with (OUT/'trace_coordinates.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['index','s_pixels','s_mm','x_level4','y_level4','z_level4','x_native_nominal','y_native_nominal','z_native_nominal','ct_bilinear'])
        for k,((x,y),arc,value) in enumerate(zip(xy,s,intensity)):
            writer.writerow([k,arc,arc*.12656,x,y,144,x*16,y*16,2304,value])
    fig,ax=plt.subplots(figsize=(10,10));ax.imshow(image,cmap='gray',vmin=0,vmax=255)
    ax.plot(xy[:,0],xy[:,1],color='red',lw=.7);ax.plot(*xy[0],'go',ms=4);ax.plot(*xy[-1],'bo',ms=4)
    ax.set_title('P0.1 candidate, Z=144 level 4; green start / blue end; NOT validated sheet identity')
    fig.tight_layout();fig.savefig(OUT/'traced_slice.png',dpi=160);plt.close(fig)
    variants={}
    fig,axes=plt.subplots(2,2,figsize=(16,10))
    for ax in axes.ravel():
        ax.imshow(image,cmap='gray',vmin=0,vmax=255,interpolation='nearest')
        ax.set_xlim(165,295);ax.set_ylim(170,100)
    axes[0,0].plot(guide[:,0],guide[:,1],'y--',lw=.8)
    axes[0,0].scatter(GUIDE[:,0],GUIDE[:,1],s=12,c='yellow')
    axes[0,0].set_title('Manual guide, before optimization')
    axes[0,1].plot(path[:,0],path[:,1],'r-',lw=.8);axes[0,1].set_title('Baseline candidate')
    axes[1,0].plot(path[:,0],path[:,1],'r-',lw=.8)
    for name,kwargs in [('sigma_0.4',{'sigma':.4}),('sigma_0.8',{'sigma':.8}),
                        ('guide_shift_-0.5',{'shift':-.5}),('guide_shift_+0.5',{'shift':.5})]:
        q,_,_=trace(image,guide,normal,**kwargs)
        np.save(OUT/(name+'_path_xy.npy'),q)
        d=np.linalg.norm(q-path,axis=1)
        variants[name]={'median_displacement_px':float(np.median(d)),
                        'p95_displacement_px':float(np.percentile(d,95)),
                        'max_displacement_px':float(d.max()),
                        'ct_same_guide_station_rmse':float(np.sqrt(np.mean((sample(image,q)-sample(image,path))**2)))}
        axes[1,0].plot(q[:,0],q[:,1],lw=.65,label=name)
    axes[1,0].legend(fontsize=7);axes[1,0].set_title('All planned sensitivity variants (not candidates selected for success)')
    axes[1,1].imshow(response,cmap='coolwarm',vmin=-40,vmax=40,alpha=.8)
    axes[1,1].set_title('Gaussian 0.6 minus Gaussian 2.0 response')
    fig.tight_layout();fig.savefig(OUT/'trace_audit.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,3));ax.plot(s,intensity,lw=.7);ax.set(xlabel='arc length (level-4 pixels)',ylabel='CT uint8 units',ylim=(0,255));fig.tight_layout();fig.savefig(OUT/'intensity_profile.png',dpi=150);plt.close(fig)
    # Independent four-corner interpolation verification, not scipy's sampler.
    x,y=xy.T;x0=np.floor(x).astype(int);y0=np.floor(y).astype(int);dx=x-x0;dy=y-y0
    independent=(image[y0,x0]*(1-dx)*(1-dy)+image[y0,x0+1]*dx*(1-dy)+image[y0+1,x0]*(1-dx)*dy+image[y0+1,x0+1]*dx*dy)
    repeat,_,_=trace(image,guide,normal)
    angle=np.unwrap(np.arctan2(path[:,1]-140,path[:,0]-240))
    metrics={'trace_points':len(xy),'optimizer_points':len(path),'arc_length_px':float(s[-1]),
             'arc_length_mm':float(s[-1]*.12656),'row_width':row.shape[1],'row_height':1,
             'row_step_px':float(s[1]),'max_optimizer_step_px':float(np.linalg.norm(np.diff(path,axis=0),axis=1).max()),
             'endpoint_distance_px':float(np.linalg.norm(path[-1]-path[0])),
             'angular_sweep_degrees_about_240_140':float(np.degrees(angle[-1]-angle[0])),
             'corridor_boundary_fraction':float(np.mean((states==0)|(states==len(OFFSETS)-1))),
             'intensity_min_mean_max':[float(intensity.min()),float(intensity.mean()),float(intensity.max())],
             'variants':variants,'deterministic_repeat_identical':bool(np.array_equal(path,repeat)),
             'independent_bilinear_max_error':float(np.abs(independent-intensity).max()),
             'source_slice_sha256':hashlib.sha256(image.tobytes()).hexdigest(),
             'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    assert metrics['deterministic_repeat_identical'] and metrics['independent_bilinear_max_error']<1e-10
    (OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':
    main()
