"""Measured frame/placement audit; does not modify any surface or input."""
from pathlib import Path
import json
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t10'

def unit(a):
    return a / np.linalg.norm(a, axis=-1, keepdims=True)

def depth_summary(a):
    mean = a.mean((1,2))
    peak = np.argmax(a.astype(float), axis=0)
    return {'mean_intensity_by_depth': mean.tolist(), 'mean_profile_peak_index': int(mean.argmax()), 'column_peak_at_depth_edge_fraction': float(((peak==0)|(peak==16)).mean()), 'column_peak_within_center2_fraction': float((abs(peak-8)<=2).mean()), 'zero_voxel_fraction':float((a==0).mean()), 'warning':'Intensity peaks are CT occupancy diagnostics, not ink labels or an oracle surface.'}

def main():
    points = np.stack([tifffile.imread(ROOT/f'data/raw/p0_t10/geometry_{axis}.tif') for axis in 'xyz'],axis=-1)
    meta = json.loads((OUT/'official_geometry_metadata.json').read_text())
    selection = json.loads((OUT/'crop_selection.json').read_text())['selected']
    y0,x0=selection['y0'],selection['x0']
    # Grid regions are approximate render-UV correspondence for frame diagnostics,
    # not a fresh CT render. Producer registration and voxel-center offsets are retained.
    gy0=int(y0*4*meta['scale'][1]);gx0=int(x0*4*meta['scale'][0])
    gy1=int((y0+512)*4*meta['scale'][1])+1;gx1=int((x0+512)*4*meta['scale'][0])+1
    region=points[gy0-1:gy1+1,gx0-1:gx1+1].astype(float)
    du=(region[1:-1,2:]-region[1:-1,:-2])/2
    dv=(region[2:,1:-1]-region[:-2,1:-1])/2
    valid=(region[1:-1,1:-1]>0).all(-1)&(np.linalg.norm(du,axis=-1)>0)&(np.linalg.norm(dv,axis=-1)>0)
    normal=unit(np.cross(du,dv))[valid];tangent=unit(du)[valid];longitudinal=unit(dv)[valid]
    representative=unit(normal.mean(0));variation=np.degrees(np.arccos(np.clip(normal@representative,-1,1)))
    frame={'grid_region_yx':[gy0,gy1,gx0,gx1], 'valid_grid_points':int(valid.sum()),'surface_normal_abs_z_quantiles':np.quantile(abs(normal[:,2]),[0,.5,.95,1]).tolist(),'U_tangent_abs_z_quantiles':np.quantile(abs(tangent[:,2]),[0,.5,.95,1]).tolist(),'V_deviation_from_scanZ_degrees_quantiles':np.quantile(np.degrees(np.arccos(np.clip(abs(longitudinal[:,2]),0,1))),[0,.5,.95,1]).tolist(),'normal_variation_about_mean_degrees_quantiles':np.quantile(variation,[0,.5,.95,1]).tolist(),'convention':'geometric diagnostic N=dP/dU cross dP/dV; delivered depth sign preserved, producer flip flag not recorded in source attrs','coordinate_caveat':'approximate grid bounds from documented scale; not an independent subpixel registration or rerender'}
    comparison={'official_surface':{'point':'producer tifxyz curved XYZ map in2.4um scan coordinates','tangent_frame':'varies over both U and V; measured finite differences of trusted grid','normal':'full3D cross(U,V); delivered render sign used unchanged','depth':'producer109planes, nominal center54; source[13,97) pooled4; selected pooled[2,19)','model_depth_center':'pooled10 averages source53..56; nominal center54.5, +1.2um relative to plane54','mapping':'XYmean pyramidlevel2 + depthmean4; metadata2.4um ->9.6um','interpolation':'consume exact producer-rendered CT; current reference renderer uses trilinear CT sampling; production interpolation/flip flags not recorded','normalization':'same official robust_mad/tifxyz_robust function as T9','patch':[17,128,128],'sampling':'contiguous curved-surface pixels at9.6um; stride64, Hann blend'},'T9':{'point':'independently traced transverse CT contour on126.56um grid, mapped native=16*level4','tangent_frame':'single fixed XYtangent per patch; other in-plane axis exactly scanZ','normal':'saved(-ty,tx,0), constant per patch; no out-of-plane normal component','depth':'17symmetric offsets -8..8; center exactly saved trace point','model_depth_center':'index8 at frozen trace; actual papyrus seating uncalibrated','mapping':'native7.91um trilinear sample step9.6/7.91; no native high-resolution averaging','interpolation':'scipy order1 float32; half-even round to uint8','normalization':'same official robust_mad/tifxyz_robust function as control','patch':[17,128,128],'sampling':'one center probability per63.28um s and126.56um scanZ; no cross-Z alignment'},'measured_official_frame':frame,'confounders':['different scan: control78keV2.4um vs T9 54keV7.91um','training-corpus control, not generalization','control locations do not coincide with T9 samples','T9 sparse raster and independently traced paths vs contiguous surface coordinates'],'causal_limit':'These differences implicate geometric representation if control is positive; they do not isolate a sole cause or establish a repair.'}
    metrics=json.loads((OUT/'official_control_metrics.json').read_text())
    if metrics['status']=='PASS':
        volume=np.load(OUT/'official_input.npy');audit=[]
        for z,start in zip([128,152,177],[128,192,256]):
            with np.load(ROOT/f'outputs/p0_t9/smoke_{z}.npz') as archive:thread=archive['raw_patch'].copy()
            control=volume[:,start:start+128,start:start+128]
            audit.append({'T9_z':z,'control_origin_yx':[start,start],'control':depth_summary(control),'T9':depth_summary(thread)})
        comparison['fixed_patch_audit']=audit
    (OUT/'control_vs_t9_geometry.json').write_text(json.dumps(comparison,indent=2,allow_nan=False)+'\n')
    fig,axs=plt.subplots(1,3,figsize=(14,4))
    axs[0].hist(abs(normal[:,2]),bins=30);axs[0].axvline(0,color='red',label='T9 fixed nz=0');axs[0].legend();axs[0].set(xlabel='Absolute normal scan-Z component',title='Trusted geometry (different locations)')
    axs[1].hist(np.degrees(np.arccos(np.clip(abs(longitudinal[:,2]),0,1))),bins=30);axs[1].axvline(0,color='red');axs[1].set(xlabel='V tilt from scan Z (degrees)',title='T9 V axis is fixed scan Z')
    if 'fixed_patch_audit' in comparison:
        for row in comparison['fixed_patch_audit']:
            axs[2].plot(range(17),row['control']['mean_intensity_by_depth'],color='black',alpha=.5)
            axs[2].plot(range(17),row['T9']['mean_intensity_by_depth'],color='tab:orange',alpha=.5)
    axs[2].axvline(8,color='gray',ls=':');axs[2].set(xlabel='Model depth index',ylabel='Mean raw CT',title='Fixed patches: control black / T9 orange')
    fig.tight_layout();fig.savefig(OUT/'control_vs_t9_sampling.png',dpi=140);plt.close(fig)
    print(json.dumps(frame,indent=2),flush=True)

if __name__=='__main__':
    main()
