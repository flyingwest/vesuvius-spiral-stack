"""Bounded predefined summaries of the exact saved T6 normal profiles."""
import json,shutil
from pathlib import Path
import numpy as np
from p0_t6_features import ROOT,sha,dump,table,ZS,DS,SPACING
OUT=ROOT/'outputs/p0_t7';T6=ROOT/'outputs/p0_t6';T5=ROOT/'outputs/p0_t5'
WIDTHS={'narrow':np.array([2,3,4]),'medium':np.array([1,2,3,4,5]),'wide':np.arange(7)}
OFF=np.array([-4,-2,-1,0,1,2,4.])
FAMILIES=['mean','median','std','range','center_minus_mean','asymmetry','gradient_abs','curvature','ridge_proxy']
FEATURES={'center_intensity':('P(0)',(0,255))}
SCALES={'center_intensity':255.}
for width,ids in WIDTHS.items():
 radius=float(max(abs(OFF[ids])))
 specs={'mean':('mean of selected offsets',(0,255),255.),'median':('median of selected offsets',(0,255),255.),
 'std':('population standard deviation',(0,127.5),127.5),'range':('max minus min',(0,255),255.),
 'center_minus_mean':('P0 minus selected profile mean',(-255,255),255.),
 'asymmetry':('absolute difference of mean positive versus negative offset values',(0,255),255.),
 'gradient_abs':('abs(P(+R)-P(-R))/(2R)',(0,255/(2*radius)),255/(2*radius)),
 'curvature':('(2P0-P(-R)-P(+R))/R^2',(-510/radius**2,510/radius**2),510/radius**2),
 'ridge_proxy':('P0 minus mean selected points with abs(offset)>=R/2',(-255,255),255.)}
 for family in FAMILIES:
  desc,limits,scale=specs[family];FEATURES[f'{family}_{width}']=(desc+f'; width={width};R={radius}px',limits);SCALES[f'{family}_{width}']=scale

def compute(profile):
 result={'center_intensity':profile[:,3]}
 for width,ids in WIDTHS.items():
  p=profile[:,ids];off=OFF[ids];radius=max(abs(off));mean=p.mean(1)
  values={'mean':mean,'median':np.median(p,axis=1),'std':p.std(1),'range':np.ptp(p,axis=1),
   'center_minus_mean':profile[:,3]-mean,'asymmetry':np.abs(p[:,off>0].mean(1)-p[:,off<0].mean(1)),
   'gradient_abs':np.abs(p[:,-1]-p[:,0])/(2*radius),'curvature':(2*profile[:,3]-p[:,0]-p[:,-1])/radius**2,
   'ridge_proxy':profile[:,3]-p[:,abs(off)>=radius/2].mean(1)}
  for name,v in values.items():result[f'{name}_{width}']=v
 return result

def initialize():
 OUT.mkdir(exist_ok=False)
 for line in (T6/'SHA256SUMS').read_text().splitlines():
  h,p=line.split('  ',1);assert sha(ROOT/p)==h,p
 prior=[p for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents]+[ROOT/'README.md',ROOT/'AGENTS.md']
 dump(OUT/'preservation.json',{str(p.relative_to(ROOT)):sha(p) for p in prior})
 dump(OUT/'plan.json',{'experiment':'P0-T7 bounded deterministic profile summaries','source':'exact outputs/p0_t6/slice_<z>_profiles.npz; no CT queries or geometry modifications',
  'slices':ZS,'ds_um':DS,'spacing_um':SPACING,'widths':{name:OFF[ids].tolist() for name,ids in WIDTHS.items()},
  'scalar_count':len(FEATURES),'features':{name:{'definition':definition,'display_range':limits,'fixed_normalization_scale':SCALES[name]} for name,(definition,limits) in FEATURES.items()},
  'evaluation':'unchanged T6 evaluator:49 zero-shift pairs,1001 windows,±64 diagnostic lags,199 global and19 local circular controls,paired block bootstrap5pairs/2000 seed2606',
  'regions':'five fixed starting-Z groups128..137,138..147,148..157,158..167,168..176; pair assigned by first Z',
  'joint_selection':'equal-weight rank percentiles of mean adjacent r,mean local zero r,local peak control-exceedance fraction,number of5regions improving BOTH adjacent/local versus mean_wide; tie by name',
  'baseline_in_ranking':'mean_wide remains eligible; center reference excluded from candidate ranking; do not claim winner if no improvement',
  'combinations':'after scalar evaluation, choose top3 distinct normalized-signal candidates by joint rank; test only their3 unordered pair means, each0.5*(A/scaleA+B/scaleB); no fitted weights or normalization from correlations',
  'broad_gain_guide':'relative to mean_wide: adjacent gain>=0.02,local gain>=0.01,at least30/49 pairs improve,at least4/5 regions improve BOTH,local peak-control fraction not worse; descriptive operational guide fixed before evidence',
  'classification':'A broad guide satisfied with reproducibility; B smaller consistent gains but guide not met; C no meaningful joint improvement beyond T6; disclose tradeoffs rather than single-metric winner',
  'limits':'28 scalar stacks +<=3 actual combinations; no iterative expansions; no learning,OCR,letters,ink claims,geometry changes'})
 dest=OUT/'scalars';dest.mkdir();mask=np.load(T6/'mask.npy');stacks={name:np.full(mask.shape,np.nan) for name in FEATURES}
 for i,z in enumerate(ZS):
  p=np.load(T6/f'slice_{z}_profiles.npz');ii=p['sample_indices']
  for name,values in compute(p['profiles']).items():stacks[name][i,ii]=values
 for name,a in stacks.items():
  assert np.array_equal(np.isnan(a),mask);np.save(dest/f'{name}_stack.npy',a)
 assert np.array_equal(stacks['center_intensity'],np.load(T5/'stack.npy'),equal_nan=True)
 assert np.array_equal(stacks['mean_wide'],np.load(T6/'profile_mean_stack.npy'),equal_nan=True)
 for folder in [OUT,dest]:
  for name in ['mask','provenance','coverage','s_um','z_indices','z_relative_um']:shutil.copyfile(T6/f'{name}.npy',folder/f'{name}.npy')
 dump(dest/'FEATURES_FINALIZED.json',{str(p.relative_to(dest)):sha(p) for p in sorted(dest.glob('*.npy'))})
 print('28 predefined scalars finalized; both baselines bit-identical; no source data or geometry changed')
if __name__=='__main__':initialize()
