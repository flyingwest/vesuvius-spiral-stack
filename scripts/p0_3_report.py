"""Present the frozen P0.3 rasters and measurements without changing their samples."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_3'
def read(name):return json.loads((OUT/name).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,value):p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    metrics=read('metrics.json');before={p.name:sha(p) for p in OUT.iterdir() if p.is_file()}
    cmap=plt.get_cmap('gray').copy();cmap.set_bad('magenta')
    # Equal physical scale along both axes; display only, no array resampling.
    fig,axes=plt.subplots(1,2,figsize=(9,8))
    for ax,fid in zip(axes,metrics):
        a=np.load(OUT/f'{fid}_stack.npy');s=np.load(OUT/f'{fid}_s_px.npy')*126.56;zs=np.load(OUT/f'{fid}_z.npy')
        zpos=(zs-zs[0])*126.56
        ax.imshow(np.ma.masked_invalid(a),cmap=cmap,vmin=0,vmax=255,interpolation='nearest',aspect='equal',
            extent=[s[0]-15.82,s[-1]+15.82,zpos[-1]+63.28,-63.28])
        ax.set(title=f'{fid}: true physical aspect\nZ={zs[0]}–{zs[-1]}',xlabel='Signed arc length (µm)',ylabel='Longitudinal offset (µm)')
    fig.suptitle('P0.3 raw CT, fixed 0–255; missing = magenta; no filled gaps')
    fig.tight_layout();fig.savefig(OUT/'physical_aspect.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    a=np.load(OUT/'F02_stack.npy');s=np.load(OUT/'F02_s_px.npy')*126.56;zs=np.load(OUT/'F02_z.npy')
    for i,z in enumerate(zs):axes[0].plot(s,a[i],label=str(z),lw=1)
    axes[0].set(xlabel='s (µm)',ylabel='Raw CT units',title='F02 rows: no normalization or fitted shifts');axes[0].legend(ncol=2,fontsize=7)
    pairs=read('F02_pairs.json');values=[p['NCC'] for p in pairs]
    axes[1].bar([str(p['z_source']) for p in pairs],values,color=['tab:red' if v<.2 else 'tab:blue' for v in values])
    axes[1].axhline(0,color='black');axes[1].set(xlabel='First Z of adjacent pair',ylabel='Pearson NCC',ylim=(-1,1),title='F02 local coherence (all 13 pairs)');axes[1].tick_params(axis='x',rotation=60)
    fig.tight_layout();fig.savefig(OUT/'F02_coherence_diagnostic.png',dpi=150);plt.close(fig)
    volume=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')
    xy=np.load(OUT/'F02_sample_xy.npy');mask=np.load(OUT/'F02_mask.npy');meta=read('F02_meta.json')
    fig,axes=plt.subplots(1,4,figsize=(14,4))
    for ax,z in zip(axes,[139,140,146,147]):
        i=list(zs).index(z);points=xy[i,mask[i]];r=meta[i];center=points.mean(axis=0)
        ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
        ax.plot(*points.T,'r.-',lw=.8,ms=3)
        ax.plot([r['start_x'],r['end_x']],[r['start_y'],r['end_y']],'co',ms=4,fillstyle='none')
        ax.set(xlim=(center[0]-4,center[0]+4),ylim=(center[1]+4,center[1]-4),title=f'F02 Z={z}\nred = retained-core samples',xlabel='X',ylabel='Y')
    fig.tight_layout();fig.savefig(OUT/'F02_low_coherence_CT.png',dpi=150);plt.close(fig)
    verdict={'P0_3':'PARTIAL_LOCAL_DEMONSTRATION','F02':metrics['F02']['status'],'F03':metrics['F03']['status'],
        'both_fragments_coherence_demonstrated':False,'full_P0_PASS_claimed':False,
        'CT_or_geometry_only':True,'general_3D_mesh':False,'missing_data_filled':False,
        'suspect_low_coherence_pairs_F02':[[p['z_source'],p['z_target']] for p in pairs if p['NCC']<.2],
        'confirmed_misregistration_count':None,'reason':'Shift probe has no eligible support; low coherence is a warning, not confirmed layer switching.',
        'next_step_not_executed':'Independent geometric-correspondence validation for low-coherence F02 pairs; further F03 support validation before enlarging it.'}
    save(OUT/'verdict.json',verdict)
    table='\n'.join(f"| {p['z_source']}→{p['z_target']} | {p['samples']} | {p['NCC']:.4f} | {p['mean_absolute_difference_CT']:.3f} |" for p in pairs)
    provenance={'source_volume':{'path':'outputs/orientation/ct_level4_zyx.npy','sha256':sha(ROOT/'outputs/orientation/ct_level4_zyx.npy')},
        'source_track_files':{name:sha(ROOT/name) for name in ['outputs/p0_2/run_02/fragments.json','outputs/p0_2/run_02/tracks.json',
            'outputs/p0_2/validation/F02_measurements.json','outputs/p0_2/validation/F03_measurements.json','docs/p0_2_validation.md']},
        'spacing_um':{'s':31.64,'z':126.56},'numeric_dtype':'float64 CT and XY; bool valid mask; integer Z',
        'mask_convention':'True / PNG255 = valid; False / PNG0 = missing',
        'PNG_CT_encoding':'uint8, nearest integer via np.rint, clipped [0,255], invalid pixels zero; use accompanying mask',
        'scripts':{str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'scripts').glob('p0_3*.py')}}
    save(OUT/'provenance.json',provenance)
    report=f'''# P0.3 — first Direct Spiral-Stack rasters

## Result

**PARTIAL LOCAL DEMONSTRATION.** Separate arc-length/longitudinal CT rasters
were produced for the strictly validated F02 and F03 cores. F02 passes the
disclosed local-coherence diagnostic: mean adjacent-row NCC **0.5495**, with
**12/13 positive pairs**, versus a row-order permutation mean **0.1087**.
This is evidence of spatially coherent local structure, not independent/random
rows. F03's allowed core is only about half a source pixel long: its raster is
validly sampled, but **coherence is NOT EVALUABLE**. Do not claim that both
fragments pass or that the full frozen P0 research hypothesis has been established.

F02 has a negative-correlation transition at Z=139→140 and weak transitions at
146→147 and147→148. These are preserved warning cases. No hidden surface
reconstruction, width normalization, intensity-driven warping or gap filling
was needed to produce the rasters. Missing F03 Z=146 remains entirely masked.

## Hypothesis, scope and exact inputs

Hypothesis: geometrically corresponding fragment rows, sampled at a fixed
physical arc interval, preserve local CT structure when ordered across Z.
This is an exploratory demonstration on previously selected tracks, not an
independent held-out validation or a search for readable text.

- PHercParis4 / Scroll1, `54keV_7.91um_Scroll1A.zarr/4`.
- Remote array reference:
  `https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4`.
- Exact source: `outputs/orientation/ct_level4_zyx.npy`, uint8 ZYX=(128,493,506),
  global Z offset128. X columns/right, Y rows/down, Z approximately longitudinal.
- Physical voxel spacing: **126.56 µm isotropic** at level4; native7.91 µm.
- F02: **T0077**, Z=137–150 inclusive (local9–22), 14 observed slices.
- F03: **T0243**, Z=140–148 inclusive (local12–20), eight observations plus
  one missing row at146. The accepted145→147 identity link does not supply a
  path or intensity observation at146.
- Track coordinates: `outputs/p0_2/run_02/fragments.json` and `tracks.json`.
  Validation inputs: `docs/p0_2_validation.md` and
  `outputs/p0_2/validation/F02_measurements.json`, `F03_measurements.json`.
- `outputs/p0_3/provenance.json` binds exact volume, track and validation hashes.
  `plan.json` was saved before computation. No new data were downloaded.

### Strict interpretation of validated support

The validation report explicitly excludes uncertain varying tails. The only
unambiguous support common to **all** observations is F02 window-u=[1,3.5]px
(six vertices) and F03 u=[1,1.5]px (two vertices). Only those portions are
sampled. These u coordinates are the earlier local tangent-window coordinates,
not arc length. No old entire F02/F03 trace, unvalidated tail or other fragment
enters a raster. This intentionally sacrifices width; readiness for an experiment
does not guarantee that its permitted support can establish coherence.

## Arc length, interpolation and row alignment

For each retained ordered polyline p_i=(x_i,y_i), compute
`l_0=0`, `l_i=sum_j ||p_j-p_(j-1)||` in level4 pixel units. Locate a common
geometric transverse station `u_anchor=(u_min+u_max)/2`, namely2.25 for F02
and1.25 for F03. Piecewise-linear interpolation along the u-monotone polyline
gives its local arc coordinate `l_anchor`. Define **s=l−l_anchor**. This is a
per-row translation, not scaling. It does not assume equal fragment lengths.
The direction of increasing s is increasing the frozen window u; F03 can
therefore run toward decreasing global X.

Use global grid `s=k*0.25px`, covering all retained row endpoints, rounded
outward to grid boundaries. Thus **Δs=31.64µm**, **Δz=126.56µm**.
The four-times-finer in-plane sampling is oversampling, **not increased CT
resolution**. For every in-support grid point interpolate XY linearly in true
polyline arc length, then sample the original CT via SciPy `map_coordinates`
with **order=1, prefilter=False**, float64 input/output. This is bilinear XY
sampling within the same observed Z plane. No smoothing, row normalization,
learned signal, ink inference or Z interpolation is applied.

Coordinates outside each row's retained endpoints and all absent observations
are NaN with mask=false. No extension or extrapolation is sampled. The alignment
CSV contains exact input node, Z/localZ, endpoints, physical length, signed
offset added to local arc length (`−l_anchor`) and its change from the adjacent
observed row. The common anchor always lies at output s=0. These offsets describe
the geometric coordinate origin, not a CT-optimized translation or proof of
material-point registration. Absolute correspondence away from the anchor is
an approximation, diagnosed below.

## Raster dimensions, lengths and missing data

| Quantity | F02 | F03 |
| --- | ---: | ---: |
| Numeric array shape (Z,s) | 14×15 | 9×5 |
| PNG dimensions (width×height) | 15×14 | 5×9 |
| s sample centers, source px | −1.75…+1.75 | −0.5…+0.5 |
| Retained arc length min–max, px | 2.637613–3.043911 | 0.500010–0.536928 |
| Retained arc length min–max, µm | 333.816–385.237 | 63.281–67.954 |
| Mean / SD row length, px | 2.837154 / 0.133249 | 0.511094 / 0.010860 |
| Valid samples per observed row | 11–12 | 3 |
| Valid / total samples | 158/210 | 24/45 |
| Total missing fraction | 24.762% | 46.667% |
| Missing-observation samples | 0 | 5 (entire Z=146 row) |
| Outside-core/padding samples | 52 | 16 |

Grid-boundary padding contributes to missing fraction; it is not papyrus damage.
Even F02's complete Z coverage has masked edges. Endpoint markers show exact
retained endpoints, which need not coincide with sampling centers. F03's three
valid samples represent only half a source pixel of width, not three independent
measurements.

The raw PNGs have one image pixel per numeric sample, with no row resizing:
fixed0–255 display, nearest-integer uint8 quantization for previews only. Invalid
PNG pixels are0; **always read the mask**, where255 means valid and0 missing.
Lossless float64 NPY arrays retain fractional intensity and NaNs. A square-pixel
viewer does not know Δz/Δs=4; `physical_aspect.png` shows equal physical scales.
Diagnostic panels are explicitly enlarged/labeled and do not modify the arrays.

## Coherence measurements and negative results

The pre-output exploratory diagnostic requires ≥9 shared samples spanning ≥2
source pixels, each row standard deviation≥1 CT unit; at least75% of observed
adjacent pairs must be evaluable. Local coherence passes if mean Pearson NCC≥0.5
and at least75% of evaluable pairs are positive. These practical thresholds are
disclosed before sampling, not trained or selected after seeing the raster.
NCC removes means only for comparison; the stored rows are unnormalized.

All13 F02 adjacent pairs are eligible. Mean NCC=0.549468, positive fraction
12/13=92.308%, mean absolute adjacent CT difference=15.330456. Per-pair values:

| Z pair | Shared samples | NCC | Mean absolute CT difference |
| --- | ---: | ---: | ---: |
{table}

Preserved2000 row-order permutations (NumPy generator seed302) have mean
NCC0.108669, central95% interval[−0.139662,0.370143]; **0/2000** equal or exceed
the observed adjacent mean. This is descriptive evidence that Z order matters,
not p=0 or an independent biological accuracy test. Permuted masks/rows retain
the same data; spatial interpolation and neighboring slices are correlated.

F03 has6 adjacent observed pairs, **0 eligible**. Its three-sample descriptive
correlations are saved but must not be interpreted as adequate spatial evidence.
Average absolute CT change is10.268725; this alone cannot establish coherence.
No correlation is computed as a valid adjacent match across missing146, and
the145→147 identity link is not substituted for two missing adjacent intensity
observations. The F03 result is insufficient support, **not proof of incoherence**.

## Alignment, shear/stretch and possible misregistration

- Maximum adjacent geometric offset change: F02 **0.254038px /32.151µm**;
  F03 **0.010213px /1.293µm**. Exact signed changes appear in the CSVs; entries
  touching the F03 missing row are blank, not zero.
- Comparing homologous u stations after origin alignment gives maximum residual
  arc-coordinate mismatch0.254038px for F02 and0.010213px for F03. Thus pure
  translations do not achieve exact material-coordinate registration.
- Ratios of corresponding small segment lengths across adjacent rows span
  **0.6871–1.3804** for F02 and **0.9620–1.0089** for F03. They reflect changing
  detected geometry and possible local shear, not an applied stretch correction
  or measured tissue strain. No scale factor or elastic warp is applied.
- A diagnostic-only shift scan at−0.5,−0.25,0,+0.25,+0.5px uses the same valid
  support at every shift to avoid overlap selection bias. **Neither fragment
  retains enough support for an eligible scan**. Consequently zero scan flags
  does not mean no misregistration: CT-based optimal offsets are unresolved.
  No candidate shift changes the raster.
- F02 Z139→140 has NCC−0.4205; Z146→147 is0.0324 and147→148 is0.1735.
  These are local texture/correspondence warnings, not silently removed rows.
  `F02_low_coherence_CT.png` shows the exact retained sampling cores for the
  first two transitions; their positions remain on the validated local ridge,
  but coarse CT variability and along-fragment misregistration cannot be separated
  without additional evidence. No new physical layer switch is established.

Visually, F02 has small persistent local intensity features and changes in mean
row brightness; it is not a smooth continuous band everywhere. F03 is a narrow
stripe with an explicit empty row. These appearances are not the pass criterion;
the measured F02 order-dependent correlation provides the quantitative evidence.

## Outputs, commands and verification

All files are under `outputs/p0_3/`. Required rasters/masks/alignments:

- `F02_stack.png`, `F02_stack_mask.png`, `F02_alignment.csv`.
- `F03_stack.png`, `F03_stack_mask.png`, `F03_alignment.csv`.

For each Fxx, raw `Fxx_stack.npy` (float64), `Fxx_mask.npy` (bool),
`Fxx_sample_xy.npy` (float64), `Fxx_s_px.npy`, `Fxx_z.npy` retain samples,
coordinates and axes. Multiply s_px by126.56 to get µm; use Z differences for
physical longitudinal distances. `Fxx_meta.json`, `Fxx_pairs.json`,
`Fxx_geometry_residuals.json`, `Fxx_shift_scan.json`, `Fxx_summary.json` preserve
all quantitative results. `Fxx_permutation_null.npy` retains permutation means
(empty for unevaluable F03).

Diagnostics: `F02_diagnostic.png`, `F03_diagnostic.png` show valid CT, explicit
mask, geometric offsets and endpoints; `physical_aspect.png`,
`F02_coherence_diagnostic.png`, `F02_low_coherence_CT.png` retain display and
failure diagnostics. `metrics.json`, `verdict.json`, `plan.json`,
`provenance.json`, `preservation_before.json`, `verification.json`,
`artifact_verification.json`, `SHA256SUMS` provide provenance and checks.

Exact commands executed from `/home/wpotech/deployments/vesuvius-spiral-stack`:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p03 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_3_first_stack.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p03 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_3_first_stack.py --verify
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p03 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_3_report.py
```

Generation refuses to overwrite `outputs/p0_3/`. Use `--verify` to recompute
without replacing evidence. Checks cover exact NPY/JSON replay; an independent
four-pixel bilinear reference; a known bent-polyline arc-length calculation;
mask/NaN agreement; consecutive Z indices; missing146; PNG/numeric consistency;
and every preserved prior-file hash. No packages were installed. No existing
README.md, AGENTS.md, prior report, CT volume, P0.2 track or validation artifact
was modified.

## Conclusion and next experiment (not performed)

F02 demonstrates the narrow Direct Spiral-Stack chain from transverse fragment
to arc-length row to a coherent local 2-D CT raster, without a surface mesh.
There is **no uniform two-fragment success claim**: F03 remains too small to
evaluate and several F02 transitions require better correspondence evidence.
This establishes a local geometric/sampling proof, not ink, text or full-scroll
reconstruction. The original broader P0.2 failures remain unchanged.

Next recommended experiment: independently test material-coordinate alignment
around F02's low-coherence transitions; separately validate additional F03
support before considering any expansion. Do not fill its missing row or enlarge
its unvalidated tails to obtain a favorable result. **Stop after P0.3.**
'''
    path=ROOT/'docs/p0_3_first_stack.md';assert not path.exists();path.write_text(report)
    # Original evidence and previews must remain unmodified by this presentation step.
    for name,digest in before.items():assert sha(OUT/name)==digest,name
    for name,digest in read('preservation_before.json').items():assert sha(ROOT/name)==digest,name
    for fid in metrics:
        a=np.load(OUT/f'{fid}_stack.npy');m=np.load(OUT/f'{fid}_mask.npy');p=np.asarray(Image.open(OUT/f'{fid}_stack.png'))
        pm=np.asarray(Image.open(OUT/f'{fid}_stack_mask.png'))
        assert p.shape==a.shape and np.array_equal(pm==255,m)
        np.testing.assert_array_equal(p[m],np.rint(a[m]).astype(np.uint8));assert np.all(p[~m]==0)
        with (OUT/f'{fid}_alignment.csv').open() as f:assert len(list(csv.DictReader(f)))==a.shape[0]
    assert np.isnan(np.load(OUT/'F03_stack.npy')[6]).all()
    pngs=list(OUT.glob('*.png'))
    for p in pngs:
        with Image.open(p) as im:im.verify()
    save(OUT/'artifact_verification.json',{'frozen_raster_artifacts_unchanged':len(before),'prior_files_unchanged':len(read('preservation_before.json')),
        'PNG_matches_numeric_data_and_mask':True,'CSV_row_counts_match':True,'F03_Z146_all_missing':True,
        'readable_PNGs':len(pngs),'report_sha256':sha(path),'replay':read('verification.json')})
    paths=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='SHA256SUMS')
    (OUT/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in paths))
    print(json.dumps(verdict,indent=2));print(json.dumps(read('artifact_verification.json'),indent=2))


if __name__=='__main__':main()
