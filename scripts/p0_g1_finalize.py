"""P0-G1 reporting and independent checks; no geometry or raster modifications."""
import csv
import json
import hashlib
from pathlib import Path
import numpy as np
from PIL import Image
from p0_g1_independent import process_slice

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_g1'
def read(p):return json.loads((OUT/p).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):(OUT/p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')


def main():
    rows=read('row_geometry.json');metrics=read('metrics.json')
    a=np.load(OUT/'lossless_stack.npy');mask=np.load(OUT/'lossless_mask.npy');xy=np.load(OUT/'sample_xy.npy')
    volume=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')
    assert a.shape==mask.shape==(50,216) and len(rows)==50
    assert np.array_equal(np.isfinite(a),mask) and np.array_equal(np.isfinite(xy).all(axis=-1),mask)
    for i in range(50):
        p=xy[i,mask[i]];x=p[:,0];y=p[:,1];x0=np.floor(x).astype(int);y0=np.floor(y).astype(int);dx=x-x0;dy=y-y0
        im=np.asarray(volume[i],float)
        independent=(1-dx)*(1-dy)*im[y0,x0]+dx*(1-dy)*im[y0,x0+1]+(1-dx)*dy*im[y0+1,x0]+dx*dy*im[y0+1,x0+1]
        np.testing.assert_allclose(a[i,mask[i]],independent,rtol=0,atol=1e-10)
        numeric=np.load(OUT/'slices'/f'z{128+i}_numeric.npz')
        record=read(f'slices/z{128+i}.json');pts=record['points']
        if pts:
            pp=np.array([p['xy'] for p in pts]);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(pp,axis=0),axis=1))]
            for gap in record['gaps']:
                interval=(numeric['s_px']>arc[gap['from_point']]+1e-12)&(numeric['s_px']<arc[gap['to_point']]-1e-12)
                assert np.isnan(numeric['row'][interval]).all()
        assert np.isnan(a[i,len(numeric['row']):]).all()
    # Reordered isolated calls: no previous-Z initialization or persistent state.
    for z in [177,146,128]:
        result=process_slice(np.asarray(volume[z-128]));old=read(f'slices/z{z}.json')
        for k in ['points','decisions','gaps','summary']:assert result[k]==old[k]
        np.testing.assert_array_equal(result['row'],np.load(OUT/'slices'/f'z{z}_numeric.npz')['row'])
    png=np.array(Image.open(OUT/'raw_stack.png'));pm=np.array(Image.open(OUT/'stack_mask.png'))
    np.testing.assert_array_equal(pm==255,mask);np.testing.assert_array_equal(png[mask],np.rint(a[mask]).astype(np.uint8));assert np.all(png[~mask]==0)
    overlay_files=sorted((OUT/'all_50_slice_geometry_overlays').glob('z*_geometry.png'));assert len(overlay_files)==50
    for p in overlay_files+[OUT/'raw_stack.png',OUT/'stack_mask.png',OUT/'geometry_diagnostics.png']:
        with Image.open(p) as im:im.verify()
    prior=read('preservation_before.json')
    for name,digest in prior.items():assert sha(ROOT/name)==digest,name
    padding=sum(a.shape[1]-r['native_row_width_samples'] for r in rows);native=sum(r['native_row_width_samples'] for r in rows)
    valid=int(mask.sum());internal=native-valid
    distributions={}
    for key in ['masked_fraction','recovered_arc_px','continuation_decisions','accepted_gap_continuations','interruptions','angular_coverage_deg']:
        v=np.array([r[key] for r in rows]);distributions[key]={'mean':float(v.mean()),'min_p10_p25_median_p75_p90_max':np.percentile(v,[0,10,25,50,75,90,100]).tolist()}
    save('distributions.json',distributions)
    outcome={'result':'RASTER_CONSTRUCTION_ACHIEVED_WITH_PARTIAL_GEOMETRIC_RECOVERY',
        'rows_processed':50,'rows_with_sampled_geometry':45,'fully_unresolved_rows':[r['z'] for r in rows if r['valid_samples']==0],
        'complete_turns_demonstrated':0,'masked_fraction':float(1-mask.mean()),
        'valid_samples':valid,'native_row_samples':native,'rectangular_padding_samples':padding,
        'inside_native_row_missing_samples':internal,'inside_native_row_missing_fraction':internal/native,
        'visible_coherent_structure':'Weak short start-aligned patches; no convincing broad coherent structure in the sparse raster.',
        'blanks':'Dominate both rectangular output and native row domains.',
        'cross_Z_persistence_used_or_required':False,'cross_Z_identity_failure_claimed':False,
        'recovered_spiral_is_partial':True,'same_slice_continuation_truth':'geometric plausibility only; no physical-layer oracle',
        'continuation_decisions_total':sum(r['continuation_decisions'] for r in rows),
        'interruptions_total':sum(r['interruptions'] for r in rows)}
    save('outcome.json',outcome)
    save('implementation_correction.json',{'interrupted_attempt':'outputs/p0_g1_initial_interrupted',
        'completed_slice_files_in_interrupted_attempt':len(list((ROOT/'outputs/p0_g1_initial_interrupted/slices').glob('*.json'))),
        'reason':'Repeated whole-image sin/cos orientation computation inside each peak loop.',
        'change':'Cache sin(2*tangent) and cos(2*tangent) once per slice; numerical formula and all geometric rules unchanged.',
        'initial_script_sha256':sha(ROOT/'outputs/p0_g1_initial_interrupted/initial_script.py'),
        'executed_script_sha256':sha(ROOT/'scripts/p0_g1_independent.py')})
    table='\n'.join(f"| {r['z']} | {r['recovered_arc_px']:.3f} | {r['angular_coverage_deg']:.1f} | {r['interruptions']} | {r['continuation_decisions']} / {r['accepted_gap_continuations']} | {100*r['masked_fraction']:.2f}% | {r['native_row_width_samples']} | {'no seed' if r['notes']=='NO_ELIGIBLE_START_RIDGE' else 'STOP; no sampled segment' if r['valid_samples']==0 else 'very short; STOP' if r['recovered_arc_px']<5 else 'partial arc; STOP'} |" for r in rows)
    report=f'''# P0-G1 — independent per-slice geometric spiral stack

## Corrected research direction

This experiment implements the user's corrected hypothesis:
**independently trace visible spiral geometry on each XY slice, leave breaks
blank, resume geometrically where credible, unwind the result to a row, and
stack the rows in Z order.** It does not require physical-fragment identity or
persistence between Z levels. P0.2–P0.4 are preserved historical experiments,
not acceptance gates for P0-G1. The current user instruction supersedes their
tracking interpretation; README.md and AGENTS.md remain unchanged as requested.

## Result

**A50-row raster was formed from independent slice processing, with explicit
masked gaps; geometric recovery is PARTIAL and sparse.** The output is50×216
(Z,s), with1336 valid samples. Forty-five slices contain sampled geometry;
five are fully unresolved. The tracer accepted435 **same-slice geometric**
continuations across gaps, and none of those gap interiors was sampled.

Blanks dominate:87.630% of the rectangular raster is masked. Median recovered
arc length is3.301px; the best slice recovers49.123px across disconnected visible
pieces. Maximum angular progress is155°, so **no complete turn or multi-turn
spiral has been demonstrated**. Short start-aligned patches are visible, but
there is no convincing broad coherent structure in the sparse stack.

Thus the requested unwind/mask/stack mechanism is exercised, but this bounded
tracer does not provide a strong full-spiral reconstruction. That limitation is
about **within-slice geometry**, not cross-Z physical identity. No row is
rejected because a fragment fails to persist into an adjacent slice.

## Input, physical coordinates and scope

- Scroll1 / PHercParis4, `54keV_7.91um_Scroll1A.zarr/4`.
- Remote source reference:
  `https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4`.
- Exact file: `outputs/orientation/ct_level4_zyx.npy`, shape(128,493,506),
  uint8, axesZYX, slab global Z offset128.
- All50 slices: global Z128–177 inclusive, local0–49; nominal native
  indices2048–2832 at stride16.
- Level4 voxel spacing126.56µm isotropic; catalog native spacing7.91µm.
  X increases right, Y down, Z approximately along the scroll.
- No new CT acquisition, labels, model, prior fragment or prior track is used
  in construction. Prior files are read only for context and preservation hashes.

The test follows **one independently selected candidate winding per slice**,
starting in a central region. It does not claim to extract every sheet in the
entire cross-section. The fixed polar reference is an approximate image-plane
convention, not a measured anatomical center. This choice is a limitation of
this first deterministic implementation, not a change to the research hypothesis.

## Reproducible same-slice geometry method

### Reference and detection

Use fixed centerC=(240,140), elliptical coordinates
`p(theta,r)=C+(r*cos(theta),0.5*r*sin(theta))`.
The start ray is +X, theta=0, with increasing theta clockwise in image pixels.
Radii8–180px are sampled every0.25px; angles every0.5°. On each ray detect peaks
of `Gaussian_0.6(CT)−Gaussian_2(CT)` with prominence≥2.5CT units,
response>3 and radial separation≥1px. Peak positions receive bounded parabolic
subpixel refinement. A Gaussian-gradient structure tensor (gradient sigma1,
product averaging sigma1.2) supplies local unoriented ridge tangent/coherence.
These operations detect CT geometry; they are not ink inference.

Among start-ray peaks of radius20–80 with tangent alignment to the upward/downward
ray-transverse direction≥0.65 and coherence≥0.2, select the peak nearest r=40,
then higher coherence, then deterministic candidate index. **No seed is supplied
from a neighboring slice or from F02/F03.** A missing eligible seed yields an
entirely masked row. Each slice has its own selected start at s=0.

### Ordered traversal and credible continuation

Follow increasing angular stations. Predict radius from the last up-to7
consecutive points, using a local linear radial slope clipped to±120px/radian;
after a gap use the observed ridge tangent to restart the slope estimate.
Enumerate all detected candidates within±6px of predicted radius and chord
distance0.03–14px from the last accepted point. Measure predictive proximity,
forward tangent agreement, turning angle, candidate ridge orientation and local
radial peak spacing. Exact per-candidate values/failures are saved.

Score components, bounded0–1:

- Proximity: `exp(-distance_to_prediction/max(1,0.4*spacing))`, weight0.40.
- Tangent: nonnegative cosine with outgoing direction, weight0.35.
- Turn consistency: `exp(-turn_angle_radians/0.7)`, weight0.15.
- Ridge orientation: absolute cosine with connection direction, weight0.10.

Require score≥0.68, winner margin≥0.06, ridge alignment≥0.65, turning angle≤65°,
orientation coherence≥0.15 and predicted normal miss≤0.45local spacing.
Spacing is the nearest neighboring radial-peak interval transformed into XY
distance along that ray. These are disclosed geometric heuristics, not physical
sheet-identity probabilities. The curvature term is a local turn proxy, not a
fitted material curvature field.

If the next station is untraceable or ambiguous, leave it unresolved and search
forward up to15° using the last accepted geometry. The first later station
with a uniquely credible candidate can resume tracing, but **the entire
intervening connector remains masked**. If no candidate qualifies, stop that
slice and leave the rest unresolved. Multiple interruptions are permitted.
Traversal has a three-turn cap and stops on near-revisiting old geometry rather
than forcing closure; neither limit was reached in the observed results.

This is not a cross-Z tracker. The numerical entry point is
`process_slice(image)`: no Z value, neighboring image, old path or persistent
track state is an input. There is no 3-D surface or mesh.

### Arc length, gap width and sampling

Within recovered runs, use cumulative Euclidean chord arc length. For a gap,
reserve the straight endpoint-to-endpoint chord length as a **lower-bound
geometric placeholder**. Do not claim it is the true length of torn/folded or
occluded papyrus. The missing path and its actual arc length are UNKNOWN.
After such a gap the row coordinate is therefore approximate, even though
sample spacing along recovered segments is fixed.

Sample at **Δs=0.5level4 px=63.28µm**, using piecewise-linear XY interpolation
along recovered segments and bilinear original uint8 CT evaluated in float64.
No gap connector is sampled, no missing Z row is interpolated, and no intensity
is used to imitate letters. Coordinates of missing samples are NaN.

Append one masked sample beyond a terminated row as an open-ended STOP sentinel.
This sentinel does not estimate the unknown remainder's physical length. Pad
all rows to the largest observed geometric extent (216samples) with masked
values. There is no row resizing, nonlinear alignment or cross-Z registration.
The output's s coordinate is recovered arc length plus approximate gap lengths,
not a globally registered physical sheet coordinate.

## Global measurements

| Quantity | Result |
| --- | --- |
| Slices processed / output rows | 50 /50 |
| Rows with sampled geometry | 45 |
| Fully unresolved rows | 145,170,171,175,177 |
| No eligible start | 4:170,171,175,177 |
| Stops after unsuccessful ≤15° continuation search | 46 |
| Output dimensions | height50,width216 |
| Valid / total samples | 1336 /10800 |
| Total missing fraction | 87.6296% |
| Recovered arc, min /median /mean /max, px | 0 /3.300722 /13.171451 /49.122820 |
| Recovered arc, median /max, mm | 0.417739 /6.216984 |
| Angular progress, median /maximum | 26.75° /155° |
| Same-slice accepted masked-gap continuations | 435 |
| Continuation episodes, including terminal STOP decisions | 481 |
| Interruptions including four no-start cases | 485 |
| Complete angular turns demonstrated | 0 |

Z145 found a seed and candidate continuation geometry but no connected sampled
segment at the fixed grid interval; its row remains entirely masked. A selected
point alone is not counted as recovered arc. This is distinct from the four
no-seed slices.

| Distribution | Minimum | P10 | P25 | Median | P75 | P90 | Maximum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Masked fraction | 0.5463 | 0.6569 | 0.7870 | 0.9653 | 0.9861 | 0.9958 | 1.0000 |
| Recovered arc, px | 0 | 0.1669 | 1.0652 | 3.3007 | 23.3378 | 36.1199 | 49.1228 |
| Accepted gap continuations | 0 | 0 | 0 | 4.5 | 16 | 24.1 | 28 |
| Continuation episodes incl STOP | 0 | 1 | 1 | 5.5 | 17 | 25.1 | 29 |

Blank accounting:7285samples are rectangular padding beyond the individual
row domains; inside the3515native-domain samples,2179 are still missing
({100*internal/native:.3f}%). These include masked connectors, STOP sentinels and
unsampled isolated-point regions. Thus blanks are not solely an artifact of
rectangular padding: they dominate even inside the traversed row domains.
The mask does not imply physical absence of papyrus.

## Per-slice results

All rows use Δs63.28µm, final width216 and the geometric start convention above.
The CSV additionally includes physical arc length, exact valid-sample count,
gap-placeholder arc length and the repeated reference string for standalone use.
“Decisions/accepted” counts gap-level continuation episodes (plus terminal STOP)
versus accepted masked-gap continuations; ordinary local angular steps are
reported separately in `row_geometry.csv`. No-start cases have zero continuation
decisions but one interruption. Native width excludes rectangular padding.

| Z | Recovered arc px | Angular progress° | Interruptions | Decisions /accepted gap continuations | Masked fraction | Native width | Notes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
{table}

## Qualitative observations and geometric failure modes

`all_50_slice_geometry_overlays/` contains all50 original/overlay image pairs.
Green is sampled local geometry; dashed magenta is a **masked diagnostic
connector**, never a recovered path; yellow is the independent start and red
the terminal STOP. The fixed ray/reference is drawn separately.

Examples: Z140,146 and154 recover recognizable lower-right/lower central
winding fragments interspersed with gaps. Z128 stops almost immediately.
Z170 has no eligible ray seed despite substantial papyrus elsewhere in the
image. These demonstrate both real per-slice geometry use and incomplete
search coverage. The method is not allowed to interpret a failed seed as an
empty CT slice.

Main limitations:

1. **Fixed ray/center and single-seed selection.** Changes in local shape can
   select different radial bands or produce no eligible seed. This is geometric
   reference instability, not failure of physical-fragment persistence.
2. **Angular-monotone traversal.** Flattened folds and cusps may require backward
   angular motion, which this implementation cannot represent. A ≤15° look-ahead
   cannot bypass every fold or occlusion.
3. **Local orientation/spacing ambiguity.** Dense ridges, partial-volume effects
   and non-radial layer spacing limit the continuation heuristic. A passed score
   is geometric plausibility, not independently established same-sheet identity.
4. **Many small gaps.** Frequent blank connectors and early STOPs leave little
   contiguous sampled geometry on many slices. Later resumption is allowed but
   does not fabricate intervening tissue.
5. **Unknown gap arc length.** Chord placeholders underestimate curved hidden
   paths and may cause downstream s shifts. No cross-Z adjustment conceals this.

The raster contains sparse horizontal pieces and weak local start-aligned
patches; a convincing extended coherent pattern is **not demonstrated**.
Cross-Z identity or correlation is not used as a success/failure gate. The
negative evidence concerns geometric coverage and ambiguity within each slice.

## Artifacts and reproduction

Required files under `outputs/p0_g1/`:

- `all_50_slice_geometry_overlays/z128_geometry.png` through `z177_geometry.png`.
- `raw_stack.png`, `stack_mask.png`, `geometry_diagnostics.png`.
- `row_geometry.csv`, `lossless_stack.npy`, `lossless_mask.npy`.

Additional lossless data: `sample_xy.npy`, `s_px.npy`, `z.npy`;
`slices/zNNN_numeric.npz` contains each independent row, sample XY and all
ridge candidates. Candidate columns are angular-station index,radius,x,y,
ridge-angle radians,orientation coherence,DoG response,nearest radial spacing.
`slices/zNNN.json` retains accepted points, every enumerated continuation score,
failure reason, gap and stop. `row_geometry.json`, `metrics.json`,
`distributions.json`, `outcome.json`, `plan.json`, `provenance.json`,
`verification.json`, `artifact_verification.json`, `SHA256SUMS` preserve results.

PNG data use fixed0–255 CT, nearest-integer uint8 preview; invalid pixels are0
and require the companion mask (255valid/0missing). Numeric CT is float64 with
NaNs, masks are bool. The raster pixel spacing is63.28µm in s versus126.56µm in Z;
display panels are enlarged for inspection, not used as physical measurements.
No enhancement, OCR, ink inference, language priors, learned reconstruction or
simulation was performed.

Commands executed from repository root:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-g1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_g1_independent.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-g1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_g1_independent.py --verify
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-g1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_g1_finalize.py
```

The first attempt was stopped before any slice completed because sin/cos
orientation arrays were recomputed for every candidate. The interrupted
directory and original script are preserved at
`outputs/p0_g1_initial_interrupted/`. The only correction caches those same
arrays once per slice; geometric equations, parameters and decisions are
unchanged. `implementation_correction.json` binds both script versions.
No thresholds were adjusted after observing results.

Verification exactly reproduces all50slice decisions and numeric arrays;
reordered isolated calls for Z177,146,128 also match, demonstrating no prior-row
state. Independent four-corner bilinear interpolation agrees within1e-10.
All gap interiors, terminal/padded domains, masks and PNG values are checked;
all50overlays decode successfully. Prior-phase hashes, README.md and AGENTS.md
remain unchanged. Generation refuses to overwrite the evidence directory.

## Conclusion and next recommendation (not executed)

The corrected independent-slice experiment has been executed, producing the
requested50-row raw raster with explicit unresolved spans. This demonstrates
the mechanics of independent geometric unwinding with masked continuation,
but **useful spiral recovery remains partial**, and blank regions dominate.
It neither rejects the original idea on cross-Z grounds nor claims a complete
spiral was recovered.

A next geometry-only experiment could evaluate a per-slice fragment graph or
local-tangent traversal that handles angular reversals, with an explicit
same-slice start-selection test and continuation audit. That would address the
observed 2-D failure modes without reverting to cross-Z tracking or a mesh.
No follow-up is implemented here. **Stop after P0-G1.**
'''
    path=ROOT/'docs/p0_g1_independent_spiral_stack.md';assert not path.exists();path.write_text(report)
    save('provenance.json',{'source_sha256':sha(ROOT/'outputs/orientation/ct_level4_zyx.npy'),
         'executed_script_sha256':sha(ROOT/'scripts/p0_g1_independent.py'),'report_script_sha256':sha(Path(__file__)),
         'report_sha256':sha(path),'source_z_inclusive':[128,177],'source_zyx_shape':list(volume.shape),
         'source_dtype':str(volume.dtype),'voxel_spacing_um':126.56,'arc_sampling_um':63.28,
         'constructed_from_prior_tracks':False,'cross_Z_geometry_inputs':False,'download_bytes':0})
    save('artifact_verification.json',{'all50slice_replay':read('verification.json'),
        'reordered_isolated_slice_replay':[177,146,128],'independent_bilinear_max_tolerance':1e-10,
        'gap_interiors_and_padding_missing':True,'mask_nan_png_consistency':True,'overlay_count':len(overlay_files),
        'prior_files_hash_preserved':len(prior),'README_AGENTS':'unchanged'})
    paths=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS')
    (OUT/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT)}\n' for p in paths))
    print(json.dumps(outcome,indent=2));print(json.dumps(read('artifact_verification.json'),indent=2))


if __name__=='__main__':main()
