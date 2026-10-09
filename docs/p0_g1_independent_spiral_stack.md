# P0-G1 — independent per-slice geometric spiral stack

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
(61.991%). These include masked connectors, STOP sentinels and
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
| 128 | 2.721 | 8.5 | 2 | 2 / 1 | 97.22% | 9 | very short; STOP |
| 129 | 0.742 | 2.0 | 1 | 1 / 0 | 99.07% | 3 | very short; STOP |
| 130 | 0.939 | 5.5 | 2 | 2 / 1 | 98.61% | 7 | very short; STOP |
| 131 | 34.871 | 126.0 | 24 | 24 / 23 | 68.06% | 170 | partial arc; STOP |
| 132 | 1.359 | 3.5 | 1 | 1 / 0 | 98.61% | 4 | very short; STOP |
| 133 | 4.577 | 44.0 | 9 | 9 / 8 | 95.37% | 54 | very short; STOP |
| 134 | 2.047 | 7.0 | 2 | 2 / 1 | 98.15% | 7 | very short; STOP |
| 135 | 35.924 | 120.0 | 17 | 17 / 16 | 65.74% | 161 | partial arc; STOP |
| 136 | 20.514 | 142.0 | 17 | 17 / 16 | 80.09% | 157 | partial arc; STOP |
| 137 | 11.827 | 64.5 | 10 | 10 / 9 | 88.89% | 75 | partial arc; STOP |
| 138 | 0.992 | 24.0 | 5 | 5 / 4 | 99.07% | 18 | very short; STOP |
| 139 | 11.340 | 62.0 | 8 | 8 / 7 | 88.89% | 65 | partial arc; STOP |
| 140 | 49.093 | 132.5 | 20 | 20 / 19 | 54.63% | 208 | partial arc; STOP |
| 141 | 37.884 | 137.0 | 26 | 26 / 25 | 65.28% | 212 | partial arc; STOP |
| 142 | 2.708 | 8.5 | 1 | 1 / 0 | 97.22% | 7 | very short; STOP |
| 143 | 2.224 | 7.0 | 1 | 1 / 0 | 97.69% | 6 | very short; STOP |
| 144 | 4.475 | 15.5 | 2 | 2 / 1 | 95.83% | 12 | very short; STOP |
| 145 | 0.000 | 0.0 | 1 | 1 / 0 | 100.00% | 2 | STOP; no sampled segment |
| 146 | 49.123 | 155.0 | 29 | 29 / 28 | 54.63% | 216 | partial arc; STOP |
| 147 | 0.197 | 0.5 | 1 | 1 / 0 | 99.54% | 2 | very short; STOP |
| 148 | 1.062 | 2.5 | 1 | 1 / 0 | 98.61% | 4 | very short; STOP |
| 149 | 33.399 | 115.5 | 21 | 21 / 20 | 68.52% | 143 | partial arc; STOP |
| 150 | 26.602 | 134.5 | 28 | 28 / 27 | 74.54% | 196 | partial arc; STOP |
| 151 | 7.246 | 46.0 | 6 | 6 / 5 | 93.06% | 42 | partial arc; STOP |
| 152 | 35.354 | 94.0 | 16 | 16 / 15 | 66.67% | 122 | partial arc; STOP |
| 153 | 21.276 | 120.0 | 17 | 17 / 16 | 80.56% | 127 | partial arc; STOP |
| 154 | 42.021 | 153.5 | 29 | 29 / 28 | 60.65% | 211 | partial arc; STOP |
| 155 | 34.563 | 144.5 | 24 | 24 / 23 | 67.59% | 184 | partial arc; STOP |
| 156 | 2.523 | 14.0 | 4 | 4 / 3 | 97.22% | 14 | very short; STOP |
| 157 | 24.025 | 135.0 | 27 | 27 / 26 | 78.24% | 208 | partial arc; STOP |
| 158 | 20.775 | 79.5 | 16 | 16 / 15 | 81.02% | 92 | partial arc; STOP |
| 159 | 12.957 | 40.0 | 9 | 9 / 8 | 87.50% | 56 | partial arc; STOP |
| 160 | 30.224 | 128.0 | 19 | 19 / 18 | 72.22% | 164 | partial arc; STOP |
| 161 | 10.047 | 87.5 | 17 | 17 / 16 | 90.28% | 102 | partial arc; STOP |
| 162 | 1.076 | 3.0 | 1 | 1 / 0 | 98.61% | 4 | very short; STOP |
| 163 | 43.734 | 129.5 | 25 | 25 / 24 | 59.72% | 202 | partial arc; STOP |
| 164 | 2.960 | 9.5 | 2 | 2 / 1 | 96.76% | 9 | very short; STOP |
| 165 | 3.013 | 10.0 | 2 | 2 / 1 | 97.22% | 8 | very short; STOP |
| 166 | 1.882 | 5.5 | 1 | 1 / 0 | 98.15% | 5 | very short; STOP |
| 167 | 1.728 | 6.0 | 2 | 2 / 1 | 98.15% | 6 | very short; STOP |
| 168 | 0.185 | 0.5 | 1 | 1 / 0 | 99.54% | 2 | very short; STOP |
| 169 | 2.664 | 29.5 | 11 | 11 / 10 | 98.15% | 52 | very short; STOP |
| 170 | 0.000 | 0.0 | 1 | 0 / 0 | 100.00% | 1 | no seed |
| 171 | 0.000 | 0.0 | 1 | 0 / 0 | 100.00% | 1 | no seed |
| 172 | 3.588 | 44.0 | 6 | 6 / 5 | 96.30% | 31 | very short; STOP |
| 173 | 0.512 | 1.5 | 1 | 1 / 0 | 99.07% | 3 | very short; STOP |
| 174 | 21.178 | 143.5 | 15 | 15 / 14 | 81.02% | 127 | partial arc; STOP |
| 175 | 0.000 | 0.0 | 1 | 0 / 0 | 100.00% | 1 | no seed |
| 176 | 0.421 | 1.5 | 1 | 1 / 0 | 99.54% | 2 | very short; STOP |
| 177 | 0.000 | 0.0 | 1 | 0 / 0 | 100.00% | 1 | no seed |

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
