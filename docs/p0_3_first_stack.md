# P0.3 — first Direct Spiral-Stack rasters

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
| 137→138 | 11 | 0.8277 | 21.770 |
| 138→139 | 11 | 0.7928 | 20.412 |
| 139→140 | 11 | -0.4205 | 16.094 |
| 140→141 | 11 | 0.8088 | 9.855 |
| 141→142 | 11 | 0.7323 | 6.314 |
| 142→143 | 11 | 0.4431 | 8.578 |
| 143→144 | 11 | 0.4617 | 12.585 |
| 144→145 | 11 | 0.5607 | 28.083 |
| 145→146 | 11 | 0.8512 | 27.654 |
| 146→147 | 11 | 0.0324 | 16.256 |
| 147→148 | 11 | 0.1735 | 13.962 |
| 148→149 | 11 | 0.9273 | 10.317 |
| 149→150 | 12 | 0.9518 | 7.417 |

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
