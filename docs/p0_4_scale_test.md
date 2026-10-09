# P0.4 — longitudinal scale test of F02

## Verdict

**P0.4 FAIL — longitudinal scale-up under the frozen method.**

Fifty consecutive Z levels were attempted. The F02 seed-connected identity
remains supported only on the original **14 slices, Z=137–150**. No new row is
supported in either direction. The other **36 rows remain entirely masked**.
The original coherent region survives bit-identically; this is **not an observed
collapse of its CT coherence**. The failure is inability to establish a longer
physical identity without accepting unsupported continuations.

This is not PARTIAL under the requested definition: no second bounded region
has an accepted connection to F02. Disconnected local tracklets cannot be
relabeled F02 merely because they look coherent. No continuation was forced,
and no mesh or revised scoring method was introduced. The result does not prove
that Direct Spiral-Stack inherently needs a 3-D surface, nor that papyrus ends
at the current track boundaries.

## Hypothesis and bounded experiment

Hypothesis: the existing F02 winding-fragment identity and fixed-scale
arc-length raster can remain supported over substantially more Z levels with
the existing deterministic method.

Input dataset: Scroll1 / PHercParis4, standardized
`54keV_7.91um_Scroll1A.zarr/4`, remote reference:
`https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4`.

- Exact CT file: `outputs/orientation/ct_level4_zyx.npy`, uint8,
  shape ZYX=(128,493,506), global Z offset128.
- Available slab: global Z128–255. Attempted range: **128–177 inclusive**,
  local indices0–49, nominal native Z2048–2832 at stride16.
- The initial 50-plane interval is asymmetric around the original track because
  the cached slab begins at128. It still provides candidate evidence in both
  directions, including immediate and one-missing-slice alternatives at both
  track ends. Neither end reached a data boundary, so this asymmetry did not
  cause the observed stop.
- Isotropic source spacing: **126.56 µm** at level4; catalog native7.91 µm.
  Axes remain ZYX; X rightward, Y downward, Z approximately longitudinal.
- Seed: `z144_control_F02_r02`, previously T0077. Track IDs created by this run
  are local graph identifiers; node lineage, not the new numeric track ID,
  establishes connection to the validated F02 seed.
- Frozen inputs include P0.2 run02 fragments/windows/tracks and P0.3 F02
  numeric rows, axes and validation evidence. No F03 work is performed.
- No downloads or dependency installation. Existing CT is read without modifying
  it. The local partial Zarr cache is not treated as a complete scroll.

`outputs/p0_4/plan.json` was written before computation. Extension beyond the
initial interval was conditional on the seed-connected track reaching within
one slice of its boundary. This did not happen. Continuing to100–200 attempts
would add disconnected hypotheses, not evidence of a reliable F02 extension;
the conditional expansion was therefore not executed.

## Frozen method and parameters

The implementation **imports** the earlier functions rather than rewriting their
rules: `extract_window` and `associate` from `scripts/p0_2_local_ridges.py`,
`registration` and `tracks_from` from `scripts/p0_2_cross_z.py`, and
`arc_geometry`, `sample_row`, `correlation` from `scripts/p0_3_first_stack.py`.
All prior source hashes are checked. No thresholds were relaxed after failure.

The fixed F02 window has origin (270.916036,132.741970), tangent
(0.776693,0.629880), u=−4…4 in0.5px steps and normal v=−8…8 in0.25px steps.
Each plane independently supplies DoG(σ0.6−σ2) peaks with prominence≥1,
response>1 and separation≥1px. Mutual-nearest adjacent-u peaks require normal
step≤0.75px and competitor margin≥0.125px. Retain ≥7 vertices and arc≥3px.
This is the earlier local detector, not a general surface extractor.

Registration uses the same fixed ROI x160:296,y100:174, Gaussianσ1 and bounded
translation search ±2px at0.5px increments. Frozen association evidence includes
spatial residual, common-u overlap, tangent angle, curvature, local spacing,
normal CT profiles and radial displacement around proxy center(240,140).
Weights: position0.20, overlap0.15, texture0.35, orientation0.10,
curvature0.05, spacing0.05, radial proxy0.10.

Acceptance still requires ≥5 common stations, normal residual≤1px,
NCC≥0.55, angle≤25°, overlap≥0.6, score≥0.70, mutual best matching,
margin≥0.06 and displacement≤0.45 local spacing. Hard rejects remain
residual>2.5px, NCC<0.1 or angle>50°. One-missing-slice associations retain
the stronger score≥0.76, NCC≥0.70 and margin≥0.10 requirements, with original
adjacent-link endpoint exclusions. There are no longer-gap links or manual
restarts. All candidates and failures are preserved in `edges.json`.

Only the connected component containing the exact F02 seed may supply rows.
Disconnected tracklets are retained as diagnostic evidence, never joined into
the F02 raster. A detected ridge is not sufficient to assign physical identity.

Rasterization remains the P0.3 core u=[1,3.5], origin at u=2.25,
**Δs=0.25px=31.64µm**, **Δz=126.56µm**. Each observed polyline is parameterized
by cumulative Euclidean arc length, translated to the fixed anchor and sampled
via piecewise-linear XY interpolation and bilinear original CT. No equal-width
rescaling, local stretch correction, CT contrast processing or Z interpolation
is performed. Unsupported endpoints/rows are NaN with false masks. The old
14-row overlap is **bit-identical to P0.3**, including fractional CT values.

## Scale metrics

| Metric | Result |
| --- | --- |
| Z levels attempted | 50, Z128–177 |
| Identity-supported / valid raster rows | 14 / 14 |
| Entirely masked rows | 36:128–136 and151–177 |
| New supported observations beyond P0.3 | 0 |
| Longest continuous supported run | 14, Z137–150 |
| Seed-connected supported regions | 1 |
| Accepted F02 adjacent / gap links | 13 / 0 |
| Raster dimensions (Z,s) | 50×15; PNG width15,height50 |
| Valid samples / total samples | 158 / 750 |
| Missing samples / fraction | 592 / 78.9333% |
| Entire-row missing samples | 540 |
| Edge padding outside retained cores | 52 |
| Evaluable adjacent-row correlations | 13 of49 attempted adjacent pairs |
| Unevaluable adjacent pairs | 36; not counted as correlation zero |
| Positive correlations | 12/13 =92.3077% of evaluable pairs |
| Low-coherence transitions, NCC<0.2 | 3 |
| Suspected layer switches within accepted F02 links | 0 apparent, per CT review; true count UNKNOWN |
| New links extending F02 | 0; no extended false-join rate can be estimated |
| Conflicting accepted two-step targets in F02 | 0 |
| Independently detected candidate fragments | 182 |
| All local tracklets / unconnected to F02 | 69 / 68 |
| Candidate associations: accepted/rejected/unknown | 113 /150 /101 |

The113 accepted associations include other window-local tracklets; **only13
belong to the F02 seed component**. Their total must not be reported as113
validated winding continuations. Zero apparent switches does not establish a
zero biological error rate: there is no registered layer oracle here, and no
new F02 links were admitted. `metrics.json` leaves the automated switch count
null; `visual_review_scope.json` records the narrower apparent count and scope.

### Distributions, not just averages

| Measurement | Min | 5th percentile | 25th | Median | Mean | 75th | 95th | Max | SD |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Adjacent NCC (13 pairs) | −0.420502 | −0.148756 | 0.443129 | 0.732341 | 0.549468 | 0.827692 | 0.937119 | 0.951823 | 0.392399 |
| Adjacent origin shift, px (13 pairs) | −0.254038 | −0.215019 | −0.086114 | −0.024144 | −0.000875 | 0.091703 | 0.209128 | 0.240790 | 0.135532 |
| Retained row length, px (14 rows) | 2.637613 | 2.667675 | 2.756490 | 2.803284 | 2.837154 | 2.954926 | 3.043554 | 3.043911 | 0.133249 |

Multiply shift/length pixels by126.56 for µm. Maximum absolute shift is32.151µm;
row lengths span333.816–385.237µm. These statistics describe only the unchanged
14 supported rows, **not successful longitudinal scale-up**. The original local
coherence criterion remains satisfied within that region but cannot compensate
for failure to extend identity.

`F02_correlation_by_z.csv` contains every adjacent pair, with explicit support,
eligibility, NCC, descriptive NCC, CT difference and low-coherence indicator.
`F02_extended_diagnostic.png` plots correlation against longitudinal Z; blanks
outside137–150 are unevaluable, not a hidden continuous correlation curve.

## Why the two directions stopped

### Lower-Z boundary:136→137 and135→137

The best136→137 candidate, `z136_control_F02_r01`, has score0.675,
CT-profile NCC0.585 and median normal residual0.936px, but fails orientation,
overall score and competitor margin. Another nearby piece r04 is closer at
0.600px but has NCC0.365, score0.636 and fails texture/competition. More distant
bands are rejected. The135→137 missing-slice alternative r03 has score0.512,
NCC0.290, residual1.855px and fails both ordinary and stronger gap gates.

The raw CT and multiple proposals in `boundary_CT_overlays.png` support a
competition/shape-change ambiguity, not a confident physical endpoint.
No candidate is outside the image and the cached lower bound128 is not reached.
We do not select the visually closest piece after observing the failure.

### Upper-Z boundary:150→151 and150→152

Z150's selected fragment is already short. At151 the nearby r01 proposal has
only3 common stations, residual1.566px, NCC0.437 and score0.510; it fails
support, displacement, texture, overlap, score and competition. Other enumerated
bands are farther away and rejected. The150→152 candidate r01 has seven common
stations and residual0.904px but NCC0.465 and score0.676, failing texture,
score and strict gap evidence. The fact that a plausible band is visible at152
does not establish an accepted connection across151.

Detection fragmentation/support instability and changing local appearance are
plausible contributors. Identity is **UNKNOWN**, not “papyrus disappeared.”
Longer look-ahead, window adaptation or relaxed thresholds would change the
tested method; none was introduced to recover a favorable scale result.

## The three P0.3 low-coherence transitions

These were audited on original CT, frozen polylines, full association evidence
and six homologous core-u vertices. The latter comparison is diagnostic only:
no raster is warped into u or corrected to increase correlation. Correlations
from six coarse, dependent samples are not a physical truth oracle.

| Transition | Arc-row NCC | Homologous-u NCC | Normal-profile association NCC | Normal residual, px | Max anchored arc residual, px | Mean CT change at u |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 139→140 | −0.420502 | −0.169957 | 0.783419 | 0.307520 | 0.240790 | +2.965511 |
| 146→147 | 0.032407 | −0.174864 | 0.655299 | 0.572168 | 0.105894 | −18.220735 |
| 147→148 | 0.173549 | 0.016927 | 0.852306 | 0.465907 | 0.112254 | +8.441670 |

For **each transition**, primary causal classification is **UNKNOWN**:

- **Legitimate CT/material change:** CT intensity/profile changes are measured,
  and the selected band remains physically plausible, but the biological cause
  is not established. Mean intensity differences alone do not explain NCC,
  which subtracts each row's mean.
- **Alignment error:** possible. Pure arc-origin translation leaves small
  homologous-coordinate residuals. Using homologous-u samples improves139→140
  somewhat but does not restore positive coherence; it worsens the other two.
  Thus this diagnostic does not support arc-coordinate mismatch as the sole
  explanation. It does not rule out material-point misregistration.
- **Detector fragmentation:** not the direct cause of these three raster
  transitions: all six core vertices are present in both observations and
  neither row has an interior detector gap. This differs from the boundary
  failures and from the earlier F03 gap.
- **Possible layer switch:** no positive evidence in the accepted CT review or
  F02 shortcut audit; physical band identity remains provisional, not proven.
- **Insufficient support:** the existing ±0.5px shift-scan diagnostic cannot
  retain its required2px span at all offsets. It cannot decisively adjudicate
  alignment; no threshold was weakened to produce a cause label.

The strongest counterexample is139→140: accepted normal-profile matching can
coexist with negative **along-fragment** intensity correlation. These are
different signals and different questions. A high ridge-identity score must
not be substituted for coherent row content.

Individual evidence is retained in `low_coherence_transition_audit.json` and
`transition_139_140.png`, `transition_146_147.png`, `transition_147_148.png`.
None of the three was deleted, realigned for improved metrics or called a
confirmed biological change. Their ambiguity is an outcome of this analysis.

## Deliverables and exact commands

All requested files are under `outputs/p0_4/`:

- `F02_extended_stack.png`, `F02_extended_stack_mask.png`.
- `F02_extended_diagnostic.png` with CT/mask/endpoints, correlation-vs-Z and offsets.
- `F02_correlation_by_z.csv`, `F02_alignment.csv`.
- `F02_extended_stack.npy` float64 (NaN means missing), `F02_extended_mask.npy`
  bool (true means valid), `F02_extended_sample_xy.npy` float64,
  `F02_extended_s_px.npy`, `F02_extended_z.npy`.

The PNG is fixed-range raw CT quantized with nearest-integer rounding; missing
pixels are0 and must be interpreted with the companion mask (255 valid,0
missing). No appearance enhancement is used. Numeric sample intervals are
31.64µm horizontally and126.56µm vertically; the PNG has one pixel per sample,
so square-pixel viewing is not an equal-physical-aspect display. The diagnostic
panels are enlarged for legibility, not used to measure geometry.

Additional reproducibility files: `nodes.json`, `edges.json`, `tracks.json`,
`selected.json`, `registration.json`, `boundary_candidates.json`,
`shortcut_audit.json`, `alignment.json`, `pairs.json`, `metrics.json`,
`plan.json`, `frozen_tracking_parameters.json`, `preservation_before.json`,
`visual_review_scope.json`, `verdict.json`, `provenance.json`,
`verification.json`, `artifact_verification.json`, `SHA256SUMS`.

Executed from `/home/wpotech/deployments/vesuvius-spiral-stack`:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p04 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_4_scale.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p04 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_4_scale.py --verify
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p04 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_4_audit.py
.venv/bin/python -B scripts/p0_4_finalize.py
```

Generation refuses to overwrite its output directory. Verification recomputes
the experiment and compares saved arrays/JSON exactly. It also checks unchanged
extraction on prior slices, bit-identical P0.3 sample overlap and all prior
hashes. Final checks verify graph direction/one-to-one membership, PNG/NPY and
CSV agreement, independent bilinear sampling, source bounds and missing masks.
README.md, AGENTS.md and all prior-phase evidence remain unchanged.

## Conclusion / next recommendation (not executed)

The earlier F02 local result survives, but the frozen extraction/matching rules
cannot establish a substantially longer F02 identity. Report **FAIL for this
scale test**, not a successful50-row reconstruction and not proof that the
physical winding terminates. Maintaining identity beyond the present boundaries
would require presently unverifiable joins under this method.

Next recommended experiment: obtain independent same-region identity evidence
at136↔137 and150↔151↔152, then separately test the specific detection or
correspondence limitation against known wrong-layer candidates. Preserve this
failure before proposing any changed method; simply loosening gates is not a
supported remedy. Do not proceed to ink, OCR, learned correction, mesh
reconstruction or Idea B. **Stop after P0.4.**
