# P0-T11 — local 3-D frame correction on frozen Direct Threads

## Result and scope

**C = FRAME_CORRECTION_NO_BENEFIT under the frozen patch-geometry gate.**
Frame estimation passes for one of two candidates, but the required geometry
improvement does not. Smoke and full model inference were therefore NOT RUN.
This is a completed, gate-limited negative experiment, not an observed failure
of T11 ink probabilities. Ink visibility is NOT ASSESSED; the requested enum
NO means not established. T9 remains FAIL and T10 remains A, unchanged.

The hypothesis was that an orientation-only correction at frozen 2-D thread
centers would improve model-input patch geometry before any learned readout.
The comparison found weak improvements in some CT proxies, not enough to open
the original inference gate. No threshold was relaxed after adding coverage.
No candidate was selected for inference.

## Frozen inputs and methods

The full target remains the50independent rows Z128..177 and130,629valid frozen
T5/T9 centers. The bounded sample comprises72fixed centers: groups128..130,
151..153,175..177, eight evenly indexed common-valid s stations per group.
Exact coordinates, row/index mappings, tangents and normals are saved in
`outputs/p0_t11/samples.json`. Selection used geometry identities, not CT
appearance or probabilities. Coordinates, s, bridge decisions, masks,
provenance and ordering were never changed.
Native coordinates retain the T9 scale-only mapping
P=(16*x_level4,16*y_level4,16*Z), with no fitted origin offset.

The original plan and scientific implementations were committed before target
evidence as `9f27070`. See `docs/p0_t11_plan.md` and
`outputs/p0_t11/PRE_EVIDENCE_FREEZE.json`. All frozen hashes still match.

Only two deterministic estimators were tested: constrained3-D structure tensors
using Gaussian derivatives sigma1 or sigma2native voxels, on a33³local cube,
with Gaussian spatial weighting sigma6. The dominant normal is constrained
perpendicular to the exact within-slice tangent. This preserves the tangent
of the frozen transverse trace while allowing the normal and longitudinal
direction to tilt in3-D. There is no point relocation, cross-Z path fitting,
general surface mesh, smoothing of the recovered path, or curvature fitting.

Normal sign aligns with the saved T9 normal. The reported orthonormal frame
satisfies t_s cross t_z = n, with t_z = n cross t_s. Array rendering uses
image Y = -t_z because the original T9 image frame is left-handed. This retains
T9's image orientation and depth sign at zero tilt; it is not a depth reversal.

Rendering remains17×128×128 at9.6µm spacing, native scale9.6/7.91, centered depth
offsets-8..8, trilinear float32 sampling and half-even uint8 rounding. Every old
control patch was independently compared to the unchanged sealed T9 renderer:
72points×2stored copies, maximum uint8 error0.

The checkpoint remains
`scrollprize/ink_9um/hybrid_3d2d-seed42/step-075000.pth`, SHA256
`e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab`.
Its existing bytes were freshly verified, but no T11 model load/forward pass
was needed because the geometry gate closed. Official normalization, precision,
batch8 and central2×2 output definitions in the conditional runner are unchanged.
No local LLM, OCR, text priors, training, checkpoint search, or target-region
trusted geometry was used.

## Missing coverage: infrastructure only

The first CT-only attempt completed141rotated patch arrays and retained144
sample/method records, then stopped because three patch reads lacked source
chunks outside T9's original rectangular sampling footprint. Its files,
failure log, blocker record and original freeze are preserved byte-identically
under `outputs/p0_t11/frames/`, `frames.log`, and `DATA_BLOCKER.json`.
The147entries in `blocked_attempt_preservation.json` include the initial
attempt and original plan/freeze/sample files. This blocker had no influence
on scientific classification.

Authoritative source, unchanged from T9:

https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/

Native level0,7.91µm,128³uint8 chunks in ZYX order. Exact half-open native CT
bounds, trilinear upper-neighbor support, URLs, expected sizes and sample/method
dependencies are recorded in `missing_ct_inventory.json`.

| Supplemental chunk Z/Y/X | Bytes |
|---|---:|
|22/17/31|817,557|
|22/17/32|830,799|
|22/17/33|833,032|
|22/18/33|817,319|
|22/19/33|819,155|
|22/19/46|717,385|
|Total6chunks|4,835,247|

Only these inventoried objects were downloaded, under
`data/raw/p0_t11/native_chunks/`. The existing5,523chunks totaling
4,443,682,806bytes were freshly size/hash verified and were not replaced or
redownloaded. Supplemental verification checked HEAD sizes, Blosc decoding
to128³uint8, and saved SHA256. The server publishes opaque ETags, not a
cryptographic content digest; local SHA256 establishes subsequent integrity,
not an independent publisher-signed digest.

REQUIRED_PATCHES=144; FULLY_COVERED_PATCHES=144. Every preregistered rotated
patch has complete rectangular CT support. No clipping or padding was used.
This coverage claim does not extend to every potential full-inference patch.

The support-only wrapper prefers original T9 chunks and falls back to the
T11 supplement. All original scientific modules are imported unchanged.
The successful retry is isolated in `complete_support_attempt/`; replay is
in its separate `frame_replay/` directory. The user explicitly authorized
acquisition while review was running, before the support wrapper hash freeze.
That chronology is disclosed, not backdated. The wrapper, inventory and
support receipts were hash-frozen before the retry. Following the later user
instruction, the support amendment and result are committed only at sealing.

## Frame sanity and comparison to T10

Reject an estimate if constrained eigenvalue confidence<0.25, energy is
degenerate, longitudinal tilt>75°, unconstrained normal/tangent incompatibility
>30°, or shifting the diagnostic neighborhood by±2native Zvoxels changes the
axial normal by>20°. The diagnostic shifts never move the rendered center.
Sanity requires≥70%accepted and median perturbation≤15°.

Distributions below include all72estimates per candidate, including rejects:

| Diagnostic | Sigma1 | Sigma2 | Trusted T10 control |
|---|---:|---:|---:|
|Accepted estimates|50/72(69.44%)|54/72(75.00%)|Not comparable|
|Sanity gate|FAIL|PASS|Not applicable|
|Median absolute normal Z|0.55133|0.41837|0.34915|
|95th-percentile absolute normal Z|0.99116|0.95461|0.72360|
|Median longitudinal tilt|33.46°|24.73°|25.03°|
|95th-percentile tilt|82.52°|72.84°|48.87°|
|Median confidence|0.65269|0.77758|Not applicable|
|Median local perturbation angle|1.44°|1.56°|Not measured|
|95th-percentile perturbation angle|5.78°|5.89°|Not measured|
|Median sparse-station normal angle|58.78°|54.55°|Not comparable|

Rejection causes overlap: sigma1 has4low-confidence,1perturbation,
11tangent-incompatibility and8excessive-tilt failures; sigma2 has0,1,13,4
respectively. Distinct rejected counts are22and18.

Sparse-station angles are between widely spaced stations, not local
along-thread smoothness. The local check measures only±2native-Z perturbation
stability. Matching T10's median tilt is diagnostic, not an acceptance criterion
or evidence that the same physical sheet is being followed. T10 uses a different
scan and locations; its geometry was never applied to the T11 target.

## Paired patch-quality result

Metrics compare old/new on the same accepted sample IDs within each candidate.
Sigma1 and sigma2 accepted subsets differ, so their medians are not a paired
between-estimator comparison. Distant model-Y wings hold out the region used
for orientation estimation. Intensity-based proxies are neither ink evidence
nor a registered surface oracle.

| Median CT proxy | Sigma1 old → new | Sigma2 old → new |
|---|---:|---:|
|Longitudinal gradient energy fraction|0.35393→0.27343|0.35001→0.29396|
|Depth-peak standard deviation(indices)|5.55634→5.40515|5.54503→5.45033|
|Normalized distance from depth center|0.61558→0.59890|0.61456→0.60489|
|Empty-volume fraction|0.15348→0.13656|0.14594→0.12287|
|Repeated-band peak count proxy|4.000→4.333|4.333→4.333|
|Center-slab occupancy|0.68561→0.72984|0.71742→0.73226|

Sigma2 reduces longitudinal gradient fraction16.01% and empty fraction15.81%
relative, but depth-peak spread improves only1.7079%, below the frozen10%
requirement. Centeredness improves only1.57% and repeated-band counts do not
improve. Every other sigma2geometry gate passes; this one failure suffices.
Sigma1 also fails sanity, majority-improvement and band-crossing gates, with
depth-spread improvement2.7211%. Neither candidate qualifies.

The fixed samples0,32,64 are displayed for both candidates with identical
uint8 range0..255. Bands crossing the patch and sizable empty regions remain.
Some local tilts change their appearance without demonstrating consistent
surface seating across the whole model patch. Orthogonal cuts use the same
aspect-auto convention for all comparisons, equally exaggerating depth.
Native XY/XZ/YZ overlays mark the unchanged center, frozen tangent(cyan),
image longitudinal axis(yellow), and depth normal(red).

## Inference and coherence

SMOKE_INFERENCE=NOT_RUN; FULL_INFERENCE=NOT_RUN. The original progression
requires geometry improvement before smoke. Adding CT coverage does not waive
that gate. No CUDA/OOM failure occurred in this phase's scientific inference,
because no model inference was attempted. The renderer-only validation ran in
the pinned scientific container without GPU passthrough; its startup CUDA
notice was not a model failure and it completed the CPU sampling check.

There is no T11 probability stack, distribution, threshold-positive fraction,
connected-probability structure, or neighboring-Z probability statistic to
report. T11 inferred0of130,629centers. The planned full raster shape remains
50×4482, but no such T11 raster was generated. The two requested probability
PNG filenames contain explicit NOT RUN status notices, not probability data.

Sealed baseline mean/median adjacent-Z correlations, unchanged:

|Readout|Mean adjacent-Z|Median adjacent-Z|Mean local-window|
|---|---:|---:|---:|
|T5raw intensity|0.051516|0.059245|0.012505|
|T6profile mean|0.108094|0.138312|0.029243|
|T9pretrained old frame|0.002067|-0.003840|-0.000613|
|T11|NOT RUN|NOT RUN|NOT RUN|

The coherence comparison figure labels T11 NOT RUN. Existing T5/T6/T9 values
must not be interpreted as a measured T11 comparison. C means insufficient
patch-geometry improvement under this bounded protocol, not proof that every
orientation-only method fails or that the frozen thread is useless.

## Review, replay, verification and next step

Online independent-context reviewer `/root/t11_preflight_review` used the same
model/provider as the author, not an independent provider. Preflight defects
in replay/missing-data paths were fixed before evidence and documented.
The reviewer independently recomputed final acceptance predicates and gates
from all144records, confirmed frozen centers, tangent inputs, handedness and
depth sign, and confirmed C with no smoke/full inference. See the three
review records in the output directory.

All144frame estimates, metrics and rotated patch hashes replayed exactly.
Old and new arrays replayed pixel-identically; the141completed arrays from
the interrupted attempt also match the complete-support retry. Original
T9 controls match its sealed renderer exactly at all72points. Smoke/full
inference replay is NOT APPLICABLE, not a claimed model replay PASS.

Final verification checks3,307predecessor entries,147blocked-attempt/freeze
entries, all original scientific-source hashes, the support-only hash freeze,
six supplemental hashes, coordinates, normal sign, coverage and replay.
SHA256SUMS covers reports, source, receipts, initial/retry/replay arrays and
supplemental data. Lossless NPZs and supplemental CT remain local under existing
Git ignores. P0-T1throughT10, initial T11 evidence and original plan are unchanged.

Not ready for full corrected-frame inference. A future separately authorized
experiment could investigate why a stable small-neighborhood orientation does
not establish patch-wide seating or reduce repeated-sheet intersections,
without moving the frozen threads. This experiment does not perform such a
redesign, create a mesh, or scale to more slices.

Primary paths:

* `outputs/p0_t11/missing_ct_inventory.json`
* `outputs/p0_t11/supplemental_ct_manifest.json`
* `outputs/p0_t11/coverage_validation.json`
* `outputs/p0_t11/frame_metrics.json`, `frame_distributions.png`, `frame_overlays/`
* `outputs/p0_t11/t9_vs_t11_patch_comparison.png`, `patch_quality_metrics.csv`
* `outputs/p0_t11/smoke_metrics.json`, `metrics.json`, `coherence_comparison.png`
* `outputs/p0_t11/frame_replay.json`, `verification.json`, `SHA256SUMS`
* `outputs/p0_t11/complete_support_attempt/frames/patch_*.npz` (lossless local)
* `data/raw/p0_t11/native_chunks/` (six supplemental chunks, local)

Reproduction uses unchanged source modules and the frozen support wrapper.
`scripts/p0_t11_support.py validate` checks coverage/hashes. A fresh replay
must use a new output path; the existing exclusive replay directory is sealed
and must not be overwritten. Verify the complete seal from repository root
with `sha256sum -c outputs/p0_t11/SHA256SUMS`.
