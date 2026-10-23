# P0-T13 control-thread prerequisite

CONTROL_THREAD_PREREQUISITE_COMPLETE = YES  
INDEPENDENT_CONTROL_THREADS = PASS  
CONTROL_SLICE_COUNT = 7  
TRUSTED_OVERLAP = PARTIAL  
LEAKAGE_AUDIT = PASS  
READY_TO_RESUME_P0_T13 = NO

The bounded experiment is complete, but its prerequisite acceptance gate **fails**.
Seven independent CT-only threads were reconstructed and replayed exactly. A total
of 503 observed samples have post-freeze trusted associations. However, no central
five-slice window satisfies the frozen overlap-support criterion. No threshold was
relaxed, no alternate path selected after evaluation, and no T13 work was executed.

## Hypothesis and exact input

Hypothesis: independent maximum-recovery transverse threads in a small crop of the
T10 native scan would include enough of the trusted sheet across five consecutive
slices to support later local-frame evaluation. This tests control-input feasibility,
not T13 correspondence, ribbon quality, ink sensitivity or text recovery.

Trusted control: PHercParis4-w02, segment `20231031143852`, unchanged from P0-T10.
T10 cached a **surface-aligned rendered volume**, not native transverse slices.
Tracing that render would inherit its geometry. This experiment instead acquired
bounded native CT from the public scan asset matching the trusted mesh/render ID:

`https://vesuvius-challenge-open-data.s3.amazonaws.com/PHercParis4/volumes/20260411134726-2.400um-0.2m-78keV-masked.zarr/`

The sealed render metadata names its source as
`/volumes/esrf/20260311/2.4um_PHerc-Paris4_masked.zarr/`; the trusted mesh is
`20231031143852-on-20260411134726-2.4um.tifxyz`. No old 7.91 µm target scan was
substituted. Native storage is ZYX, shape 75784×32693×32693, uncompressed uint8
128³ chunks. Transverse slices are native XY planes at fixed Z, with 2.4 µm
native spacing. The public multiscale metadata records index scale; physical
spacing is bound by the matching scan identity and sealed T10 metadata.

**Frozen selected interval: native Z44453..44459 inclusive.** Seven slices were
chosen as the smallest count within the requested approximate 7–15 range. Central
slices 44455,44456,44457 each have Z±2 neighbors. Their total center span is 14.4 µm;
this differs substantially from the original T13 target's 126.56 µm slice pitch.
These controls do not establish performance at that larger spacing.

Selection used only the same T10 crop's approximate documented UV/grid bounds:
grid Y[819,922), X[1049,1153). Its trusted native XYZ bounding box is
(19167.207,15823.432,43383.137) through (21296.160,16447.305,45535.910).
The interval center is rounded median trusted Z44456. The entire projected crop
XY bounding box plus 128 native pixels was rounded outward to chunk boundaries:
**X[18944,21504), Y[15616,16640)**. This is one rectangular acquisition region,
not a sheet-shaped mask. Geometry selected acquisition bounds only, as permitted
by Task 1; it supplied no seed, branch, path, bridge or reconstruction coordinate.
The inherited T10 crop was originally label-selected; this experiment read no labels.

`selected_slice_range.json` was written before reconstruction and CT acquisition.
160 HTTP206 byte ranges supplied **18,350,080 bytes**, below the predeclared128MB
ceiling. Each range is seven native planes from one uncompressed chunk. Receipts
retain URLs, Content-Range, headers, sizes and hashes. No full segment was acquired.
Each slice is reduced by exact4×4 XY mean, half-even uint8 rounding, to256×640;
Z is not pooled. Trace XY spacing is9.6µm and arc sampling4.8µm. Native mapping:

`(X,Y,Z) = (18944 + 4*x + 1.5, 15616 + 4*y + 1.5, native_slice_index)`.

## Frozen reconstruction method

`plan.json` and `PRE_RECONSTRUCTION_FREEZE.json` bound code, parameters, input
hashes and acceptance thresholds before tracing. This is an implementation-time
hash preregistration, not a separately reviewed or pre-run Git commit.

Same Direct Thread-Stack philosophy as T2/T3, adapted to a local native-control
crop rather than reusing their old scan's seed and exterior box:

1. Independent same-slice DoG(σ0.6−σ2), response>3, skeletonization, maximal
   degree-two graph chains; exclude a four-pixel crop border. No ideal spiral.
2. Archive extracted chains with arc≥3px; shorter chains are excluded. Eligible
   search chains have no
   internal geometric contacts/crossings. Open a closed cycle by removing its
   duplicate-endpoint closing edge for search; preserve the opened graph chain.
3. Start both orientations of the eight longest eligible CT chains. No trusted
   reference point or fixed surface-aligned endpoint is given to the worker.
4. Maximize observed arc minus local link costs with beam128, depth64, up to12
   lowest-cost outgoing links per state. Bridge distance≤12XYpx (115.2µm), tangent
   cosine≥−0.5, both connector cosines≥−0.35. Cost is2×distance +3×tangent turn
   +2×connector bends +4×normal miss/4px +2×unsupported-gap term. Local DoG support
   is secondary CT evidence. This bounded search is not a global maximum.
5. Forbid source-chain reuse and reject extensions introducing proper crossings,
   non-adjacent point contacts or collinear overlaps. This prevents topology
   defects during search rather than repairing paths using trusted geometry.
6. Preserve observed edges as1, inferred nonzero bridges as2, and weak-link
   observed neighborhoods as3. Code3 applies to the first5px after tangent
   cosine<0.7 or normal miss>1.8px. Sample at0.5XYpx along continuous Euclidean
   arc, including bridges; bridge CT remains NaN/masked. Low-confidence observed
   CT remains available. No later coordinate, order or mask correction occurs.

The local seed and objective intentionally need not select the annotated sheet.
Success here means independent usable contour geometry, not proven single-sheet
identity. Raw CT brightness is never treated as ink.

|Native Z|Observed arc mm|Bridge arc mm|Observed samples|Bridges|Near trusted samples|Near fraction|
|---|---:|---:|---:|---:|---:|---:|
|44453|14.445|0.867|3011|19|0|0%|
|44454|9.823|0.555|2048|10|110|5.371%|
|44455|13.044|0.895|2719|16|170|6.252%|
|44456|8.217|0.276|1712|7|115|6.717%|
|44457|11.604|0.890|2419|19|24|0.992%|
|44458|9.257|0.499|1930|13|84|4.352%|
|44459|7.774|0.330|1620|8|0|0%|

Every final thread has zero proper crossings, non-adjacent contacts and collinear
overlaps, independently checked using the sealed T3 exhaustive segment-pair auditor.
Geometric cleanliness does not establish physical sheet identity.

## Leakage audit and preserved infrastructure failure

Every worker ran in a pinned Docker image with no network, read-only root,
all capabilities dropped and no-new-privileges. The only host mounts were one
CT slice, the worker code, Python libraries and a fresh output directory. Neither
trusted geometry nor labels, probabilities or neighboring CT slices were mounted.
Python file-open auditing records one `/input.npy` read per successful worker.

A second exact replay additionally ran under `strace` to record file-access
syscalls, including interpreter/native-library accesses. It mounted only the
additional strace executable. Raw syscall logs, successful read-only file paths,
container configurations, source/library/runtime hashes and exact input hashes are
retained. All21 successful runs (seven primary, seven reverse-order replay, seven
instrumented replay) satisfy isolation and Python-open checks. Instrumented replay
is byte-identical at the array/decision level. No geometry/label reads occurred.

The first worker failed at output-directory creation before `np.load`: dropping
capabilities prevented container root from writing a host-owned directory. Its
configuration/log and original launcher remain under `attempt_01_permission_denied/`.
The infrastructure amendment runs as the host UID/GID; CT, scientific code,
parameters and isolation remain unchanged. No contaminated scientific attempt
occurred. The original launcher hash is verified against its archived copy.

Selection/evaluation are separate host stages and intentionally can read trusted
geometry. Evaluation checks `independent_threads_FROZEN.json` before its first
geometry read. No thread is changed based on evaluation. The audit claim concerns
this enforced stage separation, not absence of geometry from the entire experiment.

## Post-freeze registration and trusted association

Registration uses the known crop/pooling transform into the matching native scan's
XYZ frame; no ICP, fitted displacement or coordinate repair is performed. Maximum
forward/inverse arithmetic residual is **6.536993168992922e−13 native voxel**.
This is numerical transform consistency, **not physical registration accuracy**.
An independent producer/fiducial registration-error bound is unavailable and is
stored as null. It must not be reported as zero uncertainty.

The already trusted local grid is triangulated for evaluation only. Each slice's
Z plane intersects206 trusted triangles. For every frozen sampled point, the
nearest point on the union of these same-Z intersection segments is calculated.
This is exact distance to the piecewise-linear trusted transverse section, an
upper bound on unrestricted 3-D nearest-surface distance. It is not a new scroll
surface reconstruction. Neither the triangles nor their intersection curves are
input to tracing or a T13 correspondence algorithm.

Trusted normals come from grid finite differences, interpolated at the association.
Longitudinal direction is normalized scan-Z projected into that trusted tangent
plane; the original trusted V direction is saved separately. Frame association
requires projection norm≥0.25. Normal sign is the producer-grid cross-product
convention, not a calibrated physical ink-side sign. Later angular evaluation must
handle axial sign explicitly.

**503/15,459 observed samples (3.253768%)** lie within19.2µm with valid frames.
Within the five slices having associations, median near-surface distances are
4.150,8.968,5.642,4.057,8.254µm; maximum19.2µm threshold is not exceeded.
All-observed median distances range247.851–695.668µm, showing that most recovered
material is elsewhere. Near-point median endpoint-normal variation is1.26–1.60°;
this is a local interpolation diagnostic, not calibrated angular truth uncertainty.

Saved uncertainty includes native voxel2.4µm, XY pooling-cell half diagonal6.788µm,
trusted UV grid pitch20native parameter units, per-point normal variation, measured
surface distance, and explicitly unknown producer physical registration error.
Arc span is min-to-max near-sample s, not necessarily continuous overlap length:
Z44455's9.830mm span includes separated crossings and must not be called one long
matched sheet. Proximity itself cannot prove physical identity.

The fixed three central overlays show independently selected paths on multiple
bands, with localized trusted overlap. Cyan path selection changes substantially
between slices. The visual observation agrees with sparse measured overlap; no
acceptance rests solely on appearance. The separate provenance overview shows
observed, inferred and low-confidence parts without evaluation geometry.

## Frozen acceptance gate and control input

The qualitative meaningful-overlap requirement was operationalized **before
reconstruction** as ≥25 observed points within19.2µm on each of five consecutive
slices, and≥100µm central near-point arc span, with finite trusted frames. Other
conditions require≥5 independent threads, a central Z±2 window and leakage PASS.
These are feasibility thresholds, not calibrated confidence or power estimates.

All requirements except meaningful five-slice overlap pass. Near-point counts by
central window are:

- Z44455: [0,110,170,115,24] — FAIL.
- Z44456: [110,170,115,24,84] — FAIL (24<25).
- Z44457: [170,115,24,84,0] — FAIL.

The middle window is a threshold near-miss; do not interpret this as evidence that
24 points are physically unusable while25 would prove identity. Nevertheless,
changing25 to24 after seeing results would weaken the frozen gate and was not done.
No window has been certified to support later correspondence evaluation.

`control_input_manifest.json` retains the partial central subsets (170,115,24
points) and references every frozen neighboring independent thread, with complete
XYZ coordinates, trusted positions/normals/longitudinal/V directions and uncertainty
in `registered_evaluation/`. It explicitly says gate FAIL and readiness NO.
No cross-slice correspondence was chosen, no local ribbon was built, and no model
was loaded. Partial inputs are evidence, not an accepted T13 prerequisite.

## Replay, preservation, publication and next recommendation

All seven graph/search/decision/thread results replay exactly twice, including
reverse slice order and syscall instrumentation. Native HTTP-range payloads
reassemble the inputs and independent pooling matches exactly. Independent
four-corner bilinear CT checks have maximum error4.2633e−14. Registration arrays,
unit frames, provenance, masks and topology all verify.

**4,002 predecessor entries remain byte-identical**, including T1–T12, the blocked
T13 checkpoint and every `8fdf4cd` evidence file. Their manifests were not changed.
Execution verification PASS is separate from scientific prerequisite-gate FAIL.

Artifacts are under `outputs/p0_t13_control/`: all requested selection, metrics,
overlap/leakage audits, manifest, verification and SHA256SUMS; primary threads;
CT-only inputs; two replay directories; raw access logs; frozen plans and source
hashes; acquisition receipts; and diagnostic PNGs. Native payloads remain under
`data/raw/p0_t13_control/`, checksummed. Small CT-only input arrays, primary thread/
graph arrays and registered associations are explicitly published despite generic
array ignore rules. Redundant replay arrays and downloaded range payloads remain
local and sealed; `git_artifact_inventory.json` records the split.

Entry points: `scripts/p0_t13_control_acquire.py`, `_reconstruct.py`, `_run.py`,
`_audit_run.py`, `_evaluate.py`, `_verify.py`, `_seal.py`. Generation stages refuse
existing destinations. Do not rerun them over sealed evidence; use new isolated
output paths for a future replay. Read-only seal check from repository root:

```sh
sha256sum -c outputs/p0_t13_control/SHA256SUMS
```

The separately published commit contains this completed negative/partial
prerequisite milestone, not a T13 result. Commit/remote identity is reported in
the final response, with no force push.

Next recommended decision, **not executed**: externally review whether a new,
prospectively specified CT-only selection method or overlap-support criterion is
appropriate given the sparse/unstable path coverage and short physical Z baseline.
Do not select a better branch using these trusted associations. Do not automatically
resume T13. Stop after this prerequisite.
