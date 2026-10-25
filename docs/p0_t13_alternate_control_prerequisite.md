# P0-T13 alternate trusted-control prerequisite

```text
ALT_CONTROL_COMPLETE = YES
ALT_CONTROL_RESULT = FAIL
SELECTED_REGION = native X[17280,19584), Y[17664,19072), Z51447..51453
CONTROL_SLICE_COUNT = 7
BEST_5_SLICE_WINDOW = 51449..51453
BEST_WINDOW_COUNTS = [0, 0, 0, 0, 133]
TRUSTED_OVERLAP = PARTIAL
LEAKAGE_AUDIT = PASS
READY_TO_RESUME_P0_T13 = NO
```

This single alternate-region experiment is complete. 0 of
3 five-slice windows pass the unchanged count gate. No second
region was reconstructed or selected after overlap evaluation. No P0-T13
correspondence, ribbon, frame-error experiment or model inference was executed.
Prior threads, thresholds and all sealed evidence remain unchanged.

## Hypothesis and prospective neutral selection

Hypothesis: an alternate interior region of the same trusted segment might yield
independently recovered CT threads with sufficient local trusted overlap on every
slice in one five-slice window. This is control-input feasibility, not proof of
physical sheet identity, cross-slice correspondence, coherent reconstruction or ink.

The selection policy was written and hash-frozen **before reading candidate geometry
for selection**, acquiring CT intensities or reconstructing threads. Candidate
centers were fixed by fractions of the trusted grid dimensions: quarter-height /
mid-width, three-quarter-height / mid-width, mid-height / quarter-width, mid-height /
three-quarter-width. Take the first neutral-eligible candidate, never rank by overlap.
Only the first candidate was needed; all later candidates remained unexamined.
The policy and selected-region/source hashes were frozen before reconstruction;
this is local hash preregistration, not a pre-run Git publication.

Use PHercParis4-w02 / trusted segment20231031143852 and the same sealed T10
native-XYZ geometry. Selected trusted grid bounds are
**Y[376,479), X[1048,1153)**, with a one-grid-point derivative halo. Selection gates:
all positive coordinates; adjacent grid-edge lengths0.01..200 native voxels;
maximum neighboring normal angle≤45°; one connected transverse intersection of
arc≥500 native voxels on each of seven consecutive slices; section-to-acquisition-
crop margin≥64 native voxels; separation≥1024 native Z slices from the old range;
and every required CT chunk present at expected size.

Measured geometry: edge lengths **15.912377..21.152258**
native voxels; maximum neighboring normal angle **18.461978°**.
All seven intersections are single connected polylines with208 segments, arc
2042.734535..2043.204250
native voxels (approximately4.9mm), and minimum CT crop margin
165.966457 native voxels (approximately398µm).
These are continuity/data-selection checks, not an independent calibration of
physical registration accuracy.

The slice center is rounded median native Z of the candidate grid,51450;
**Z51447..51453** gives the smallest seven-slice range in the allowed7–15 bound.
The XY region is the candidate's projected bounding box plus128 native voxels,
rounded outward to128-voxel boundaries. It is a rectangular crop, not a sheet mask.
Trusted geometry determines region bounds only: it supplies no worker seed, path,
branch, bridge, topology correction or coordinate adjustment. Neither ink labels,
model outputs nor thread-overlap counts are used in selection. Prior artifact
hashing is integrity checking, not a region-selection input.

## Native CT and unchanged reconstruction

Native source:
`https://vesuvius-challenge-open-data.s3.amazonaws.com/PHercParis4/volumes/20260411134726-2.400um-0.2m-78keV-masked.zarr/`

This is the same scan asset20260411134726 as the trusted T10 surface, rather than
the old7.91µm target scan or T10's geometry-resampled surface-volume planes.
Native storage is ZYX, transverse slices XY, spacing2.4µm. All198 needed
uncompressed128³uint8 chunks passed HEAD size2097152. HTTP206 byte ranges acquired
only the seven selected native planes: **22,708,224 bytes**, below128MB.
URL/range/header/hash receipts are retained. No full segment acquisition occurred.

Unchanged exact4×4 XY mean and half-even uint8 rounding produce9.6µm XY sampling;
Z remains native2.4µm, seven-center span14.4µm. Arc sampling remains4.8µm. This is a
shorter Z baseline than the original T13 target grid; no transfer of longitudinal
performance to126.56µm target spacing is established.

The worker is byte-identical to commit576b562:
`scripts/p0_t13_control_reconstruct.py`, SHA256
`5bba58e70be384ccd7d5830ae3870d57607a235a7a4df73ea8a0ab18c58c354e`.
All414 prior host code/input/library hashes match. The pinned container image,
CPU-only execution and deterministic thread settings are unchanged.

No ridge, graph, search, bridge, topology or provenance parameter changed:
DoG(.6,2)>3; maximal skeleton chains with arc≥3px; both orientations of eight
longest eligible chains; beam128/depth64/outgoing12 maximum-recovery search;
≤12px bridges with the original tangent/chord/cost rules; reuse/intersection
guards; provenance1 observed,2 inferred bridge,3 low-confidence; and0.5XYpx arc
sampling. Bridges retain geometry/arc but NaN/masked CT. Search is bounded rather
than globally optimal, and geometric cleanliness does not prove physical identity.

## Independence, freeze and replay

Each worker runs with no network, read-only root, all capabilities dropped and
no-new-privileges. Only its single CT slice, frozen source, libraries and fresh
output are mounted. Trusted coordinates, normals, labels, prior overlap scores,
probabilities and neighboring CT are absent. Python-open audits cover every
successful primary/replay process; syscall file-access auditing covers a second
exact replay, including native-library reads. File paths and hashes, mount
configurations and image identity are saved. Leakage audit PASS.

All seven threads were hash-frozen before post-reconstruction trusted association.
Both complete replays reproduce every graph, path, decision, mask/provenance and
numeric array exactly. Native range reassembly and independent pooling are exact.
Independent exhaustive segment-pair audits find zero proper crossings,
non-adjacent contacts and collinear overlaps on every thread. Independent
four-corner bilinear and arc/coordinate checks pass. The worker has never been
retuned or rerun with selected alternative parameters.

## Registration, uncertainty and all window results

After freeze, use the unchanged trusted-section and nearest-point helpers on the
selected grid. Native mapping is
`XYZ=(17280+4*x+1.5,17664+4*y+1.5,nativeZ)`; no fitted transform, ICP, corrected
thread or alignment based on overlap. Distance is to the exact nearest point on
the piecewise-linear trusted surface's same-Z intersection. This is an upper bound
on unrestricted3D nearest-surface distance. The trusted grid is tessellated only
for evaluation, never used as a reconstructed thread.

A qualifying point is observed/unmasked, within **19.2µm**, and passes the unchanged
trusted-frame eligibility (scan-Z projected tangent norm≥0.25). Distance-only counts
are also saved. Trusted normal, projected longitudinal direction, producer V
direction, surface position and local interpolation variation are retained only
for later evaluation. No cross-slice correspondences or local ribbons are chosen.

Maximum numerical coordinate roundtrip residual:
**4.5474735088646412e-13 native voxel**.
Physical registration uncertainty remains **uncalibrated/null**. Numerical residual
is not physical accuracy. The original pooling, mesh-interpolation and producer
registration uncertainty caveats remain; no new calibration project was attempted.

|Native Z|Observed samples|Qualifying samples|Qualifying fraction|
|---|---:|---:|---:|
|51447|1255|0|0.0000%|
|51448|2912|24|0.8242%|
|51449|1014|0|0.0000%|
|51450|1777|0|0.0000%|
|51451|769|0|0.0000%|
|51452|755|0|0.0000%|
|51453|2116|133|6.2854%|

Total descriptive overlap is157/10598 observed samples
(1.4814%). This pooled fraction is not an acceptance criterion.
Near-point s span is min-to-max and may join separated overlap islands; it is not
necessarily continuous single-sheet support.

|Five-slice window|Counts|Minimum|Summed deficit|Pass|
|---|---|---:|---:|---|
|51447..51451|[0, 24, 0, 0, 0]|0|101|False|
|51448..51452|[24, 0, 0, 0, 0]|0|101|False|
|51449..51453|[0, 0, 0, 0, 133]|0|100|False|

Every consecutive five-slice window was tested. PASS requires≥25 qualifying
points on **each** of its five slices, without averaging. The inherited central
100µm span safeguard is recorded;25 unique4.8µm samples already imply≥115.2µm span,
so it cannot turn a count PASS into FAIL. Best-window selection was frozen:
earliest PASS; otherwise largest minimum count, smallest summed deficit, largest
total and earliest start. Best counts are **[0, 0, 0, 0, 133]**, with
4 failing slices and summed deficit100.
The gate was not weakened or reinterpreted.

Qualitative inspection: The first and middle displayed independent paths follow CT bands spatially separated from the yellow trusted section. The last displayed path has several localized green overlap islands, separated along its arc. The images do not establish a persistent same-sheet thread across the seven slices; no ink or text interpretation is made.
These observations do not identify ink and do not override numeric acceptance.

## Preservation, artifacts and stopping point

**5009 predecessor entries are unchanged**, including
P0-T1 throughT12, blocked T13, the original control prerequisite, its seven threads,
and the bounded extension. All15 prior root seals
verify. Execution verification PASS is separate from scientific prerequisite result.

Outputs: `outputs/p0_t13_control_alt/`, including selection policy/audit and hash
freezes, `selected_region.json`, acquisition receipts, `independent_threads/`,
`thread_metrics.csv`, `per_slice_overlap.csv`, `five_slice_windows.csv`,
`trusted_overlap_audit.json`, `leakage_audit.json`, `registered_evaluation/`,
`verification.json`, both replays, file-access hashes/logs, visualizations and
`SHA256SUMS`. A `control_input_manifest.json` is emitted only on PASS and contains
frozen thread references plus evaluation-only trusted frames for the selected
window; it does not authorize or execute T13.

Small CT-only inputs, primary arrays and registered associations are published;
raw range payloads and redundant replay NPZs remain local but sealed. The Git
artifact inventory records that distinction. Code entry points are
`scripts/p0_t13_control_alt_select.py`, `p0_t13_control_alt.py`,
`p0_t13_control_alt_evaluate.py`, `p0_t13_control_alt_verify.py`, and
`p0_t13_control_alt_finish.py`. Stage creation refuses existing output directories;
use fresh isolated destinations for any later replay, not these sealed paths.
Read-only verification from repository root:

```sh
sha256sum -c outputs/p0_t13_control_alt/SHA256SUMS
```

The result is committed and pushed separately to origin/main, without force.
The exact commit/remote identity is reported in the execution response. Next
recommendation: let the user and ChatGPT assess this sealed feasibility result
before authorizing any next experiment. Do not automatically search another region
or execute P0-T13, regardless of prerequisite status. Stop here.
