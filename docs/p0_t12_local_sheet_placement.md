# P0-T12 — bounded local sheet placement

## Outcome

**P0_T12_RESULT = D (operational LOCAL_SHEET_IDENTITY_AMBIGUOUS).**
LOCAL_SHEET_FIT=FAIL; PATCH_GEOMETRY_GATE=FAIL; SMOKE_INFERENCE=NOT_RUN;
FULL_INFERENCE=NOT_RUN; INK_STRUCTURE_VISIBLE=NOT_ASSESSED.

This bounded detector/placement method did not establish reliable local sheet
identity:0/72target fits and0/3trusted-control fits were accepted. It does NOT
prove intrinsic CT ambiguity, inevitable layer switching, or a need for global
tracking. The control failures mean this placement procedure is not validated
on known-good geometry. The validated T10 ink model is not implicated.

The experiment is complete at its preregistered early-stop gate. No model was
loaded, no ink inference was performed, no new surface-following target patch
was rendered, and no displacement was applied. Rejected placements and every
candidate record are preserved. No gate or threshold was relaxed.

## Hypothesis, immutable inputs and design

Hypothesis: small local centering and plane/low-curvature fitting around frozen
thread anchors could improve model-input geometry beyond orientation alone.
The exact72T11targets and three T10control patch locations were retained.
Target groups are Z128..130,151..153,175..177, eight common-valid s stations per
group, with native P=(16*x_level4,16*y_level4,16*Z). Exact sample identities,
coordinates, tangents, normals and indices are in `outputs/p0_t12/samples.json`.
All primary thread coordinates, paths, masks, bridges, provenance, arc length
and ordering remain frozen. No neighboring-Z evidence alters a thread.

The reviewed plan/source/input contract was committed as `7b16ff9` before
target fitting; all83hash-bound entries match. Original T11sigma2orientation
is used for its54accepted targets; the other18retain the original T9frame.
Saved raw17×128×128patches supply the CT evidence. Native source remains
`54keV_7.91um_Scroll1A.zarr`,7.91µm, with9.6µm model-space sampling, original
depth/image signs, trilinear float32 interpolation and half-even uint8 rounding.

The intended scientific checkpoint remains
`scrollprize/ink_9um/hybrid_3d2d-seed42/step-075000.pth`, with unchanged official
preprocessing. It was not loaded in this stopped experiment. No checkpoint
search, training, OCR, language priors, global mesh, long-range tracking, or
local LLM was used.

## Preregistered local-placement methods

Central6×6depth profiles are Gaussian-smoothed(sigma0.75model voxel). Local
intensity peaks require prominence≥8, with bounded parabolic peak refinement.
The nearest peak to depth0must be within±2model voxels(19.2µm) and uniquely
closer by≥1depth unit than a competitor. No strongest-distant-peak substitution.
All candidate peaks, heights, prominences and reasons are retained.

At9×9fixed local stations(X,Y=15,27,...111), the same detector uses3×3profile
averages. Require≥57/81unique assignments near the central anchor and≥20sites
in each checkerboard train/holdout split. A local plane is fit by three fixed
weighted least-squares rounds. Displacement, residual, competing-ridge,
neighbor-continuity and slope guards are evaluated without model feedback.

A bounded quadratic is eligible only if unique sheet evidence is adequate
and the plane's heldout median residual exceeds0.15or95th exceeds0.35depth
units. It retains all admissibility guards and must improve heldout median
residual≥20%. Max surface/center displacement remains19.2µm; normalized
curvature is bounded. These sharper adequacy triggers were fixed using only
synthetic curved-sheet tests before target evidence, after reviewer feedback;
the user's10%geometry gate never changed. Synthetic tests exercised an
accepted curved fit and three critical guard mutants were correctly rejected.

The renderer represents only an independent local height graph, using local
normals for depth sampling.9.6µm parameter-grid/depth spacing is retained;
this does not claim exact isometric arclength on a curved graph. No target
fit reached rendering, so no such geometric distortion was applied here.

The two policies were plane-only and plane-then-eligible-quadratic. Neither
selected a method. Rejected fits retain the T11-policy patch for all-target
metrics; they are not omitted. No target quadratic fit was attempted because
the prerequisite local identity evidence failed first. Consequently this is
not empirical evidence against all possible low-curvature surface following.

## Fit rejection and coverage

| Target outcome | Count |
|---|---:|
|No candidate within19.2µm displacement bound|31|
|No detected central peak|18|
|Insufficient unique local ridge evidence|23|
|Accepted local placements|0|
|Reliable identity guards passed|0|
|Quadratic fits attempted|0|

The23targets with a central candidate failed the required57/81local unique
ridge assignments, before coefficient fitting. Thus accepted surface residuals
and fitted-center displacement distributions are undefined, not zero.
Median/max APPLIED displacement are both0µm because every fit was rejected.
CSV fallback centers equal the original anchors and are explicitly marked
accepted=False; they are not fitted-sheet estimates.

Coverage was checked for every accepted placement: zero accepted placements,
zero new support requirements, zero missing chunks. This is a vacuous accepted-
placement coverage pass, NOT a claim that72new fitted sheets were rendered or
that future full-run coverage is complete. No data were acquired, clipped,
padded or replaced. Existing T9/T11input patches and their source/coordinate
contracts are hash-bound. All72original and orientation-policy patches remain
available and are used in comparison and replay.

## Exact geometry metrics and gate

The unchanged T11quality function computes depth-argmax spread and center
distance on fixed Ywings; empty fraction; center-plane repeated-band peak
counts; longitudinal/depth gradient-energy ratio; and center-slab occupancy.
These intensity proxies are not physical sheet identity or ink labels.

All medians below use the entire unchanged72target set, including rejected
fit fallbacks. This differs from T11's quoted54accepted-estimate subset.

|Metric|T9|T11policy|T12(both policies)|T12−T9|
|---|---:|---:|---:|---:|
|Depth-peak spread(indices)|5.621225894654443|5.575175322509516|5.575175322509516|-0.046050572144927|
|Empty-volume fraction|0.132884306066176|0.122274959788603|0.122274959788603|-0.010609346277574|
|Repeated-band crossings|4.333333333333333|4.333333333333333|4.333333333333333|0|
|Normalized center-distance proxy|0.61728515625|0.61439208984375|0.61439208984375|-0.002893066406250|
|Longitudinal gradient fraction|0.368776169269553|0.320013873695020|0.320013873695020|-0.048762295574532|
|Center-slab occupancy|0.720922851562500|0.732568359375000|0.732568359375000|+0.011645507812500|

Depth-spread reduction=0.8192264998408862%, below the required10%.
Empty fraction reduces7.983897114449354%relative; band crossings do not
increase. The local-persistence proxy(longitudinal gradient fraction) reduces
13.222734991558005%, but does not replace the frozen primary criterion.
Every T12−T11delta is EXACTLY0: all improvements over T9 are inherited T11
fallback performance, not local-placement success.

Both policies fail:

* Reliable identity:0/72, required≥70%.
* Accepted fits:0/72, required≥70%.
* Depth-spread improvement:0.8192265%, required≥10%.

Band-crossing and empty-fraction nondegradation checks pass. The fraction
passing all local fit/quality checks is0. No alternate control-derived gate
was introduced. With reliable identity below70%, the unchanged classification
rule returns operational D rather than C. No accepted movement exists on
which to establish or refute physical layer jumps.

## Trusted-control comparison

The same detector, fit functions and CT-quality metrics were applied to the
three frozen T10patches at local origins(128,128),(192,192),(256,256) in the
trusted PHercParis4-w02segment20231031143852crop. No control geometry was
transferred to the target and no thresholds were calibrated from these results.

All three control fits rejected: one central candidate outside the displacement
bound; two insufficient local ridge support(33/81and32/81). Their raw CT
mean-profile peak lies at model depth4(relative−4), which is not the same
quantity as the annotated ink plane or a physical papyrus-surface oracle.

|Control metric|Three values|Median|
|---|---|---:|
|Depth-peak spread|5.078053740,4.621231605,4.727582793|4.727582793|
|Empty fraction|0,0,0|0|
|Band crossings|2,2.666666667,2|2|
|Normalized center distance|0.581530762,0.634997559,0.595190430|0.595190430|
|Longitudinal gradient fraction|0.163870031,0.131957980,0.207914891|0.163870031|
|Center-slab occupancy|0.990417480,0.994116211,0.985668945|0.990417480|
|Accepted center/surface-fit residual|Not available for all3|Not available|

Known-good rendered/model-positive inputs can fail this local brightness-ridge
detector. Thus the placement method itself remains unvalidated. Control
depth-argmax spread is also nonzero: it cannot be treated as a true surface
residual. These observations limit interpretation; they do not permit a
post-hoc change of the detector, displacement bound, sample set or gate.

## Visual inspection and inference stop

Fixed targets0,32,64are shown as T9, T11and explicitly labeled T12fallback
at identical0..255display scaling. T11and T12arrays are byte-identical.
Repeated bands and empty areas remain; there is no fitted T12surface image
to judge. Control plots show individual values(n=3), not a fabricated dense
reference distribution. Model XYand depth-Xcuts use identical aspect-auto
conventions across methods.

No pretrained model was loaded. No smoke, probability map, neighboring-Z ink
metric, full stack, or model replay exists. No probability images were
fabricated. Ink visibility is NOT_ASSESSED. The geometry gate blocks all
conditional model stages, even though T10previously validated the model.

## Replay, review and preservation

All72target detector/fit records and three control records replay exactly.
All72×2policy patch constructions, full metric records and gate replay exactly
into distinct replay directories. Each patch equals its T11fallback and
retains its exact frozen coordinate identity. Rejected candidates and null
fit residuals remain available in `sheet_candidate_audit.json` and the
lossless local NPZs; replay did not manufacture accepted fits.

Independent-context reviewer `/root/t12_preflight_review` used the same online
model/provider as author, not an independent provider. It approved the
prospectively corrected CT-only plan and independently confirmed the final
counts, exact frozen hashes, failed gate and operational-D qualification.
The reviewer explicitly required disclosure that all controls also rejected.

Verification checks3,807predecessor entries byte-identical and83pre-evidence
hashes unchanged. SHA256SUMS includes report/source, frozen inputs, candidates,
all target/control/replay records and locally retained lossless patches.
No T1–T11artifact, original T11blocked attempt, or thread geometry was altered.

Not ready for ink inference or full-scale placement. A future separately
authorized experiment would first need a local sheet-identity/placement
measure that succeeds on trusted control geometry before interpreting target
rejection as a property of the CT. No such redesign or threshold change was
performed here.

## Artifact paths

Report: `docs/p0_t12_local_sheet_placement.md`.
Under `outputs/p0_t12/`:

* `local_fit_metrics.csv`, `centering_offsets.csv`
* `sheet_candidate_audit.json`, `control_metric_distributions.png`
* `t9_t11_t12_patch_comparison.png`, `patch_quality_metrics.json`
* `geometry_gate.json`, `decision.json`, `coverage.json`
* `fit_replay.json`, `patch_replay.json`, `verification.json`, `SHA256SUMS`
* `fit/`, `fit_replay/`, `geometry/`, `geometry_replay/` retain complete evidence

No smoke/full artifacts are expected when those stages are NOT_RUN. Existing
replay directories are exclusive and sealed; do not overwrite them. Verify
the seal from repository root with `sha256sum -c outputs/p0_t12/SHA256SUMS`.
