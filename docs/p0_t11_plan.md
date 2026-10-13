# P0-T11 frozen-point local 3-D orientation experiment

Hypothesis: orientation-only correction at frozen T5/T9 centers improves
surface-following patch geometry, then pretrained ink coherence. No paths,
arc-length coordinates, bridges, source masks, provenance, row ordering,
source CT, checkpoint, normalization or depth offsets change. No surface mesh,
target-region trusted geometry, OCR, learned geometry, or local LLM.

## Fixed design before target evidence

72 geometry-only sample identities in samples.json: groups Z128..130,
151..153,175..177, eight evenly indexed common valid s stations per group.
No CT/model output was used to select them. All masks remain frozen.

Only two candidates: Gaussian-derivative structure tensors with native-voxel
sigma1 and sigma2, in that order. Both use a33³ cube centered at the exact
native point, Gaussian spatial weights sigma6. Tensor gradients are XYZ.
Find the dominant normal within the plane perpendicular to the exact saved
XY tangent. This is a constrained local plane orientation, not a surface fit
or relocation. No smoothing along paths or across Z. Sign aligns with saved
T9 normal. Let t_s = saved tangent, t_z = n cross t_s. The geometric frame is
right-handed. Rendering uses image Y = -t_z because T9's stored array frame
is left-handed: t_xy cross scanZ = -old_normal. Thus zero-tilt rendering is
identical to T9, without an unintended image or depth reversal.

Confidence=(largest-smallest constrained eigenvalue)/largest. Reject if
confidence<.25, tensor energy<=1e-6, tilt>75degrees, unconstrained normal's
angle out of the plane perpendicular to tangent>30degrees, or a +/-2native
voxel Z displacement changes estimated axial normal by>20degrees. These
displacements test estimator stability only; final P never moves.
Sanity requires>=70%accepted samples and median perturbation<=15degrees.
Record all candidates, rejects and finite distributions, compare to T10
diagnostically only. Sparse-station normal angles are not local derivatives.

## Geometry gate (no ink-model evidence)

Render original and new 17×128×128 patches with9.6/7.91 native step,
same trilinear interpolation, half-even uint8 rounding, centered offsets.
For held-out patch wings Y0..39 and88..127 (outside the estimation cube),
compute smoothed raw-CT longitudinal gradient energy divided by longitudinal
plus depth energy; depth argmax standard deviation and mean |argmax-8|/8.
Other proxies: zero-voxel fraction, mean number of center-plane X peaks with
prominence20 at Y16,64,111, and occupancy>32 in depths6..10.
These are geometry/occupancy proxies, NOT ink or a registered surface oracle.

For each estimator on its accepted samples, require ALL:
sanity; median longitudinal gradient fraction reduction>=10%; median depth
peak spread reduction>=10%; >=60%points improve longitudinal gradient;
median empty fraction increase<=.03; median normalized center-distance
increase<=.05; median band crossings do not increase. Pair old/new on exactly
the same accepted sample IDs and report rejected samples separately.
Require all72patch pairs available before opening inference. Choose first
candidate passing all gates, not the best probability/image.
Both candidates run even if the first passes. Missing CT is explicit, never
filled with invented voxels or treated as absence of anatomy.

If neither passes sanity: D, stop model inference. If a candidate passes
sanity but neither passes geometry: C, stop model inference. In either case,
ink/probability changes are NOT ASSESSED. Requested stack PNGs are explicitly
marked NOT RUN placeholders, not fabricated predictions. Preserve diagnostics
and seal the completed gate-limited experiment.

## Conditional smoke and full inference

If geometry passes, use exact seed42 step-075000 checkpoint hash
e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab,
official preprocessing, float16 autocast, batch8, mean central2×2 output.
Run all72old/new patches, using unchanged old-frame fallback for rejected
estimates consistently with full inference. Require finite [0,1] probabilities and nonconstant
new centers; report distributions, six neighboring-Z correlations (eight
fixed stations each), and patch-map persistence. No method selection by ink.
Only then run all valid centers Z128..177. Rejects use their unchanged T9
frame as an explicitly recorded fallback; never silently remove centers or
change the source mask. Report correction and fallback fractions separately.
No independently stretched rows, probability-driven geometry, or bridge inference.

Full comparison uses frozen T5/T6/T9 grids, mean/median neighboring-row
Pearson correlation, 64column local-window correlations with>=16paired
finite samples, fraction p>.5,4-connected components and vertical spans.
Coherence success requires T11 mean adjacent correlation>=T9+.05 AND local
correlation>=T9+.03, plus visibly persistent nonrandom structure, not isolated
high values. Full result A if geometry and coherence improve; B if geometry
passes but coherence does not. Report valid-center and padded-raster fractions.

## Replay and preservation

Before evidence, verify all T10-manifest and predecessor hashes; commit the
reviewed design/source. Repeat all72×2 frame/patch estimates into a distinct
replay directory; require identical numeric frame records and raw patch hashes.
If inference runs, replay every smoke patch with maximum probability error
<=.002; for full inference also verify every saved coordinate/mask and replay
one deterministic first accepted center per row. These are bounded model
replays, not a second full run. Hash raw arrays locally even if ignored by Git.
Record all gates, limits, resource state, review disposition and exact scope.
Commit/push the sealed T11 result to origin/main without force.
