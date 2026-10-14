# P0-T12 preregistration: local sheet placement

Hypothesis: bounded local placement/surface following improves patch geometry
beyond frozen T11 orientation. T11 results are disclosed design input, not new
T12 evidence. Preserve every T1–T11 byte, initial blocked T11 attempt, and the
exact72samples. Original P and primary XY thread never change; a separate
model-patch center may shift along its normal. No global mesh, long-range
tracking, OCR, training, local LLM, or model-driven selection.

## Input and bounded methods

Use saved T11 sigma2frame where accepted(54/72); otherwise original T9 frame.
Use the corresponding saved17×128×128 raw CT patch for ridge evidence. Retain
image/depth signs and9.6/7.91 native mapping. T9 and T11 raw patches are exact
saved controls, not rerun with new preprocessing. Native source remains
54keV_7.91um_Scroll1A level0. No sample expansion.

Central6×6profile: smooth in depth with Gaussian sigma0.75model voxel, find
local intensity peaks with prominence≥8uint8units, refine each peak by a
parabolic offset clamped±0.5. Choose the uniquely closest peak to depth0,
requiring displacement≤2model voxels(19.2µm) and distance margin≥1over the
next candidate. No peak, excessive displacement or competing near-equal
ridge means ambiguity/rejection, never switch to the strongest distant ridge.
CT brightness is a sheet-evidence proxy, not ink evidence.

At a fixed9×9grid(X,Y=15,27,...111), use3×3profile averages and the same peak
detector. Choose nearest to the central anchor, within2depth units and with
margin≥1. Require≥57/81unique assignments and≥20each in a checkerboard
train/holdout split. All competitors, rejections and selected candidates saved.

Method1: local plane height c+a*u+b*v, u/v=(pixel-63.5)/63.5. Fit only training
sites by three deterministic weighted least-squares rounds, prominence sqrt
weights and residual weight min(1,1/abs(residual)). No random sampling.
Require: abs(c)≤2; every rendered height≤2; heldout median residual≤1 and
95th≤2; competitor conflict≤5%(fitted point must favor assigned ridge by
≥0.5); maximum neighboring assigned-ridge step≤2; max slope≤0.1; center stays
within1depth of its chosen anchor and favors that anchor over any rival by≥0.5.

Method2 is a bounded plane-then-quadratic policy. Only attempt quadratic
c+a*u+b*v+d*u²+e*u*v+f*v² when plane heldout median residual>.15 OR95th>.35,
AND its
center/competitor/neighbor identity guards pass. Never use curvature to rescue
ambiguous sheet assignments. These adequacy thresholds are stricter than
plane admissibility(1/2), fixed using synthetic curved-sheet checks before
target evidence. Same fit/checks plus max(|2d|,|e|,|2f|)≤1 and≥20%heldout-median
residual improvement over plane. If quadratic qualifies, policy2uses it;
otherwise retain an accepted plane. If neither passes, retain the original T11-policy
patch, flag the fit rejected, and assign zero applied displacement. Do not
discard such points from aggregate evaluation. Both policies are evaluated;
first passing policy wins, never the best model probability or visual output.

Fitted height is a graph over the existing local frame; sample depth along
the graph's normalized local normal, retaining its signed orientation. XY
parameter-grid spacing and depth spacing remain9.6µm. Surface slopes produce
a bounded graph-metric deviation from9.6µm along-surface arclength; record it,
do not claim exact isometric flattening. Interpolation remains trilinear
float32, half-even uint8. Maximum surface/center translation≤19.2µm.

## Control and geometry gate

Apply the identical ridge detector, local fits and quality metrics to the
three preexisting T10control patches at local origins(128,128),(192,192),
(256,256). Do not refit their thresholds or substitute control geometry on
the target. Report raw depth-position proxy, depth spread, empty fraction,
crossings, center residual, holdout surface residual, and ambiguous fits.
These are reference diagnostics, not forced target distributions. No alternate
control-calibrated gate will replace the explicit10%criterion after evidence.

Use unchanged T11 CT-quality function: Gaussian-smoothed depth argmax spread
on heldout Ywings0..39,88..127; empty fraction; center-plane Xpeak crossings
(prominence20 at Y16,64,111); normalized center distance; longitudinal/depth
gradient fraction and central occupancy. Fitting uses a checkerboard holdout;
quality pixels overlap the local neighborhood, so these are internal geometry
diagnostics, not independent generalization evidence.

Each policy's gate requires all:≥70%reliable sheet identity;≥70%accepted fits;
median depth-peak spread over ALL72points improves≥10%relative to T9; median
band crossings do not increase; median empty fraction does not increase; all
accepted fits satisfy signed bounded-offset/competitor/continuity guards.
Report per-sample all-quality-pass fraction, accepted/rejected stratification,
median/max applied offset, and all medians for T9/T11/T12. No masking of failed
fits or sample selection from output. These conservative anti-jump checks are
not an oracle proof of same-sheet identity.

If reliable sheet identity<70%: D. Otherwise if no policy passes: C.
D is operational: these bounded local detector/fit guards could not establish
reliable identity. It is not proof that CT intrinsically lacks sheet identity
or that non-local tracking is unavoidable.
No ink inference in either branch. Missing CT is infrastructure BLOCKED,
not C/D; inventory exact new support before reading, acquire missing-only from
the same source with separate receipts and≤20MBinitial cap, never pad/clip.

## Conditional model stages

Only after geometry PASS, run T9/T11/selected-T12on all72points using exact
seed42step-075000, official preprocessing,17×128×128,batch8,float16 autocast,
center2×2mean. Same checkpoint SHA256e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab.
Require finite[0,1]maps and nonconstant centers. Positive smoke evidence also
requires mean adjacent-Z correlation over the six fixed8station row pairs
improves≥.05over matched T9,≥4/6pairs improve, and median within-patch nearest-
neighbor probability correlation improves≥.03. No shifting/stretching rows.
If geometry passes but smoke is weak: B, no full run. If positive, run all
valid frozen centers Z128..177, preserve bridges/masks, with rejected fits
falling back to frozen T11policy(T9when its orientation fails). No sample/path
selection by model. Full A requires adjacent-Z mean improvement≥.05 and
local64column-window correlation improvement≥.03over matched T9 plus visible
persistent nonrandom structure; otherwise B. Conditional model/full runners
must receive implementation review before executing; thresholds cannot change.

## Reproducibility

Synthetic sheet/competition/curvature tests before evidence; independent
preflight review, then hash-freeze and commit plan/source. Deterministic
OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1 enforced at entry.
Replay all72local fits, three control fits and every constructed patch from
saved immutable inputs into fresh directories; require exact numeric records
and raw hashes. Saved input patches must match their frozen coordinates and
T11 manifests. If model runs, replay all smoke maps(maxerror≤.002); if full
runs, verify all coordinate/mask rows and replay one fixed center per row.
Preserve lossless arrays locally, checksum, report exact stopped stage and
negative results honestly, commit/push sealed T12to origin/main without force.
