"""Generate the signal-only T6 report from preserved measurements."""
import json,csv
import numpy as np
from p0_t6_features import ROOT,OUT,FEATURES,DS,SPACING,dump

def main():
 rows=json.loads((OUT/'feature_metrics.json').read_text());by={r['feature']:r for r in rows};b=by['center_intensity'];best=by['profile_mean'];v=json.loads((OUT/'verification.json').read_text());assert v['status']=='PASS'
 sampling=list(csv.DictReader((OUT/'sampling_metrics.csv').open()));display=json.loads((OUT/'display.json').read_text())
 definitions='\n'.join(f'| `{name}` | {definition} | {limits[0]}..{limits[1]} |' for name,(definition,limits) in FEATURES.items())
 main_table='\n'.join(f"| {r['feature']} | {r['mean_r']:.5f} | {r['median_r']:.5f} | {100*r['positive_fraction']:.2f} | {r['delta_mean_vs_center']:+.5f} | {r['pairs_better_than_center']}/49 | {r['best_lag_mean_r']:.5f} |" for r in rows)
 local_table='\n'.join(f"| {r['feature']} | {r['local_mean_r']:.5f} | {r['local_median_r']:.5f} | {100*r['local_positive_fraction']:.2f} | {r['local_best_mean_r']:.5f} | {r['local_peak_above_control95']}/1001 |" for r in rows)
 contrast='\n'.join(f"| {name} | {d['contrast_range'][0]:.8f} | {d['contrast_range'][1]:.8f} |" for name,d in display.items())
 dump(OUT/'assessment.json',{'P0_T6_COMPLETE':'YES','P0_T6_RESULT':'PASS','BEST_FEATURE':'profile_mean','READY_FOR_NEXT_STEP':'YES',
   'meaning':'PASS only for reproducible adjacent-Z CT feature repeatability improvement; not ink sensitivity or broad reconstruction',
   'decision':'A: bounded deterministic feature research; D: insufficient evidence for ink probability or choosing learned models; geometry not changed or disproved',
   'visual_inspection':['profile_examples.png','normal_sampling_diagnostic.png','feature_comparison.png','feature_stack_overview.png'],
   'visual_findings':'Profile mean emphasizes broad light/dark patches and smoother texture, not persistent manuscript-like fine structure. Other feature stacks remain largely irregular texture. Mask geometry is identical; no visible pattern is called text.'})
 text=f'''# P0-T6 — local CT features on frozen threads

## Result and explicit answers

**P0_T6_COMPLETE = YES**  
**P0_T6_RESULT = PASS**  
**BEST_FEATURE = profile_mean**  
**READY_FOR_NEXT_STEP = YES — bounded deterministic-feature investigation only.**

This is a **narrow CT repeatability PASS**, not ink detection, readable text, broad
2-D reconstruction coherence, or vindication/rejection of every geometric choice.
The user criterion allows clear reproducible improvement in adjacent-Z **and/or**
local-window coherence. Profile mean provides the adjacent-Z improvement; robust
fine-scale/local improvement has **not** been established.

1. **Does a feature materially outperform center intensity?** Yes for whole-row
   zero-shift repeatability: mean r **{b['mean_r']:.5f} → {best['mean_r']:.5f}**,
   gain **{best['delta_mean_vs_center']:.5f}**, with {best['pairs_better_than_center']}/49
   pairs improving. The fixed moving-block bootstrap gives a descriptive95%
   paired-gain interval [{best['paired_block_ci_low']:.5f},
   {best['paired_block_ci_high']:.5f}]. This is a modest absolute but consistent
   gain, not proof of an ink-specific observable.
2. **Best adjacent-Z feature:** `profile_mean`, also best median r
   ({best['median_r']:.5f} versus {b['median_r']:.5f}).
3. **Best local-window feature:** `profile_mean` by unshifted mean local r
   ({best['local_mean_r']:.5f} versus {b['local_mean_r']:.5f}), but the paired
   per-Z local-mean improvement interval includes zero. Optimized local peaks
   alone are not used to claim a robust improvement.
4. **Broad or localized?** Whole-row gains occur across a majority of pairs, mainly
   in broad intensity variation. Fine-scale gains remain weak/localized. Positive
   pair fraction actually declines from {100*b['positive_fraction']:.2f}% to
   {100*best['positive_fraction']:.2f}%, and the worst pair becomes more negative.
   No broad manuscript-like structure is demonstrated.
5. **Direction:** support **A**, a bounded deterministic-feature follow-up with
   independent validation and controls for ordinary averaging/background texture.
   **D** applies to ink sensitivity: no labels, known ink reference or detection
   accuracy is available. These results do not justify **B**, a learned
   ink-probability model, or establish that **C**, changing geometry, is necessary.
   Geometry remains a possible confound, not a demonstrated failure. No follow-up
   experiment was started.

## Objective and immutable inputs

Test whether a local CT neighborhood supplies a more repeatable scalar row than
one center sample on exactly the same frozen independent threads. All geometry
comes from `outputs/p0_t5/slice_<z>/path_and_samples.npz`, **Z128..177 inclusive**.
The authoritative baseline is `outputs/p0_t5/stack.npy` and its mask/provenance.
T5 remains unchanged, including its five flagged topology-defective rows; T6 neither
repairs nor removes them. No thread search is invoked.

CT source: `outputs/orientation/ct_level4_zyx.npy`, PHercParis4 / Scroll1,
`54keV_7.91um_Scroll1A.zarr`, level4, standardized uint8, saved shape(128,493,506),
ZYX with globalZ origin128. T6 accesses only the same-Z plane for each row.
Original acquisition references are in `docs/p0_orientation.md` and
`outputs/orientation/manifest.json`; no download was needed.
Voxel spacing is **126.56µm** in X/Y/Z. Arc sampling remains **63.28µm**.
The 50-row stack shape is **50×4482**, first-to-last Z center span6.20144mm,
common last s sample283.55768mm. s=0 remains each T5 seed; no alignment or resizing.

`plan.json` and `preservation.json` were written before new feature evidence;
`FEATURES_FINALIZED.json` locks all extracted numeric features before evaluation.
Feature definitions were not adjusted after inspecting scores or plots. The plan
specified indicative materiality guides (adjacent mean gain around0.03 or local
mean gain around0.05 plus consistency), not text appearance or a new ink-success gate.
This is an exploratory eight-feature comparison, not a preregistered clinical-style
trial or held-out generalization study.

## Tangent, normal and profile sampling

For each source sample with provenance1 or3 and mask=false:

- Locate its containing edge and the maximal contiguous edge run labeled observed
  or low-confidence (1/3), excluding inferred edges2.
- Interpolate the frozen polyline at s−2px and s+2px, clipping both to that observed
  run. Normalize their chord as tangent t. The normal is n=(−t_y,t_x).
- A chord norm below1e−8 would fall back to the containing nonzero edge direction;
  **{sum(int(r['fallback_tangents']) for r in sampling)} fallbacks** occurred.
  **{sum(int(r['clipped_tangent_windows']) for r in sampling)}** windows were clipped
  at observed-run boundaries. This estimates orientation without modifying the path
  or using geometry from another Z.
- Sample bilinearly at XY + u·n for
  **u = −4, −2, −1, 0, +1, +2, +4 pixels**. The physical offsets are
  −506.24, −253.12, −126.56, 0, +126.56, +253.12, +506.24µm.
  No wider profile was used; full normal span is1.01248mm.
- All queries were inside the CT image; bounds were asserted before interpolation.
  No extrapolation or cross-Z sampling was required. The seven samples are equally
  weighted points at **nonuniform distances**, not a spatial integral.

There are **130629** valid center locations and
**{sum(int(r['profile_queries']) for r in sampling)}** CT profile queries.
Profiles are never computed at inferred or unresolved centers. Off-center queries
may include air, neighboring papyrus layers or other local material; they are actual
CT neighborhood values, **not all certified measurements of one papyrus sheet**.
At the current coarse scale this mixture is a central limitation.

## Fixed feature set

Let P(u) be the sampled normal profile at offset u pixels. Reversing the tangent
flips normal sign but leaves all eight selected scalar definitions unchanged.
No feature uses expected letters or ink labels.

| Feature | Definition | Fixed display range |
| --- | --- | --- |
{definitions}

Feature values remain float64 without quantization. Gradient and curvature units
above use pixels; conversion to per-µm units would multiply by a fixed constant
and would not change Pearson correlations. These are finite normal derivatives,
not a general 2-D Laplacian or a learned surface detector.

## Masks, provenance and lossless representation

Every feature has the **same exact mask** as T5, and the same provenance:
0 unresolved/padding;1 observed;2 inferred bridge;3 low-confidence observed.
Codes0/2 are NaN in every feature. Code3 remains measured data. No bridge is filled,
no observation is dropped, and no extra s range is introduced. Valid fraction is
**{100*b['valid_fraction']:.5f}%** for every stack. Shared `mask.npy`,
`provenance.npy`, `coverage.npy`, `s_um.npy`, `z_indices.npy` and
`z_relative_um.npy` are byte-identical T5 copies and apply to all feature arrays.

Lossless `<feature>_stack.npy` files contain each numeric raster. Each
`slice_<z>_profiles.npz` stores valid sample indices, frozen XY, tangent/normal,
actual queried XY, seven CT values, clipped tangent s bounds and fallback flags.
The center-intensity numeric stack is **bit-identical to T5**, not merely correlated.

## Adjacent-Z evaluation and baseline comparison

All **49 consecutive pairs** are evaluable for each of8 features, using the
identical jointly valid s overlap. Pearson r is computed without shifts. First
sample differences exclude gaps. Diagnostic lags scan−64..+64 samples
(±4.04992mm) with at least128 joint values; **no shift is applied to any raster**.
All individual pair measurements and tested lags are saved.

| Feature | Mean r | Median r | Positive % | Mean gain vs center | Pairs improved | Diagnostic best-lag mean r |
| --- | --- | --- | --- | --- | --- | --- |
{main_table}

The best feature's range is [{best['min_r']:.5f}, {best['max_r']:.5f}], compared
with center [{b['min_r']:.5f}, {b['max_r']:.5f}]. Thus a higher mean does not mean
uniform improvement; worsening anticorrelations are retained. Standard deviation,
range, derivatives, ridge contrast and asymmetry do **not** outperform the baseline
on mean adjacent r. Their negative results are preserved rather than hidden.

Paired uncertainty uses2000 circular moving-block resamples of the49 adjacent-pair
sequence, block length5, fixed RNG seed2606, saved indices. The same indices are
used for every feature. Blocks acknowledge neighboring-pair dependence but do not
prove independent samples. Intervals are descriptive and **not corrected for
selection among eight related features**. No confirmatory significance claim is made.
The profile-mean interval is positive and improvement occurs in35/49 pairs,
supporting the limited repeatability result under the stated T6 criterion.

For additional context,199 deterministic circular shifts of each second row within
its shared native prefix (offsets128..N−128, masks moved with values) give descriptive
whole-row controls. Profile mean exceeds its own95th control percentile in
{best['pairs_above_global_shift95']}/49 pairs, versus {b['pairs_above_global_shift95']}/49
for center. This is limited additional evidence, not a formal independent p-value.

## Local/windowed coherence and diagnostic lags

Use256-sample windows (16.13640mm between first/last centers; nominal width16.19968mm),
stride128 samples (8.09984mm), at least64 joint values. Each feature has1001 matching
windows. Compute both zero-shift r and the diagnostic local maximum over±64 samples.
Repeat the same local peak search for19 circular-shift controls per window, retaining
all peak controls. Maxima are selection-biased and are not achieved raster alignment.

| Feature | Mean local zero r | Median local zero r | Positive local % | Mean local diagnostic peak r | Peaks above matched control95 |
| --- | --- | --- | --- | --- | --- |
{local_table}

`profile_mean` is the best feature by mean **unshifted** local correlation, but its
improvement is small: {b['local_mean_r']:.5f}→{best['local_mean_r']:.5f}.
Averaging windows per adjacent pair before the paired block comparison gives gain
{best['delta_local_pair_mean']:.5f}, interval
[{best['local_block_ci_low']:.5f}, {best['local_block_ci_high']:.5f}], which includes
zero. This equally weighted pair statistic differs from the pooled-window mean
because pairs have different numbers of windows.

Local optimized mean increases {b['local_best_mean_r']:.5f}→{best['local_best_mean_r']:.5f},
but control-exceeding windows decrease **{b['local_peak_above_control95']}→{best['local_peak_above_control95']}**.
Thus optimized local peaks must not be presented as robust extra fine-scale recovery.
First-difference adjacent mean also changes from {b['difference_mean_r']:.5f} to
{best['difference_mean_r']:.5f}, not an improvement. The evidence is consistent with
normal-neighborhood averaging exposing broader CT variation, rather than revealing
stable ink-like high-frequency detail. That explanation is an inference, not a
proven causal decomposition.

## Visual inspection

Inspected `profile_examples.png`, `normal_sampling_diagnostic.png`,
`feature_comparison.png` and `feature_stack_overview.png`. Examples were selected
before scores by valid-sample fractions0.1/0.5 in Z128,152,177, not by attractive
appearance. Profiles have peaks, troughs and asymmetry; the normal lines visibly
span multiple CT structures. These do not identify which structures contain ink.

Profile mean emphasizes broader light/dark patches and smoother-looking texture.
It still does not establish long coherent manuscript-like fine structure. Other
feature stacks remain predominantly irregular texture. The same structured support
and missing-tail pattern appears in every feature because masks were copied exactly.
No visible dark/light region is interpreted as text; readability was not a gate.

`<feature>_stack.png` uses the fixed theoretical display range in the definition
table, clipped and rounded to0–255. `<feature>_stack_contrast.png` uses one global
pooled valid1st–99th percentile range per feature. No row is enhanced independently
and no numeric array is modified by display. Missing pixels are black in native
PNGs; masks distinguish them from real low values. The enlarged overview uses
magenta for missing data and is not a true-physical-aspect rendering.

| Feature | Global display low | Global display high |
| --- | --- | --- |
{contrast}

## Replay, preservation and reproducibility

`verification.json` records:

- All50 same-Z source profiles replayed with independent four-corner bilinear CT
  interpolation; maximum error **{v['max_bilinear_error']:.3g}**, tolerance1e−10.
- Tangent bounds never cross an inferred edge; unit tangents/perpendicular normals,
  frozen coordinates and all profile locations reproduced.
- All eight feature formulas independently recomputed, within1e−10.
- Center stack bit-identical to T5; masks/provenance/coordinates byte-identical.
- All392 adjacent correlations recomputed exactly; numeric masks and image mappings
  checked, all PNGs decoded.
- Feature-finalization hashes and all earlier preservation hashes pass.

`metric_verification.json` adds paired-bootstrap replay and scalar zero/selected-lag
checks for all8008 local windows. It does not claim an independent reimplementation
of every null-control search. T5 geometry search was never rerun. README/AGENTS,
T5 paths and authoritative stacks, and all prior recorded evidence remain unchanged.
All new artifacts are covered by `outputs/p0_t6/SHA256SUMS`, verified immediately.

Exact commands from repository root (extraction/evaluation refuse overwriting an
existing run; verification is from saved geometry, not reconstruction):

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t6 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t6_features.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t6 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t6_evaluate.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t6 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t6_verify.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t6 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t6_audit_metrics.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t6 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t6_report.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t6 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t6_seal.py
sha256sum -c outputs/p0_t6/SHA256SUMS
```

## Outputs and limitations

Under `outputs/p0_t6/`, all8 named features have `_stack.png`,
`_stack_contrast.png`, `_coherence.png` and lossless `_stack.npy` files.
Required summary files are `feature_comparison.png`, `feature_metrics.csv`,
`feature_metrics.json`, `profile_examples.png`, `normal_sampling_diagnostic.png`
and `verification.json`. Shared mask/provenance arrays, per-slice profiles, per-feature
adjacent/window/lag/control CSVs, display mapping, bootstrap indices, rankings,
plan, preservation and assessment files complete the evidence package.
Lossless arrays are kept locally even when excluded by Git rules; a Git clone alone
is not a complete replay dataset. `git_artifact_inventory.json` records that split.

Limitations: one coarse scan/region, no independent ink truth, no held-out region,
related features and autocorrelated samples, inferred geometry and imperfect T5
layer choices, and normals that may cross nearby material. Averaging can improve
repeatability without improving sensitivity to ink. Five T5 topology-defective rows
remain present and identified upstream; neither their removal nor any geometry repair
was used to improve these metrics. Failure of several simple derivatives is not
proof that geometry failed or that a learned ink model would succeed.

**Conclusion:** the local profile mean clearly improves the narrow adjacent-Z CT
repeatability metric on frozen paths. This satisfies P0-T6's OR criterion, while
local/fine-scale and ink-specific evidence remain insufficient. Continue only with
a separately authorized bounded deterministic-feature test, ideally distinguishing
ordinary averaging/background effects from a target-specific signal. No learned
model, OCR, language prior, mesh, cross-Z geometry or path correction was used.
**Stop after P0-T6.**
'''
 (ROOT/'docs/p0_t6_local_ct_features.md').write_text(text)
 print('Report written: PASS (CT repeatability only); best profile_mean; next direction A')
if __name__=='__main__':main()
