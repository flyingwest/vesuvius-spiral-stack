"""Report the completed fixed-budget scalar and combination experiment."""
import json,csv
from p0_t7_family import ROOT,OUT,FEATURES,WIDTHS,OFF,dump

def main():
 decision=json.loads((OUT/'decision.json').read_text());metrics=json.loads((OUT/'feature_metrics.json').read_text());combos=json.loads((OUT/'combination_metrics.json').read_text());plan=json.loads((OUT/'combination_plan.json').read_text());verification=json.loads((OUT/'verification.json').read_text());assert verification['status']=='PASS'
 by={r['feature']:r for r in metrics};base=by['center_intensity'];mean=by['mean_wide'];best=by[decision['BEST_INDIVIDUAL_FEATURE']];bc=decision['best_tested_combination_metrics'];comparison=decision['combination_vs_individual']
 scalar_table='\n'.join(f"| {r['feature']} | {r['mean_r']:.5f} | {r['median_r']:.5f} | {100*r['positive_fraction']:.2f} | {r['local_mean_r']:.5f} | {100*r['local_peak_above_control95']/1001:.2f} | {r['regions_both_improved']} | {r['best_lag_mean_r']:.5f} |" for r in metrics)
 combo_table='\n'.join(f"| {r['feature']} | {r['mean_r']:.5f} | {r['median_r']:.5f} | {r['local_mean_r']:.5f} | {100*r['local_peak_above_control95']/1001:.2f} | {r['delta_vs_T6_mean']:+.5f} | {r['regions_both_improved']} |" for r in combos)
 regionrows=list(csv.DictReader((OUT/'selected_region_comparison.csv').open()));region_table='\n'.join(f"| {r['feature']} | {r['z_start']}–{r['z_end_pair_start']} | {float(r['delta_r_vs_profile_mean']):+.5f} | {float(r['delta_local_vs_profile_mean']):+.5f} | {r['both_improved']} |" for r in regionrows if r['feature'] in [best['feature'],bc['feature']])
 dump(OUT/'visual_assessment.json',{'inspected':['feature_comparison.png','coherence_comparison.png','per_z_comparison.png'],
  'observation':'Mean-medium, wide-mean and their mixture emphasize similar broad intensity patches with unchanged jagged coverage and gaps. Short local correlation patches persist but positive/negative regions alternate; no new long vertical/diagonal fine structure is established. Ranking and decisions used metrics, never letter-like appearance.',
  'no_text_or_ink_claim':True})
 text=f'''# P0-T7 — bounded deterministic profile-feature search

## Decision and direct answers

P0_T7_COMPLETE = YES

P0_T7_RESULT = **C — SIGNAL_INSUFFICIENT**

BEST_INDIVIDUAL_FEATURE = **mean_medium**, under the frozen joint-ranking policy

BEST_COMBINATION = **NONE accepted as an improvement**

READY_FOR_LEARNED_MODEL = **NO**

The strongest **tested** combination is `combo_1_3`, an equal-weight normalized
medium/wide mean mixture. Its image is retained even though it is not accepted as
an upgrade. `mean_medium` is a screening selection, **not a claim that it beats
T6**. T6 `mean_wide` remains the better reference for mean adjacent and pooled local
zero-shift correlation. No scalar or combination satisfies the prespecified
broad-improvement guide; the conclusion also follows from the actual metric trade-offs,
not just a numerical threshold.

1. **Does an individual clearly beat profile mean? No.** The joint-ranked candidate
   has adjacent mean r={best['mean_r']:.5f} versus {mean['mean_r']:.5f}; local mean
   r={best['local_mean_r']:.5f} versus {mean['local_mean_r']:.5f}. It has more matched
   control-exceeding local windows ({best['local_peak_above_control95']} versus
   {mean['local_peak_above_control95']} of1001), but this does not compensate for
   weaker zero-shift metrics across most Z regions.
2. **Does a combination clearly beat the selected individual? Not jointly.**
   `combo_1_3` improves its adjacent mean by{comparison['mean_adjacent_gain']:.5f},
   but local improvement is small/uncertain and control exceedances fall to
   {bc['local_peak_above_control95']}/1001. It still does not beat the T6 baseline.
3. **Broad or localized gains? Localized/trade-offs.** Both selected candidate and
   mixture improve adjacent and local means together in only1 of5 fixed Z regions.
4. **Persistent non-random structure?** Bounded correspondence already visible in
   T5/T6 remains; these features do not establish stronger broad fine-scale patterns.
   Images are still dominated by irregular texture, broad patches and identical gaps.
5. **Strong enough to scale? No new support for scaling.** More Z with the same
   bounded scalar summaries is not justified by this result.
6. **Does this justify a learned ink-sensitive model now? No.** A negative bounded
   deterministic search does not establish available ink information, suitable labels,
   resolution, or a valid independent ink benchmark. A separately authorized feasibility
   investigation could address those prerequisites; no model was trained or proposed
   as a proven solution here.

C means **no meaningful multi-criterion improvement beyond T6 in this predefined
family**. It does not mean all CT signals are useless, that geometry has failed, or
that the narrow T6 repeatability PASS was incorrect. B would imply credible small
joint gains beyond T6; here the gains in local-control counts trade off against
weaker adjacent/local correlations, so that stronger conclusion is not made.

## Immutable inputs and scope

Use the exact saved seven-point normal profiles from
`outputs/p0_t6/slice_<z>_profiles.npz`, Z128..177. Their lineage is the immutable T5
paths `outputs/p0_t5/slice_<z>/path_and_samples.npz`, with CT source
`outputs/orientation/ct_level4_zyx.npy`: PHercParis4 / Scroll1,
`54keV_7.91um_Scroll1A.zarr` level4, standardized uint8, ZYX, localZ origin128.
The original array shape is(128,493,506); the unchanged T5 raster is50×4482.

The grid remains63.28µm in s and126.56µm in Z, span6.20144mm between first and last
Z centers, common last sampled s283.55768mm. Every stack has130629 valid samples,
**58.29050%** valid and the identical T5 mask/provenance. Codes1/3 remain real
observations;0/2 remain masked NaN, including inferred bridges. Frozen s coordinates,
missing tails and uncertain observations are not resized, shifted or filled.

T7 extraction does **not** query new CT or estimate a new tangent. It summarizes
T6 profiles at offsets−4,−2,−1,0,+1,+2,+4px along T6's saved normal. Verification
alone recomputes bilinear CT at those saved coordinates to confirm lineage. It
never runs thread search or changes a path. T5's five topology-defective rows remain
flagged upstream; they are neither repaired nor excluded to improve results.

## Predefined scalar family:28 stacks, no open-ended generation

`plan.json` freezes formulas, widths, ranking, controls and combination policy before
new feature evaluation. `scalars/FEATURES_FINALIZED.json` freezes the generated
arrays before scoring. The finite code loop implements only the explicitly listed
formula/width Cartesian product; it is not automated expression discovery.

One shared center baseline plus nine families at three widths gives **28** named
scalar stacks, below the30-stack cap:

| Width | Selected saved offsets px | Radius µm | Points |
| --- | --- | --- | --- |
| narrow | −1,0,+1 | 126.56 | 3 |
| medium | −2,−1,0,+1,+2 | 253.12 | 5 |
| wide | −4,−2,−1,0,+1,+2,+4 | 506.24 | 7 |

Let P0 be the center, P±R the selected endpoints and P the selected point set:

| Family | Fixed definition |
| --- | --- |
| center_intensity | P0, shared across widths |
| mean | Arithmetic mean(P), equally weighted samples |
| median | Median(P) |
| std | Population standard deviation(P) |
| range | max(P)−min(P) |
| center_minus_mean | P0−mean(P) |
| asymmetry | abs(mean(positive offsets)−mean(negative offsets)) |
| gradient_abs | abs(P+R−P−R)/(2R), CT/pixel |
| curvature | (2P0−P−R−P+R)/R², CT/pixel² |
| ridge_proxy | P0−mean(selected samples with abs(offset)≥R/2) |

Gradient/curvature use each width's endpoints. They therefore differ from T6's
central ±1 derivatives except at narrow width. `mean_wide` is **bit-identical** to
T6 profile mean. `center_intensity` is **bit-identical** to the T5 raw stack.
No signed-normal ambiguity occurs: all chosen definitions are invariant to normal
reversal. Some narrow responses are affine-equivalent (for example ridge/curvature
variants);28 named stacks are not28 independent hypotheses. All are retained rather
than silently increasing or shrinking the search after scores were seen.

## Identical evaluation for every stack

Reuse T6's numerical evaluator, adapting only plot height for the larger candidate
list. All feature/width candidates and combinations use:

- All49 adjacent-Z zero-shift Pearson correlations on the same jointly valid samples.
- Mean, median, range, positive fraction and pairwise improvements over raw center
  and T6 wide mean; every pair is saved, not just averages.
- Diagnostic global lag−64..+64 samples (±4.04992mm), minimum128 valid samples.
  No chosen lag is written into an authoritative stack.
- Local256-sample windows, stride128, minimum64 valid;1001 windows per stack.
- The same local small-lag search, and19 deterministic circular-shift controls per
  window subjected to the same peak search. Report the **fraction/count of windows
  exceeding matched peak controls**, not just optimistic peak correlations.
- 199 whole-row circular-shift controls per pair, masks moved with values.
- Fixed2000 circular moving-block bootstrap resamples, length5 adjacent pairs,
  RNG2606, saved indices; paired intervals versus each baseline.

Controls and confidence intervals are descriptive. Overlapping windows, reused
pairs, related features, same-data ranking and combination selection mean there is
no held-out or multiplicity-adjusted confirmation. The local-control fraction is not
an ink detection rate. Scalar-phase copies of the inherited single-metric
`rankings.json` are **diagnostics only**; the root joint selection described below
is authoritative.

### Joint ranking and regional robustness

Five fixed groups use the starting Z of each adjacent pair:
128–137,138–147,148–157,158–167,168–176 (10,10,10,10,9 pairs).
A region improves only when both its mean adjacent r and its mean within-pair local
r exceed T6 wide mean. This demands more than one attractive row or local maximum.

Rank candidates using the equally weighted average of four rank percentiles:
mean adjacent r; pooled mean local zero-shift r; local peak-control exceedance
fraction; and number of regions improving **both** metrics. Ties break by feature
name. Center is a reference, not a candidate; T6 wide mean remains eligible.
An ordinal screen may select a trade-off rather than dominate the baseline. Hence
selection is explicitly separate from acceptance.

The frozen broad-gain guide requires, versus T6 mean: adjacent mean gain≥0.02,
local pooled gain≥0.01,≥30/49 improved adjacent pairs,≥4/5 jointly improved regions,
and no decrease in matched local-control fraction. No candidate meets it.
The actual candidates also fail joint superiority without that guide.

## Scalar results (all28 preserved)

Local-control % is the proportion of1001 windows whose optimized peak exceeds its
matched19-shift peak-control95th percentile. Region count compares both zero-shift
metrics to T6 mean. Full precision, gains to both baselines, intervals and all49
pair measurements are in the JSON/CSV evidence.

| Scalar | Mean adjacent r | Median r | Positive % | Mean local zero r | Local-control % | Regions /5 | Best-lag mean r |
| --- | --- | --- | --- | --- | --- | --- | --- |
{scalar_table}

The fixed joint leader is **{best['feature']}**, not a new best zero-shift average.
Its mean gain relative to T6 is{best['delta_vs_T6_mean']:+.5f}, with paired95%
interval[{best['T6_mean_gain_ci_low']:.5f},{best['T6_mean_gain_ci_high']:.5f}].
Only{best['pairs_better_than_T6_mean']}/49 adjacent pairs improve. Local-control
exceedances rise104→129, but mean local r drops0.02924→0.02784. This is a trade-off,
not a broad repeatability upgrade. The selected feature beats raw center in mean
r, as T6 already did; that is not the novel T7 criterion.

## Exactly three fixed deterministic combinations

After scalar evaluation, the frozen rule selects the top three non-affine-equivalent
signals: **mean_medium, median_wide, mean_wide**. Affine-duplicate screening was
specified as ensuring distinct signals; no candidates needed skipping among these
three. The scalar ranking and selected list are preserved in `combination_plan.json`.

All three are divided by the fixed physical intensity range255, then combined with
weights0.5/0.5. No learned normalization, regression, optimization of weights or
iterative combination proposals occurs.

| Combination | Formula |
| --- | --- |
| combo_1_2 | 0.5×(mean_medium/255 + median_wide/255) |
| combo_1_3 | 0.5×(mean_medium/255 + mean_wide/255) |
| combo_2_3 | 0.5×(median_wide/255 + mean_wide/255) |

Two baselines are reevaluated in the combination phase for consistency; they do
not count as extra tested combinations. All three combined arrays are hash-locked
before their evaluation. The same four-criterion rank selects the strongest tested
combination, but no combination is accepted as an upgrade beyond T6.

| Combination | Mean adjacent r | Median r | Mean local zero r | Local-control % | Mean gain vs T6 | Regions /5 |
| --- | --- | --- | --- | --- | --- | --- |
{combo_table}

The strongest joint-ranked mixture is **{bc['feature']}**. Relative to the selected
medium mean, its adjacent gain is{comparison['mean_adjacent_gain']:.5f}, paired95%
interval{comparison['paired95']}. Its mean within-pair local gain is
{comparison['mean_pair_local_gain']:.5f}, interval{comparison['paired_local95']},
which includes zero. Local-control counts fall129→121. Relative to T6 mean, the
mixture's adjacent gain is{bc['delta_vs_T6_mean']:+.5f}, interval
[{bc['T6_mean_gain_ci_low']:.5f},{bc['T6_mean_gain_ci_high']:.5f}], and its pooled
local mean remains slightly lower. Thus **BEST_COMBINATION=NONE accepted**;
`best_combined_feature_stack.png` deliberately shows the best tested negative
result rather than hiding it.

## Where gains occur

Differences below are against T6 mean, with regions indexed by pair-start Z:

| Candidate | Pair-start Z region | Adjacent mean gain | Mean local-pair gain | Both improve |
| --- | --- | --- | --- | --- |
{region_table}

Only the first region improves both means for the selected individual and mixture.
The other four regions have lower adjacent means, with mixed local changes.
The mean of per-pair local means differs slightly from pooled local-window mean
because path widths give different numbers of windows. Both statistics are retained.
`per_z_comparison.png` displays every pair so these discontinuities are not concealed
by one aggregate score.

## Visual evaluation

Inspected `feature_comparison.png`, `coherence_comparison.png`, and
`per_z_comparison.png`. Medium/wide means and their mixture emphasize similar broad
patches; jagged support and missing tails remain identical. Short local positive
correlation patches occur, interspersed with negative regions and discontinuities.
There is no visually established new persistent fine-scale structure across most
Z levels. Selection was numerical, not based on letter-like appearance. No patch
is interpreted as ink or text.

Required root stack PNGs use each signal's **single global1st–99th percentile**
display from its phase `display.json`; they are native4482×50 images. This is
visualization only, never row-wise contrast or numeric rescaling. Missing is black
in native images and magenta in the enlarged comparison. Lossless arrays retain
original CT/derived units, and per-feature fixed theoretical-range images are also
preserved under each phase. No transformed/lag-aligned raster is produced.

## Replay and verification

`verification.json` confirms:

- Every saved T6 profile agrees with independently interpolated original CT at
  frozen T5-derived XY/normal positions, tolerance1e−10; no reconstruction rerun.
- All28 scalar formulas replay from saved profiles with independent arithmetic,
  including sorted-value medians and explicit derivative formulas.
- Both authoritative baseline stacks are bit-identical; shared masks/provenance
  and physical coordinate arrays are byte-identical T5/T6 copies.
- All3 fixed-weight combinations replay exactly.
- **{verification['adjacent_correlations_replayed']} adjacent correlations** replay
  exactly (28 scalar plus5 combination-phase stacks, including2 repeated baselines).
- **{verification['local_windows_replayed']} local windows** have exact zero-shift
  replay and scalar checks of the selected diagnostic peak, maximum discrepancy
  {verification['local_selected_peak_max_error']:.3g}. This does not claim an independent
  reimplementation of every local null-control maximization.
- Phase array seals, earlier evidence preservation and PNG decoding pass.

No finalized path, source row, profile, T5 stack, T6 evidence, README or AGENTS was
modified. Negative feature/combination results and all scores remain available.
`SHA256SUMS` covers all new evidence, scripts and this report; every entry is verified
immediately after generation. Large numeric arrays/logs remain local if excluded by
Git, but are included in the full checksum seal and Git artifact inventory.

## Reproduction and artifact inventory

Commands from repository root; extraction, combination creation and evaluation
refuse overwriting their run outputs:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t7 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t7_family.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t7 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t7_run_metrics.py scalars
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t7 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t7_select.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t7 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t7_run_metrics.py combinations
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t7 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t7_summary.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t7 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t7_verify.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t7 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t7_report.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t7 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t7_seal.py
sha256sum -c outputs/p0_t7/SHA256SUMS
```

Required root outputs: `baseline_center_stack.png`, `baseline_profile_mean_stack.png`,
`best_individual_feature_stack.png`, `best_combined_feature_stack.png`,
`feature_comparison.png`, `coherence_comparison.png`, `per_z_comparison.png`,
`feature_metrics.csv`, `feature_metrics.json`, `combination_metrics.csv`,
`verification.json`. Root `mask.npy`, `provenance.npy`, coverage and coordinate arrays
apply unchanged to all numeric stacks in `scalars/` and `combinations/`.

Each phase retains every feature's NPY, fixed/global-contrast PNGs, coherence plots,
adjacent/window/lag/circular-control CSVs, detailed metrics, bootstrap indices and
array seal. `scalar_region_metrics.csv`, `combination_region_metrics.csv`,
`selected_region_comparison.csv`, `combination_plan.json`, `decision.json`,
`visual_assessment.json`, plan and preservation manifests document decisions and limits.
The Git milestone excludes numeric arrays under repository rules; a clone alone is
not a complete input/replay package. Originals remain preserved locally.

## Limitations and stopping point

One coarse scan/region, no ink ground truth, related deterministic summaries,
within-dataset selection, irregular independent path ordering and uncertain layer
identity limit inference. ±4px profiles can span neighboring material; smooth broad
CT repeatability is not ink sensitivity. These constraints apply equally to baselines
and new features and are not evidence that geometry failed.

**C: no meaningful joint improvement beyond T6 in this bounded search.** Retain T6
profile mean as the reference; do not scale or train an ink model solely on these
scores. A separately authorized target-specific data/label and resolution feasibility
study could determine whether a learned model has an evaluable task. No such study,
OCR, mesh, language prior, learned model, cross-Z reconstruction or new feature search
was performed. **Stop after P0-T7.**
'''
 (ROOT/'docs/p0_t7_deterministic_feature_search.md').write_text(text)
 print('Report written: C SIGNAL_INSUFFICIENT beyond T6; no accepted combination; learned readiness NO')
if __name__=='__main__':main()
