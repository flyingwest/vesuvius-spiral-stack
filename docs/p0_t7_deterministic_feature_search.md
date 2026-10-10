# P0-T7 — bounded deterministic profile-feature search

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
   has adjacent mean r=0.09396 versus 0.10809; local mean
   r=0.02784 versus 0.02924. It has more matched
   control-exceeding local windows (129 versus
   104 of1001), but this does not compensate for
   weaker zero-shift metrics across most Z regions.
2. **Does a combination clearly beat the selected individual? Not jointly.**
   `combo_1_3` improves its adjacent mean by0.00914,
   but local improvement is small/uncertain and control exceedances fall to
   121/1001. It still does not beat the T6 baseline.
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
| center_intensity | 0.05152 | 0.05924 | 77.55 | 0.01251 | 10.79 | 0 | 0.12852 |
| mean_narrow | 0.08193 | 0.09860 | 69.39 | 0.02544 | 10.59 | 0 | 0.17330 |
| median_narrow | 0.06956 | 0.07743 | 71.43 | 0.02101 | 10.39 | 0 | 0.15320 |
| std_narrow | 0.00485 | 0.01075 | 55.10 | -0.00150 | 9.19 | 0 | 0.09635 |
| range_narrow | 0.00443 | 0.00984 | 55.10 | -0.00168 | 8.99 | 0 | 0.09614 |
| center_minus_mean_narrow | 0.00084 | -0.00361 | 46.94 | -0.00378 | 10.99 | 0 | 0.09283 |
| asymmetry_narrow | 0.00052 | -0.00009 | 48.98 | 0.00098 | 8.29 | 0 | 0.07529 |
| gradient_abs_narrow | 0.00052 | -0.00009 | 48.98 | 0.00098 | 8.29 | 0 | 0.07529 |
| curvature_narrow | 0.00084 | -0.00361 | 46.94 | -0.00378 | 10.99 | 0 | 0.09283 |
| ridge_proxy_narrow | 0.00084 | -0.00361 | 46.94 | -0.00378 | 10.99 | 0 | 0.09283 |
| mean_medium | 0.09396 | 0.11780 | 73.47 | 0.02784 | 12.89 | 1 | 0.19250 |
| median_medium | 0.08070 | 0.09177 | 67.35 | 0.01847 | 10.89 | 0 | 0.17258 |
| std_medium | 0.00045 | 0.01774 | 55.10 | -0.00492 | 9.89 | 0 | 0.09797 |
| range_medium | 0.00142 | 0.01313 | 57.14 | -0.00305 | 9.99 | 0 | 0.09623 |
| center_minus_mean_medium | 0.00194 | 0.00703 | 53.06 | -0.00601 | 11.79 | 0 | 0.09572 |
| asymmetry_medium | 0.00081 | -0.00023 | 48.98 | 0.00281 | 9.19 | 0 | 0.07967 |
| gradient_abs_medium | -0.00063 | -0.01214 | 44.90 | -0.00234 | 9.39 | 0 | 0.08928 |
| curvature_medium | 0.00079 | -0.00072 | 48.98 | -0.00634 | 10.99 | 0 | 0.09536 |
| ridge_proxy_medium | 0.00194 | 0.00703 | 53.06 | -0.00601 | 11.79 | 0 | 0.09572 |
| mean_wide | 0.10809 | 0.13831 | 73.47 | 0.02924 | 10.39 | 0 | 0.20604 |
| median_wide | 0.09717 | 0.11077 | 73.47 | 0.01947 | 10.99 | 0 | 0.18610 |
| std_wide | -0.00123 | -0.00601 | 46.94 | -0.01042 | 10.89 | 0 | 0.10074 |
| range_wide | 0.00021 | 0.00203 | 51.02 | -0.00741 | 11.09 | 0 | 0.09348 |
| center_minus_mean_wide | -0.00017 | 0.00332 | 51.02 | -0.00698 | 9.29 | 0 | 0.09220 |
| asymmetry_wide | 0.00729 | 0.01267 | 59.18 | -0.00308 | 8.79 | 0 | 0.09652 |
| gradient_abs_wide | 0.02155 | 0.01447 | 65.31 | 0.01233 | 8.19 | 0 | 0.09943 |
| curvature_wide | 0.00458 | 0.00055 | 53.06 | 0.00138 | 8.99 | 0 | 0.08986 |
| ridge_proxy_wide | 0.00093 | -0.00396 | 46.94 | -0.00714 | 9.19 | 0 | 0.09724 |

The fixed joint leader is **mean_medium**, not a new best zero-shift average.
Its mean gain relative to T6 is-0.01414, with paired95%
interval[-0.02657,-0.00186].
Only17/49 adjacent pairs improve. Local-control
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
| combo_1_2 | 0.10160 | 0.12932 | 0.02580 | 11.39 | -0.00649 | 1 |
| combo_1_3 | 0.10310 | 0.13135 | 0.02916 | 12.09 | -0.00499 | 1 |
| combo_2_3 | 0.10618 | 0.12614 | 0.02548 | 10.69 | -0.00191 | 0 |

The strongest joint-ranked mixture is **combo_1_3**. Relative to the selected
medium mean, its adjacent gain is0.00914, paired95%
interval[0.002714854905942418, 0.015474502145210407]. Its mean within-pair local gain is
0.00187, interval[-0.003516202159640981, 0.007869771991162894],
which includes zero. Local-control counts fall129→121. Relative to T6 mean, the
mixture's adjacent gain is-0.00499, interval
[-0.01136,0.00136], and its pooled
local mean remains slightly lower. Thus **BEST_COMBINATION=NONE accepted**;
`best_combined_feature_stack.png` deliberately shows the best tested negative
result rather than hiding it.

## Where gains occur

Differences below are against T6 mean, with regions indexed by pair-start Z:

| Candidate | Pair-start Z region | Adjacent mean gain | Mean local-pair gain | Both improve |
| --- | --- | --- | --- | --- |
| mean_medium | 128–137 | +0.00472 | +0.01422 | True |
| mean_medium | 138–147 | -0.04019 | -0.02260 | False |
| mean_medium | 148–157 | -0.01460 | -0.01163 | False |
| mean_medium | 158–167 | -0.01111 | +0.00618 | False |
| mean_medium | 168–176 | -0.00898 | +0.00289 | False |
| combo_1_3 | 128–137 | +0.00652 | +0.00985 | True |
| combo_1_3 | 138–147 | -0.01745 | -0.01002 | False |
| combo_1_3 | 148–157 | -0.00767 | -0.00647 | False |
| combo_1_3 | 158–167 | -0.00202 | +0.00344 | False |
| combo_1_3 | 168–176 | -0.00428 | +0.00128 | False |

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
- **1617 adjacent correlations** replay
  exactly (28 scalar plus5 combination-phase stacks, including2 repeated baselines).
- **33033 local windows** have exact zero-shift
  replay and scalar checks of the selected diagnostic peak, maximum discrepancy
  3.21e-14. This does not claim an independent
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
