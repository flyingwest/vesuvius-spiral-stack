"""Summarize frozen T5 artifacts. Run only after evaluation and visual inspection."""
import json
from pathlib import Path
import numpy as np
from p0_t5_consecutive import ROOT,OUT,ZS,dump

def make_report():
    m=json.loads((OUT/'stack_metrics.json').read_text());pairs=json.loads((OUT/'adjacent_correlations.json').read_text())
    rows=[json.loads((OUT/f'slice_{z}/metrics.json').read_text()) for z in ZS]
    display=json.loads((OUT/'display.json').read_text());ver=json.loads((OUT/'verification.json').read_text());geom=json.loads((OUT/'geometry_replay_validation.json').read_text())
    assessment=json.loads((OUT/'visual_assessment.json').read_text())
    lag=json.loads((OUT/'lag_metrics.json').read_text())
    assert ver['status']==geom['status']=='PASS'
    a=np.load(OUT/'stack.npy');mask=np.load(OUT/'mask.npy');prov=np.load(OUT/'provenance.npy');coverage=np.load(OUT/'coverage.npy')
    path_table='\n'.join(f"| {r['z']} | {r['total_length_mm']:.3f} | {r['observed_length_mm']:.3f} | {r['bridge_length_mm']:.3f} | {100*r['observed_utilization_fraction']:.2f} | {r['bridges']} | {100*r['low_confidence_fraction_samples']:.2f} | {r['corrections']} | {r['final_defect_count']} | {'YES' if r['exterior_reached'] else 'NO'} | {r['last_sample_s_um']/1000:.3f} |" for r in rows)
    pair_table='\n'.join(f"| {r['z_a']}–{r['z_b']} | {r['joint_valid']} | {r['r']:.5f} | {r['difference_r']:.5f} | {r['best_lag_samples']} |" for r in pairs)
    total_before={k:sum(json.loads((OUT/f'slice_{z}/topology_before.json').read_text())[k] for z in ZS) for k in ['proper_crossings','nonadjacent_point_contacts','collinear_overlaps']}
    total_after={k:sum(json.loads((OUT/f'slice_{z}/validation.json').read_text())['topology'][k] for z in ZS) for k in total_before}
    residual=[{'z':r['z'],'proper':r['proper_crossings'],'contacts':r['point_contacts'],'overlap':r['collinear_overlaps'],'reason':r['cleanup_termination']} for r in rows if r['final_defect_count']]
    common_valid_count=np.sum(~mask,axis=0)
    extra={'total_topology_before':total_before,'total_topology_after':total_after,'residual_defects':residual,
           'total_corrections':sum(r['corrections'] for r in rows),'zero_edges_removed':sum(r['zero_length_edges_removed'] for r in rows),
           'columns_valid_in_at_least25rows':int((common_valid_count>=25).sum()),'columns_valid_in_at_least40rows':int((common_valid_count>=40).sum()),
           'median_abs_adjacent_length_change_mm':float(np.median([abs(r['path_length_change_mm']) for r in pairs])),
           'max_abs_adjacent_length_change_mm':float(max(abs(r['path_length_change_mm']) for r in pairs)),
           'median_abs_diagnostic_lag_um':float(np.median([abs(r['best_lag_um']) for r in pairs])),
           'observed_utilization_minmax':[min(r['observed_utilization_fraction'] for r in rows),max(r['observed_utilization_fraction'] for r in rows)]}
    dump(OUT/'report_measurements.json',extra)
    text=f'''# P0-T5 — consecutive-slice Direct Thread-Stack

## Result and scope

**P0_T5_COMPLETE = YES. P0_T5_RESULT = {assessment['result']}.**

{assessment['conclusion']}

All 50 independent reconstructions produced numeric rows. {m['exterior_rows']}/50
({m['exterior_rows']*2}%) reached the inherited exterior target; {m['clean_exterior_rows']}/50
are both exterior-reaching and topology-clean. The stack shape is **50 × {m['shape'][1]}**.
Residual topology defects are disclosed below, not hidden or declared corrected.

P0-T4's five widely separated slices did not adequately test consecutive-slice
coherence. This experiment uses **every Z from 128 through 177 inclusive**.
No physical-fragment persistence requirement is used to judge it. The evaluated
hypothesis is independent 2-D maximum-recovery reconstruction → physical-arc row →
consecutive-Z stack, with visible/inferred/low-confidence provenance.

## Exact input and physical coordinates

Source: `outputs/orientation/ct_level4_zyx.npy`, uint8, shape (128,493,506), ZYX.
The dataset is PHercParis4 / Scroll 1, standardized
`54keV_7.91um_Scroll1A.zarr`, pyramid level 4. Original acquisition and source URLs
are recorded in `docs/p0_orientation.md` and `outputs/orientation/manifest.json`.
T5 uses local planes [0:50,:,:], global level-4 Z=128..177; XY is neither transposed
nor reversed. No data were downloaded. The source SHA-256 is in `plan.json`.

Voxel spacing is isotropic **126.56 µm**. Rows use **Δs=63.28 µm**, exactly the T3/T4
sampling grid. Consecutive Z center spacing is **126.56 µm**, giving
**{m['z_center_span_mm']:.5f} mm** between first and last centers (50 voxel cells
span {50*.12656:.5f} mm). Use `z_indices.npy`, `z_relative_um.npy` and `s_um.npy` to
interpret numeric axes. These coordinates are physical distances, not calibrated
material registration.

## Fixed independent reconstruction method

`scripts/p0_t5_consecutive.py` reuses the existing T1 ridge graph, corrected T2
papyrus-domain builder/search, and T3 run04 inclusive-contact exchange. No search
or cleanup parameters were changed in response to neighboring-row correlations.
The existing five slices are prior observed design input; the other 45 extend the
same fixed method. This is exploratory scale testing, not a new preregistered trial.

1. Difference-of-Gaussians response: σ=0.6 minus σ=2 pixels, threshold >3;
   skeletonize within X[25,475), Y[25,325), and extract maximal graph chains.
2. Use the unchanged deterministic T1 seed rule: nearest point to (246,112)
   on a chain of arc ≥20 px. T2 explores both orientations of that seed chain.
   No seed or branch is selected from another Z level.
3. Use the fixed polygon in `p0_t2_papyrus_domain.py` to exclude holder geometry.
   The inherited exterior target is X=430..460, Y=150..210: a peripheral papyrus
   region, not a newly certified outermost boundary.
4. Search oriented chains with virtual endpoint bridges ≤12 px; tangent cosine
   ≥−0.5 and endpoint/chord cosines ≥−0.35. The unchanged penalty is
   2×distance + 3×turn + 2×bend + 4×normal-miss/spacing + 3×crossing-proxy +
   2×unsupported-gap term. Reward observed length plus soft 1.5×target progress.
   No ideal spiral or turn-count objective is present.
5. Beam width 256, maximum depth 180, outgoing cap 12, spatial/directional beam
   diversity; forbid repeated observed chains/pixels and bridge/history crossings
   under the inherited sampled proxy. This is bounded search, not a global optimum.
   Retain complete-candidate archives, top alternatives, graph scores and logs.
6. Apply up to 20 T3 two-edge exchanges, admitting proper and nonparallel point
   contacts. Require strictly fewer detector intersection events, then minimize
   observed loss and new bridge length. Reverse the intervening edge sequence,
   retain its observed geometry, mark replacement edges inferred, and mark nearby
   observed edges low-confidence within the existing 5 px arc neighborhood.
7. Preserve stalled/capped corrections as negative evidence. No new topology
   optimizer or cross-Z rescue is substituted. Zero-length edges introduced by an
   exchange have no arc or signal; remove only those before sampling to avoid 0/0.
   {extra['zero_edges_removed']} such edges were removed; no positive-length observed
   edge is deleted by this numerical handling.
8. Parameterize the complete polyline, including all bridges, by cumulative arc.
   Linear XY interpolation at 0.5 px steps; bilinear original CT only for labels
   1/3 whose interpolated ridge response exceeds 3. Bridge label 2 is always NaN.
   Retain the historical local Z146/Z160 uncertainty marks from T2; these neither
   select geometry nor use another slice's observations.

### Independence barrier

Each worker receives its own XY plane and fixed constants. Its newly extracted
ridge graph lives under `slice_<z>/graph_input/outputs/p0_t1/`; this compatibility
root lets the unchanged T2 loader read only that slice's graph. The T1 conservative
trace is retained for lineage/seed selection; its STOP decision does not truncate
T2 maximum recovery. Two workers perform independent shell/Python work, not
cross-Z tracking.

Each slice has a `FINALIZED.json` hash receipt. `ROWS_FINALIZED.json` is created
only after all 50 receipts exist. `p0_t5_evaluate.py` refuses to compare or stack
before that barrier and verifies all row hashes. All correlations are post-finalization
evaluation. There was no retuning, row shifting, warping, rescaling, or inferred
CT filling after viewing coherence.

## Per-slice results

Observed length includes low-confidence detector-supported contour edges; it is not
proof of physical layer identity. Utilization is observed geometric length / total
geometric length. Low-confidence percentage below is of native row samples.
Bridge count counts final inferred edges; corrections can split/reverse bridge runs.
Full precision, fragment counts, seed/end coordinates, retained fraction and usable
sample counts are in `path_metrics.csv` and each slice's `metrics.json`.

| Z | Total mm | Observed mm | Bridge mm | Utilization % | Bridges | Low-conf % | Corrections | Defects | Exterior | Last sampled s mm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
{path_table}

Observed-contour retention relative to each slice's precleanup T2 best is
**{100*m['observed_retention_range'][0]:.3f}–{100*m['observed_retention_range'][1]:.3f}%**.
Observed/total geometric utilization spans
{100*extra['observed_utilization_minmax'][0]:.3f}–{100*extra['observed_utilization_minmax'][1]:.3f}%.
Path lengths [minimum, Q1, median, Q3, maximum] mm are
{[round(x,3) for x in m['length_mm_quantiles']]}.

### Topology limitations and negative evidence

The independent exhaustive segment-pair auditor, distinct from the KD-tree cleanup
detector, reports aggregate before → after counts:

- Proper crossings: {total_before['proper_crossings']} → {total_after['proper_crossings']}.
- Non-adjacent point-contact segment pairs: {total_before['nonadjacent_point_contacts']} → {total_after['nonadjacent_point_contacts']}.
- Non-adjacent collinear overlaps: {total_before['collinear_overlaps']} → {total_after['collinear_overlaps']}.

There were {extra['total_corrections']} accepted exchanges. Residual defects occur at
Z={m['topology_defective_z']}. Detailed typed counts and termination reasons:

```json
{json.dumps(residual,indent=2)}
```

These rows remain in the requested complete 50-row diagnostic raster, clearly
flagged in metrics. They are not certified single-thread topology successes.
The run04 correction heuristic does not guarantee cleanup on arbitrary new graphs;
its contact detector excludes collinear intersections, and a strict-decrease local
exchange can stall. Final exhaustive counts, rather than the termination string,
determine cleanliness. No unsupported claim that remaining defects are unfixable is
made. This is a limitation of the unchanged method tested here.

## Stacking and provenance

Set s=0 at each independent seed, retain the saved traversal order and **zero offsets**.
Copy each numeric row verbatim into the common rectangle; pad only absent tails
with NaN/mask=true/provenance=0. Inferred bridges keep their geometric s width but
no invented signal. `coverage.npy` separates native domain from padded tails.

- Code 0: no usable observation (or padding); masked NaN.
- Code 1: observed CT-supported contour; original sampled CT.
- Code 2: inferred bridge; masked NaN.
- Code 3: low-confidence observed region; real sampled CT, not discarded.

`stack.npy` is float64, `mask.npy` boolean (true=missing), `provenance.npy` uint8.
Rows are strictly Z128..177. The common last sample is **{m['common_last_s_mm']:.5f} mm**;
all native domains cover s=0 through **{m['all_native_overlap_last_s_mm']:.5f} mm**.
There are {m['all50_joint_valid_columns']} columns valid in all 50 rows,
{extra['columns_valid_in_at_least40rows']} valid in at least 40, and
{extra['columns_valid_in_at_least25rows']} in at least 25. This is numeric domain
coverage, not proof equal s indexes equal material.

Valid CT fills **{100*m['valid_fraction']:.3f}%** of the rectangle
({m['valid_samples']} measurements); **{100*m['masked_fraction']:.3f}%** is masked.
Of the rectangle, {100*m['bridge_fraction']:.3f}% is inferred bridge,
{100*m['native_unresolved_fraction']:.3f}% native unresolved, and
{100*m['padding_fraction']:.3f}% out-of-domain padding. Low-confidence measurements
occupy {100*m['low_confidence_fraction']:.3f}%. Within native row domains,
{100*m['valid_fraction_native']:.3f}% are usable CT samples.

## Consecutive-Z coherence

Every pair has ΔZ=126.56 µm. Pearson r excludes gaps and uses shared native s
coverage. First differences use only consecutive valid samples, never span a masked
gap, and are evaluation-only. The lag column below is a diagnostic maximum within
±64 samples, never an applied alignment. Positive lag compares row z at k with
row z+1 at k+lag. Full valid-overlap counts and controls are retained in CSV.

| Adjacent Z | Joint valid | Raw r | First-difference r | Best diagnostic lag samples |
| --- | --- | --- | --- | --- |
{pair_table}

Mean adjacent r = **{m['mean_adjacent_r']:.5f}**; median = **{m['median_adjacent_r']:.5f}**;
positive fraction = **{100*m['positive_adjacent_fraction']:.2f}%**.
Raw-r [min,P10,Q1,median,Q3,P90,max] =
{[round(v,5) for v in m['raw_r_quantiles']]}.
Mean first-difference r = {m['mean_difference_r']:.5f}.
The {m['clean_pair_count']} pairs with both rows topology-clean have mean raw
r={m['clean_pair_mean_r']:.5f}; topology-defective rows are not silently removed
from the primary all-49 analysis.

Each pair is also compared to 199 deterministic circular shifts (128..N−128)
of its second row within the common native prefix, moving masks with values.
{m['pairs_above_shift95']}/49 raw correlations exceed their pair's shift-control
95th percentile. These are descriptive controls, not independent replicates or
formal p-values: nonstationary CT, masks, circular boundaries and multiple tests
limit interpretation. They prevent equating any small positive correlation with
successful structure recovery.

Local windows span 256 samples (~16.20 mm), stride128 (~8.10 mm), minimum64 valid
pairs. Of {m['local_window_count']} windows, {100*m['local_positive_fraction']:.2f}%
have positive r. The longest fixed-window run of consecutive adjacent pairs with
r≥0.2 is {m['max_local_pair_run_r_ge_0_2']} pairs; this threshold is a descriptive
summary, not a success gate. `local_persistence.json` records exact locations.
`coherence_diagnostic.png` shows correlations versus Z, local windows versus s/Z,
and diagnostic lag sensitivity, rather than only a global average.

### Drift, shear and discontinuities

Median absolute adjacent path-length change is
{extra['median_abs_adjacent_length_change_mm']:.3f} mm; maximum is
{extra['max_abs_adjacent_length_change_mm']:.3f} mm. The diagnostic best-lag sample
[min,Q1,median,Q3,max] is {m['lag_samples_quantiles']}, with median absolute lag
{extra['median_abs_diagnostic_lag_um']/1000:.3f} mm. Weak/multiple maxima cannot be
interpreted as measured physical displacement. No physical shear field is inferred
and no lag is used to improve the stack. Independent branch/path ordering, different
seed locations, inferred-distance accumulation and topology exchanges can cause
abrupt s misregistration. `adjacent_correlations.csv` records length changes and
seed displacement for every pair; these are post hoc diagnostics only.

### Diagnostic best-lag analysis (not an adjusted reconstruction)

The authoritative zero-shift raster is unchanged. Whole-row diagnostic maxima
within ±64 samples (±4.04992 mm) have mean **{lag['global_best_mean']:.5f}**, median
**{lag['global_best_median']:.5f}**, minimum {lag['global_best_min']:.5f} and maximum
{lag['global_best_max']:.5f}. {lag['global_best_positive_pairs']}/49 are positive;
{lag['global_best_boundary_pairs']} maxima occur at the search boundary. These are
optimistically selected maxima, not an achieved reconstruction correlation.
The modest improvement over raw mean {m['mean_adjacent_r']:.5f} does not establish
broad alignment by one horizontal displacement. `lag_correlations.csv` gives the
zero/best comparison for all49 pairs; `diagnostic_lags.csv` retains every tested lag.
`lag_diagnostic.png` plots both versus Z and their histograms.

The zero-shift distribution [min,max] is [{lag['global_raw_min']:.5f},
{lag['global_raw_max']:.5f}]; {sum(p['r']>0 for p in pairs)}/49 are positive.
Exact 0.1-wide histogram bins and counts are in `lag_metrics.json` alongside the
quartile/decile summary above. All49 pairs are evaluable; no pair is omitted.

For each of the **{lag['local_windows']}** local 256-sample windows, independently
scan the same ±64 sample lag range, with at least64 valid overlapping values.
Mean local zero-shift r is {lag['local_zero_mean']:.5f}; mean local peak r is
**{lag['local_best_mean']:.5f}**, median **{lag['local_best_median']:.5f}**.
Local first-difference correlation at the selected raw-signal lag averages
{lag['difference_at_local_raw_peak_mean']:.5f}; it shares the selection bias.
{lag['local_peak_boundary_count']} local maxima hit the search boundary.
No local shift, stretch or interpolated registration is written into a row.

To distinguish displaced structure from trivial peak-search inflation, repeat the
identical local lag maximization for **19 circular-shift controls per window**
(deterministic offsets128..N−128, masks carried with data). Only
**{lag['local_peak_above_control95_count']}/{lag['local_windows']}** observed peaks
exceed their matched peak-control95th percentile. This supports localized
correspondence, not broad fine-scale coherence; overlapping windows and coarse
controls are exploratory, not calibrated independent hypothesis tests.
The longest consecutive-pair run in one fixed window exceeding these controls
and r≥0.3 is **{max(v['longest_adjacent_pair_run'] for v in lag['persistent_control_exceeding_runs'])} pairs**.
The two controls (199 whole-row shifts and19 locally maximized shifts) answer
different questions and must not be conflated.

`local_window_coherence.png` separately shows zero-shift r, best-lag r, selected
local lag and control exceedances. `local_lag_windows.csv` preserves every
window's lag, overlap, r, difference-r and control comparison.
`local_lag_examples.png` shows explicitly selected strong examples, not typical
performance: Z135–136 at s=8.10..24.24mm improves approximately0.612→0.685 with
lag−15 samples; Z167–168 at s=72.90..89.03mm improves0.083→0.639 with lag+44.
These examples contain shared broad intensity variation, with unresolved fine
detail and gaps. Their selection is disclosed; they do not demonstrate coherent
structure across a substantial majority of the stack.

The local lag field changes sign irregularly across s and Z rather than forming
a stable smooth shear. Major length discontinuities include Z146→147
(−114.64mm), Z167→168 (+99.98mm), and Z140→141 (−98.46mm). Full locations and
seed displacements are saved in `lag_metrics.json`. Distortion/ordering instability
can explain some weak zero-shift correspondence, but a simple small-lag model
does not explain the whole stack. The data also do not justify saying every
adjacent row is unrelated.

## Visual evaluation and image mappings

{assessment['visual_observations']}

Inspected: {', '.join(assessment['inspected'])}.

`raw_stack.png` is the native {m['shape'][1]}×50 image, fixed original 0–255 CT
mapping rounded to uint8; missing values display black and are distinguished by
`mask_stack.png` (255 missing, 0 valid). `raw_stack_contrast.png` applies one global
pooled 1st–99th percentile window [{display['low']:.12f}, {display['high']:.12f}]
for display only. No row-wise enhancement, smoothing or text-oriented processing
is used. Provenance colors: black0, cyan1, magenta2, orange3.

`physical_aspect_stack.png` repeats each original display row twice vertically
because ΔZ/Δs=2, without horizontal resizing; it is {m['shape'][1]}×100 pixels.
Underlying numeric data remain 50 rows. This gives correct physical pixel aspect,
not an enlarged square illustration. `stack_inspection.png` intentionally enlarges
vertical display for inspection and is not a physical-aspect image. Original
within-row sampling has never been resized.

## Reproduction and verification

Commands from repository root (generation refuses to overwrite the run):

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t5 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t5_consecutive.py --run
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t5 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t5_evaluate.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t5 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t5_lag_diagnostics.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t5 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t5_validate.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t5 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t5_evaluate.py --verify
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t5 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t5_validate_diagnostics.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t5 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t5_report.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t5 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t5_seal.py
sha256sum -c outputs/p0_t5/SHA256SUMS
```

The report consumes the recorded human-readable visual assessment; it does not
select a verdict automatically from correlations. `dependency_versions.json` and
`plan.json` record environment and all inherited parameters. `run.log` preserves
execution output. `geometry_replay_validation.json` confirms all 50 source planes,
exact regenerated 2-D ridge graphs, saved correction replay, continuous geometry
and cumulative arc, independently interpolated XY and four-corner bilinear CT
within1e-10, masks/provenance and exhaustive topology-count reproduction.
This replay does not rerun the expensive beam search or change any path.

`verification.json` confirms exact rectangular-stack replay, masks/provenance,
all49 correlation recomputations, PNG decoding and preservation of prior evidence.
The five historical Z rows reproduce sealed run04 row/mask/provenance arrays:
**{'ALL EXACT' if all(v[k] for v in ver['historical_T3_replay'] for k in ['row','mask','provenance']) else 'see verification.json for differences'}**.
Validation PASS means recorded data and counts are reproducible, not that all
rows have clean topology or broad coherence. README/AGENTS and T1–T4 evidence
remain byte-identical. A final SHA256SUMS covers this report, T5 scripts and all
T5 evidence; immediate verification seals the result.

## Artifact inventory and stopping point

All outputs are under `outputs/p0_t5/`:

- Required PNGs: raw_stack, raw_stack_contrast, provenance_stack, mask_stack,
  physical_aspect_stack, coherence_diagnostic, lag_diagnostic, local_window_coherence.
- Required tables: row_extents.csv, adjacent_correlations.csv, lag_correlations.csv,
  path_metrics.csv.
- Lossless stack.npy, mask.npy, provenance.npy, plus coverage and coordinate arrays.
- Per-slice CT, candidate graph/response/skeleton, virtual-edge scores, all complete
  and top paths, precleanup path, corrections, final path-and-samples NPZ,
  path/sample CSVs, row/mask/provenance NPYs, geometry PNG and validation/seal.
- Post-finalization diagnostics: local_windows.csv, diagnostic_lags.csv,
  shift_controls.csv, local_persistence.json, stack_metrics.json and report_measurements.json.

{assessment['limitations']}

**READY_FOR_NEXT_STEP = NO**, interpreted as readiness for an unchanged larger-Z
scale-up. The completed experiment supports further bounded diagnosis, not a
larger acquisition or longer identical run yet. This is a scientific recommendation,
not an unfinished execution task or a requirement for cross-Z physical tracking.

The Git milestone excludes lossless NPY/NPZ arrays, logs and full candidate-path
archives under existing/local ignore rules. They remain on disk, unmodified and
covered by the full SHA256SUMS. A Git clone alone is therefore not the complete
replay package; it needs these preserved local artifacts. `git_artifact_inventory.json`
records included/excluded paths and sizes. The commit/push outcome is reported
in the final status; no force push is permitted.

`diagnostic_validation.json` independently replays all1001 local windows across
129 lags with scalar Pearson calculations, validates histogram totals and checks
contrast/mask/provenance images. No authoritative sample changes. The validation
is numerical correctness, not a claim that the selected peaks establish physical
layer correspondence.

No OCR, ink detector, language prior, mathematical-spiral fit, learned correction,
mesh, 3-D surface or cross-Z tracking was used. No readable-letter criterion was
applied. {assessment['next_proposed']}

**Stop after P0-T5. No next experiment executed.**
'''
    (ROOT/'docs/p0_t5_consecutive_thread_stack.md').write_text(text)
    print('Report written:',assessment['result'])

if __name__=='__main__':make_report()
