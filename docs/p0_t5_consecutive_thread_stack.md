# P0-T5 — consecutive-slice Direct Thread-Stack

## Result and scope

**P0_T5_COMPLETE = YES. P0_T5_RESULT = PARTIAL.**

A valid 50-row consecutive-Z raster preserves maximum-recovery contour data, inferred gaps and uncertainty. Some localized CT correspondence survives small-lag diagnostics, but substantial path-length/order instability and irregular misregistration prevent broad fine-scale coherence. This is neither a PASS nor evidence that every neighboring row is unrelated.

All 50 independent reconstructions produced numeric rows. 50/50
(100%) reached the inherited exterior target; 45/50
are both exterior-reaching and topology-clean. The stack shape is **50 × 4482**.
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
**6.20144 mm** between first and last centers (50 voxel cells
span 6.32800 mm). Use `z_indices.npy`, `z_relative_um.npy` and `s_um.npy` to
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
   10 such edges were removed; no positive-length observed
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
| 128 | 159.726 | 139.303 | 20.423 | 87.21 | 38 | 17.23 | 1 | 0 | YES | 159.719 |
| 129 | 198.373 | 172.580 | 25.793 | 87.00 | 53 | 22.30 | 3 | 1 | YES | 198.320 |
| 130 | 249.602 | 214.233 | 35.369 | 85.83 | 65 | 21.88 | 3 | 0 | YES | 249.576 |
| 131 | 199.995 | 172.957 | 27.037 | 86.48 | 50 | 19.20 | 1 | 0 | YES | 199.965 |
| 132 | 195.671 | 163.672 | 31.999 | 83.65 | 54 | 21.63 | 2 | 0 | YES | 195.662 |
| 133 | 219.879 | 190.032 | 29.847 | 86.43 | 63 | 20.98 | 5 | 0 | YES | 219.835 |
| 134 | 249.176 | 215.486 | 33.690 | 86.48 | 62 | 20.62 | 1 | 0 | YES | 249.133 |
| 135 | 206.072 | 179.441 | 26.631 | 87.08 | 54 | 23.03 | 2 | 1 | YES | 206.040 |
| 136 | 189.945 | 162.790 | 27.155 | 85.70 | 57 | 23.92 | 2 | 0 | YES | 189.903 |
| 137 | 225.569 | 194.203 | 31.365 | 86.09 | 60 | 20.14 | 2 | 0 | YES | 225.530 |
| 138 | 197.857 | 171.119 | 26.738 | 86.49 | 54 | 17.97 | 1 | 0 | YES | 197.813 |
| 139 | 198.253 | 175.012 | 23.242 | 88.28 | 51 | 18.86 | 0 | 0 | YES | 198.193 |
| 140 | 231.259 | 199.988 | 31.271 | 86.48 | 66 | 21.78 | 3 | 0 | YES | 231.225 |
| 141 | 132.796 | 115.188 | 17.608 | 86.74 | 32 | 24.49 | 0 | 0 | YES | 132.761 |
| 142 | 120.989 | 105.489 | 15.500 | 87.19 | 28 | 19.14 | 1 | 0 | YES | 120.928 |
| 143 | 197.292 | 171.535 | 25.757 | 86.94 | 51 | 18.22 | 2 | 0 | YES | 197.244 |
| 144 | 238.662 | 204.841 | 33.821 | 85.83 | 62 | 20.07 | 0 | 0 | YES | 238.629 |
| 145 | 283.569 | 243.164 | 40.404 | 85.75 | 61 | 14.06 | 0 | 1 | YES | 283.558 |
| 146 | 224.136 | 193.802 | 30.334 | 86.47 | 54 | 18.27 | 3 | 0 | YES | 224.074 |
| 147 | 109.496 | 96.743 | 12.753 | 88.35 | 26 | 15.19 | 0 | 0 | YES | 109.474 |
| 148 | 192.714 | 166.602 | 26.111 | 86.45 | 47 | 12.80 | 1 | 0 | YES | 192.688 |
| 149 | 143.051 | 126.414 | 16.637 | 88.37 | 29 | 15.70 | 1 | 0 | YES | 143.013 |
| 150 | 155.360 | 135.392 | 19.968 | 87.15 | 41 | 20.11 | 1 | 0 | YES | 155.352 |
| 151 | 206.018 | 179.190 | 26.828 | 86.98 | 50 | 19.32 | 3 | 1 | YES | 205.976 |
| 152 | 239.981 | 203.836 | 36.145 | 84.94 | 58 | 20.35 | 1 | 0 | YES | 239.958 |
| 153 | 208.464 | 179.983 | 28.481 | 86.34 | 57 | 23.70 | 3 | 0 | YES | 208.444 |
| 154 | 202.030 | 176.058 | 25.972 | 87.14 | 58 | 19.98 | 1 | 0 | YES | 201.990 |
| 155 | 241.137 | 206.759 | 34.378 | 85.74 | 59 | 20.07 | 3 | 0 | YES | 241.097 |
| 156 | 213.633 | 182.816 | 30.817 | 85.57 | 59 | 18.95 | 1 | 0 | YES | 213.633 |
| 157 | 222.233 | 188.837 | 33.396 | 84.97 | 67 | 21.24 | 2 | 0 | YES | 222.176 |
| 158 | 192.037 | 165.901 | 26.136 | 86.39 | 52 | 20.43 | 0 | 0 | YES | 191.992 |
| 159 | 181.220 | 152.904 | 28.316 | 84.37 | 50 | 20.08 | 1 | 0 | YES | 181.171 |
| 160 | 167.533 | 147.661 | 19.872 | 88.14 | 41 | 19.26 | 1 | 0 | YES | 167.502 |
| 161 | 202.335 | 174.189 | 28.146 | 86.09 | 62 | 26.02 | 2 | 0 | YES | 202.306 |
| 162 | 224.561 | 194.625 | 29.936 | 86.67 | 55 | 19.86 | 1 | 0 | YES | 224.517 |
| 163 | 198.356 | 175.450 | 22.907 | 88.45 | 47 | 23.48 | 1 | 0 | YES | 198.320 |
| 164 | 239.239 | 205.407 | 33.832 | 85.86 | 58 | 19.31 | 2 | 1 | YES | 239.198 |
| 165 | 225.419 | 193.647 | 31.772 | 85.91 | 59 | 20.71 | 2 | 0 | YES | 225.403 |
| 166 | 135.072 | 120.542 | 14.530 | 89.24 | 30 | 18.50 | 0 | 0 | YES | 135.040 |
| 167 | 131.270 | 116.564 | 14.705 | 88.80 | 28 | 19.08 | 0 | 0 | YES | 131.243 |
| 168 | 231.247 | 197.823 | 33.424 | 85.55 | 63 | 20.33 | 1 | 0 | YES | 231.225 |
| 169 | 136.328 | 118.009 | 18.319 | 86.56 | 30 | 17.68 | 1 | 0 | YES | 136.305 |
| 170 | 195.909 | 168.854 | 27.055 | 86.19 | 51 | 24.06 | 0 | 0 | YES | 195.852 |
| 171 | 179.361 | 156.358 | 23.003 | 87.17 | 43 | 20.28 | 0 | 0 | YES | 179.336 |
| 172 | 221.790 | 190.277 | 31.513 | 85.79 | 60 | 20.20 | 0 | 0 | YES | 221.733 |
| 173 | 180.259 | 158.589 | 21.670 | 87.98 | 42 | 20.71 | 1 | 0 | YES | 180.221 |
| 174 | 209.662 | 185.184 | 24.478 | 88.32 | 48 | 20.10 | 1 | 0 | YES | 209.647 |
| 175 | 176.328 | 153.790 | 22.537 | 87.22 | 42 | 17.87 | 1 | 0 | YES | 176.298 |
| 176 | 226.855 | 196.348 | 30.507 | 86.55 | 56 | 21.79 | 2 | 0 | YES | 226.796 |
| 177 | 145.880 | 127.616 | 18.263 | 87.48 | 32 | 20.08 | 0 | 0 | YES | 145.860 |

Observed-contour retention relative to each slice's precleanup T2 best is
**99.690–100.000%**.
Observed/total geometric utilization spans
83.646–89.242%.
Path lengths [minimum, Q1, median, Q3, maximum] mm are
[109.496, 179.585, 199.184, 224.455, 283.569].

### Topology limitations and negative evidence

The independent exhaustive segment-pair auditor, distinct from the KD-tree cleanup
detector, reports aggregate before → after counts:

- Proper crossings: 51 → 4.
- Non-adjacent point-contact segment pairs: 48 → 1.
- Non-adjacent collinear overlaps: 1 → 0.

There were 66 accepted exchanges. Residual defects occur at
Z=[129, 135, 145, 151, 164]. Detailed typed counts and termination reasons:

```json
[
  {
    "z": 129,
    "proper": 1,
    "contacts": 0,
    "overlap": 0,
    "reason": "no crossing-decreasing exchange"
  },
  {
    "z": 135,
    "proper": 1,
    "contacts": 0,
    "overlap": 0,
    "reason": "no crossing-decreasing exchange"
  },
  {
    "z": 145,
    "proper": 1,
    "contacts": 0,
    "overlap": 0,
    "reason": "no crossing-decreasing exchange"
  },
  {
    "z": 151,
    "proper": 0,
    "contacts": 1,
    "overlap": 0,
    "reason": "no crossing-decreasing exchange"
  },
  {
    "z": 164,
    "proper": 1,
    "contacts": 0,
    "overlap": 0,
    "reason": "no crossing-decreasing exchange"
  }
]
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
Rows are strictly Z128..177. The common last sample is **283.55768 mm**;
all native domains cover s=0 through **109.47440 mm**.
There are 1 columns valid in all 50 rows,
1716 valid in at least 40, and
3106 in at least 25. This is numeric domain
coverage, not proof equal s indexes equal material.

Valid CT fills **58.290%** of the rectangle
(130629 measurements); **41.710%** is masked.
Of the rectangle, 9.353% is inferred bridge,
1.851% native unresolved, and
30.506% out-of-domain padding. Low-confidence measurements
occupy 13.963%. Within native row domains,
83.878% are usable CT samples.

## Consecutive-Z coherence

Every pair has ΔZ=126.56 µm. Pearson r excludes gaps and uses shared native s
coverage. First differences use only consecutive valid samples, never span a masked
gap, and are evaluation-only. The lag column below is a diagnostic maximum within
±64 samples, never an applied alignment. Positive lag compares row z at k with
row z+1 at k+lag. Full valid-overlap counts and controls are retained in CSV.

| Adjacent Z | Joint valid | Raw r | First-difference r | Best diagnostic lag samples |
| --- | --- | --- | --- | --- |
| 128–129 | 1785 | 0.12190 | -0.05798 | -44 |
| 129–130 | 2162 | 0.11563 | -0.01239 | -64 |
| 130–131 | 2187 | 0.15079 | 0.05517 | -36 |
| 131–132 | 2090 | 0.03050 | -0.03613 | -16 |
| 132–133 | 2112 | -0.18063 | 0.02210 | 64 |
| 133–134 | 2435 | 0.08900 | 0.00191 | 62 |
| 134–135 | 2311 | 0.02336 | -0.01030 | -5 |
| 135–136 | 2085 | 0.15396 | 0.01730 | 49 |
| 136–137 | 2028 | 0.05198 | -0.03369 | -15 |
| 137–138 | 2172 | 0.00694 | 0.01104 | -64 |
| 138–139 | 2182 | 0.05924 | 0.00332 | 13 |
| 139–140 | 2196 | 0.03552 | -0.02128 | 43 |
| 140–141 | 1457 | 0.14249 | -0.03315 | 5 |
| 141–142 | 1315 | 0.08784 | 0.00559 | -50 |
| 142–143 | 1362 | -0.00402 | -0.01603 | 62 |
| 143–144 | 2132 | 0.07794 | 0.03870 | -35 |
| 144–145 | 2627 | 0.11278 | -0.00734 | -51 |
| 145–146 | 2537 | -0.05141 | -0.04674 | -14 |
| 146–147 | 1257 | 0.04856 | -0.01253 | -31 |
| 147–148 | 1245 | 0.06758 | 0.02338 | 34 |
| 148–149 | 1623 | 0.11114 | 0.01933 | -19 |
| 149–150 | 1618 | 0.15406 | 0.02625 | 27 |
| 150–151 | 1721 | 0.02775 | 0.01333 | -64 |
| 151–152 | 2271 | -0.07999 | 0.01273 | -39 |
| 152–153 | 2262 | -0.11981 | -0.04962 | -35 |
| 153–154 | 2216 | 0.12987 | 0.00579 | -47 |
| 154–155 | 2211 | -0.13641 | 0.02120 | -23 |
| 155–156 | 2289 | 0.09960 | -0.01027 | -2 |
| 156–157 | 2299 | 0.02381 | 0.05393 | 22 |
| 157–158 | 2104 | 0.12090 | 0.00564 | 21 |
| 158–159 | 1929 | 0.01221 | -0.01415 | 35 |
| 159–160 | 1831 | 0.05785 | 0.01030 | -34 |
| 160–161 | 1908 | 0.17390 | -0.00624 | 63 |
| 161–162 | 2238 | 0.15562 | 0.01650 | 13 |
| 162–163 | 2248 | 0.07276 | -0.00401 | 56 |
| 163–164 | 2251 | 0.11040 | -0.03896 | 30 |
| 164–165 | 2458 | 0.02291 | 0.01614 | -64 |
| 165–166 | 1574 | -0.00643 | -0.00277 | 54 |
| 166–167 | 1552 | 0.20058 | 0.03557 | -33 |
| 167–168 | 1471 | -0.13526 | 0.00671 | 44 |
| 168–169 | 1515 | -0.05659 | -0.02630 | -64 |
| 169–170 | 1521 | 0.03248 | -0.00505 | -54 |
| 170–171 | 2010 | 0.15237 | 0.01521 | -11 |
| 171–172 | 2029 | 0.13773 | 0.04366 | -46 |
| 172–173 | 2110 | 0.15993 | 0.02763 | 49 |
| 173–174 | 2116 | -0.11543 | -0.01166 | -21 |
| 174–175 | 2045 | -0.03384 | 0.01520 | 62 |
| 175–176 | 1957 | 0.04161 | 0.01938 | 21 |
| 176–177 | 1680 | 0.07065 | -0.02196 | 2 |

Mean adjacent r = **0.05152**; median = **0.05924**;
positive fraction = **77.55%**.
Raw-r [min,P10,Q1,median,Q3,P90,max] =
[-0.18063, -0.08708, 0.01221, 0.05924, 0.1209, 0.15398, 0.20058].
Mean first-difference r = 0.00131.
The 39 pairs with both rows topology-clean have mean raw
r=0.05044; topology-defective rows are not silently removed
from the primary all-49 analysis.

Each pair is also compared to 199 deterministic circular shifts (128..N−128)
of its second row within the common native prefix, moving masks with values.
9/49 raw correlations exceed their pair's shift-control
95th percentile. These are descriptive controls, not independent replicates or
formal p-values: nonstationary CT, masks, circular boundaries and multiple tests
limit interpretation. They prevent equating any small positive correlation with
successful structure recovery.

Local windows span 256 samples (~16.20 mm), stride128 (~8.10 mm), minimum64 valid
pairs. Of 1001 windows, 51.55%
have positive r. The longest fixed-window run of consecutive adjacent pairs with
r≥0.2 is 3 pairs; this threshold is a descriptive
summary, not a success gate. `local_persistence.json` records exact locations.
`coherence_diagnostic.png` shows correlations versus Z, local windows versus s/Z,
and diagnostic lag sensitivity, rather than only a global average.

### Drift, shear and discontinuities

Median absolute adjacent path-length change is
34.802 mm; maximum is
114.640 mm. The diagnostic best-lag sample
[min,Q1,median,Q3,max] is [-64.0, -36.0, -11.0, 34.0, 64.0], with median absolute lag
2.215 mm. Weak/multiple maxima cannot be
interpreted as measured physical displacement. No physical shear field is inferred
and no lag is used to improve the stack. Independent branch/path ordering, different
seed locations, inferred-distance accumulation and topology exchanges can cause
abrupt s misregistration. `adjacent_correlations.csv` records length changes and
seed displacement for every pair; these are post hoc diagnostics only.

### Diagnostic best-lag analysis (not an adjusted reconstruction)

The authoritative zero-shift raster is unchanged. Whole-row diagnostic maxima
within ±64 samples (±4.04992 mm) have mean **0.12852**, median
**0.13492**, minimum -0.09284 and maximum
0.28978. 45/49 are positive;
6 maxima occur at the search boundary. These are
optimistically selected maxima, not an achieved reconstruction correlation.
The modest improvement over raw mean 0.05152 does not establish
broad alignment by one horizontal displacement. `lag_correlations.csv` gives the
zero/best comparison for all49 pairs; `diagnostic_lags.csv` retains every tested lag.
`lag_diagnostic.png` plots both versus Z and their histograms.

The zero-shift distribution [min,max] is [-0.18063,
0.20058]; 38/49 are positive.
Exact 0.1-wide histogram bins and counts are in `lag_metrics.json` alongside the
quartile/decile summary above. All49 pairs are evaluable; no pair is omitted.

For each of the **1001** local 256-sample windows, independently
scan the same ±64 sample lag range, with at least64 valid overlapping values.
Mean local zero-shift r is 0.01251; mean local peak r is
**0.30084**, median **0.29722**.
Local first-difference correlation at the selected raw-signal lag averages
0.13145; it shares the selection bias.
56 local maxima hit the search boundary.
No local shift, stretch or interpolated registration is written into a row.

To distinguish displaced structure from trivial peak-search inflation, repeat the
identical local lag maximization for **19 circular-shift controls per window**
(deterministic offsets128..N−128, masks carried with data). Only
**108/1001** observed peaks
exceed their matched peak-control95th percentile. This supports localized
correspondence, not broad fine-scale coherence; overlapping windows and coarse
controls are exploratory, not calibrated independent hypothesis tests.
The longest consecutive-pair run in one fixed window exceeding these controls
and r≥0.3 is **2 pairs**.
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

The authoritative raw and globally contrasted overviews show short patches of approximately vertical or diagonal intensity agreement, but predominantly interrupted, row-specific texture rather than long persistent fine structures. Right-tail coverage changes abruptly between rows. The zero-shift local map alternates positive and negative regions. Local best-lag plots show stronger selected correspondence, including shared broad valleys/peaks in the disclosed examples, but the lag map is irregular and matched-control exceedances are sparse. No stable global shear is apparent. Missing CT remains visibly masked; neither bridge filling nor a registered raster was created.

Inspected: stack_inspection.png, coherence_diagnostic.png, lag_diagnostic.png, local_window_coherence.png, local_lag_examples.png, slice_129/geometry.png.

`raw_stack.png` is the native 4482×50 image, fixed original 0–255 CT
mapping rounded to uint8; missing values display black and are distinguished by
`mask_stack.png` (255 missing, 0 valid). `raw_stack_contrast.png` applies one global
pooled 1st–99th percentile window [54.199009993845, 140.526781965985]
for display only. No row-wise enhancement, smoothing or text-oriented processing
is used. Provenance colors: black0, cyan1, magenta2, orange3.

`physical_aspect_stack.png` repeats each original display row twice vertically
because ΔZ/Δs=2, without horizontal resizing; it is 4482×100 pixels.
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
**ALL EXACT**.
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

The fixed coarse ridge graph can fuse or switch apparent layers; high observed-ridge utilization is not certification of physical sheet identity. Forty-five rows are topology-clean; five retain one residual defect and are flagged rather than repaired post hoc. Independent seed positions, branch choices, different path orders, and inferred distances can all disturb s correspondence. Diagnostics cannot uniquely attribute every discontinuity to one cause. Strong windows were selected after evaluation and are not typical performance. The learned/3-D alternatives and text readability were not tested.

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
applied. A separately authorized bounded investigation of per-slice seed/order and topology robustness would be more informative than immediately extending the unchanged method to a larger Z span. No such follow-up was performed.

**Stop after P0-T5. No next experiment executed.**
