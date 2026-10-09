# P0-T4 — first maximum-recovery Direct Thread-Stack

## Result

**P0_T4_COMPLETE = YES. P0_T4_RESULT = PARTIAL.**

The five unchanged, independently unwound P0-T3 rows form a valid **5 × 3772**
numeric raster on a common physical arc-distance grid. Masks, inferred gaps and
low-confidence observations survive exactly. No row was shifted, resized, reversed,
warped or retraced. The representation succeeds; consistent fine-scale spatial
coherence is **not established**. This is not a failure based on letter readability
or absence of proven cross-Z fragment identity.

## Hypothesis, inputs and scope

Test whether independently reconstructed center-to-exterior threads can be stacked
by Z on their existing physical arc-length coordinates while preserving CT-supported
signal, and whether neighboring available rows show non-random correspondence.

Authoritative inputs are `docs/p0_t3_unwindable_threads.md` and
`outputs/p0_t3/run_04/`. For each **Z = 128, 144, 146, 160, 170**, load
`slice_<z>/row.npy`, `mask.npy`, `provenance.npy` and `path_and_samples.npz`;
run-level `metrics.json` supplies complete polyline lengths. These are PHercParis4 /
Scroll 1 transverse XY slices from the preserved
`outputs/orientation/ct_level4_zyx.npy`, level 4 CT with isotropic **126.56 µm**
voxel spacing. T3 sampled the original CT bilinearly on its finalized polylines;
T4 copies those numeric measurements, without sampling CT again.

These five slices are **not consecutive or uniformly spaced**. Their center-to-center
Z span is **5.31552 mm**; successive available-row separations are 2.02496,
0.25312, 1.77184 and 1.26560 mm. There are 38 unobserved integer Z levels within
128–170. No observations are invented at those levels. The primary five-row PNG
uses row index vertically; `z_indices.npy` and `z_relative_um.npy` supply actual
positions. `row_extent_diagnostic.png` also shows the 43-level physical-Z grid
with unobserved levels blank, not interpolated.

## Coordinate and representation

Use **s = 0 at each unchanged T3 interior seed**, increasing in its saved traversal
order. Every row uses **Δs = 63.28 µm = 0.5 CT pixel**. All offsets are zero.
This is a common seed-relative distance coordinate, **not a registered material
coordinate**. Independent seeds, inferred bridge lengths and path ordering can
give the same s different anatomical meanings. T3 topology exchanges sometimes
reversed long subsequences; T4 does not undo or adjust those finalized paths.

The union grid has 3772 sample centers from 0 to **238.62888 mm**. Its nominal
3772-cell width is 238.69216 mm; the longest continuous input path is 238.661696 mm.
These are different quantities, not resizing. All five native sample domains
overlap from 0 through **159.71872 mm**. There are **1027** columns with valid
measurements simultaneously in all five rows.

`stack.npy` is float64, shape (5, 3772), with NaN wherever missing. `mask.npy`
is boolean, true for missing. `provenance.npy` is uint8 with unchanged codes:

| Code | Meaning | Numeric signal |
| --- | --- | --- |
| 0 | Unresolved/no usable data, or padded tail | NaN |
| 1 | Observed contour | Original T3 CT value |
| 2 | Inferred bridge | NaN; bridge distance remains in s |
| 3 | Low-confidence observed contour | Original T3 CT value |

`coverage.npy` distinguishes native unresolved code 0 from out-of-domain padding.
No inferred bridge receives CT intensity. Low-confidence samples remain valid.
The fixed s grid and actual Z coordinates are saved losslessly alongside the stack.

## Row extents and data accounting

Fractions below refer to sampled data, not T3 geometric observed-arc percentages.
Some observed edges have unresolved samples when T3 ridge support was inadequate.

| Z | Full path length mm | Native width | Valid (1+3) | Bridge (2) | Native unresolved (0) | Low confidence (3) | Padded tail | Native valid % |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 128 | 159.726081 | 2525 | 2130 | 320 | 75 | 435 | 1247 | 84.356 |
| 144 | 238.661696 | 3772 | 3124 | 540 | 108 | 757 | 0 | 82.821 |
| 146 | 224.136168 | 3542 | 2960 | 484 | 98 | 647 | 230 | 83.569 |
| 160 | 167.532660 | 2648 | 2269 | 315 | 64 | 510 | 1124 | 85.687 |
| 170 | 195.909455 | 3096 | 2587 | 424 | 85 | 745 | 676 | 83.559 |

Of the 18,860 cells in the rectangular stack, **13,070 (69.3001%) are valid**,
and **30.6999% are masked**. Missing cells comprise 2083 inferred bridge samples
(11.0445% of the rectangle), 430 native unresolved samples, and 3277 padded cells.
There are 3094 low-confidence measurements: 16.4051% of the rectangle and 23.6725%
of valid data. Across native row domains, 83.8735% of samples are valid.
The exact per-row lengths, last sample coordinates and fractions are in
`row_extents.csv`; no row-width normalization has been applied.

## Coherence measurements

Pearson correlations use only jointly finite samples within each pair's shared
native domain. “Neighbor” means adjacent **available** rows, not adjacent CT levels.
First-difference correlations use only valid consecutive sample differences and
do not differentiate across masked gaps. Local correlations use 256-sample windows,
stride 128, requiring at least 64 joint observations. No intensity filtering is
used in the raster; differencing is diagnostic only.

| Pair Z | Shared last s mm | Joint valid samples | Joint valid / shared % | Raw Pearson r | First-difference r | Median local r | Positive local windows % |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 128–144 | 159.71872 | 1778 | 70.416 | 0.10433 | −0.00143 | −0.11678 | 27.778 |
| 144–146 | 224.07448 | 2452 | 69.226 | −0.06004 | 0.04211 | 0.01367 | 53.846 |
| 146–160 | 167.50216 | 1878 | 70.921 | 0.05190 | 0.03242 | 0.08434 | 68.421 |
| 160–170 | 167.50216 | 1910 | 72.130 | 0.20925 | −0.01968 | 0.06270 | 73.684 |

To avoid equating any positive r with meaningful alignment, compare against 199
deterministic circular shifts of the second row within the shared native prefix,
including its mask. Shift offsets are unique integers linearly spaced from 128
through N−128 samples. Raw-control 95th percentiles are respectively 0.13788,
0.16903, 0.16911 and 0.11030. Only 160–170 exceeds its raw-control 95th percentile;
none of its 199 controls matches its raw r. But its first-difference r is negative,
so this is evidence for broad intensity correspondence, not convincing fine detail.
The closest pair, 144–146, has negative raw correlation; its small difference r
exceeds its difference-control 95th percentile (0.03992) only narrowly.

These exploratory shifted controls are not independent replicates or a formal
significance test. Nonstationary CT, mask layout, circular boundaries and multiple
diagnostics limit inference. No zero-tail-frequency result is called a zero p-value.
All individual controls and measurements are preserved in CSVs.

### Lateral displacement and shear

A **diagnostic-only** lag sweep of −64 to +64 samples (±4.04992 mm), with at least
128 joint values, finds pairwise maximum r at −44, −52, −11 and −49 samples
(−2.78432, −3.29056, −0.69608 and −3.10072 mm). Their maximum r values are
0.21897, −0.00737, 0.10595 and 0.24544. Positive lag compares the first row at k
with the second at k+lag. **No such lag is applied to any output row.**

Different weak maxima do not establish a measurable physical shear or justify
registration. Even a shift does not yield positive raw r for the closest pair in
this bounded range. Path-order differences and accumulated distance through
inferred bridges can produce nonuniform misregistration that one offset cannot
describe. The stack has zero imposed stretch/shear; an anatomical shear estimate
is unsupported. No post-hoc transforms were introduced to improve metrics.

## Visualization and inspection

`raw_stack.png` is exactly 3772 × 5 pixels: original CT rounded to uint8, fixed
0–255 mapping, masked NaNs displayed black. The separate mask distinguishes
missing black from real dark measurements. `raw_stack_contrast.png` uses a single
pooled valid-data 1st–99th percentile window, **49.96173966396096 to
141.13257486883427**, clipped and rounded to 0–255. This visualization-only mapping
does not modify numeric arrays or enhance rows separately.

`mask_stack.png` is white (255) for missing, black (0) for valid.
`provenance_stack.png` uses black for 0, cyan for 1, magenta for 2 and orange for 3.
`stack_inspection.png` enlarges all five rows with nearest-neighbor display, shows
fixed/global-contrast/provenance panels, and colors missing CT magenta.

Visual inspection of that overview, `coherence_diagnostic.png`, and
`row_extent_diagnostic.png` shows intact within-row intensity variation, localized
internal gaps, and substantial right-tail blank regions from unequal lengths.
Some broad light/dark regions approximately correspond, but narrow structures
do not form consistently aligned features across the five rows. Local correlations
change sign along s. No stable shear direction is visually identifiable. At true
Z spacing the stack is extremely sparse. These observations support PARTIAL,
not a claim that the independent-slice hypothesis has been disproved.

## Reproduction, artifacts and validation

Implementation: `scripts/p0_t4_stack.py`. Parameters and preservation baseline
are saved in `outputs/p0_t4/plan.json` and `preservation.json`. The generation
command, run from the repository root, was:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t4 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t4_stack.py
```

Generation requires a nonexistent output directory and refuses to overwrite this
run. The read-only numeric replay command is:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t4 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t4_stack.py --verify
sha256sum -c outputs/p0_t4/SHA256SUMS
```

Replay passed: T3 seal, exact reconstruction of stack/mask/provenance/coverage,
raw and mask PNG agreement, recomputed pair diagnostics, and all previously
recorded file hashes including README/AGENTS. `replay_validation.txt` records the
result. Supplemental representation checks are in `representation_validation.json`.
The T4 checksum manifest covers this report, implementation and all T4 evidence
other than the manifest itself. It was verified immediately after generation.
No tracing experiment was replayed and no previous evidence was modified.

Required files under `outputs/p0_t4/`:

- `raw_stack.png`, `raw_stack_contrast.png`, `mask_stack.png`, `provenance_stack.png`.
- `row_extent_diagnostic.png`, `row_extents.csv`.
- `stack.npy`, `mask.npy`, `provenance.npy` (lossless).

Additional evidence: `coverage.npy`, `s_um.npy`, `z_indices.npy`,
`z_relative_um.npy`, `stack_inspection.png`, `coherence_diagnostic.png`,
`metrics.json`, `display.json`, `neighbor_correlations.json`, and CSVs
`neighbor_correlations.csv`, `window_correlations.csv`, `lag_correlations.csv`,
`shift_null.csv`. All diagnostics and negative evidence are retained.

## Limitations and conclusion

The seed-relative physical coordinate is well-defined without arbitrary transforms,
but has not been shown to align equivalent structures between independent paths.
Five unequally spaced observations, independent seed/path choices, coarse CT,
inferred bridge lengths and retained uncertain connections limit coherence.
Masks preserve provenance rather than erase downstream recovery. T3 topology
cleanliness is inherited, not newly claimed as proof of a physical sheet.

**PARTIAL:** a valid maximum-recovery 2-D raster exists with exact source data,
spacing, gaps and provenance preserved. Broad correspondence occurs in one pair,
but the available evidence does not establish reliable fine-scale coherence.
Readable letters were not evaluated. No cross-Z tracking, 3-D surface, OCR, ink
inference, learned correction or path modification was used.

A future separately authorized experiment could test a reproducible per-slice
seed/order convention on more closely spaced independent slices; this report does
not execute that experiment or prescribe physical-fragment tracking. **Stop at P0-T4.**
