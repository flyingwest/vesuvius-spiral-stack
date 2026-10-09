# P0.1 — one XY slice, candidate trace and grayscale row

## Result: physical success criterion NOT MET

A continuous open candidate curve was traced, parameterized by arc length,
and sampled into a **556 × 1 pixel grayscale row**. Numerical replay is exact.
However, **one physically plausible winding without obvious layer jumps was
not established**. The candidate has suspected layer changes in the lower arc,
ambiguous contacts in the upper arc, and sensitivity to a small guide shift.
Treat `traced_slice.png`, `trace_coordinates.csv`, and `unwrapped_row.png` as
**unvalidated candidate artifacts**, not a successful papyrus unwrap.

This negative result concerns this simple method on this coarse slice. It does
not falsify Direct Spiral-Stack on native-resolution CT. No second slice was
traced or stacked; no OCR, language priors, learned detector, or 3-D mesh method
was used. The row contains CT intensity, **not ink probability**.

## Hypothesis and exact input

Exploratory hypothesis: an open near-turn around a visible central winding can
be followed by a deterministic, manually guided intensity-ridge optimizer
without losing physical layer identity. This is a falsifiable local pilot,
not a formally preregistered acceptance experiment. The research-experiments
skill's preservation principle was applied; no commits or independent-review
claims are made.

- Dataset: **Scroll 1 / PHercParis4**, `54keV_7.91um_Scroll1A.zarr`, 54 keV.
- Exact remote array:
  `https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4`.
- Exact local source: `outputs/orientation/ct_level4_zyx.npy`, array shape
  `(128,493,506)`, stored ZYX, dtype uint8.
- Selected slice: **local index 16, global level-4 Z=144**, nominal native
  Z=2304. Read exactly `volume[16,:,:]`, shape `(493,506)`; no transpose.
- Level-4 spacing: **126.56 µm/pixel** in both XY directions, from 7.91 µm
  native spacing ×16. This is published downsampled CT, not a native Z=2304
  plane or a claim of exact cross-level voxel-center registration.
- Uncompressed source-slice bytes SHA-256:
  `93f77b3f25e3d33736dfe838b6a00b92766df326efa43042c9c359bc19cdcef2`.
- Full source-object URLs, byte sizes, hashes and acquisition provenance remain
  in `outputs/orientation/manifest.json`; see `docs/p0_orientation.md` for the
  empirical orientation evidence. No new CT was downloaded for P0.1.

Selection rationale: this was the already inspected representative XY plane
matching the scout's axial neighborhood. Its central nested winding geometry
is visible, and it avoids choosing a new slice after seeing tracing outcomes.
The manual guide was selected from this slice before optimization; its shape
is a strong human-supplied prior, not an independently inferred spiral.

## Preprocessing and deterministic method

The source PNG preserves uint8 values exactly. No local histogram rescaling,
thresholding, denoising of the sampled signal, or learned processing was used.
For the tracing cost only, compute a difference of Gaussian-smoothed images:
`R = Gaussian(I, sigma=0.6 px) - Gaussian(I, sigma=2.0 px)` using SciPy defaults
(reflecting boundaries, truncate=4). This emphasizes bright narrow features;
it is a structural heuristic, not a material or ink classifier.

Exact implementation: `scripts/p0_1_single_slice.py`.
`outputs/p0_1/trace_plan.json` was saved before running the optimizer and records
all parameters and the planned sensitivity variants. The 14 manual XY guide
vertices, in global level-4 pixels, are:

```text
(246,111), (258,120), (271,131), (281,141), (281,147),
(270,155), (245,159), (210,158), (183,154), (175,150),
(187,141), (211,132), (229,122), (241,114)
```

1. Parameterize the guide vertices by their cumulative chord lengths. Fit a
   natural cubic spline in XY and evaluate 477 stations with guide-parameter
   spacing at most 0.5 px. Calculate normals from its derivative.
2. At each station form 25 candidate points at normal offsets −3 to +3 px in
   0.25 px steps. Bilinearly interpolate the ridge response at those points.
3. Dynamic programming minimizes the summed unary costs
   `-R/255 + 0.001*offset²` plus transition cost
   `0.01*(offset_previous-offset_current)²`. Forbid offset changes above 0.5 px
   or Euclidean station-to-station moves above 1.25 px. Initial and final
   offsets are free; ties follow NumPy's first-index `argmin`.
4. Backtrack the minimum-cost **open** path. Do not impose a closed loop, fill
   gaps, join endpoints, repair contacts, or claim mathematical continuity
   proves material continuity.
5. Sum Euclidean segment lengths along the resulting polyline. Interpolate XY
   linearly at equally spaced arc-length stations, including both endpoints.

The guide's upper/left ambiguity was known visually; the optimizer was still
run to test whether this simple bounded approach could maintain identity.
No baseline parameter or guide was retuned after seeing its result.

## Arc length and sampling

| Quantity | Measured result |
| --- | ---: |
| Optimizer path vertices | 477 |
| Arc-length-resampled trace points | **556** |
| Total polyline arc length | **277.1273271765 level-4 pixels** |
| Physical length | **35.0732345275 mm** |
| Uniform arc-length interval | 0.4993285174 px |
| Maximum optimizer step | 0.8396006440 px |
| Start/end separation (not joined) | 4.6098760208 px |
| Angular sweep about descriptive center (240,140) | 350.683648° |
| Stations on the corridor boundary | 4.612159% |
| CT sample min / mean / max | 24.387204 / 98.876390 / 134.226512 |
| Output row width × height | **556 × 1**, grayscale uint8 |

The sweep describes an almost complete turn of the candidate; the center is
not an estimated umbilicus, and the sweep is not proof of one physical winding.
Length is that of the unsmoothed optimized polyline and includes its small
zigzags. The physical length inherits the catalog spacing and coarse resolution.

Sample **the original unsmoothed source slice**, not the ridge-response image,
at the arc-length XY coordinates using `scipy.ndimage.map_coordinates`,
`order=1`, `prefilter=False`, `mode='nearest'`. Coordinates stay strictly inside
the image, so boundary mode has no effect. Bilinear interpolation uses the four
neighboring pixel centers. Save float64 samples in `sampling_values.npy` and
the CSV; create the row PNG with `clip(rint(value),0,255).astype(uint8)`.
There is no cross-path averaging, signal smoothing, or intensity normalization.

CSV columns include point index, `s_pixels`, `s_mm`, global level-4 XY/Z,
nominal native XYZ (scale ×16), and floating CT intensity. Increasing columns
of the row follow the CSV order. The 48-pixel-high preview merely repeats the
same row for visibility; **it is not a multi-slice stack**.

## Physical audit, gaps, and negative evidence

Inspect `trace_audit.png`, comparing the unoccluded source/guide panel with the
baseline and sensitivity overlays. The right-hand arc follows a visibly curved
band for part of its course. The lower section, approximately x=185–245,
y=153–160, contains diagonal moves between nearly parallel bands and low-contrast
contacts; physical layer jumps are suspected. The upper-left/apex region,
approximately x=220–250, y=110–125, does not uniquely constrain one sheet.
These are visual audit regions, not measured ground-truth jump locations.

There are no omitted path stations or explicit gap-fill operations, but the
continuity constraint forces the optimizer through ambiguous evidence. Thus a
continuous CSV can conceal unresolved physical gaps or contacts. The intentional
open endpoint gap is retained. A bright ridge can also be a phase fringe or
merged structure rather than one papyrus sheet.

Four predeclared variants were run and preserved, with no best-variant selection:

| Variant | Median displacement (px) | 95th percentile (px) | Maximum (px) | CT RMSE at matched guide stations |
| --- | ---: | ---: | ---: | ---: |
| Gaussian sigma 0.4 | 0 | 0.50 | 1.00 | 3.543080 |
| Gaussian sigma 0.8 | 0 | 0.75 | 2.00 | 6.212085 |
| Guide shifted −0.5 px along normal | 0 | 3.55 | **6.50** | 7.545548 |
| Guide shifted +0.5 px along normal | 0 | 0.50 | 2.25 | 3.045723 |

These distances compare corresponding guide stations, not nearest points on
curves. The RMSE is not an equal-arc-length row comparison. The large local
displacement under a half-pixel guide shift exposes competing routes in the
upper arc. Reproducibility of a fixed path does not imply anatomical stability.

## Artifacts, commands, and verification

Required artifacts are all in `outputs/p0_1/`:

- `source_slice.png`: exact 506 × 493 source pixels.
- `traced_slice.png`: full-slice candidate overlay; green start, blue end.
- `unwrapped_row.png`: exact 556 × 1 sampled candidate row.
- `trace_coordinates.csv`: arc-length coordinates and floating sample values.

Additional evidence retained: `core_inspection.png`, `core_grid.png`,
`trace_audit.png`, `intensity_profile.png`, `unwrapped_row_preview.png`,
`ridge_response.npy`, `optimizer_path_xy.npy`, all four variant path NPYs,
`sampling_values.npy`, `trace_plan.json`, `metrics.json`,
`preservation_before.json`, `verification.json`, and `SHA256SUMS`.

Executed from the repository root:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p01 .venv/bin/python -B scripts/p0_1_single_slice.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p01 .venv/bin/python -B scripts/verify_p0_1.py
```

The tracing command refuses to overwrite completed metrics. Reproduce the
numerical result offline with the verifier, which reruns the deterministic
trace without replacing candidate evidence. Initial inspection PNGs were
generated from `volume[16]`; `core_inspection.png` is crop `[70:220,100:390]`
enlarged 4× with nearest-neighbor interpolation, and `core_grid.png` is a
coordinate-grid view of x=120–370, y=75–205. These plots were for guide selection.

Verification checks exact deterministic path replay, saved CSV coordinates,
arc length, NPY samples, raw source pixels, and rounded row pixels. An independent
four-corner bilinear formula agreed with SciPy to maximum absolute error **0**.
Every pre-existing file in the preservation inventory is hash-checked,
including orientation data/evidence and README/AGENTS. No prior result is edited.

## Conclusion and next recommended experiment

**Geometric/sampling machinery works; physical P0.1 success is not demonstrated.**
The baseline row is numerically stable under exact replay, but plausible
single-winding identity and stability to guide ambiguity fail this inspection.
No claim of text or ink is made, and no attractive overlay is used as a pass.

Recommended next experiment, **not run**: inspect a bounded higher-resolution
XY crop of the same axial neighborhood, focusing on the flagged contacts, before
attempting another single-slice trace. Higher resolution may resolve competing
bands; if it does not, preserve that negative result as well. Do not proceed
to multi-slice stacking on this candidate. Stop at P0.1.
