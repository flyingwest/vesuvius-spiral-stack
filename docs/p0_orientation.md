# P0 data acquisition and orientation verification

## Conclusion and scope

**Orientation established for this region: use stored Z (array axis 0) as the
approximate longitudinal `z` in Direct Spiral-Stack. The transverse plane is
XY, accessed as `ct[k, :, :]`.** No transpose or axis reversal is required.

This is an exploratory acquisition/orientation check, not P0 reconstruction,
not a preregistered tracing experiment, and not a P0 PASS. No spiral tracing,
surface fitting, mesh processing, unwrapping, OCR, or ink inference was done.
`README.md`, `AGENTS.md`, and the earlier scout report remain unchanged.

## Hypothesis, selection, and footprint

Hypothesis before downloading voxels: the catalog's stored Z direction is
approximately longitudinal; XY should display the scroll's rolled cross-section.
The acquisition plan was saved before voxel GET requests in
[`acquisition_plan.json`](../outputs/orientation/acquisition_plan.json).

The source is the Scroll 1 / PHercParis4 recommendation in
[`data_scout.md`](data_scout.md). A peripheral 128-cubed crop cannot show the
whole transverse geometry, so the first useful input here is a coarse slab
with full XY coverage. The initial level-3 (8×) HEAD survey required 95,994,379
bytes for available chunks; one corner key `3/2/7/7` returned 404. No level-3
voxel chunks were downloaded. A smaller level-4 (16×) slab had all 16 required
chunks available and was selected before viewing voxels. This is a bounded,
useful first survey, not a proof of the absolute minimum possible byte count.

Downloaded CT chunks total **22,039,544 bytes**; three metadata files add
**3,818 bytes**, giving **22,043,362 bytes (21.0222 MiB)** of retained HTTP GET
payload. HEAD responses and exploratory metadata requests are additional small
network overhead, not included in that retained-file total. HTTP headers and
transport overhead are not counted. No entire scroll or paired NRRD cube was
downloaded. The paired cube was unnecessary once orientation was established.

The decoded array has 31,930,624 bytes; the saved `.npy` has 31,930,752 bytes.
Compressed source objects are retained separately for offline reproducibility.
Before this report, the orientation output directory including images,
manifests, metrics, and checksums occupied 34,262,116 logical file bytes;
filesystem allocation overhead is not included.

## Dataset identity and coordinates

| Field | Value |
| --- | --- |
| Scroll / object | Scroll 1 / PHercParis4 |
| Scan/store ID | `54keV_7.91um_Scroll1A.zarr` |
| Energy | 54 keV, from the scan catalog |
| Source representation | Published standardized uint8 Zarr v2, not original scanner uint16 |
| Array dimension order | Z, Y, X; C order |
| Full level-0 shape | `(14376, 7888, 8096)` |
| Native voxel spacing | Isotropic 7.91 µm, catalog value |
| Acquired pyramid level | 4; scale `(16,16,16)` relative to level 0 |
| Level-4 full-array shape | `(899,493,506)` |
| Acquired bounds, half-open, level-4 ZYX | `[128:256, 0:493, 0:506]` |
| Corresponding native index bounds | `[2048:4096, 0:7888, 0:8096]` |
| Saved array shape | `(128,493,506)` |
| Saved dtype | `uint8` |
| Acquired voxel spacing | Isotropic **126.56 µm** |
| Acquired field of view, ZYX | `(16.19968,62.39408,64.03936)` mm |
| Chunk dimensions | `(128,128,128)` |
| Compression | Blosc / Zstd, level 3, shuffle 0 |

The scan URL comes from Villa's
[`scrolls.yaml:5–11`](../external/villa/vesuvius/src/vesuvius/install/configs/scrolls.yaml).
Villa revision inspected: `fb8c2c4c2705746eefb42e9bb359560906f391e2`.
The remote [root metadata](https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/.zattrs)
declares ZYX and pyramid factors, not physical units. Physical spacing above
combines that factor with the catalog spacing. Both
[native metadata](https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/0/.zarray)
and [level-4 metadata](https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4/.zarray)
were saved verbatim. Exact subvoxel center registration across pyramid levels
is not independently calibrated; native bounds are scale-based index bounds.

For saved local index `(k,j,i)`, the level-4 global index is `(128+k,j,i)`.
Images display increasing horizontal coordinate rightward and vertical coordinate
downward, with equal physical aspect ratio. Plot axes are scan-coordinate mm,
not anatomical labels. The XZ and YZ displays retain the true short axial extent;
they are not stretched to squares.

## Exact downloaded files

Every key below is appended literally to this base URL:

```text
https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/
```

Every local filename is `data/raw/p0_orientation/scroll1.zarr/` plus that key.
The 19-file inventory, exact per-file URL, byte count, SHA-256, and response
headers including server dates are recorded in
[`manifest.json`](../outputs/orientation/manifest.json).

| Key | Bytes |
| --- | ---: |
| `.zattrs` | 3018 |
| `0/.zarray` | 402 |
| `4/.zarray` | 398 |
| `4/1/0/0` | 1633272 |
| `4/1/0/1` | 1729614 |
| `4/1/0/2` | 1806018 |
| `4/1/0/3` | 1092425 |
| `4/1/1/0` | 1789403 |
| `4/1/1/1` | 1823708 |
| `4/1/1/2` | 1835501 |
| `4/1/1/3` | 1768688 |
| `4/1/2/0` | 1299239 |
| `4/1/2/1` | 1708927 |
| `4/1/2/2` | 1624696 |
| `4/1/2/3` | 1000541 |
| `4/1/3/0` | 678836 |
| `4/1/3/1` | 1143457 |
| `4/1/3/2` | 707962 |
| `4/1/3/3` | 397257 |

**The local Zarr is a partial cache, not a full volume.** Only the specified
bounds have been acquired. Opening unacquired regions could return Zarr fill
values; those would not be observations. All chunks intersecting the actual
orientation slab returned HTTP 200 and were decoded successfully.

## Commands and environment

Executed from `/home/wpotech/deployments/vesuvius-spiral-stack`:

```bash
.venv/bin/python -B scripts/p0_orientation.py acquire
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-orientation .venv/bin/python -B scripts/p0_orientation.py render
.venv/bin/python -B scripts/p0_orientation.py verify
```

[`scripts/p0_orientation.py`](../scripts/p0_orientation.py) is the exact
acquisition/inspection implementation. It checks all object sizes using HEAD,
enforces a 24 MiB payload cap, fetches only the explicit inventory with four
workers, checks HEAD/GET length agreement, and records byte hashes. It refuses
to overwrite an existing acquisition. On this checkout, use only `render` and
`verify` to reproduce the results offline; acquisition is for a fresh data path.

The existing environment was used without installation: Python 3.12.3,
NumPy 2.5.3, Zarr 3.4.0, numcodecs 0.17.0, Matplotlib 3.11.2, Pillow 12.3.0.
No Villa package import was necessary: bounded HTTP object downloads followed
by a standard local Zarr read avoid its Python-3.14 installation requirement.
No intensity normalization, filtering, learned processing, or resampling was
applied locally. The downsampling is the publisher's existing pyramid.

## Representative orthogonal views and observations

All files below are under `outputs/orientation/`. Each annotated PNG has a
matching `_raw.png` containing the exact uint8 slice pixels. Display limits
are fixed at 0–255 for every image, using grayscale and nearest-neighbor display.

| Filename | Exact saved-array view | Observation |
| --- | --- | --- |
| `XY_z0144.png` | `ct[16,:,:]` | Full flattened/oval scroll outline with nested winding bands around a central region; native nominal z=2304, matching the scout neighborhood. |
| `XY_z0192.png` | `ct[64,:,:]` | Same enclosing rolled geometry at nominal native z=3072; strongly deformed interior and open gaps. |
| `XY_z0240.png` | `ct[112,:,:]` | Rolled cross-section persists at nominal native z=3840 despite substantial internal changes. |
| `XZ_y0246.png` | `ct[:,246,:]` | Longitudinal section through folded layers; local curls and bends, not the enclosing transverse roll. |
| `YZ_x0253.png` | `ct[:,:,253]` | Longitudinal layer courses and changing gaps; no enclosing whole-roll cross-section. |

These images were inspected, not inferred solely from metadata. The repeated
nested winding geometry in XY, and its persistence across different stored Z
positions, establish the useful transverse orientation. Smooth surrounding
support-like material is not interpreted as papyrus. XZ/YZ also contain curves:
curvature alone does not establish a winding plane. This scan is substantially
deformed, not an ideal circular roll with straight axial layers.

## Quantitative checks

[`metrics.json`](../outputs/orientation/metrics.json) records exact statistics.
Nonzero voxel fraction is 0.928550; intensity percentiles 0/1/50/99/100 are
0/0/62/142/248. These are data-integrity descriptors, not papyrus or ink labels.

A fixed central cuboid, global level-4 ZYX `[128:256,80:420,80:420]`, was used
for a descriptive adjacent-voxel continuity comparison:

| Direction | Mean absolute intensity difference | Pearson correlation | Pairs |
| --- | ---: | ---: | ---: |
| Z | 16.730864 | 0.831424 | 14681200 |
| Y | 20.693668 | 0.746178 | 14753280 |
| X | 17.020012 | 0.824193 | 14753280 |

Z is the smoothest of these three directions, but X is close. This metric is
only corroborative: spatial correlation, support material, folds, and pyramid
smoothing confound it. It is not an independent statistical axis classifier or
a measured tilt angle. The primary evidence is the orthogonal image geometry,
with metadata fixing the array-axis names.

The verifier passed source size/hash checks, equality between a fresh offline
Zarr read and the saved NPY, pixel equality for all five raw orthogonal views,
and unchanged README/AGENTS hashes. Output checksums are in
[`SHA256SUMS`](../outputs/orientation/SHA256SUMS). The acquisition script's hash
and dependency versions are in the manifest.

## Uncertainties and stopping point

- Z is approximately longitudinal in this slab, not a fitted local umbilicus
  tangent or proof of a straight axis throughout the full scroll.
- 126.56 µm sampling resolves gross winding structure here, but can merge fine
  sheets and gaps. No individual winding identity or traceability is established.
- Anatomical top/bottom, handedness, recto/verso, and reading direction were not
  empirically established. Preserve stored axis signs; do not infer these from
  the plotted winding appearance.
- No paired labels were acquired, so CT-to-label registration remains untested.
- This standardized intensity is not raw calibrated density or an ink signal.
- The 16× slab extends beyond the scout's narrow z range to use a single axial
  chunk layer economically. It does not expand acquisition to the full scroll.

Conclusion: orientation-only goal met; **stop here**. The next recommended
experiment, not executed, is a separately authorized native-resolution local
data/label inspection before choosing any direct spiral-tracing experiment.
