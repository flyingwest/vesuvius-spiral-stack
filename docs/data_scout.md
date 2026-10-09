# Real CT data scout

Scouted 2026-10-06. Villa checkout:
`fb8c2c4c2705746eefb42e9bb359560906f391e2`.
Read the repository's `README.md` and `AGENTS.md` completely, and Villa's
`AGENTS.md` and `SCROLL_ANATOMY.md`. This is a source/metadata investigation,
not a reconstruction experiment. No software was installed, no complete CT
chunk or label volume was downloaded, and no spiral-stack algorithm was implemented.
Remote checks were limited to text/metadata, HEAD requests, and two 1 KiB
NRRD header ranges (HTTP 206).

An existing untracked version of this report was present at session start.
It was preserved and checked against source and live endpoints. Recheck evidence
is saved in [`data_scout_evidence.json`](data_scout_evidence.json): request
URLs/methods, HTTP statuses, response headers, metadata bodies, and NRRD header
prefixes. The mask prefix is not its complete NRRD header.

## Recommendation and frozen hypothesis

Use **Scroll 1 / PHercParis4, 54 keV, 7.91 µm**, through the canonical
standardized Zarr. Begin with the native-resolution, chunk-aligned crop
`[z=2304:2432, y=2560:2688, x=4864:4992]` (half-open bounds).
It is wholly inside the published manually annotated cube's bounding box
`02256_02512_04816`. One CT chunk costs **719,102 bytes**, verified by HEAD.[4]
For a self-contained comparison input, download that cube's paired CT and
harmonized instance labels instead: **15,972,223 bytes combined**.[5][6]

This is the recommended **local data/correspondence preflight**, not a substitute
for P0. A 128³ or 256³ peripheral crop does not contain the full transverse
spiral. To investigate P0 itself, retain full transverse context using the
same scan and axial neighborhood; see the staged slab plan below.

The frozen method remains:

```text
real transverse CT slice -> trace papyrus spiral -> arc length s
                        -> sample grayscale / ink-related signal -> one 1-D row
stack neighboring axial rows -> I(s,z)
```

PASS requires a spatially coherent 2-D reconstruction; readable letters are
not required. Stretch, shear, or curvature can be acceptable. Stop/reconsider
if reliable winding correspondence requires effectively reconstructing the
full 3-D surface manifold. Existing meshes and labels are comparison evidence,
not the reconstruction method. Raw CT darkness is not an ink label. No OCR,
Greek-language priors, language-model repair, or synthetic CT is introduced.

## 1. Simplest supported access

Villa's documented interface is `vesuvius.Volume`, with NumPy-like
`volume[z_slice, y_slice, x_slice, level]` indexing. It reads intersected remote
Zarr chunks on demand; do not recursively download the store. A direct array
URL also works and avoids ambiguity about pyramid levels.

Sources in the inspected checkout:

- [`vesuvius/docs/accessing_data.md`](../external/villa/vesuvius/docs/accessing_data.md):
  `Working with Volume`, indexing, normalization, and `Cube` examples.
- [`vesuvius/src/vesuvius/data/volume.py`](../external/villa/vesuvius/src/vesuvius/data/volume.py):
  constructor at line 146, direct Zarr initialization at 393, indexing at
  approximately 850–980, and `Cube` loading at 1194–1375.
- [`vesuvius/src/vesuvius/data/utils.py`](../external/villa/vesuvius/src/vesuvius/data/utils.py):
  `open_zarr` at 186; HTTP read-only handling at 276; optional compressed-chunk
  memory cache at 329. This is not a persistent disk cache.
- [`scrolls.yaml`](../external/villa/vesuvius/src/vesuvius/install/configs/scrolls.yaml):
  canonical Scroll 1 URL at line 11.

### Commands for a later acquisition, not executed

From the repository root, minimal supported package setup for this checkout is:

```bash
python3.14 -m venv .venv-data
.venv-data/bin/python -m pip install './external/villa/vesuvius[volume-only]'
# Review the package's data terms; run the following only if accepting them:
.venv-data/bin/vesuvius.accept_terms --yes
```

[`pyproject.toml`](../external/villa/vesuvius/pyproject.toml) requires Python
`>=3.14,<3.15` in this checkout. The core/volume-only installation does not
require the models, GUI, or rendering extras. Dependency installation size is
not measured and is separate from the data budget. Do not assume an older
PyPI release exactly matches this checkout. The terms command writes an
agreement file; it was not run during scouting.

No installation is necessary for the paired-NRRD download commands below;
they only require `curl`. Reading those files later requires a NRRD reader.
Unlike `Volume`, `Cube` downloads both complete NRRDs eagerly. Its
`cache=False` branch still writes temporary files with `delete=False`
(`volume.py:1330`); it is not a zero-disk streaming option. Prefer explicitly
named downloads for a bounded, auditable paired-cube acquisition.

Exact native-resolution crop acquisition:

```bash
.venv-data/bin/python - <<'PY'
from pathlib import Path
import hashlib
import json
import numpy as np
from vesuvius import Volume

url = ('https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/'
       'volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/0')
v = Volume(type='zarr', path=url, normalization_scheme='none',
           return_as_type='none', return_as_tensor=False, cache=False)
a = v[2304:2432, 2560:2688, 4864:4992]
assert a.shape == (128, 128, 128) and a.dtype == np.uint8
out = Path('data/scout/scroll1_native_chunk')
out.mkdir(parents=True, exist_ok=True)
np.save(out / 'ct_zyx.npy', a)
(out / 'manifest.json').write_text(json.dumps({
    'url': url, 'axes': 'zyx', 'level': 0,
    'bounds_half_open_zyx': [[2304, 2432], [2560, 2688], [4864, 4992]],
    'voxel_size_um': 7.91, 'shape': list(a.shape), 'dtype': str(a.dtype),
    'normalization': 'none',
    'villa_commit': 'fb8c2c4c2705746eefb42e9bb359560906f391e2',
    'array_sha256': hashlib.sha256(a.tobytes(order='C')).hexdigest(),
}, indent=2) + '\n')
PY
```

This snippet only acquires data; it contains no tracing or reconstruction.
It was checked against the source API but not executed against CT data.

## 2. Coordinates, axes, and intensity

Live [root metadata](https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/.zattrs)
explicitly declares axes **Z, Y, X**, with pyramid scales 1, 2, 4, 8, 16, 32
and no translation. These scales are relative voxel factors, not micrometres;
the physical base spacing comes from the scan catalog (7.91 µm).[1]

Live [level-0 metadata](https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/0/.zarray):

```json
{"shape":[14376,7888,8096],"chunks":[128,128,128],"dtype":"|u1",
 "order":"C","dimension_separator":"/","fill_value":0,"zarr_format":2,
 "compressor":{"id":"blosc","cname":"zstd","clevel":3,"shuffle":0,"blocksize":0}}
```

The metadata above was re-fetched successfully during this scout.[2]

- A transverse plane is `a[z, :, :]`; rows are y, columns are x. Use image
  display origin `upper` to show increasing x rightward and increasing y
  downward, without an accidental transpose or vertical flip.
- z is the scan-stack direction, approximately the scroll's longitudinal
  direction; the actual umbilicus can bend/drift. Do not equate it to a fitted
  straight physical axis.
- Villa's [FAQ, Scroll 1 orientation](../external/villa/scrollprize.org/docs/09_faq.md)
  (lines 420–435) describes the original TIFF stack: increasing slice number
  points toward the viewer, from the middle toward the top of the scanned
  scroll, with a counterclockwise central spiral. This is a historical
  physical-orientation reference, not an independently measured orientation
  of this Zarr. Its metadata confirms ZYX and positive pyramid scales but
  does not establish anatomical handedness. Preserve stored axes and verify
  landmarks against the original/reference images before a physical claim.
- Mesh vertices use **XYZ**, so sample a volume with `(vertex.z, vertex.y,
  vertex.x)`. For a local crop subtract its global ZYX origin after reordering.
  Pyramid coordinates must also account for the scale; avoid assuming exact
  subvoxel center registration from a factor alone.
- The selected standardized CT is **uint8**, not original uint16 TIFF data.
  `normalization_scheme='none'` preserves the published values, not the original
  scanner dynamic range. Do not infer ink or calibrated density from intensity.

The selected NRRD headers declare `sizes: 256 256 256`, identity space
directions, and `space origin: (2256,2512,4816)`. The CT type is uint8; the mask
type is little-endian uint16; both are gzip encoded. Although the header also
says `space: left-posterior-superior`, its origin matches the documented
**ZYX cube filename**, not a justified anatomical frame for the scroll.
Villa's `Cube` uses `nrrd.read` with its default ordering and indexes the
returned array as ZYX. Use `nrrd.read(..., index_order='F')` explicitly to match
that implementation; do not independently reverse axes because of the LPS
label. Confirm alignment visually/numerically before treating it as verified.

## 3. Dataset choice and region size

| Candidate | Suitability for this task |
| --- | --- |
| Canonical Scroll 1, 54 keV / 7.91 µm Zarr | Preferred: full spiral context is available lazily, native-resolution local labels exist, and historical unwrapping references are plentiful. Villa's VC tutorial calls Scroll 1 the easiest of its listed scrolls to segment. |
| Scroll 1 harmonized instance cube | Best small paired preflight: actual CT plus manual sheet-instance labels, approximately 2.025 mm per side. Cannot establish full-spiral P0 alone. |
| Grand Prize banner region | Useful later mesh comparison, but the volume masks non-GP regions. Missing spiral context makes it a poor primary input for the frozen P0. |
| Scrolls 2–5 | Canonical URLs exist in `scrolls.yaml`, but no equally direct paired local reference was established in this scout. No claim that these scans are intrinsically unusable. |
| New spiral-input/surface-label bundles | Potential later references; unnecessary initial footprint/registration work. No spiral-fitting pipeline is adopted. |

Sources: [VC tutorial](../external/villa/scrollprize.org/docs/06_tutorial_VC.md)
line 175; [local dataset catalog](../external/villa/scrollprize.org/docs/02_data_datasets.md)
lines 140–176; [current public catalog](https://scrollprize.org/data_datasets).
The current catalog marks the older cube and Grand Prize releases as archived;
the specific objects above remain reachable as of the scout date.[10]

The local `spiral-input-PHercParis4-README.md` lists approximately 49.6 GB of
annotations and inputs, including manually verified and unverified patches;
those two categories must not be conflated. The local `surface-labels-README.md`
still has a TODO for the actual file layout. These are concrete reasons to
prefer the directly addressable paired cube for the first acquisition, not
claims that the newer releases lack useful reference data. Both READMEs are in
`external/villa/scrollprize.org/static/data/datasets/`.

Suggested staged inputs, all from the same canonical scan:

1. **Local preflight:** level 0, ZYX `[2304:2432,2560:2688,4864:4992]`.
   128 consecutive native slices over approximately 1.012 mm. Verify actual
   papyrus, local labels, axis ordering, and continuity before extending scope.
2. **Full transverse survey:** level 2, ZYX `[576:608,0:1972,0:2024]`.
   This spans the same axial neighborhood at 31.64 µm spacing. The live
   [level-2 metadata](https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/2/.zarray)
   reports shape `[3594,1972,2024]`, uint8, and 128³ chunks.[3] Downsampling may
   merge windings; this is a survey, not evidence that native data fail P0.
3. **Native full-context candidate:** level 0, ZYX
   `[2304:2336,0:7888,0:8096]`. Preserve all transverse windings and 32
   consecutive slices, then select tighter x/y bounds only after inspecting
   the survey. This is a proposal, not a confirmed tractable spiral region.

## 4. Trusted comparison overlap

**Yes, at the bounding-box level: a published manually annotated segmentation
covers the chosen crop's nominal region.** Nonzero label coverage and exact
registration have not been measured. The parent cube bounds are ZYX
`[2256:2512,2512:2768,4816:5072]`; containment is exact. The server's
[instance-label README](https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumetric-instance-labels/README.txt)
describes manual sheet-instance labels and recommends the harmonized variant,
whose IDs are consistent across neighboring cubes. Small labeled components
were removed during harmonization. These are human reference annotations,
not infallible ground truth, and an instance ID is not automatically a winding
number or recto/verso identity.[7]

[Verified directory](https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumetric-instance-labels/instance-labels-harmonized/02256_02512_04816/)
contains both the CT and harmonized mask. Villa's
[`cubes.yaml`](../external/villa/vesuvius/src/vesuvius/install/configs/cubes.yaml)
line 11 and `Cube` example select the **unharmonized** directory. To obtain the
recommended harmonized pair, use these exact commands later:

```bash
mkdir -p data/scout/02256_02512_04816
cube_url='https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumetric-instance-labels/instance-labels-harmonized/02256_02512_04816'
curl --fail --location "$cube_url/02256_02512_04816_volume.nrrd" -o data/scout/02256_02512_04816/volume.nrrd
curl --fail --location "$cube_url/02256_02512_04816_mask.nrrd" -o data/scout/02256_02512_04816/mask.nrrd
sha256sum data/scout/02256_02512_04816/*.nrrd
```

The selected 128³ crop corresponds to local cube slices
`[48:176,48:176,48:176]`. The paired NRRD CT is the safest immediate reference
for the mask. Exact intensity equality and voxel-center registration between
that CT and the standardized Zarr remain **unverified**; dimensions/origins
and catalog identity establish nominal overlap, not a measured registration.
Before a Zarr-versus-label comparison, inspect overlays and measure any offset.

**A flattened unwrap at these exact coordinates was not verified.** The
[Grand Prize README](https://dl.ash2txt.org/datasets/grand-prize-banner-region/README.txt)
documents semi-manual OBJ meshes (XYZ), associated scan data (ZYX), and
uncropped TIFF dimensions that preserve alignment. The full mesh archive is
listed as 288 MB. No mesh bounds were downloaded, so neither its intersection
with this particular cube nor pixel-to-volume mapping for a flattened image
is asserted. Establish those bounds before acquiring an unwrap for comparison.[8]
ML medial-surface predictions are not a replacement for the manual reference.

One important negative finding: that README describes `gp_volume.zarr` as
multiresolution OME-Zarr, but the live root has a **single-array `.zarray`**,
shape `[14376,7888,8096]`, uint8/128³ chunks, and root `.zattrs` returned 404.
The exact array root is
`https://dl.ash2txt.org/datasets/grand-prize-banner-region/volumes/gp_volume.zarr`.[9]
Its numbered directories are chunk coordinates, not pyramid levels; do not
open `/0` as a resolution level. Requests omitting `/volumes/` also returned
404; those failed exploratory URLs are preserved in the evidence log.
The README's historical `foundation/scroll2zarr`
source link also names a directory absent from the inspected checkout.

## 5. Minimum footprint and reproducible checks

Data bytes below exclude HTTP headers, retries, Python/runtime dependencies,
and filesystem overhead. MiB means 2²⁰ bytes. Compressed transfer size must not
be confused with decoded chunk working-set size.

| Acquisition | Network payload | Local data footprint |
| --- | --- | --- |
| Selected 128³ native crop | One object: **719,102 B**, plus a few KiB of metadata/probes | 2 MiB decoded; `.npy` adds a small header. No persistent remote cache needed. |
| Paired 256³ NRRD CT | **14,799,078 B** | Same compressed file; 16 MiB decoded |
| Paired harmonized NRRD labels | **1,173,145 B** | Same compressed file; 32 MiB decoded |
| Both NRRDs | **15,972,223 B ≈ 15.23 MiB** | Approximately 15.23 MiB disk if retaining only originals; 48 MiB decoded arrays |
| One Zarr crop plus both NRRDs | **16,691,325 B**, plus metadata | Approximately 17.23 MiB disk for NRRDs plus decoded crop |
| Level-2 full transverse survey, 32 planes | 256 compressed chunk objects; total transfer **not measured** | About 121.8 MiB output; 512 MiB aggregate decoded chunk payload |
| Level-0 full transverse slab, 32 planes | 3,968 compressed chunk objects; total transfer **not measured** | About 1.90 GiB output; 7.75 GiB aggregate decoded chunk payload |

For a crop `[a:b]`, the intersected chunk count on an axis is
`floor((b-1)/128) - floor(a/128) + 1`; multiply across axes. Even one plane
requires whole compressed 3-D chunks. Aggregate decoded payload is not
necessarily peak RAM; sequential reads can limit residency. Compression varies
spatially: do not extrapolate the selected 719,102-byte object to a whole slab.
Obtain HEAD sizes for a proposed slab's exact chunk keys before downloading it.

Minimal CT-only acquisition is therefore under 1 MiB of network data for the
selected crop. Minimal convenient **paired comparison** is the 15.23 MiB NRRD
pair. A full-context native P0 input is substantially larger; no claim is made
that the one-chunk preflight suffices for P0 PASS. Avoid the entire 77 GB GP
Zarr, 389.9 GB TIFF archive, 1.4 GB label archive, or recursive whole-scroll copy.

Reproduce the small metadata/HEAD checks without fetching CT chunks:

```bash
store_url='https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr'
cube_url='https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumetric-instance-labels/instance-labels-harmonized/02256_02512_04816'
curl --fail --silent --show-error "$store_url/.zattrs"
curl --fail --silent --show-error "$store_url/0/.zarray"
curl --fail --silent --show-error "$store_url/2/.zarray"
curl --fail --head "$store_url/0/18/20/38"
curl --fail --head "$cube_url/02256_02512_04816_volume.nrrd"
curl --fail --head "$cube_url/02256_02512_04816_mask.nrrd"
curl --fail --range 0-1023 --max-filesize 1024 "$cube_url/02256_02512_04816_volume.nrrd" -o /home/wpotech/.hermes/cache/scratch/scroll1-volume-header.bin
```

Scouting used equivalent bounded Python `urllib.request` requests. The chunk's
HEAD returned HTTP 200, ETag `"673c9205-af8fe"`, and last-modified
`2024-11-19 13:26:29 GMT`. These are provenance markers, not content hashes.
The NRRD range requests returned HTTP 206 and 1,024 bytes each. No decoded CT
statistics, image assessment, registration metrics, or P0 results exist yet.

Verification also checked the crop containment/chunk indices and recomputed
decoded sizes from metadata. The acquisition Python snippet was syntax-checked,
not run; neither dependency availability nor voxel decoding was tested.
`README.md` and `AGENTS.md` were left unchanged. The frozen `AGENTS.md` SHA-256 is
`6859793a0a690642dc93fc34d4a69cd712a13703592583d6320b081931851741`.

## Next experiment record (proposed, not run)

- **Hypothesis:** native CT within the specified cube contains locally separable
  papyrus paths whose identity persists across neighboring transverse planes.
- **Input:** paired cube above, plus optionally the exact one-chunk Zarr crop;
  preserve origins, byte hashes, dtype, and spacing. Do not normalize per slice.
- **Parameters:** 128³ local crop at the stated global bounds; stored ZYX
  orientation; no learned processing or tracing implementation in this scout.
- **Artifacts:** source manifest, acquired CT/labels, axial contact sheet,
  label overlays, and a written acquisition/registration report.
- **Metrics:** shape/dtype, nonzero fraction and intensity quantiles, reference
  coverage, measured registration offset, and subsequently adjacent-slice
  displacement and winding-identity switches if tracing is authorized.
- **Qualitative observations:** not available until voxels are inspected; record
  contact/fusion, tears, ambiguous winding identity, and downsampling losses.
- **Conclusion:** source inspection supports a small real-data preflight with
  nominal manual-reference overlap. It does not establish traceability or P0.
- **Next:** acquire the paired cube, verify coordinates, then budget and inspect
  the full transverse survey. Only after that select a native slab for the
  frozen direct spiral-stack experiment. Preserve failures as well as successes.

## Sources

[1] https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/.zattrs
[2] https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/0/.zarray
[3] https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/2/.zarray
[4] https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/0/18/20/38
[5] https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumetric-instance-labels/instance-labels-harmonized/02256_02512_04816/02256_02512_04816_volume.nrrd
[6] https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumetric-instance-labels/instance-labels-harmonized/02256_02512_04816/02256_02512_04816_mask.nrrd
[7] https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumetric-instance-labels/README.txt
[8] https://dl.ash2txt.org/datasets/grand-prize-banner-region/README.txt
[9] https://dl.ash2txt.org/datasets/grand-prize-banner-region/volumes/gp_volume.zarr/.zarray
[10] https://scrollprize.org/data_datasets
