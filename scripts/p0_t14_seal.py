"""Publish T14 receipts/report/artifacts only after verification. Refuses a second seal."""
from pathlib import Path
import json,hashlib,shutil
R=Path(__file__).resolve().parents[1];O=R/'outputs/p0_t14'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 assert not (O/'SHA256SUMS').exists()
 v=json.loads((O/'verification.json').read_text());m=json.loads((O/'stack/metrics.json').read_text());g=json.loads((O/'geometry/geometry.json').read_text());assert v['status']=='PASS' and all(m['metric_gate_frozen_before_run'].values())
 for p in (O/'stack').iterdir():
  if p.is_dir():shutil.copytree(p,O/p.name)
  else:shutil.copy2(p,O/p.name)
 m.update({'P0_T14_COMPLETE':'YES','P0_T14_RESULT':'A','TRUSTED_GEOMETRY_USED':'YES','INK_MODEL_VALID':'YES','COHERENT_INK_STRUCTURE':'YES','CORE_IDEA_VIABLE':'YES','geometry':g,'visual_observation':'Continuous localized bright stroke follows labeled ink in the raw and physical-aspect s/Z stack. No row-to-row scrambling. Crop shows part of a stroke, not whole text lines or a page; no character identification. Native CT overlays show trusted intersections along the selected local papyrus band.','conclusion':'A local oracle upper bound supports the Direct Thread-Stack representation. Automatic geometry recovery and correspondence remain unvalidated; no claim that geometry is the only automatic-system bottleneck.','readiness_for_T13':'NOT_AUTHORIZED; no T13 execution','model_patch_geometry':'Exact sealed official producer-rendered curved trusted-frame CT; freshly inferred UV patches sampled along newly intersected physical arc-length rows. No new native CT renderer and no old probability map reused as output.'})
 (O/'metrics.json').write_text(json.dumps(m,indent=2)+'\n')
 report=f'''# P0-T14 — oracle upper bound for Direct Thread-Stack

```text
P0_T14_COMPLETE = YES
P0_T14_RESULT = A
TRUSTED_GEOMETRY_USED = YES
INK_MODEL_VALID = YES
COHERENT_INK_STRUCTURE = YES
CORE_IDEA_VIABLE = YES
```

Trusted geometry plus physical arc-length transverse rows retains coherent ink in
this local control. The output shows a continuous portion of an ink stroke matching
registered annotations, rather than row-to-row scrambling. **A applies to a local
representation upper bound**, not readable text, manuscript-wide layout, automatic
geometry success, or a competition result. No OCR or language prior was used.

## Hypothesis, frozen scope and inputs

Hypothesis: exact transverse intersections of the trusted surface, independently
parameterized from one common geometric boundary and stacked in native Z order,
can retain meaningful ink structure. T13 prerequisite work was stopped before any
metric-audit artifact was written. No new control region was searched and no
automatic thread detector was run.

Use the already validated P0-T10 PHercParis4-w02 / trusted segment20231031143852.
The source is `PHercParis4/segments/20231031143852/surface-volumes/2.4um-0.22m-78keV-volume-20260411134726.zarr`.
Native geometry is the sealed x/y/z TIFF set in `data/raw/p0_t10/`, registered to
scan20260411134726 at2.4µm. The fixed T10 level2 crop is
Y[4096,4608), X[5248,5760), 512×512 at9.6µm. Existing registered training-corpus
labels and supervision come from the same T10 control. This is not held out.

The plan uses only the crop interior[64,448) in both axes, keeping model-context
margins. It selects **all integer native Z planes crossing every U column** between
the interior's top and bottom geometry boundaries. Selection uses geometry alone,
not new ink predictions or labels. This yields **Z44000..44960 inclusive,
961 consecutive slices**, center-to-center height **2.304mm**. Row lengths range
**3.701165..3.831495mm** (median3.785124mm). No acquisition was needed.
The prior T10 crop itself was originally label-selected, disclosed here.

`plan.json` and `PRE_INFERENCE_FREEZE.json` were locally hash-frozen before model
inference, not published as a pre-run Git commit. A requires verified geometry and
model replay, AUC≥0.75, AP≥prevalence+0.15, Dice≥0.5 and visible coherent organization.
These are local upper-bound decision criteria, not guarantees of page readability.
No decision threshold or scientific parameter changed after results.

## Slice → oracle thread → physical arc → row → stack

Trusted stored XYZ is bilinearly evaluated at the level2 mean-pooling centers.
The current producer source uses a centered full-resolution canvas; its offset is
approximately0.5 native UV pixel, combined with1.5 pixels for4×4 pooling. Exact
float32-scale values and resulting offsets are in `geometry/geometry.json`.
The source convention was inspected in Villa's `QuadSurface.cpp`,
`vc_render_tifxyz.cpp` and tifxyz modules. Historical producer interpolation/flip
flags remain unavailable, so this is not a new claim of exact subpixel physical
registration. The delivered CT depth ordering is preserved.

The dense trusted grid is tessellated consistently into two triangles per cell,
then intersected with each **actual XY plane Z=k**. Strictly decreasing Z along V
at every selected U proves that each intersection is one U-monotone component.
Every row reaches both fixed U boundaries; no component is selected by ink quality.
The native XYZ polylines lie exactly in their requested transverse planes.
This triangulation evaluates trusted geometry only; it is not automatic reconstruction.

Physical arc length is accumulated from successive XYZ segment lengths×2.4µm.
Each row starts at the fixed left U boundary (global U5312) and is sampled at
**s=0,9.6,19.2,…µm**. There is no row shift, UV-width normalization, equal-s
correspondence fitting, image alignment or interpolation across missing threads.
All961 rows are stacked in ascending native Z at2.4µm intervals. Actual gaps:0;
the verified selected surface is continuous. Variable row-end padding remains
masked, never bridged. Output is961×400 with378,131 valid samples and1.63085%
padding. Provenance code1=oracle trusted sample;0=invalid/padding. These are not
OBSERVED automatic-detector provenance or inferred bridges.

## Trusted-frame CT rendering and fresh ink inference

Use the **already producer-rendered curved trusted-frame CT patches**, exactly as
validated by T10, rather than introducing an unvalidated native CT renderer or
choosing a normal sign. Thus patch rendering is reused from sealed source data;
inference and the transverse-thread output grid are new. Model context remains
in the official trusted surface UV frame. Predict the scalar ink field there,
then query it at each new physical-arc thread sample through its carried UV
coordinate. This scheduling shares49 overlapping patches across961 rows; it is
not961 independent model calls. No previously predicted T10 probability map was
used as the output input; it was compared only during verification.

Checkpoint `scrollprize/ink_9um/hybrid_3d2d-seed42/step-075000.pth`, SHA256
`e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab`.
Official source revision is the sealed T10 Villa revision
`fb8c2c4c2705746eefb42e9bb359560906f391e2`. Same pinned scientific container,
17×128×128 patches,9.6µm nominal model sampling, source depth[13,97) four-plane
float32 means/half-even uint8 rounding, central pooled layers[2,19), forward
depth order, robust-MAD `tifxyz_robust` normalization, stride64, batch8,
float16 CUDA autocast, Hann blending, seed2910, no TTA. Pooling independently
recomputes exactly and official complete inference API encoding agrees exactly.
No training, parameter search, model change or producer-label inference input.

## Ink result, visual evidence and limits

|Metric|Result|
|---|---:|
|Evaluated stack samples|{m['evaluated_samples']:,}|
|Unique rounded source label pixels|{m['unique_rounded_label_pixels']:,}|
|Ink prevalence|{m['prevalence']:.6f}|
|ROC AUC|{m['ROC_AUC']:.9f}|
|Average precision|{m['average_precision']:.9f}|
|IoU at0.5|{m['IoU_at_05']:.9f}|
|Dice at0.5|{m['Dice_at_05']:.9f}|
|Mean ink / nonink probability|{m['mean_ink_probability']:.6f} / {m['mean_nonink_probability']:.6f}|
|Adjacent native-Z correlation|{m['adjacent_Z_correlation']:.9f}|
|Four-plane (9.6µm) correlation|{m['four_Z_9p6um_correlation']:.9f}|

Labels are nearest-neighbor sampled through the same oracle UV association;
probabilities use bilinear interpolation. Evaluation retains the T10 supervision,
64-pixel border exclusion and >2-pixel distance from opposite label class.
The label scores therefore test retained ink discrimination in the new raster,
not independent physical registration. A common incorrect warp could retain label
scores: visual inspection and independent geometry checks are also necessary.
Native Z oversamples the9.6µm model field;378,131 samples are not independent
observations. High adjacent-row correlation is expected and is not alone success.
No confidence interval treats these duplicated samples as independent.

`oracle_ink_stack_raw.png` is16-bit linear probability with invalid pixels encoded0
and separately masked. Display clips only the fixed T10 range[0.25,0.75]. The raw
961-row image looks vertically stretched because Z pitch is one quarter of s pitch.
`physical_aspect_stack.png` uses box display resampling to400×241, preserving
nominal s/Z physical aspect. Z height is axial height, not geodesic page height;
local tilt/shear remains. The comparison shows one coherent partial stroke with
low-probability surrounding material. It does not show complete text lines or
establish broad manuscript organization; no character was identified.

Five evenly spaced rows have trusted-render CT and native-XY thread plots. Three
additional overlays use the already available native transverse CT atZ44449,
44456,44463. These show the oracle section along the selected local papyrus band.
No extra CT or region search was needed. Physical registration uncertainty remains
**UNCALIBRATED**; numerical consistency must not be called physical accuracy.

## Verification, preserved failures and conclusion

Both geometry/stack analyses replay exactly, including arrays, JSON and images.
All49 model patches replay with exact raw patch hashes and **maximum probability
error0**. Fresh inference also equals sealed T10 probability output exactly.
Independent rank-sum AUC agrees; PNG encoding and masks verify. All sampled points
have native Z error0. Neighbor sample chords range9.566537..9.600000µm; accumulated
polyline arc step is9.6µm. A separate four-corner interpolation check differs from
tessellated/resampled geometry by at most**1.533369µm** on every101st valid sample,
below one native voxel. This is a numerical approximation check, not physical truth.

**{v['predecessor_files_unchanged']} predecessor files and all{v['prior_seals_verified']}
existing seals verify unchanged**, including P0-T1–T12, blocked T13 and every control
experiment. T13 and its metric audit were not executed. A setup reader initially
assumed all legacy seal paths were repository-relative; it was corrected before
the successful preservation baseline. A NumPy Boolean JSON serialization failure
occurred after stack arrays/images were computed; original code and both failed
outputs remain in `serialization_amendment/`. The only amendment converts NumPy
scalars to JSON primitives. Verification proves code changed only there and every
pre-amendment array/image matches the final result. No scientific retuning occurred.

Conclusion **A: ORACLE_DIRECT_THREAD_STACK_WORKS**, locally. Given trusted geometry,
a consistent boundary and validated surface-frame ink inference, physical-arc rows
in true Z order can retain meaningful ink. The representation itself is therefore
viable in this control. Automatic sheet recovery, cross-slice correspondence,
centering and frame construction remain unresolved practical bottlenecks. This
experiment neither validates them nor proves they are the only possible failure
sources. Generalization to more extensive, multiply intersecting or discontinuous
surfaces is untested. No further experiment is executed; research direction returns
to the user and ChatGPT.

## Artifacts and reproduction

Report: `docs/p0_t14_oracle_direct_thread_stack.md`.
Artifacts: `outputs/p0_t14/`, including requested root-level images, overlays,
metrics, verification and SHA256SUMS; `geometry/threads.npz` retains XYZ/UV/s/Z/mask;
`inference/` retains new model outputs/patch hashes; `oracle_stack.npz` retains
lossless probabilities, labels and masks; replay and failure evidence are separate.
All modest scientific arrays are tracked; redundant replay/source-input copies
remain locally checksummed and are reproducible from existing sealed inputs.

Entry points: `scripts/p0_t14_oracle.py`, `p0_t14_infer.py`, `p0_t14_verify.py`,
`p0_t14_seal.py`. Runtime commands are in `outputs/p0_t14/reproduction.md`.
Stage creation refuses existing directories. For a new replay use isolated fresh
output destinations/checkout; never overwrite this seal or predecessor files.
Read-only checksum verification: `sha256sum -c outputs/p0_t14/SHA256SUMS` from the
repository root. Commit/push identity is reported separately to avoid embedding a
commit's own hash. Publish to origin/main without force and stop after T14.
'''
 (R/'docs/p0_t14_oracle_direct_thread_stack.md').write_text(report)
 files=sorted([p for p in O.rglob('*') if p.is_file() and 'runtime_tmp' not in p.parts]+list((R/'scripts').glob('p0_t14*.py'))+[R/'docs/p0_t14_oracle_direct_thread_stack.md'])
 (O/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(R)}\n' for p in files))
 for line in (O/'SHA256SUMS').read_text().splitlines():h,n=line.split('  ',1);assert sha(R/n)==h
 print('sealed',len(files),'files')
if __name__=='__main__':main()
