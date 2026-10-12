# P0-T10 — official-surface positive control

## Conclusion

**A = MODEL_PIPELINE_VALID_GEOMETRY_MISMATCH.** The unchanged pretrained
checkpoint is strongly ink-sensitive on the registered official control.
Measured input-representation differences implicate T9 geometry, but do not
establish it as the sole cause of T9 failure or establish a successful repair.
P0-T9 remains FAIL, unchanged. No full T9 inference was repeated.

Hypothesis: the T9 checkpoint and official preprocessing discriminate ink on
a trusted producer-rendered surface. This bounded experiment tests integration,
not held-out generalization, text recognition, or a replacement reconstruction
method. No OCR, language priors, training, parameter search, or local LLM was used.

## Preflight and review

PREFLIGHT = CLEAR_WITH_NOTED_DIFFERENCES.

The independent-context online reviewer `/root/t10_preflight_review` approved
evidence collection after checking source/label registration, depth selection,
preprocessing, and frozen gates. The reviewer uses the same model/provider as
the author; this is not independent-provider review. Pre-evidence commits:
`9b8c7e25ce2542715b1900dec047c19f3030d818` and `61f3f4b`.
`outputs/p0_t10/REVIEW_RECORD.md` records findings and their dispositions.

The major requirement was measured geometry evidence rather than an automatic
A from different methods. The post-control audit and reviewer independently
verified the frame measurements below and closed that finding. Minor fixes for
degenerate bootstrap draws and pass-only auditing were committed before
inference. The official-API comparison shares model/helpers and establishes
conformance, not independent model correctness; registered labels provide the
ink-sensitivity evidence.

## GPU memory retry

T9_MEMORY_RETRY = PASS. The exact sealed T9 smoke function was called with only
its output directory redirected into T10. It verified all 5,523 existing CT
chunk hashes, checkpoint hash, real CT reads, 7.91 -> 9.6 µm sampling,
official preprocessing, `model.cuda()`, and finite probabilities in [0,1].
The three original probes Z128/152/177 had independent interpolation error 0.

Before load, CUDA reported 102,990,659,584 bytes free (95.92 GiB) and system
MemAvailable was 118,783,820 KiB (113.28 GiB). Host resource capture showed no
resident scientific competitor/local LLM. GPU and unified-memory accounting
are not interchangeable; GB10 nvidia-smi aggregate memory is unsupported.
The old OOM did not recur. This establishes current compatibility, not proof
that Nemotron was the sole historical cause. Existing CT was not redownloaded.

Execution used the local shell outside the failing sandbox and the pinned
scientific Docker image below. No Nemotron service was started.

## Exact control and provenance

Dataset: `phercparis4-w02`; trusted segment `20231031143852`.
Producer volume:

https://vesuvius-challenge-open-data.s3.amazonaws.com/PHercParis4/segments/20231031143852/surface-volumes/2.4um-0.22m-78keV-volume-20260411134726.zarr

Trusted mesh: segment's
`mesh/20231031143852-on-20260411134726-2.4um.tifxyz/{meta.json,x.tif,y.tif,z.tif}`.
The producer-rendered volume is consumed byte-for-byte, not rerendered with
new geometry. It is the smallest rendered area among the eight Paris4 entries
in the released source table, not the smallest possible official dataset.

Official corpus documentation:
https://huggingface.co/buckets/scrollprize/datasets/tree/ink_9um
The downloaded README and source metadata are retained. Registered ink labels
and supervision masks are in
`ink_9um/labels/aligned-scrollprizeorg-21slices/phercparis4-w02/`.
Their metadata binds the exact source URL, XY level2, source depth [13,97),
9.6 µm output, and annotation depth10. This is a training-corpus positive
control, not a held-out test. No precomputed independent official prediction
was needed; the official inference API was executed for conformance comparison.

Fixed crop: level2 Y[4096,4608), X[5248,5760), 512×512 pixels. The first qualifying
crop in a saved center-out list was selected using only label coverage before
CT/model evidence: 145,587 supervised pixels, 43,032 ink. No selection from
visible predictions. Acquired 16 original CT chunks and four geometry objects
after HEAD size checks, plus 32 label/mask chunks. Receipts preserve URLs,
byte counts, ETags where provided, and SHA256. Full source depth is 109 planes.

## Checkpoint and exact input conventions

Checkpoint: `scrollprize/ink_9um/hybrid_3d2d-seed42/step-075000.pth`.
SHA256: `e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab`.
Official source: Villa revision `fb8c2c4c2705746eefb42e9bb359560906f391e2`.

* Native producer scan: 2.4 µm, 78 keV. XY level2 is mean-downsampled to
  9.6 µm; original rendered depth remains 2.4 µm.
* Official `prepare_isotropic_input` selects source [13,97), averages successive
  groups of four in float32, rounds half-even to uint8, producing 21 planes.
  An independent calculation matches exactly.
* Forward central planes [2,19) give D×Y×X = 17×128×128 per model patch,
  isotropic 9.6 µm. Label plane10 becomes model depth8. No reversal or offset.
* Producer geometry defines curved points P(U,V), full 3-D tangents and a
  normal perpendicular to them. Nominal source center is plane54; pooled
  center10 averages source53..56, mean54.5, a +1.2 µm pooling-center displacement.
* Current reference renderer uses trilinear CT interpolation and supports
  normal reversal and different geometry interpolation paths. Historical
  production flip/interpolation flags are not recorded in the delivered attrs.
  We preserve the registered delivered depth ordering unchanged; we do not
  claim an independently reproduced producer rerender or known physical sign.
* Checkpoint-selected `tifxyz_robust` / robust-MAD normalization: clip at
  percentiles1/99, subtract median, divide by 1.4826 MAD, with official fallback
  behavior for degenerate scale. The same official helper is used by T9.
* Forty-nine patches, stride64, batch8, float16 CUDA autocast, no TTA,
  official Hann overlap blending, fixed seed2910. No post-evidence tuning.

## Quantitative and visual result

Probability map: 512×512, float32 retained losslessly. All pixels finite and
in [0,1]. PNG probability is 16-bit; display uses frozen [0.25,0.75] clipping.
Evaluation excludes a 64-pixel crop border and points within two pixels of the
opposite label class, and requires producer supervision.

| Metric | Result |
|---|---:|
| Evaluated pixels | 101,188 |
| Ink / nonink | 29,245 / 71,943 |
| ROC AUC | 0.9921295 |
| PR AUC (average precision) | 0.9826163 |
| Ink prevalence | 0.2890165 |
| Spatial-bootstrap 95% AUC interval | [0.9829574, 0.9982971] |
| Spatial-bootstrap 95% AP interval | [0.9567838, 0.9954991] |
| Mean ink / nonink probability | 0.7320911 / 0.2840592 |
| Standardized difference, equal-class variance average | 4.870369 |
| IoU / Dice at p>0.5 | 0.8297872 / 0.9069767 |
| TP / FP / FN / TN | 28,470 / 5,065 / 775 / 66,878 |

All preregistered gates pass. Bootstrap: 1,000 resamples of 64×64 spatial
blocks, no degenerate draws. These intervals describe this crop, not
generalization across scrolls. Independent rank-sum AUC and score-bin AP
checks agree; rank sums require float64 accumulation with this SciPy version.

The registered comparison image shows spatially localized high probabilities
following labeled ink regions and low probabilities over most labeled nonink.
No letter identification or OCR is needed. The complete official flat-inference
API gives identical orientation/dimensions and encoded probabilities
(maximum uint8 difference 0); this comparison shares the configured model.

## T9 comparison and representative patch audit

T9 uses independently traced transverse contours, with native XY = 16×saved
level4 coordinates and Z = 16×slice index. Each patch fixes an XY tangent,
normal (-ty,tx,0), and scan-Z longitudinal axis over the entire patch. It is
a straight local plane, not a curved surface-following patch. Centering is at
the frozen contour with no calibrated ink-bearing surface registration.
T9's depth index8 is exactly that center, with symmetric -8..8 offsets.

The trusted geometry audit measures 10,712 grid points near the control crop:
median |normal Z| = 0.349152 (95th percentile0.723601), median longitudinal
tilt from scan Z = 25.0344° (95th48.8687°). T9 imposes zero for both.
Normals vary over the trusted patch (median deviation from regional mean
28.2168°); T9 holds one normal constant. These are different locations, not
measurements of the true normals at T9 centers. Bounds use approximate
documented UV-grid correspondence, not a new subpixel registration.

The three preregistered T9 smoke patches Z128/152/177 show repeated cross-cut
bands and appreciable empty regions compared with the producer surface
patches. Their zero-voxel fractions are 0.1625/0.2472/0.0995, versus 0 for
the fixed control patches. Mean depth-intensity peaks occur at T9 indices
16/0/16 versus control4/4/4. These are occupancy diagnostics, not definitive
papyrus seating or ink localization; control intensity itself need not peak
at the annotated surface. Orthogonal cuts use the same 0..255 display and
aspect-auto convention for both sets, exaggerating depth equally.

No preprocessing or model patch-size mismatch was found. Both use 17×128×128,
9.6 µm target spacing, the same normalization and checkpoint. T9 resamples
native7.91 µm using scipy trilinear order1, float32, half-even uint8 rounding,
without additional anti-alias filtering. Control instead consumes trilinearly
rendered high-resolution CT with XY/depth mean pooling. No sign error was
demonstrated; a producer flip flag is unavailable. Depth placement/centering
is demonstrably defined differently but not independently calibrated at T9.

Additional confounders: 2.4 µm/78 keV versus7.91 µm/54 keV scan; training-corpus
control; different physical locations; contiguous 9.6 µm surface inference
versus T9 center estimates spaced63.28 µm along s and126.56 µm across scan Z,
without cross-Z alignment. T9's sealed correlations remain T5=0.0515,
T6=0.1081, T9=0.0021, T9 local-window=-0.0006; no re-evaluation changes them.

These measured representation differences support A and justify a separately
preregistered bounded correction experiment. They do not prove that changing
geometry alone will recover ink. Next experiment recommendation: isolate
surface following, depth seating, and longitudinal frame on a few fixed
locations, controlling scan and sampling effects where possible. No such
correction is performed here, and the Direct Thread-Stack paths remain frozen.

## Replay, preservation and artifacts

All 49 raw patches replayed with exact hashes, independently recomputed pooling
was exact, and maximum float probability error was 0 (frozen tolerance0.002).
Final verification checks 3,165 predecessor files byte-identical, all frozen
pre-evidence code/plan hashes, acquisition hashes, independent metrics, PNG
dimensions, and finite probabilities. `verification.json` records actual checks.
`SHA256SUMS` covers local raw/lossless evidence as well as tracked deliverables.
Large arrays, TIFFs, and Zarr payloads remain locally under existing Git ignores.

Scientific runtime image:
`nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a`.
Replay command from the repository root (no download):

```sh
docker run --rm --gpus all --ipc=host \
  -v "$PWD:/work:ro" -v "$PWD/outputs/p0_t10:/work/outputs/p0_t10" \
  -w /work -e PYTHONPATH=/work/data/cache/p0_t9/python_site:/work/external/villa/vesuvius/src \
  -e TMPDIR=/work/outputs/p0_t10/runtime_tmp \
  nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a \
  python -B scripts/p0_t10_infer.py --replay
.venv/bin/python -B scripts/p0_t10_verify.py
sha256sum -c outputs/p0_t10/SHA256SUMS
```

Primary artifacts are in `outputs/p0_t10/`: `gpu_memory_retry.json`,
`official_control_probability.png`, `official_control_display.png`,
`official_control_metrics.json`, `control_vs_t9_sampling.png`,
`control_vs_t9_geometry.json`, `representative_patch_comparison.png`,
`registered_control_comparison.png`, `replay.json`, `verification.json`, and
`SHA256SUMS`. Lossless `official_control_probability.npy`, `official_input.npy`,
patch outputs and original chunks remain locally available for replay.
