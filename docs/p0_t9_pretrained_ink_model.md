# P0-T9 — official pretrained readout on frozen Direct Thread-Stack geometry

Final status: **P0_T9_COMPLETE=YES; P0_T9_RESULT=FAIL**. Full inference and replay are complete. The pretrained readout does not produce a coherent 2-D ink-like reconstruction on the frozen geometry. Runtime/reproducibility checks PASS; scale-up is not justified. Earlier blocker entries below are preserved historical records; the final results section supersedes their pending status.

## Final results (2026-10-12 UTC)

### Hypothesis, input and unchanged settings

The official pretrained ink model, sampled on the frozen Direct Thread-Stack,
would reveal coherent localized responses persisting across neighboring Z rows,
with stronger adjacent/local correspondence than raw CT and profile mean.
This hypothesis was not supported under the frozen configuration.

Input is the existing native `54keV_7.91um_Scroll1A.zarr` CT acquisition listed
in `acquisition.json`: 5,523 chunks, 4,443,682,806 bytes. The saved inventory has
chunk-coordinate bounding ranges Z=15..22, Y=3..37, X=16..58; it is a sparse
inventory, not every chunk in that bounding box. Geometry is exactly T5 Z=128..177
at 126.56 µm spacing, with its saved T6 tangents/normals. Native center Z indices
are 2048..2832. Actual neighborhoods and all source chunks are specified by the
frozen coordinates, patch function and acquisition manifest.

The checkpoint, 7.91→9.6 µm mapping, 17×128×128 patches, batch size 8,
float16 autocast, depth/orientation, four-center-pixel readout and masks are
unchanged. Official checkpoint `robust_mad` maps to runtime `tifxyz_robust` in
the official implementation; this is not a preprocessing substitution.
All bridge and padding exclusions remain masked. No OCR, language priors,
training, geometry repair, row alignment or signal-driven setting selection
was performed.

### Completed inference and quantitative comparison

The lossless float32 probability raster is **50 rows × 4,482 columns**.
All **130,629 valid frozen centers** were inferred (100% of valid T5 centers;
58.290495% of the padded raster). The remaining 41.709505% is masked/NaN,
not a failed prediction. Finite probabilities range from 0.203806 to 0.817268;
median 0.272446. The physical-aspect PNG is 4,482×100 pixels because Z spacing
is twice the s spacing; the raw/display/mask/provenance PNGs are 4,482×50.

| Metric | T5 raw CT | T6 profile mean | T9 probability |
|---|---:|---:|---:|
| Mean adjacent-Z zero-shift Pearson r (49 pairs) | 0.051516 | 0.108094 | 0.002067 |
| Mean local zero-shift r (1,001 windows) | 0.012505 | 0.029243 | -0.000613 |
| Adjacent pairs above row-shift control 95th percentile | 9 | 12 | 3 |
| Local best-lag peaks above corresponding control 95th percentile | 108 | 104 | 104 |

T9 minus T5 mean adjacent r is -0.049449, with the frozen paired block
bootstrap interval [-0.069650, -0.028303]. The analogous local pair-mean
difference is -0.012874, interval [-0.023109, -0.001531]. Diagnostic lags are
not applied to the reconstruction. T6 values are the preserved T6 evaluation;
T5 is recomputed by the same frozen evaluator. These results do not support
improved cross-Z coherence.

At raw p>0.5, 6,019 samples form 3,971 four-connected components. Fourteen
components span at least three rows (99 samples), versus control-95 values
10 and 71.1. Their empirical upper-tail p values are both 0.005 against the
199 fixed row-shift controls. This is a limited localized persistence signal,
not evidence that those components are ink. The largest vertical span is only
four rows, equal to control-95 (upper-tail p=0.23); the largest component has
20 samples. These related diagnostics are not multiplicity-corrected, and
overlapping model neighborhoods induce dependence. The small component excess
does not override failed broad coherence and visual gates.

### Qualitative assessment and conclusion

Inspected `comparison.png` and `fixed_window_inspection.png` (four evenly
spaced 512-sample windows, chosen without response-based selection). T9 shows
scattered bright responses and short within-row streaks, not a coherent
cross-row pattern. White missing areas in these plots are masks, not signal.
No convincing localized/persistent ink-like structure or letters are visible.
The result is **FAIL for coherent reconstruction**, while execution and
reproducibility PASS. `INK_STRUCTURE_VISIBLE=NO; READY_FOR_SCALE_UP=NO`.
This does not establish that the CT contains no ink or that every possible
Direct Spiral-Stack implementation fails.

The frozen 96-probe sensitivity check was not used to change the primary
output. Normal reversal changes probability by mean absolute 0.076550
(maximum 0.412325); depth offsets -2 and +2 give means 0.048131 and 0.038125.
This reinforces the unresolved physical placement/orientation limitation.
Other limitations remain: local straight tangent-plane approximation,
uncalibrated native voxel-center registration, inherited row correspondence
instability, and inclusion of Paris4 in the model's training corpus.

Next recommended experiment (not executed): a separately specified bounded
surface-placement/orientation validation against independently registered
native-CT reference locations, before any scale-up. Preserve this negative
stack; do not tune against these responses or apparent characters.

### Execution recovery and replay

The sandbox RTM_NEWADDR failure was bypassed with the approved local
`exec_command` shell (`require_escalated`). Acquisition was already complete;
no downloads or data modifications occurred. Fresh checks verified all chunk
sizes/hashes, the checkpoint and original frozen code, and 2,892 predecessor
files. The first fresh GPU smoke failed while Nemotron was resident. With
explicit approval, the exact existing container was stopped cleanly (exit 0),
removing its reported 91,740 MiB GPU allocation. CUDA then reported
100,837,851,136 free bytes (93.91 GiB). Aggregate GPU memory is unsupported on
this unified-memory GB10; the removed process allocation is not a directly
measured aggregate delta.

The exact fresh real-CT smoke passed after the stop. Full inference then
exposed a shared-NPZ concurrent decompression error before any row was written.
A focused test reproduced it 20/20 times; eager serial loading of all 50
archives, followed by parallel array access, preserved exact array values.
The sole runner repair materializes the saved arrays before worker access.
No numerical setting or input changed. The original source and hash-bound
infrastructure amendment are retained in `nemotron_recovery/`; the original
`FROZEN_RUN_PLAN.json` was not rewritten. The post-repair smoke also passed.

All 50 rows finished in the repaired run (logged elapsed 1,467.4 seconds).
The Z=139 progress update was stale: the workflow failed to advance after
inference completion. Fresh observations at 2026-10-12 03:44:47.793152 and
03:46:12.437788 UTC, separated by an explicit 60-second wait (84.644636 seconds
overall), both found Z=177, 50 rows, no inference process, GPU utilization 0%,
and the completion log. The newest inference artifact was unchanged at
2026-10-11 22:16:43.890004 UTC. This was completed inference with stalled
orchestration, not an active/hung inference job. CPU and disk-I/O snapshots
are retained. No finished slice was rerun during this recovery.

Replay verifies all 50 row assemblies and 130,629 coordinate mappings,
byte-identical mask/provenance/physical-coordinate arrays, all 2,892 prior
preservation entries, and 96 bit-exact native-CT patches. All 384 fixed primary
and sensitivity predictions replay within the frozen 0.002 tolerance;
maximum error is 0.000747442. This is a bounded inference replay, not a second
full 130,629-center inference. See `verification.json`.

### Artifacts and reproduction

All requested products are under `outputs/p0_t9/`: raw/display probability
PNGs, mask/provenance PNGs, physical-aspect PNG, `metrics.json`,
`model_metadata.json`, `inference_parameters.json`, `verification.json`,
and `SHA256SUMS`. Lossless `ink_probability.npy`, all 50 row NPZ files,
96 probe NPZ files and smoke arrays remain locally preserved and Git-ignored.
`comparison.png`, fixed-window inspection, sensitivity/component summaries,
evaluation CSVs/controls and recovery logs supply additional evidence.

Reproduction uses the pinned container shown in the historical commands below
and the preserved local data; no acquisition step is needed. The inference
runner recognizes completed row files. Run `p0_t9_verify.py` in that container
for GPU replay. `p0_t9_evaluate.py` requires a fresh evaluation output directory;
do not overwrite the sealed evaluation. `p0_t9_inspect.py` reproduces the
fixed visual windows; `p0_t9_seal.py` validates required artifacts, scientific
verdict, original/fixed source binding, predecessor hashes and restoration,
then seals every T9 file including ignored arrays. Validate the final manifest
from the repository root with `sha256sum -c outputs/p0_t9/SHA256SUMS`.

### Nemotron restoration and publication

The exact original container `nemotron-nano-agent`
(`461b46345448e70ea9de9fdb196a1ed9d957a21f8cea8446b6b06b990ad1825d`),
image `vllm/vllm-openai:v0.20.0` / image ID
`sha256:bf05c8a0d56f7ed2fe279e0c0c876dd2bb813dd34702d975fa6ab640860ce825`,
was restarted with `docker start`, not recreated. Recorded command, image,
mounts, ports, restart policy, GPU request and other selected configuration
fields match the pre-stop record exactly. No model files or configuration
were removed or changed. The container is running, `/health` returns HTTP 200,
`/v1/models` responds, and a short `/v1/chat/completions` request returned `OK`.
Operational restoration health is PASS.

The requested name `nemotron_3_nano_omni` was **not** the pre-existing served
model ID: the endpoint returned `nemotron-nano-agent` both before stopping and
after restoration. The exact-name check is therefore FAIL, separately from
operational health. It was not renamed because restoration was explicitly
required to preserve the existing configuration. See
`nemotron_recovery/before.json` and `restoration.json` for readback evidence.

The final T9 commit includes this negative result, scoped source, compact
artifacts, recovery logs, and the checksum manifest. Large/lossless arrays and
native CT remain local and Git-ignored, with checksums preserved. Publication
targets `origin/main` without force; final commit/remote identity is reported
in the task response rather than embedded in its own hashed report.

## Frozen scope

P0-T5 geometry Z128..177, 50 independent transverse slices; no tracing, topology, bridge, coordinate, mask or provenance changes. No training, fine-tuning, OCR, text priors, cross-Z path correction, or new surface reconstruction.

Selected model: official `scrollprize/ink_9um`, `hybrid_3d2d-seed42/step-075000.pth`, revision `7109667e2607db1b90c37c8b09cb876ea7fe7bb1`. Checkpoint SHA256 `e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab`. Strict weight loading and a zero-input runtime-only GPU test passed; this is not real-CT evidence.

## Important scale correction

The stored T5 geometry and coarse CT use126.56µm level4 spacing. The original source has7.91µm native spacing. Published Zarr metadata maps level4 indices to native by a factor16 with no translation. T9 uses this mapping without fitting an offset; exact voxel-center registration remains uncalibrated. No coarse-volume upsampling is substituted for native CT.

The official model expects a17×128×128 surface-oriented input. The frozen plan uses9.6µm model spacing, trilinear native CT sampling at9.6/7.91 native-index steps, then half-even quantization to uint8 and official robust-MAD preprocessing. Model depth follows the saved T6 XY normal, model X follows its XY tangent, model Y follows native scan Z. This is a local straight tangent-plane approximation, not a recovered 3D sheet. Its validity over a1.22mm neighborhood is uncertain. The probability at the geometric center is the mean of the four central output pixels. Every inferred bridge remains masked.

Only the saved T6 frame is reused; no adjacent-slice geometry is consulted. A fixed96-probe sign/depth sensitivity check compares reverse normal and depth offsets±2model voxels, without selecting settings by apparent letters. The primary output remains unchanged.

## Recovery audit

The earlier12-worker downloader was making progress but was interrupted during an intended concurrency change. The subsequent restart did not execute. On recovery inspection no T9 downloader was alive; two measurements found417files /341,180,981bytes and unchanged mtimes. Every cached chunk passed size and Blosc-decoding checks. Recovery resumed the same list and destination with48workers, preserving cache and the4,500,000,000-byte CT cap. Measured complete CT payload is4,443,682,806bytes in5523chunks; no full-scroll download is requested.

See `outputs/p0_t9/recovery_observation_1.json`, `recovery_observation_2.json`, `download_recovery.log`, `native_chunk_budget.json`, and `FROZEN_RUN_PLAN.json`.

## Verified recovery and exact blocker

A60-second wait separated fresh observations (recorded snapshot interval81.2846seconds including command overhead). Both found4,443,682,806bytes and5523complete chunks. All expected filenames, sizes and saved SHA256hashes passed; no partial files remain. Unchanging counters indicate completed acquisition, not a hung transfer. See `stall_check_before.json`, `stall_check_after.json`, and `acquisition.json` for per-file records. No active downloader remained.

The queued real-CT smoke test passed after validating all5523files and the checkpoint hash. Three frozen probes atZ128,152,177 verified17×128×128neighborhood extraction at9.6/7.91native-index spacing, official preprocessing, exact independent interpolation checks at selected voxels, model loading and finite probabilities. Probe center probabilities were0.3057094,0.2728068,0.3312626. These are runtime checks, not evidence of ink structure. Full details and saved probe arrays are in `smoke_test.json` and `smoke_<z>.npz`.

The subsequent full-inference container exited before producing any finalized row: `torch.AcceleratorError: CUDA error: out of memory` at `model.cuda()`. A separate `VLLM::EngineCore` process PID2659531 reported91740MiB of GPU memory. Competing resource use is a plausible cause; it is not proven to be the sole cause. No unrelated service was stopped. Scientific parameters, geometry, model, checkpoint and byte cap remain unchanged. See `full_inference.log` and `recovery_status.json`.

All50row inference, evaluation, visual assessment, replay, final seal and Git publication remain pending. This is an infrastructure blocker, not a scientificFAIL. No T9 worker remains running.

## Exact resume command

After the owner makes sufficient GPU memory available, run from the repository root (no redownload):

```sh
docker run --rm --gpus all -v /home/wpotech/deployments/vesuvius-spiral-stack:/work -w /work -e PYTHONPATH=/work/data/cache/p0_t9/python_site:/work/external/villa/vesuvius/src nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a python -B scripts/p0_t9_infer.py
```

Do not interpret this incomplete report as PASS/PARTIAL/FAIL. No completed-milestone commit or push was performed.

## Official references

- Model and display convention: https://huggingface.co/scrollprize/ink_9um
- Recommended checkpoint: https://scrollprize.substack.com/p/from-ct-scan-to-ancient-text-a-first
- Official implementation: `external/villa/vesuvius/src/vesuvius/ink_detection/`, Villa revision `fb8c2c4c2705746eefb42e9bb359560906f391e2`.

Native CT payload, model responses and final classification will be recorded once execution completes. Earlier T1–T8 evidence remains unchanged.

## Fresh execution-backend recovery attempt (2026-10-11)

The local `exec_command` shell with `sandbox_permissions=require_escalated` successfully bypassed the RTM_NEWADDR sandbox startup failure. No acquisition command was run and no scientific input or parameter was changed.

Fresh validation matched all 5,523 chunk filenames and sizes against the saved acquisition and budget manifests (4,443,682,806 bytes), checked 33 representative hashes initially, then independently rechecked all 5,523 chunk hashes. The checkpoint and frozen code hashes match. All 2,892 predecessor preservation entries match. Manifest digests and receipts are under `outputs/p0_t9/execution_recovery/`.

The unchanged prepared smoke test was rerun in the pinned container. It exited 1 at `model.cuda()` with `torch.AcceleratorError: CUDA error: out of memory`, before new CT-patch inference. The previous smoke PASS is preserved and is not being represented as a fresh PASS. Full inference was not started in this recovery attempt. No probability stack, visual verdict, or scale-up readiness is established.

The resident `nemotron-nano-agent` server is healthy (HTTP 200), has a local model mount, and its VLLM process reports 91,740 MiB GPU allocation. It remains running: stopping another service requires explicit approval. The resource contention is a plausible cause, not a proven sole cause. No commit or push of a completed milestone was performed.

Fresh failing command (repository root):

```sh
docker run --rm --gpus all -v /home/wpotech/deployments/vesuvius-spiral-stack:/work -w /work -e PYTHONPATH=/work/data/cache/p0_t9/python_site:/work/external/villa/vesuvius/src nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a python -B scripts/p0_t9_smoke.py > outputs/p0_t9/execution_recovery/fresh_smoke.log 2>&1
```
