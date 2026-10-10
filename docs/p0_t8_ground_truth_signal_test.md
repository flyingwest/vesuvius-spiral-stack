# P0-T8 — Ground-truth signal test

## Final status

P0_T8_COMPLETE = YES

P0_T8_RESULT = C — NO_USEFUL_DETERMINISTIC_INK_SIGNAL

BEST_FEATURE = NONE

LEARNED_INK_MODEL_JUSTIFIED = NO

READY_FOR_NEXT_STEP = YES — the evidence is complete for deciding a bounded follow-up, not an authorization or justification to train a model.

No frozen feature established reproducible ink discrimination in this test. The highest confirmation point estimate was **range_wide**, ROC AUC 0.5254, 95% spatial-bootstrap CI [0.4814, 0.5708]. This is not evidence that the underlying CT contains no ink information: the experiment is small, single-fragment, fixed-depth, and uncertainty includes weak effects.

## Status audit and scope

The feature computation and plots were already complete on October 9. At the final October 10 audit, no T8 analysis worker was running; the pending saved-data verifier had completed successfully. Reporting, consolidated metrics, sealing and Git publication remained unfinished. This was an incomplete handoff/reporting stage, not an observed stalled CT computation. No claim is made about why the interaction paused. Finalization reused all completed metrics without rerunning feature evaluation.

This remains the deterministic labeled-data experiment. It was **not superseded by an official pretrained model experiment**. No model was trained or invoked, no OCR or language prior was used, and no Direct Thread-Stack geometry was modified. Existing official surface-aligned data were consumed as a registered test fixture; no new mesh or surface reconstruction was built.

## Hypothesis and dataset selection

Test whether the existing T6/T7 profile statistics differ between registered known-ink and labeled non-ink locations. Repeatability across slices is not an ink discriminator by itself.

Data used in this experiment were obtained from the **EduceLab-Scrolls** dataset: Frag3 / PHercParis1Fr34, 54 keV exposed surface, surface-volume UUID `20230301213423`. The official package supplies CT surface layers, infrared photography and ink labels. The dataset paper describes spectral-photograph-to-CT alignment; it does not provide a per-pixel registration-error bound for this crop. [Parsons et al., 2023](https://arxiv.org/abs/2304.02084).

Source directory: https://dl.ash2txt.org/fragments/Frag3/PHercParis1Fr34.volpkg/working/54keV_exposed_surface/

The four inspected legacy labeled fragments had layer-32 TIFF sizes 103,571,722 (Frag1), 281,948,145 (Frag2), 79,847,973 (Frag3), and 118,745,433 bytes (Frag4). Frag3 was the smallest of this comparison, not a claimed global minimum over every official dataset. No suitable local ink-label fixture was identified; prior manual segmentation references are not ink truth.

Full surface volume: 65 × 7606 × 5249, axes depth/Y/X, uint16. Metadata records 3.24 µm voxel spacing. Crop: global Y `[2784,4832)`, X `[1600,3648)`, depth layers **28,30,31,32,33,34,36**. Seven layers form an array of shape 7 × 2048 × 2048. Offsets around nominal layer32 are −4,−2,−1,0,+1,+2,+4, corresponding to ±12.96 µm at the outer samples. The nominal XY crop extent is 6.63552 mm; actual surface metric can vary with flattening.

The [official surface-volume tutorial](https://scrollprize.org/tutorial3) documents the normal-sampled layer convention, nominal center layer32, and PPM mapping `pixel(y,x) → original XYZ + surface normal`. These are producer-supplied coordinates, not a thread reconstructed here. Exact ink depth is unverified; layer32 was fixed before test metrics, not optimized.

Acquired payload: **171,297,546 bytes**, below the frozen220,000,000-byte cap. HTTP206 byte ranges fetched just the required uncompressed CT row spans and sampled PPM rows. Full PNG labels/mask/IR and small metadata were acquired. `acquisition.json` records URLs, half-open byte ranges, response headers, payload hashes and exact local files. Sparse TIFF metadata files have large logical sizes but only headers/IFDs stored; acquisition payload and local logical file size are different quantities.

## Registration verification

All label/mask/IR images and TIFF headers agree on7606×5249 XY dimensions. Label and mask encoding is **0/1**, not0/255; the initial validator's encoding assumption failed and was corrected without changing pixels, sampling parameters or feature definitions. Labels are not model predictions. Pixel locations map directly into the official surface volume; no registration shift, rotation or scale was fitted against ink statistics.

The overlay was inspected: supplied ink boundaries follow dark registered-IR strokes; CT and IR share material features but have different contrast. Original PPM normals at usable grid points have unit length to numerical precision. All4,194,304 crop pixels lie within the supplied papyrus mask. Globally48 positive label pixels lie outside the papyrus mask; these are not in the evaluated crop. CT nonzero fraction within the crop mask is 0.99999465. IR medians are134 at positives and152 at negatives, a registration sanity check only, never a CT feature.

Registration status: **documented producer registration plus consistent dimensions, coordinate-map checks and visual overlay; suitable for this bounded evaluation, not independently certified subpixel registration**. A two-pixel class-boundary exclusion reduces boundary ambiguity but cannot prove perfect alignment. Label0 means supplied no-ink annotation, not independently chemically verified absence; faint or omitted strokes may contaminate negatives.

## Frozen sampling and amendment

`plan.json` and `PLAN_FROZEN.json` were saved before feature-label metrics. The initial top/bottom split was amended after examining the registration overlay, which revealed two separated ink bands, but **before any CT discrimination metrics**. `analysis_plan_amendment.json` and `ANALYSIS_PLAN_FROZEN.json` preserve the reason and code hashes. This is disclosed label-informed sampling design, not a fully label-blind preregistration. The crop, depths, features, control calipers and inferential rules were unchanged.

Final discovery region: local X<1024. Confirmation region: local X≥1024. Four512×512 confirmation corner windows check regional consistency. The middle1024 rows remain included in overall confirmation but not these four regional checks. `final_sampling_diagnostic.png` is authoritative for the final split. The original `registration_diagnostic.png` deliberately retains its earlier top/bottom annotation as pre-amendment evidence.

Samples use a16-pixel grid starting at local(8,8). Crop-mask distance must exceed16 pixels and the distance to the opposite label class must exceed2 pixels. The distance transform is calculated within the crop; all-true mask edge handling can exclude an extra crop-edge sample, and is reproduced unchanged. Positive counts: discovery 1475, confirmation 1160. Negative counts: discovery 6597, confirmation 6947. Total **16179** saved profiles. Grid pitch is51.84 µm nominally. No neighboring T5 slice or T5 path is used.

Each feature's score direction was chosen from discovery AUC only and then frozen; all ten happened to select larger values for ink. No formula, depth, threshold or weight was fitted. Confirmation prevalence is 0.143086; that is the uninformative precision reference, not0.5. Matched cohorts are balanced by construction and have AP reference0.5.

## Frozen features

Native uint16 intensities are promoted to float64, with no normalization learned from labels. Let P(t) be the CT value at offset t. Definitions are imported unchanged from `scripts/p0_t7_family.py`:

- center_intensity: P(0).
- mean_wide: equal-weight mean over all seven offsets (T6 profile mean).
- median_wide: median over the seven values.
- std_wide: population standard deviation over seven values.
- range_wide: maximum minus minimum.
- center_minus_mean_wide: P(0) minus seven-point mean.
- asymmetry_wide: absolute difference between three positive-offset and three negative-offset means.
- gradient_abs_narrow: |P(+1)−P(−1)|/2 (T6 normal gradient).
- curvature_narrow: 2P(0)−P(−1)−P(+1) (T6 negative second derivative).
- ridge_proxy_wide: P(0)−mean(P(−4),P(−2),P(+2),P(+4)) (T6 ridge contrast).

The **pixel-domain formulas** are unchanged, but the physical profile width is not: T5/T6 used126.56 µm voxels and here3.24 µm. Nor is this surface-normal profile the same geometry as the original transverse2D normal. Therefore this test cannot establish ink sensitivity of the coarse T5 raster by transfer. It tests these summaries on native-resolution registered fragment neighborhoods.

## Statistical methods

ROC AUC uses tied-rank probability; PR AUC is **average precision**, the step-integral convention, not trapezoidal interpolation. Effect size d is the ink-minus-nonink mean difference divided by the square root of the mean within-class population variances. Tables give95% percentile spatial-bootstrap intervals.

Bootstrap:2000 draws, seed2808, resampling256×256 spatial tiles with replacement. Each selected tile carries all its samples, preserving local dependence; pixels are not independently bootstrapped. The same procedure gives99.5% intervals for a ten-feature Bonferroni-style familywise screen. Their extreme tails are only approximately estimated with2000 draws. Remaining dependence between adjacent blocks may make intervals optimistic; a single crop is not independent-fragment replication.

Frozen A/B evidence requires confirmation99.5% AUC lower bound above0.5, at least3/4 regions above chance, and controlled paired-concordance95% lower bound above0.5. A requires AUC≥0.60; B permits smaller effects. No feature qualifies, even before multiple-comparison correction on the main confirmation AUC. This is lack of established useful signal, **not an equivalence test proving exactly zero effect**.

## Confirmation metrics for every feature

| Feature | Ink / non-ink | ROC AUC [95% CI] | PR AUC/AP [95% CI] | Effect d [95% CI] |
|---|---:|---|---|---|
| center_intensity | 1160 / 6947 | 0.5008 [0.4506, 0.5439] | 0.1509 [0.0852, 0.2239] | 0.0305 [-0.1255, 0.1603] |
| gradient_abs_narrow | 1160 / 6947 | 0.5016 [0.4639, 0.5391] | 0.1407 [0.0821, 0.2068] | -0.0150 [-0.1424, 0.1094] |
| curvature_narrow | 1160 / 6947 | 0.5033 [0.4741, 0.5287] | 0.1488 [0.0867, 0.2164] | 0.0292 [-0.0575, 0.1024] |
| mean_wide | 1160 / 6947 | 0.5038 [0.4505, 0.5499] | 0.1511 [0.0856, 0.2243] | 0.0342 [-0.1258, 0.1681] |
| median_wide | 1160 / 6947 | 0.5000 [0.4491, 0.5434] | 0.1507 [0.0850, 0.2239] | 0.0303 [-0.1274, 0.1617] |
| std_wide | 1160 / 6947 | 0.5250 [0.4820, 0.5703] | 0.1491 [0.0884, 0.2196] | 0.0674 [-0.0771, 0.2177] |
| range_wide | 1160 / 6947 | 0.5254 [0.4814, 0.5708] | 0.1510 [0.0901, 0.2227] | 0.0763 [-0.0676, 0.2259] |
| center_minus_mean_wide | 1160 / 6947 | 0.5024 [0.4643, 0.5347] | 0.1514 [0.0881, 0.2205] | 0.0170 [-0.0988, 0.1175] |
| asymmetry_wide | 1160 / 6947 | 0.5034 [0.4624, 0.5418] | 0.1423 [0.0830, 0.2094] | 0.0005 [-0.1350, 0.1293] |
| ridge_proxy_wide | 1160 / 6947 | 0.5022 [0.4634, 0.5347] | 0.1515 [0.0883, 0.2207] | 0.0157 [-0.1026, 0.1184] |

## Distribution summaries

Values are in native CT-derived units; medians and interquartile ranges are shown. Full5/25/50/75/95% quantiles are in `feature_metrics.json` and `metrics.json`.

| Feature | Ink median [Q1,Q3] | Non-ink median [Q1,Q3] |
|---|---|---|
| center_intensity | 24142.000 [19642.500, 32751.000] | 24587.000 [20167.500, 31495.000] |
| gradient_abs_narrow | 2225.750 [1088.875, 3587.250] | 2167.500 [1022.750, 3667.000] |
| curvature_narrow | -238.500 [-1369.000, 1092.000] | -211.000 [-1250.000, 960.000] |
| mean_wide | 25666.500 [22045.679, 30306.893] | 25692.143 [22255.214, 29999.429] |
| median_wide | 24222.000 [20861.250, 31557.000] | 24541.000 [21018.500, 30802.500] |
| std_wide | 5268.657 [3555.603, 6937.657] | 4972.392 [3297.631, 6835.535] |
| range_wide | 15356.000 [10421.250, 20145.250] | 14499.000 [9722.000, 19669.000] |
| center_minus_mean_wide | -445.214 [-2668.036, 2180.500] | -571.286 [-2415.571, 1774.571] |
| asymmetry_wide | 7994.833 [4099.500, 12055.583] | 7675.333 [3797.167, 12212.000] |
| ridge_proxy_wide | -732.375 [-4417.875, 3461.500] | -933.000 [-3939.125, 2857.250] |

`feature_distributions.png` shows strongly overlapping distributions. Plot bins span pooled0.5–99.5% quantiles for readability only; **all values** enter metrics. No conclusion is drawn from apparent writing or plot appearance.

## Nearby and nuisance-matched controls

Greedy deterministic nearest-negative matching without replacement: same discovery/confirmation half and same256×256 tile, separation≤128 pixels. This yields1050 confirmation pairs. The stricter match additionally requires surface-normal dot product≥0.98, Gaussian-sigma32 local CT-baseline difference≤0.25 pooled SD, and Gaussian-sigma2 spatial-gradient-strength difference≤0.5 pooled SD, yielding561 pairs. These nuisance descriptors are controls, not extra candidate ink features. Matching changes the target population and may remove genuine material-mediated effects as well as confounding; unadjusted results remain primary.

Median separation is64 pixels for nearby matching and65.96969 for nuisance matching. Minimum normal dot improves from0.83879 to0.98004. Mean absolute baseline difference falls from0.76747 to0.12950 SD; mean absolute gradient-strength difference from0.93899 to0.23761 SD. Full balance is saved in `matching_balance.csv`.

| Feature | Nearby-pool AUC [95% CI] | Nuisance-pool AUC [95% CI] | Within-pair concordance [95% CI] | Regions AUC>0.5 |
|---|---|---|---|---:|
| center_intensity | 0.5201 [0.4694, 0.5658] | 0.5033 [0.4787, 0.5287] | 0.4938 [0.4525, 0.5340] | 2/4 |
| gradient_abs_narrow | 0.5159 [0.4832, 0.5426] | 0.5038 [0.4794, 0.5269] | 0.4973 [0.4617, 0.5368] | 3/4 |
| curvature_narrow | 0.5119 [0.4770, 0.5420] | 0.5009 [0.4704, 0.5350] | 0.4920 [0.4568, 0.5358] | 3/4 |
| mean_wide | 0.5294 [0.4749, 0.5769] | 0.5105 [0.4875, 0.5339] | 0.5187 [0.4858, 0.5504] | 2/4 |
| median_wide | 0.5204 [0.4682, 0.5671] | 0.5046 [0.4808, 0.5283] | 0.4831 [0.4376, 0.5186] | 2/4 |
| std_wide | 0.5453 [0.5074, 0.5757] | 0.5412 [0.5079, 0.5696] | 0.5294 [0.4958, 0.5589] | 3/4 |
| range_wide | 0.5441 [0.5031, 0.5777] | 0.5394 [0.5029, 0.5709] | 0.5223 [0.4818, 0.5584] | 3/4 |
| center_minus_mean_wide | 0.5113 [0.4716, 0.5489] | 0.5005 [0.4728, 0.5323] | 0.5045 [0.4715, 0.5394] | 2/4 |
| asymmetry_wide | 0.5194 [0.4834, 0.5517] | 0.5111 [0.4869, 0.5338] | 0.5116 [0.4842, 0.5357] | 3/4 |
| ridge_proxy_wide | 0.5111 [0.4708, 0.5494] | 0.5001 [0.4725, 0.5315] | 0.5027 [0.4640, 0.5416] | 2/4 |

Some **nominal** control findings are positive and must not be hidden: standard deviation nuisance-pool AUC0.5412 (95%CI0.5079–0.5696), range0.5394 (0.5029–0.5709). Their99.5% intervals include0.5; their within-pair concordance intervals include0.5; and main confirmation intervals include0.5. These are weak exploratory hints, not replicated controlled discrimination sufficient for categoryB. Discovery AUC for both was about0.556, attenuating to about0.525 in confirmation.

## Replay, preservation and evidence

`verification.json` reports PASS:

- 2857 prior files unchanged, including frozen geometry, rows and predecessor evidence.
- 163 acquisition files hash-verified.
- CT crop reconstructed exactly from saved HTTP-range payloads; all profiles and label coordinates replayed exactly.
- All ten feature formulas independently recomputed exactly.
- All40 feature/cohort ROC, AP and effect-size summaries independently checked.
- 400 bootstrap draws independently regenerated by literal block replication (10 of2000 per feature/cohort); quantiles of all saved bootstrap distributions verified. This is not a claim that every bootstrap draw was independently regenerated.
- Matching constraints, discovery-only directions, frozen plan hashes and PNG decoding checked.

`SHA256SUMS` covers this report, all T8 outputs except the manifest itself, T8 scripts, and the saved acquisition payloads/arrays. It is generated and every entry verified by the finalization command. Numeric arrays and raw data remain locally available even where Git ignores them; the Git milestone is not a self-contained raw-data bundle.

Preserved negative/debug evidence: initial documentation HTTP403 metadata, 0/255 validator assumption corrected to accept official0/1, original split plus explicit amendment, all non-winning feature metrics. No failed scientific run was replaced by a pretrained model run.

## Reproduction commands

From repository root, using the existing environment:

```sh
# Acquisition and analysis already completed; do NOT rerun merely to finalize.
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t8 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t8_acquire.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t8 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t8_registration.py
# Apply/read the preserved analysis_plan_amendment.json before evaluation.
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t8 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t8_evaluate.py
# Saved-data replay only, no downloads or path search:
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t8 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t8_verify.py
# Report, metrics consolidation and seal only:
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t8 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t8_finalize.py
sha256sum -c outputs/p0_t8/SHA256SUMS
```

The evaluator recreates analysis outputs if explicitly run; use the verifier for non-mutating numeric replay checks (it only rewrites the verification receipt). Preserve the manifest and original evidence before any later replay. Exact feature/analysis hashes and package files are recorded in the frozen manifests; acquisition URLs and byte ranges are in `acquisition.json`.

## Conclusion and limitations

**C — NO_USEFUL_DETERMINISTIC_INK_SIGNAL established under this bounded test.** T6 profile mean improves repeatability but does not demonstrate ink discrimination here. Its confirmation AUC is0.5038 versus center intensity0.5008. No tested feature clears the confirmation/control criteria. BEST_FEATURE=NONE; range_wide is only the highest observed confirmation AUC, not a validated winner.

Small signals remain compatible with the intervals. Limits include one fragment/crop, non-exhaustive negative labels, imperfect original registration, uncertain ink depth, different physical scale from T5, spatial dependence, conservative matching, and coarse sampling of thin strokes. No inference that Direct Thread-Stack geometry failed is warranted.

**LEARNED_INK_MODEL_JUSTIFIED=NO from this experiment's evidence.** Do not assume training will solve the problem. A future, separately authorized study could test richer depth/lateral CT neighborhoods, additional registered fragments, or an official pretrained model with explicit transfer controls. None of those experiments was run here. Stop after T8.
