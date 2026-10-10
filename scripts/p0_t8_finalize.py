"""Document and seal existing T8 evidence without rerunning analysis."""
from p0_t8_acquire import *
import csv
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 m=json.loads((OUT/'feature_metrics.json').read_text());v=json.loads((OUT/'verification.json').read_text());a=json.loads((OUT/'acquisition.json').read_text());d=json.loads((OUT/'decision.json').read_text());r=json.loads((OUT/'registration.json').read_text())
 assert v['status']=='PASS'
 sample=np.load(OUT/'sampled_features.npz');reg=np.load(OUT/'registration_arrays.npz')
 counts={}
 for name,key in [('discovery','discovery'),('confirmation','confirmation')]:
  y=sample['labels'][sample[key]];counts[name]={'positive':int(y.sum()),'negative':int((~y).sum()),'total':len(y)}
 empirical=max(m,key=lambda k:m[k]['confirmation']['auc'])
 summary={'P0_T8_COMPLETE':'YES','P0_T8_RESULT':d['result'],'BEST_FEATURE':d['best_feature'],'LEARNED_INK_MODEL_JUSTIFIED':d['learned_model_justified'],'READY_FOR_NEXT_STEP':'YES','readiness_scope':'bounded follow-up research may be proposed; does not authorize training or a new experiment','counts':counts,'empirical_highest_auc_feature':empirical,'empirical_highest_is_not_established_signal':True,'features':m,'decision':d,'verification':v,'network_payload_bytes':a['network_payload_bytes']}
 dump(OUT/'metrics.json',summary)
 # Preserve original pre-amendment registration diagnostic; add corrected split annotation.
 fig,ax=plt.subplots(figsize=(9,8));ax.imshow(reg['ir'],cmap='gray');ax.axvline(1024,color='cyan',lw=2);ax.text(20,100,'DISCOVERY',color='cyan');ax.text(1100,100,'CONFIRMATION',color='cyan')
 for i,(x,y) in enumerate([(1024,0),(1536,0),(1024,1536),(1536,1536)]):
  ax.add_patch(plt.Rectangle((x,y),512,512,fill=False,edgecolor='yellow'));ax.text(x+20,y+200,f'Region{i}',color='yellow')
 ax.set_title('Final frozen left/right split; yellow=regional checks\nIR shown only for registration, never used as a CT feature');fig.tight_layout();fig.savefig(OUT/'final_sampling_diagnostic.png',dpi=140);plt.close(fig)
 ci=lambda x:'['+', '.join(f'{z:.4f}' for z in x)+']'
 rows=[];controls=[];dist=[]
 for k,z in m.items():
  t=z['confirmation'];c=z['nuisance'];n=z['nearby']
  rows.append(f"| {k} | {t['positive']} / {t['negative']} | {t['auc']:.4f} {ci(t['auc_ci95'])} | {t['average_precision']:.4f} {ci(t['average_precision_ci95'])} | {t['effect_d']:.4f} {ci(t['effect_d_ci95'])} |")
  controls.append(f"| {k} | {n['auc']:.4f} {ci(n['auc_ci95'])} | {c['auc']:.4f} {ci(c['auc_ci95'])} | {c['paired_concordance']:.4f} {ci(c['paired_concordance_ci95'])} | {z['regions_above_chance']}/4 |")
  pq=t['positive_quantiles'];nq=t['negative_quantiles'];dist.append(f"| {k} | {pq[2]:.3f} [{pq[1]:.3f}, {pq[3]:.3f}] | {nq[2]:.3f} [{nq[1]:.3f}, {nq[3]:.3f}] |")
 b=m[empirical]['confirmation'];mean=m['mean_wide']['confirmation'];base=m['center_intensity']['confirmation']
 text=f'''# P0-T8 — Ground-truth signal test

## Final status

P0_T8_COMPLETE = YES

P0_T8_RESULT = C — NO_USEFUL_DETERMINISTIC_INK_SIGNAL

BEST_FEATURE = NONE

LEARNED_INK_MODEL_JUSTIFIED = NO

READY_FOR_NEXT_STEP = YES — the evidence is complete for deciding a bounded follow-up, not an authorization or justification to train a model.

No frozen feature established reproducible ink discrimination in this test. The highest confirmation point estimate was **{empirical}**, ROC AUC {b['auc']:.4f}, 95% spatial-bootstrap CI {ci(b['auc_ci95'])}. This is not evidence that the underlying CT contains no ink information: the experiment is small, single-fragment, fixed-depth, and uncertainty includes weak effects.

## Status audit and scope

The feature computation and plots were already complete on October 9. At the final October 10 audit, no T8 analysis worker was running; the pending saved-data verifier had completed successfully. Reporting, consolidated metrics, sealing and Git publication remained unfinished. This was an incomplete handoff/reporting stage, not an observed stalled CT computation. No claim is made about why the interaction paused. Finalization reused all completed metrics without rerunning feature evaluation.

This remains the deterministic labeled-data experiment. It was **not superseded by an official pretrained model experiment**. No model was trained or invoked, no OCR or language prior was used, and no Direct Thread-Stack geometry was modified. Existing official surface-aligned data were consumed as a registered test fixture; no new mesh or surface reconstruction was built.

## Hypothesis and dataset selection

Test whether the existing T6/T7 profile statistics differ between registered known-ink and labeled non-ink locations. Repeatability across slices is not an ink discriminator by itself.

Data used in this experiment were obtained from the **EduceLab-Scrolls** dataset: Frag3 / PHercParis1Fr34, 54 keV exposed surface, surface-volume UUID `20230301213423`. The official package supplies CT surface layers, infrared photography and ink labels. The dataset paper describes spectral-photograph-to-CT alignment; it does not provide a per-pixel registration-error bound for this crop. [Parsons et al., 2023](https://arxiv.org/abs/2304.02084).

Source directory: {BASE}

The four inspected legacy labeled fragments had layer-32 TIFF sizes 103,571,722 (Frag1), 281,948,145 (Frag2), 79,847,973 (Frag3), and 118,745,433 bytes (Frag4). Frag3 was the smallest of this comparison, not a claimed global minimum over every official dataset. No suitable local ink-label fixture was identified; prior manual segmentation references are not ink truth.

Full surface volume: 65 × 7606 × 5249, axes depth/Y/X, uint16. Metadata records 3.24 µm voxel spacing. Crop: global Y `[2784,4832)`, X `[1600,3648)`, depth layers **28,30,31,32,33,34,36**. Seven layers form an array of shape 7 × 2048 × 2048. Offsets around nominal layer32 are −4,−2,−1,0,+1,+2,+4, corresponding to ±12.96 µm at the outer samples. The nominal XY crop extent is 6.63552 mm; actual surface metric can vary with flattening.

The [official surface-volume tutorial](https://scrollprize.org/tutorial3) documents the normal-sampled layer convention, nominal center layer32, and PPM mapping `pixel(y,x) → original XYZ + surface normal`. These are producer-supplied coordinates, not a thread reconstructed here. Exact ink depth is unverified; layer32 was fixed before test metrics, not optimized.

Acquired payload: **{a['network_payload_bytes']:,} bytes**, below the frozen220,000,000-byte cap. HTTP206 byte ranges fetched just the required uncompressed CT row spans and sampled PPM rows. Full PNG labels/mask/IR and small metadata were acquired. `acquisition.json` records URLs, half-open byte ranges, response headers, payload hashes and exact local files. Sparse TIFF metadata files have large logical sizes but only headers/IFDs stored; acquisition payload and local logical file size are different quantities.

## Registration verification

All label/mask/IR images and TIFF headers agree on7606×5249 XY dimensions. Label and mask encoding is **0/1**, not0/255; the initial validator's encoding assumption failed and was corrected without changing pixels, sampling parameters or feature definitions. Labels are not model predictions. Pixel locations map directly into the official surface volume; no registration shift, rotation or scale was fitted against ink statistics.

The overlay was inspected: supplied ink boundaries follow dark registered-IR strokes; CT and IR share material features but have different contrast. Original PPM normals at usable grid points have unit length to numerical precision. All4,194,304 crop pixels lie within the supplied papyrus mask. Globally48 positive label pixels lie outside the papyrus mask; these are not in the evaluated crop. CT nonzero fraction within the crop mask is {r['ct_nonzero_fraction_inside_mask']:.8f}. IR medians are134 at positives and152 at negatives, a registration sanity check only, never a CT feature.

Registration status: **documented producer registration plus consistent dimensions, coordinate-map checks and visual overlay; suitable for this bounded evaluation, not independently certified subpixel registration**. A two-pixel class-boundary exclusion reduces boundary ambiguity but cannot prove perfect alignment. Label0 means supplied no-ink annotation, not independently chemically verified absence; faint or omitted strokes may contaminate negatives.

## Frozen sampling and amendment

`plan.json` and `PLAN_FROZEN.json` were saved before feature-label metrics. The initial top/bottom split was amended after examining the registration overlay, which revealed two separated ink bands, but **before any CT discrimination metrics**. `analysis_plan_amendment.json` and `ANALYSIS_PLAN_FROZEN.json` preserve the reason and code hashes. This is disclosed label-informed sampling design, not a fully label-blind preregistration. The crop, depths, features, control calipers and inferential rules were unchanged.

Final discovery region: local X<1024. Confirmation region: local X≥1024. Four512×512 confirmation corner windows check regional consistency. The middle1024 rows remain included in overall confirmation but not these four regional checks. `final_sampling_diagnostic.png` is authoritative for the final split. The original `registration_diagnostic.png` deliberately retains its earlier top/bottom annotation as pre-amendment evidence.

Samples use a16-pixel grid starting at local(8,8). Crop-mask distance must exceed16 pixels and the distance to the opposite label class must exceed2 pixels. The distance transform is calculated within the crop; all-true mask edge handling can exclude an extra crop-edge sample, and is reproduced unchanged. Positive counts: discovery {counts['discovery']['positive']}, confirmation {counts['confirmation']['positive']}. Negative counts: discovery {counts['discovery']['negative']}, confirmation {counts['confirmation']['negative']}. Total **{len(sample['labels'])}** saved profiles. Grid pitch is51.84 µm nominally. No neighboring T5 slice or T5 path is used.

Each feature's score direction was chosen from discovery AUC only and then frozen; all ten happened to select larger values for ink. No formula, depth, threshold or weight was fitted. Confirmation prevalence is {base['prevalence']:.6f}; that is the uninformative precision reference, not0.5. Matched cohorts are balanced by construction and have AP reference0.5.

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
{chr(10).join(rows)}

## Distribution summaries

Values are in native CT-derived units; medians and interquartile ranges are shown. Full5/25/50/75/95% quantiles are in `feature_metrics.json` and `metrics.json`.

| Feature | Ink median [Q1,Q3] | Non-ink median [Q1,Q3] |
|---|---|---|
{chr(10).join(dist)}

`feature_distributions.png` shows strongly overlapping distributions. Plot bins span pooled0.5–99.5% quantiles for readability only; **all values** enter metrics. No conclusion is drawn from apparent writing or plot appearance.

## Nearby and nuisance-matched controls

Greedy deterministic nearest-negative matching without replacement: same discovery/confirmation half and same256×256 tile, separation≤128 pixels. This yields1050 confirmation pairs. The stricter match additionally requires surface-normal dot product≥0.98, Gaussian-sigma32 local CT-baseline difference≤0.25 pooled SD, and Gaussian-sigma2 spatial-gradient-strength difference≤0.5 pooled SD, yielding561 pairs. These nuisance descriptors are controls, not extra candidate ink features. Matching changes the target population and may remove genuine material-mediated effects as well as confounding; unadjusted results remain primary.

Median separation is64 pixels for nearby matching and65.96969 for nuisance matching. Minimum normal dot improves from0.83879 to0.98004. Mean absolute baseline difference falls from0.76747 to0.12950 SD; mean absolute gradient-strength difference from0.93899 to0.23761 SD. Full balance is saved in `matching_balance.csv`.

| Feature | Nearby-pool AUC [95% CI] | Nuisance-pool AUC [95% CI] | Within-pair concordance [95% CI] | Regions AUC>0.5 |
|---|---|---|---|---:|
{chr(10).join(controls)}

Some **nominal** control findings are positive and must not be hidden: standard deviation nuisance-pool AUC0.5412 (95%CI0.5079–0.5696), range0.5394 (0.5029–0.5709). Their99.5% intervals include0.5; their within-pair concordance intervals include0.5; and main confirmation intervals include0.5. These are weak exploratory hints, not replicated controlled discrimination sufficient for categoryB. Discovery AUC for both was about0.556, attenuating to about0.525 in confirmation.

## Replay, preservation and evidence

`verification.json` reports PASS:

- {v['prior_files_unchanged']} prior files unchanged, including frozen geometry, rows and predecessor evidence.
- {v['source_files_hash_verified']} acquisition files hash-verified.
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

**C — NO_USEFUL_DETERMINISTIC_INK_SIGNAL established under this bounded test.** T6 profile mean improves repeatability but does not demonstrate ink discrimination here. Its confirmation AUC is{mean['auc']:.4f} versus center intensity{base['auc']:.4f}. No tested feature clears the confirmation/control criteria. BEST_FEATURE=NONE; {empirical} is only the highest observed confirmation AUC, not a validated winner.

Small signals remain compatible with the intervals. Limits include one fragment/crop, non-exhaustive negative labels, imperfect original registration, uncertain ink depth, different physical scale from T5, spatial dependence, conservative matching, and coarse sampling of thin strokes. No inference that Direct Thread-Stack geometry failed is warranted.

**LEARNED_INK_MODEL_JUSTIFIED=NO from this experiment's evidence.** Do not assume training will solve the problem. A future, separately authorized study could test richer depth/lateral CT neighborhoods, additional registered fragments, or an official pretrained model with explicit transfer controls. None of those experiments was run here. Stop after T8.
'''
 report=ROOT/'docs/p0_t8_ground_truth_signal_test.md';report.write_text(text)
 dump(OUT/'finalization.json',{'pre_finalization_state':'B: numerical analysis and replay complete; reporting/sealing incomplete','no_active_T8_worker_at_audit':True,'analysis_rerun':False,'pretrained_model_experiment_run':False,'visual_review':'Distributions strongly overlap; registration overlay aligns official label strokes to IR; final split shown separately.','empirical_highest_feature':empirical,'result':'C','ready_for_next_step':'YES, research decision only','debug_notes':['Initial documentation HTTP403 retained.','Binary labels are0/1; initial0/255 assertion corrected before evaluation.','Original top/bottom annotation retained; final left/right amendment frozen before feature metrics.']})
 # All prior evidence must still be byte-identical before sealing.
 assert all(sha(ROOT/p)==h for p,h in json.loads((OUT/'preservation.json').read_text()).items())
 files=[report]+[p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS']+list((ROOT/'scripts').glob('p0_t8_*.py'))+[p for p in DATA.rglob('*') if p.is_file()]
 entries=[(sha(p),str(p.relative_to(ROOT))) for p in sorted(set(files))]
 (OUT/'SHA256SUMS').write_text(''.join(h+'  '+p+'\n' for h,p in entries))
 assert all(sha(ROOT/p)==h for h,p in entries)
 print(json.dumps({'result':'C','manifest_entries':len(entries),'manifest_verified':True,'counts':counts,'best_point_estimate':empirical,'report':str(report)},indent=2),flush=True)
if __name__=='__main__':main()
