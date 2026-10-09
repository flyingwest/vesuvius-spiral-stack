# P0.2 validation — F02/F03 identity and control failures

## Verdict

**READY_FOR_P0_3 = YES** — for a bounded experiment using the locally supported
F02/F03 tracks, with explicit missing-data masks and uncertain tails excluded.
No stacking was performed here.

F02 and F03 are **physically consistent local winding-fragment hypotheses at
the available CT resolution**. Reviewing every tracked observation, competing
ridges, fixed-station displacement, order checks and direct two-step matches
found **zero suspected cross-Z layer-switch events** in these two tracks.
This is an apparent-consistency result, not a measured zero true-switch rate:
there is no registered physical-layer oracle in this region. It does not certify
the entire original F03 fragment, a whole winding, or all P0.2 tracks.

F03's missing observation at Z=146 is demonstrably a **detector-fragmentation
artifact**, not evidence of absent papyrus. Its pre/post-gap identity is supported
by local CT geometry and texture. The missing geometry remains unobserved in the
frozen track; this validation does not fill it.

The original broader P0.2 result remains **NOT MET**: only 2/5 interior controls
pass its unchanged coverage gate, and 0/14 P0.1b endpoints are newly resolved.
The narrower readiness judgment above answers the requested F02/F03 validation;
it neither changes that experiment's gate nor upgrades its endpoint results.

## Scope, source and method

Hypothesis: the frozen control tracks retain local layer identity, including
F03's skip over one missing detection. This is retrospective validation of the
already inspected run, not an independent held-out experiment or new scorer.

- Dataset: Scroll 1 / PHercParis4, `54keV_7.91um_Scroll1A.zarr/4`.
- Source: `outputs/orientation/ct_level4_zyx.npy`, uint8, ZYX=(128,493,506).
- Remote array: `https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4`.
- Existing P0.2 range: global level-4 Z=137–151 inclusive, local slab indices
  9–23; central slice Z=144/local16. Nominal native indices 2192–2416 at stride16.
- Spacing: isotropic 126.56 µm at this level; native catalog spacing 7.91 µm.
  X is column/right, Y row/down, Z approximately longitudinal. No new downloads.
- Track inputs: `outputs/p0_2/run_02/{fragments,associations,tracks,controls,windows}.json`.
  The run and `control_tracks.png` are unchanged.
- CT displays use the same fixed uint8 0–255 range, nearest-neighbor display,
  with separate unobstructed images. No generative imaging or contrast tuning.
- Numerical validation uses the frozen ridge response, Gaussian difference
  σ=0.6 minus σ=2 pixels, and bilinear sampling. Profile peaks require prominence
  ≥1, response >1, separation ≥1 pixel; detection thresholds are not relaxed.

Radial displacement uses the **same stations in every observation**, taking the
median distance to the earlier proxy center (240,140). It is not a fitted scroll
center or a physical winding index. Using fragment centroids would confound
motion with changing fragment length. F02 has six common stations, u=1–3.5 px;
F03 has only two, u=1–1.5 px. This narrow F03 common support limits the claim.
Pairwise associations have larger overlapping portions, but variable tails are
not thereby validated. Nearest-neighbor separation is the minimum over other
retained window-local proposals of median absolute normal separation at ≥3
shared stations; it is not minimum tissue thickness or guaranteed layer spacing.

## Quantitative findings

| Control / track | Observed Z | Coverage within observed span | Accepted adjacent / gap links | Radial range / net, px | Normal range, px | Nearest retained competitor min / median, px |
| --- | --- | --- | --- | --- | --- | --- |
| F02 / T0077 | 137–150 | 14/14 | 13 / 0 | 2.1253 / +1.9382 | 2.8663 | 1.3476 / 3.1023 |
| F03 / T0243 | 140–145,147–148 | 8/9 | 6 / 1 | 0.7232 / −0.6378 | 1.0606 | 2.5178 / 3.2753 |
| F04 / T0467 | 144–145 | 2/2, but only 2/15 study slices | 1 / 0 | 0.0233 / −0.0233 | 0.0012 | 1.9840 / 2.5724 |
| F05 / T0470 | 144 only | 1/1, no longitudinal test | 0 / 0 | 0 / 0, uninformative | 0, uninformative | 3.3193 / 3.3193 |

F02 radial range is 268.97 µm, normal range 362.75 µm; F03 radial range is
91.53 µm, normal range 134.23 µm. These are displacement measures, not error
estimates. F02's nearest retained competitor minimum is 170.56 µm; F03's is
318.65 µm. Fixed-station radial steps reach 0.8697 px for F02 and 0.6215 px for
F03. Accepted pairwise median normal residuals reach 0.9768 and 0.8519 px,
respectively. Different station sets explain why these are not identical metrics.

| Diagnostic | F02 | F03 |
| --- | ---: | ---: |
| Accepted-edge minimum score (not probability) | 0.7058 | 0.7592 |
| Accepted-edge minimum CT-profile NCC | 0.6531 | 0.6126 |
| Direct two-step shortcuts: accepted / unknown | 9 / 3 | 2 / 2 |
| Conflicting accepted shortcut targets | 0 | 0 |
| Evaluable paired-layer ordering checks / reversals | 18 / 0 | 5 / 0 |
| Suspected cross-Z layer switches from this review | 0 | 0 |
| Gap recoveries in frozen track | 0 | 1 |

Unknown shortcut matches are retained as uncertainty, not counted as agreement.
Full-run audit found one conflicting accepted shortcut elsewhere; the F02/F03
result must not be generalized to all tracks. A lack of detected order reversal
alone cannot prove physical identity.

## F02: all-level inspection and competing layers

`F02_overlays.png` and `F02_raw.png` contain every Z=137–150 observation, plus
the unaccepted Z=151 boundary. The selected ridge occupies a consistent diagonal
band. Nearby bands remain distinct; no selected inter-slice link visibly crosses
a trough to a different retained band.

| Z | Visual interpretation |
| --- | --- |
| 137 | Long diagonal seed ridge between distinct neighboring bands. |
| 138–140 | Upper tail shortens; shared lower portion remains on the same band. |
| 141–144 | Longer support returns; bends and changing tails are not identity switches. |
| 145–146 | Largest outward displacement; weaker accepted match, inspected specifically. |
| 147–148 | Band bends more steeply; shared support remains separate from neighboring candidates. |
| 149–150 | Only lower tail retained; narrowing support raises boundary uncertainty. |
| 151 | No accepted continuation; leave UNKNOWN rather than extending toward another ridge. |

Specific rejection evidence (node suffixes are slice-local, not layer IDs):

- From `z138_control_F02_r02` to Z=139: selected r03 has score 0.8201,
  normal residual 0.3931 px, NCC 0.7436. Competing r00/r01 have residuals
  2.7894/3.0142 px, scores 0.4185/0.4877, and are REJECTED by distance.
  Their positions correspond to visibly different nearby ridge bands. The
  rejection is consistent with CT, not oracle-proven biological correctness.
- From `z145_control_F02_r01` to Z=146: selected r01 has score 0.7058,
  residual 0.9768 px, NCC 0.6531. Upper competitor r00 is REJECTED at 3.7090 px;
  lower r02 remains UNKNOWN at 2.2509 px, failing distance, orientation,
  score, ordering-risk and competition gates. It is correctly *not joined*
  according to the visual ridge distinction. The accepted transition is near
  the score/distance cutoffs and has lower confidence than the middle run.
- Z=150's nearby r00 is only 1.3476 px away on three shared stations. This
  short-support neighbor and the Z=151 failures prevent a confident extension.
  Cumulative radial drift across the track must not be compared to this local
  spacing as if all motion happened in one step.

Confidence: **moderate local physical consistency**, weaker at Z=145→146 and
the final shrinking tail. No apparent jump was identified; global sheet identity
and varying tail geometry remain unverified.

## F03: gap is an extraction failure, not an observed tear

`F03_overlays.png` / `F03_raw.png` show all eight observations, Z=146 and
the outer context Z=139/149. Before and after the gap the selected ridge lies
near Y≈156, between other visible bands. No switch to the retained bands near
Y≈153 or Y≈160 is supported by the accepted link.

The frozen skip `z145_control_F03_r01` → `z147_control_F03_r03` has:

- 7 common stations, overlap 0.875;
- median normal residual 0.13117 px, radial residual 0.08872 px;
- CT-profile NCC 0.71572, orientation mismatch 23.913°;
- score 0.80923; estimated spacing 4.15309 px.

Its large reported margin is partly a singleton-candidate effect after hard
rejection, not calibrated confidence. Orientation is close to the 25° cutoff;
the strong location agreement is not permission to invent missing geometry.

To diagnose Z=146, use the midpoint of matched Z=145/147 coordinates **only as
a probe**, never as a recovered trace. All **7/7** probe stations have an eligible
CT ridge peak within 1 pixel at Z=146; median DoG response at the probes is
**16.0884 CT units**. The corresponding responses at Z=145/147 are
19.2955/19.8807. Thus there is no loss of all local ridge evidence.

Replaying the unchanged detector before length filtering exposes near-ridge
components with **1, 4, 2 and 6 points**, the longest only **2.60243 px**.
All fail its ≥7-point/≥3-pixel requirements. Normal steps of 0.82415,
0.90666 and 0.92134 px exceed the 0.75-pixel connection bound and break the
chain. `F03_gap_discarded_components.png` shows those actual discarded peaks;
it does not connect them. The diagnostic reconstruction exactly reproduces all
retained detector fragments at Z=145,146,147 to numerical tolerance 1e-12.

Conclusion: **same local fragment before/after the gap is physically plausible
and supported; the missing observation is a tracker/extractor artifact**.
No physical tear has been demonstrated. Confidence is moderate for identity,
but the exact in-plane path at Z=146 and changing outer tails remain UNKNOWN.
One identity recovery is supported; zero missing-slice paths were recovered.

## F04/F05: short duration is not evidence of physical termination

Both retain CT structure in neighboring slices. Their durations should not be
interpreted as two-slice and one-slice pieces of papyrus.

**F04:** Z=143→144 offers a nearby detected ridge at residual 1.083 px and
NCC 0.533, narrowly failing distance/texture/score. Z=145→146 offers one at
0.839 px but NCC 0.278 and score 0.600, failing texture/score. Fixed central
geometry probes find eligible peaks within 1 px at 14/17 stations in Z=143 and
13/17 in Z=146. This is primarily unstable geometry/texture association with
changing support, not convincing anatomical disappearance. The neighboring
candidate identity remains uncertain; none is retroactively accepted.

**F05:** Z=143→144 has residual 0.813 px, NCC 0.483, score 0.660; texture,
score and competitor margin fail. At Z=145 two candidates score 0.645 and
0.673, with different curvatures/portions of the ridge. A 15/17 near-peak count
at fixed Z=144 geometry in Z=145 shows continued CT evidence, while split
proposals and weak discrimination prevent a unique match. By Z=146/147 the
fixed geometry is poorly supported, which may reflect displacement and altered
ridge appearance, not necessarily physical loss. Detection fragmentation and
association ambiguity dominate; a true continuation is **unresolved**.

`F04_overlays.png`, `F05_overlays.png` and their raw counterparts cover
Z=141–148. Do not fix these failures by simply loosening all thresholds.

## Unresolved-case accounting and limitations

- **0 suspected switch events** in the 20 accepted F02/F03 links reviewed;
  the true count/rate is UNKNOWN, not “0%”. No registered oracle was used.
- **1 supported gap identity recovery**, F03 145→147; **0 recovered missing
  paths**. The missing Z=146 path remains unavailable for sampling.
- **2 unresolved short-control identities**, F04 and F05 outside their accepted
  support. F06 is outside this focused review and remains unattached in P0.2.
- **3 unaccepted outer-boundary cases** shown in F02/F03 review: F02→151,
  F03←139 and F03→149. These are not physical endings.
- **14 P0.1b endpoints remain UNKNOWN**, unchanged. They are a different
  denominator from control failures and boundary cases; do not sum them as
  independent biological events.
- Common support is short, especially F03; apparent ridges can include
  partial-volume structure or phase effects at 126.56 µm sampling. Neither CT
  brightness nor any score is ink probability or a calibrated identity probability.

## Artifacts, exact commands and verification

All validation outputs are under `outputs/p0_2/validation/`:

- `F02_overlays.png`, `F02_raw.png`: all tracked levels plus stopping boundary.
- `F03_overlays.png`, `F03_raw.png`, `F03_gap_discarded_components.png`.
- `F04_overlays.png`, `F04_raw.png`, `F05_overlays.png`, `F05_raw.png`.
- `F02_profile_diagnostics.png` through `F05_profile_diagnostics.png`:
  fixed-geometry normal profiles, not a stacked reconstruction.
- `metrics.json`, `Fxx_measurements.json`, `observations.csv`:
  station-consistent drift, separation, link diagnostics and track coverage.
- `control_candidate_associations.json`: all touching candidates and failures.
- `fixed_geometry_probes.json`, `F03_gap_detector_diagnosis.json`:
  actual peaks, discarded components and exact detector break causes.
- `validation_verdict.json`, `plan.json`, `preservation_before.json`,
  `verification.json`, `final_verification.json`, `SHA256SUMS`.

Commands executed from the repository root:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_local_ridges.py --verify
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_cross_z.py --verify
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_validate.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_gap_diagnosis.py
.venv/bin/python -B scripts/p0_2_finalize.py
```

Both frozen runs replayed exactly and passed their preservation checks.
Validation initially hashed 124 prior/source files; final verification checks
them again. The validation generator refuses to overwrite its output directory.
Reproduction of figures requires a fresh output location in a copied checkout;
do not delete or regenerate frozen evidence in place. No source volumes,
README.md, AGENTS.md, prior traces or accepted tracks were modified.

Recommended next experiment, **not executed**: a bounded F02/F03 P0.3 test with
explicit fragment/station support and missing-data masks. Do not pretend that
equal row lengths imply equal physical arc coordinates; do not fill F03 Z=146
from intensity interpolation or extend uncertain tails. The readiness verdict
permits that narrow research question, not a claim of P0 PASS. Stop here.
