# P0.2 — bounded cross-Z fragment tracking

## Result and hypothesis

A local fragment identity can persist across neighboring XY slices without
building a surface mesh. However, **the original P0.2 success criterion is NOT
MET**: the preserved first run passed 0/5 easy controls; the engineering revision
passed 2/5, below its unchanged requirement of at least 3/5 observed in ≥7 slices.
No P0.1b endpoint is newly resolved. This is a partial tracking demonstration,
not conservative success or P0 PASS. The subsequent focused CT validation is in
[p0_2_validation.md](p0_2_validation.md); its narrower readiness judgment does
not change this result.

Hypothesis: cross-Z fragment identity supplies discriminating evidence for
same-slice ambiguities. Each XY is extracted independently. Track identity
does not itself join any two fragments within one slice.

## Exact input and scope

Scroll 1 / PHercParis4, standardized `54keV_7.91um_Scroll1A.zarr/4`:
`https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4`.
Local file `outputs/orientation/ct_level4_zyx.npy`: shape (128,493,506), uint8,
axes ZYX, slab global Z offset128. Use **Z=137–151 inclusive**, local9–23,
15 planes centered on Z=144. Nominal native Z=2192–2416 at stride16; these are
downsampled planes, not separately registered native slices. Spacing126.56 µm
isotropic. Support thickness1.8984 mm, first-to-last center distance1.77184 mm.
XY columns increase rightward; rows downward; Z is approximately longitudinal.
Only the existing downloaded slab was used; no network acquisition occurred.

Source-plane hashes, byte counts, shapes and nominal native indices are in
`outputs/p0_2/audit/slice_manifest.json`. Z=144 pixels retain SHA256
`93f77b3f25e3d33736dfe838b6a00b92766df326efa43042c9c359bc19cdcef2`.
The scout manual-label cube does not overlap these controls; no trusted/manual
labels construct or evaluate the tracks. Physical false-association rates
remain UNKNOWN. Local brightness is not ink probability.

## Methods, parameters and disclosed revision

`run_01/plan.json` and `scripts/p0_2_cross_z.py` preserve the first method:
DoG σ0.6−σ2, response threshold3, skeletonization in x160:296,y100:174,
remove branch nodes, retain ≥7 vertices and ≥6px arc. No gap closing within XY.
The easy-control failures were examined rather than called success: pruning
removed support near three anchors, and changing fragment length plus
inconsistent texture-profile alignment suppressed even interior matches.

`run_02/plan.json` and `scripts/p0_2_local_ridges.py` preserve an **adaptive
engineering revision after seeing run_01**, not an independent confirmation:
19 fixed windows (5 prior fragment midpoints, 14 endpoints), each independently
sampled on every plane. Window location/orientation use P0.1b; neighboring-Z
geometry never initializes extraction. Thus this is prior-guided inspection,
not blind full-slice discovery. Control u=−4:0.5:4, endpoint u=−4:0.5:10;
normal v=−8:0.25:8 px. Same DoG response; peaks prominence≥1, response>1,
spacing≥1 px. Adjacent-u peaks join only when mutually nearest with ≤0.75px
normal change and ≥0.125px competitor margin. Retain ≥7 points and ≥3px arc.
Do not join ambiguous or missing in-plane evidence. Window overlap means counts
are **window-local proposals, not distinct anatomical layers**.

Inter-slice translation comes from bounded NCC search over ±2px in 0.5px steps,
using Gaussian σ1 CT. All 14 estimated translations are (0,0); this is a coarse
ROI registration, not proof of no local motion. Matching uses common u stations,
normal displacement, overlap, tangent angle, curvature, spacing, normal-profile
CT texture and radial displacement relative to proxy center(240,140).

Run02 weighted score: position0.20, overlap0.15, texture0.35, orientation0.10,
curvature0.05, spacing0.05, radial proxy0.10; components and exact formulas are
saved in associations/code. Hard reject displacement>2.5px, NCC<0.1 or angle>50°.
Acceptance requires ≥5 common stations, displacement≤1px, NCC≥0.55, angle≤25°,
overlap≥0.6, score≥0.70, reciprocal best match and competitor margin≥0.06.
Displacement must also be ≤0.45 local spacing (or ≤0.5px if spacing missing).
Non-hard failures remain UNKNOWN. These are heuristics, not probabilities.

One-missing-slice links only consider nodes without an accepted outgoing/incoming
adjacent link and require score≥0.76, NCC≥0.70 and margin≥0.10. Tracks are
one-to-one node chains, not a general 3-D papyrus representation or triangle mesh.
No final intensity stacking occurs.

Endpoint reconsideration requires a uniquely attached source track and a candidate
landing point supported on the same independently detected fragment across ≥3
consecutive neighboring slices, ≥0.6 supported fraction and ≥2-vote lead.
Central-slice-only support is distinguished from new cross-Z support. No decision
threshold was lowered after this revision failed the control gate.

## Counts and tracks

| Z | Run 01 fragments | Run 02 local proposals | Control-window / endpoint-window proposals |
| --- | ---: | ---: | ---: |
| 137 | 88 | 94 | 18 / 76 |
| 138 | 88 | 88 | 19 / 69 |
| 139 | 73 | 87 | 22 / 65 |
| 140 | 72 | 89 | 16 / 73 |
| 141 | 82 | 98 | 18 / 80 |
| 142 | 79 | 89 | 19 / 70 |
| 143 | 80 | 106 | 19 / 87 |
| 144 | 71 | 89 | 19 / 70 |
| 145 | 79 | 98 | 18 / 80 |
| 146 | 76 | 92 | 19 / 73 |
| 147 | 93 | 94 | 18 / 76 |
| 148 | 81 | 102 | 21 / 81 |
| 149 | 82 | 104 | 16 / 88 |
| 150 | 84 | 94 | 18 / 76 |
| 151 | 69 | 93 | 18 / 75 |

| Run / candidate class | Accepted | Rejected | Unknown |
| --- | ---: | ---: | ---: |
| Run01 adjacent | 111 | 1727 | 419 |
| Run01 one missing slice | 2 | 1350 | 320 |
| Run02 adjacent | 524 | 818 | 832 |
| Run02 one missing slice | 27 | 242 | 446 |

Run01:1197 fragments,3929 candidates,1084 tracks.
Run02:1417 window-local fragments,2889 candidates,866 tracks.
Run02 track-length histogram (observed slices:count):
`{"1": 574, "14": 1, "2": 155, "3": 74, "4": 32, "5": 18, "6": 7, "7": 3, "8": 2}`.
There are25 run02 tracks containing a missing slice and27 accepted gap links.
Lengths count observations, not inclusive span.

Controls: F02/T0077 is observed137–150 (14 slices); F03/T0243 at
140–145,147–148 (8 observations); F04/T0467 at144–145 (2); F05/T0470
at144 (1); F06 has no unique center attachment (0). The ≥3/5 coverage gate
fails at2/5. F04 fails neighboring texture/displacement gates; F05 suffers
split detections and near-competing matches; F06 has competing near-anchor
attachments. Full follow-up diagnosis is in the validation report.

## Endpoint results

| Endpoint | Strongest support Z (ties with no votes are not winners) | Decision |
| --- | --- | --- |
| F01_start | F01_start_C01: [] | UNKNOWN |
| F01_end | F01_end_C01: [] | UNKNOWN |
| F02_start | F02_start_C03: [143, 144] | UNKNOWN |
| F02_end | F02_end_C02: [143, 144] | UNKNOWN |
| F03_start | F03_start_C01: [] | UNKNOWN |
| F03_end | F03_end_C01: [] | UNKNOWN |
| F04_start | F04_start_C01: [] | UNKNOWN |
| F04_end | F04_end_C01: [] | UNKNOWN |
| F05_start | F05_start_C01: [] | UNKNOWN |
| F05_end | F05_end_C02: [146] | UNKNOWN |
| F06_start | F06_start_C01: [] | UNKNOWN |
| F06_end | F06_end_C03: [146] | UNKNOWN |
| F07_start | F07_start_C01: [] | UNKNOWN |
| F07_end | F07_end_C01: [] | UNKNOWN |

All14 remain UNKNOWN. A neighboring long control track does not prove an
endpoint continuation. At most two supporting slices occur for the leading
candidates at F02_start/end; no candidate meets the discriminating persistence
gate. Oracle correctness remains unavailable.

## Disappearance/reappearance and false-association indicators

- F03 T0243 skips146 between145/147: accepted score0.80923, NCC0.71572,
  median normal residual0.13117px. The initial audit found no close retained
  middle detection; validation subsequently shows actual ridge peaks discarded
  by the minimum-component-length criterion. This is not a proven physical tear.
- T0206, F01_start, Z140→142 skips141 despite a retained middle node
  `z141_F01_start_r01` only0.32852px from the halfway prediction. Score0.81970,
  NCC0.72649, residual0.20393px. This is an association gap, not absence of CT.
  See `audit/association_gap_example.png` and its JSON. It is not an endpoint
  continuation resolution.
- Across27 gap links:22 lack a close retained intermediate fragment,5 have one.
  Those categories describe the detector/association, not physical disappearance.
-251 evaluable paired-order checks,0 reversals. Of238 two-hop chains,
  direct shortcut statuses are131 accepted,97 unknown,9 not enumerated,1 rejected.
  **One shortcut accepts a different target**: preserved as a false-association
  warning, not hidden by the order-reversal count. Detailed node IDs and gates
  are in `audit/two_step_consistency.json`.
- True false-association count/rate is UNKNOWN without registered identity truth.
  Reusing acceptance-gate measurements is not independent evidence of accuracy.

## Artifacts and exact commands

Run directories each contain plan, source identity, preservation inventory,
fragments, associations, tracks, controls, endpoint reconsideration and metrics
JSON; full parameters and every accepted/rejected/unknown candidate are retained.
`run_02/overlay_z137.png`, `overlay_z140.png`, `overlay_z144.png`,
`overlay_z148.png`, `overlay_z151.png`, and `control_tracks.png` show representative
planes/identities. `audit/` holds additional consistency diagnostics and closeups.
`validation/` contains the later all-level F02/F03 raw/overlay review.

Executed from repository root (first two generation commands were executed once,
and refuse to overwrite their run directory):

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_cross_z.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_local_ridges.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_audit.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_cross_z.py --verify
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p02 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_2_local_ridges.py --verify
.venv/bin/python -B scripts/p0_2_finalize.py
```

Both runs exactly replayed and passed preservation checks. Existing environment:
Python3.12, NumPy2.5.3, SciPy1.18.1, scikit-image0.26.0, Matplotlib3.11.2.
No new dependencies, learned models, OCR, language priors, ink inference, full
downloads, mesh reconstruction or final stacked image were used. README.md and
AGENTS.md are unchanged. Next recommendation at original run completion was
focused CT validation of successful controls before any stacking; that requested
validation is now documented separately. Stop after P0.2/its validation.
