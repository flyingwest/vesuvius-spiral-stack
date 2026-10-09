# P0.1b — fragment identity and conservative continuation

## Revised interpretation and result

**Interruption is allowed. An unsupported physical layer switch is the error.**
P0.1b does not require a closed winding, a complete turn, or a continuous trace
across the XY plane. The previous P0.1 report and evidence are preserved rather
than rewritten; this report supersedes any reading of that result in which a
gap alone would be a failure of Direct Spiral-Stack.

Re-examination of the 556-point candidate produced **seven provisional fragments**.
Five are provisionally useful at moderate confidence; two remain low-confidence
hypotheses and are withheld from the useful set. All seven and their evidence
are retained. Across the 14 endpoints, **48 bounded ridge candidates** were scored.
Every endpoint returned **STOP / identity UNKNOWN**. No connecting bridge was
appended, and no layer identity was silently assigned across an excluded interval.

This demonstrates useful fragment extraction and explicit abstention, **not a
validated low false-switch rate**. The same policy stopped at all seven interior
control probes, so it has not demonstrated selective recognition of safe versus
unsafe continuation. No overlapping registered oracle was available; true jumps,
biological interruption detection, and missed continuations remain unmeasured.
These limitations do not turn legitimate physical interruptions into failures.

## Scope, hypothesis, and exact source

Hypothesis: parts of the previous candidate can be retained as locally plausible
spiral fragments while uncertain bridges are left unresolved. This is a
**retrospective, semi-manual exploratory analysis**, not a new independent
preregistered test. Visual splitting uses the already-seen CT and P0.1 evidence;
the construction annotations are not an evaluation oracle.

- Dataset: PHercParis4 / Scroll 1, `54keV_7.91um_Scroll1A.zarr`, 54 keV.
- Remote source array:
  `https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4`.
- Exact CT input: `outputs/p0_1/source_slice.png`, 506 × 493, uint8.
- Z: global level-4 **144**, previous local slab index **16**, nominal native
  index **2304**. This is a downsampled plane, not an independently registered
  native-resolution plane.
- Axes: image columns X, rows Y; X increases rightward, Y downward. Coordinates
  below are global level-4 **(x,y)** pixel-center indices. Spacing is
  **126.56 µm/pixel**, isotropic; native catalog spacing is 7.91 µm.
- Trace input: all 556 entries of `outputs/p0_1/trace_coordinates.csv`, in the
  original order, with no coordinate changes inside retained fragments.
- Same-slice sensitivity inputs: the four preserved P0.1 variant path NPYs and
  `optimizer_path_xy.npy`. No variant was selected to replace the original.
- Provenance: `docs/p0_orientation.md`, `docs/p0_1_single_slice.md`, and
  `outputs/orientation/manifest.json`. The source pixel hash remains
  `93f77b3f25e3d33736dfe838b6a00b92766df326efa43042c9c359bc19cdcef2`.

Only this XY image is analyzed. No CT downloads, neighboring-Z analysis,
multi-slice stacking, OCR, ink inference, language priors, or 3-D meshing occur.
Prior files are read for preservation hashing, not for cross-Z evidence.

## Re-examination: interruptions, ambiguity, and likely switches

An index-labeled view was saved as `outputs/p0_1b/trace_index_audit.png`. The
following **inclusive original CSV index intervals** were excluded before
candidate scoring. Exact endpoint coordinates are in `review_intervals.csv`.
Categories mean visible evidence quality, not proven physical causes.

| Indices | Classification | Observation / reason for splitting |
| --- | --- | --- |
| 0–12 | Untraceable | Open apex with merged ridges; unique continuation is unresolved. |
| 44–63 | Likely layer switch | Abrupt cross-band movement on the upper-right diagonal. |
| 89–104 | Multiple candidates | Right shoulder bend permits competing ridge continuations. |
| 209–262 | Likely layer switch | Lower arc moves between nearly parallel bands. |
| 304–357 | Likely layer switch | Lower-left diagonal crosses a dark trough toward competing bands. |
| 408–426 | Untraceable | Left-upper contact does not resolve ridge identity. |
| 456–479 | Multiple candidates | Upper-left zigzag has several plausible continuations. |
| 509–555 | Untraceable | Upper apex merges and showed strong prior guide-shift sensitivity. |

**No genuine physical tear is conclusively identified from this coarse image.**
The untraceable intervals could reflect a tear, fold, contact, partial-volume
mixing, or insufficient resolution. A dark trough is not proof of missing
papyrus. Likewise, the three likely-switch regions are not three oracle-proven
jump events; each could contain zero, one, or several true switches.

Splitting is analyst-assisted: the implementation applies these disclosed
review intervals and retains the complementary runs of at least 15 points.
It does not claim to have automatically detected eight interruptions. There
is no interpolation across discarded intervals or forced endpoint closure.

## Fragment inventory and confidence

`fragments_overview.png` shows the fragments separately over the original CT.
Each `Fxx_coordinates.csv` contains original indices, XY, Z=144, and cumulative
fragment arc length in pixels and mm. Arc lengths are recomputed from successive
retained CSV coordinates, not inferred from the omitted intervals. Their sum
need not exactly equal the corresponding P0.1 optimizer-polyline lengths because
that polyline had additional vertices between the resampled CSV points.

| Fragment | Original indices | Start (x,y) | End (x,y) | Arc px | Arc mm | Confidence / disposition |
| --- | --- | --- | --- | ---: | ---: | --- |
| F01 | 13–43 | (248.743,116.345) | (260.585,123.130) | 14.612510 | 1.849359 | LOW; withheld hypothesis |
| F02 | 64–88 | (266.394,129.319) | (275.087,136.260) | 11.781709 | 1.491093 | MODERATE; provisionally useful |
| F03 | 105–208 | (281.489,139.908) | (245.936,157.970) | 50.851943 | 6.435822 | MODERATE; provisionally useful |
| F04 | 263–303 | (221.803,157.314) | (202.477,158.033) | 19.858558 | 2.513299 | MODERATE; provisionally useful |
| F05 | 358–407 | (180.551,150.984) | (184.340,140.816) | 24.287275 | 3.073797 | MODERATE; provisionally useful |
| F06 | 427–455 | (193.122,139.048) | (203.526,132.015) | 13.620121 | 1.723762 | MODERATE; provisionally useful |
| F07 | 480–508 | (212.774,129.305) | (224.058,123.903) | 13.734334 | 1.738217 | LOW; withheld hypothesis |

Seven provisional pieces retain 309/556 points and 148.746449 px (18.825351 mm),
approximately 53.67% of the original candidate's polyline length. The five
provisionally useful pieces span 15.237774 mm in total; this is not a joined
length or a claim that they all belong to one physical winding.

| Fragment | Why its start terminates in the backward direction | Why its end terminates forward |
| --- | --- | --- |
| F01 | Untraceable open apex, 0–12 | Suspected cross-band movement, 44–63 |
| F02 | Suspected cross-band movement, 44–63 | Competing shoulder continuations, 89–104 |
| F03 | Competing shoulder continuations, 89–104 | Suspected lower-band switch, 209–262 |
| F04 | Suspected lower-band switch, 209–262 | Suspected lower-left switch, 304–357 |
| F05 | Suspected lower-left switch, 304–357 | Untraceable upper-left contact, 408–426 |
| F06 | Untraceable upper-left contact, 408–426 | Multiple upper-left continuations, 456–479 |
| F07 | Multiple upper-left continuations, 456–479 | Untraceable apex, 509–555 |

These are construction/termination reasons. Endpoint scoring is a separate
check and its rejection reasons are preserved in `endpoints.json` and
`candidate_scores.csv`. Both endpoint decisions for every fragment are STOP.
Implementation-policy correctness was verified; **biological correctness and
oracle-correct continuation are UNKNOWN for both ends of every fragment**.

### Confidence evidence, not confidence probabilities

Confidence combines source-image review with the following same-slice support
diagnostics (`fragment_support.json`). Ridge support means fraction of points
with difference-of-Gaussians response >1 CT unit. Orientation agreement is the
mean absolute cosine between local path tangent and structure-tensor ridge
orientation. Variant displacement uses corresponding original optimizer
stations interpolated to the baseline CSV arc-length parameter; it is not a
nearest-point comparison or a fresh trace.

| Fragment | Ridge support | Orientation agreement | Variant displacement p95 / max, px |
| --- | ---: | ---: | --- |
| F01 | 0.968 | 0.933 | 4.960 / 5.000 |
| F02 | 1.000 | 0.953 | 0.969 / 1.564 |
| F03 | 0.971 | 0.888 | 0.498 / 1.000 |
| F04 | 1.000 | 0.878 | 0.500 / 0.717 |
| F05 | 0.920 | 0.873 | 1.247 / 1.250 |
| F06 | 1.000 | 0.949 | 0.476 / 0.500 |
| F07 | 1.000 | 0.817 | 2.168 / 2.250 |

The >2 px p95 sensitivity warning that withholds F01/F07 is a **post-analysis
exploratory disposition**, not a preregistered pass threshold. The original
seven-fragment extraction, all scores, and all decisions remain unchanged.
MODERATE does not mean oracle-verified, and stable ridge following can still
follow the wrong physical layer or a phase fringe.

## Bounded deterministic continuation method

The analysis plan is saved in `analysis_plan.json`; exact code is
`scripts/p0_1b_fragments.py`. Thresholds were fixed for this scoring run after
the P0.1 evidence was already known. They are heuristics, not learned weights,
not calibrated probabilities, and not a tuned validation result.

### Candidate enumeration at each endpoint

For the start endpoint reverse its fragment so the outgoing tangent points
backward along the original candidate. For the end use its forward ordering.
Fit a quadratic XY curve against cumulative chord length through the last 13
points, with parameter zero at the endpoint. Let `p` be the actual endpoint,
`t` its unit fitted tangent, `n=(-t_y,t_x)`, and signed curvature
`k=(x' y'' - y' x'') / ||p'||^3`.

Predict a location 6 px ahead:
`p_hat = p + 6*t + 18*clip(k,-0.15,0.15)*n`.
Evaluate `R=G_0.6(I)-G_2.0(I)` at `p_hat + o*n`, with
`o=-8,-7.75,...,+8` px. Enumerate **all interior profile peaks** with prominence
at least 1.5, response >1, and peak separation at least 1 px. Candidate IDs
increase with normal offset; score sorting does not change IDs.

The enumeration is bounded, not exhaustive: a true continuation outside this
cross-section, outside ±8 px, at a search boundary, or below the peak criterion
is not proposed. The furthest possible proposal is bounded by the fixed search
and curvature cap. These candidates are local landing ridges, not assertions
that another retained fragment is their continuation. No gap is bridged merely
because it has nearby geometry.

### Evidence components and score

Every component lies in [0,1]. Their equal mean is the total score. All seven
components, raw measurements, and explicit rejection reasons for **every**
candidate are in `candidate_scores.csv` / `endpoints.json`.

| Required evidence | Deterministic computation |
| --- | --- |
| Position / proximity | Euclidean endpoint-to-candidate distance `d`; score `exp(-d/12)`. |
| Tangent consistency | Candidate ridge tangent from the smaller-eigenvalue eigenvector of a 2-D structure tensor; sign aligned to endpoint tangent. For angle `a`, score `(1+cos(a))/2`. |
| Curvature consistency | Candidate curvature from oriented ridge tangents at ±2 px along its tangent; wrapped angle difference divided by 4. Score `exp(-abs(k-k_candidate)/0.08)`. |
| Radial / winding ordering | Relative to proxy center (240,140), compare candidate radius to predicted radius, `delta_r`. With local spacing `L`, score `exp(-0.5*(delta_r/(0.35*L))^2)`. Normal offset is also gated against spacing. This is a radial-order risk proxy, not a measured winding ID. |
| Local layer spacing | On endpoint/candidate normal profiles, ±8 px in 0.25 steps, find positive ridge peaks with prominence 2 and separation 1 px. Use median successive-peak distance. Score `exp(-abs(log(L_candidate/L_endpoint)))`; missing estimates score zero and require STOP. |
| Gradient / ridge orientation | Gaussian sigma 1 image gradients; tensor products smoothed at sigma 1.2. Coherence is `(lambda_max-lambda_min)/(lambda_max+lambda_min+1e-12)`. Score is mean absolute ridge-tangent agreement on a diagnostic connector × sqrt(endpoint coherence × candidate coherence) × positive-ridge support fraction. |
| Local CT texture/intensity | Compare 13-sample cross-normal CT profiles spanning ±3 px. Score `max(0,correlation)*exp(-abs(mean difference)/25)*exp(-abs(std difference)/20)`. This is texture/intensity evidence only, not ink. |

The diagnostic connector is a cubic Hermite curve from endpoint to candidate,
with endpoint tangents scaled by their Euclidean separation and 33 evaluation
stations. Its ridge support uses `R>1`. It is **only an evidence probe**: its
coordinates are saved, but no connector is appended to any fragment.

Only select CONTINUE if the highest-scoring candidate meets **all** gates:
score ≥0.80; score margin ≥0.12 if a runner-up exists; both coherences ≥0.35;
tangent mismatch ≤25°; curvature mismatch ≤0.08 rad/px; endpoint |curvature|
≤0.15; |normal offset| ≤0.4 endpoint spacing; |radial deviation| ≤0.5 spacing;
texture correlation ≥0.6; mean CT difference ≤25; connector ridge support
≥0.90; connector orientation agreement ≥0.85; both spacings available.
Otherwise choose STOP/UNKNOWN. The source image is never modified.

A real gap is allowed to terminate a fragment. This policy's ridge-support gate
makes it especially reluctant to bridge one; it does not classify such a gap
as a reconstruction failure or prove that no same-layer continuation exists.
The center, spacing, gradient and texture proxies can all be confounded by
folds, adjacent layers, phase fringes, and coarse sampling. Geometry alone is
not assigned physical identity.

## All endpoint candidates

All bounded candidates are enumerated here, including poor matches. Scores are
heuristic evidence, not probabilities; the top-ranked candidate is not an
accepted link. IDs are stable within each endpoint. Both endpoints of each
fragment are evaluated independently, including the two withheld fragments.

| Endpoint | Candidate XY and total score (descending) | Decision | Oracle-correct / biological correctness |
| --- | --- | --- | --- |
| F01_start | C02 (245.68, 111.01)=0.564; C01 (241.81, 114.56)=0.388; C03 (250.11, 106.96)=0.381 | STOP | UNKNOWN / UNKNOWN |
| F01_end | C02 (266.29, 125.06)=0.602; C03 (264.98, 128.03)=0.514; C01 (267.71, 121.86)=0.484; C04 (263.26, 131.92)=0.397 | STOP | UNKNOWN / UNKNOWN |
| F02_start | C03 (261.99, 124.82)=0.501; C04 (263.17, 122.62)=0.500; C01 (257.28, 133.64)=0.463; C02 (259.05, 130.33)=0.441 | STOP | UNKNOWN / UNKNOWN |
| F02_end | C02 (277.00, 142.43)=0.615; C01 (281.08, 139.13)=0.408; C03 (272.53, 146.05)=0.316 | STOP | UNKNOWN / UNKNOWN |
| F03_start | C03 (276.02, 137.43)=0.652; C04 (277.29, 134.71)=0.540; C02 (273.38, 143.09)=0.420; C01 (272.01, 146.04)=0.380 | STOP | UNKNOWN / UNKNOWN |
| F03_end | C01 (239.61, 163.03)=0.518; C02 (240.00, 156.79)=0.514; C03 (240.29, 152.05)=0.487 | STOP | UNKNOWN / UNKNOWN |
| F04_start | C02 (227.47, 159.99)=0.541; C03 (226.84, 164.19)=0.411; C01 (228.36, 154.05)=0.351 | STOP | UNKNOWN / UNKNOWN |
| F04_end | C02 (196.79, 160.04)=0.526; C01 (198.02, 165.14)=0.510; C03 (195.55, 154.93)=0.420 | STOP | UNKNOWN / UNKNOWN |
| F05_start | C01 (182.97, 158.07)=0.402; C02 (180.28, 161.98)=0.282 | STOP | UNKNOWN / UNKNOWN |
| F05_end | C01 (188.00, 135.01)=0.648; C03 (190.76, 140.06)=0.587; C04 (192.32, 142.91)=0.556; C02 (189.08, 136.99)=0.535 | STOP | UNKNOWN / UNKNOWN |
| F06_start | C02 (189.70, 145.00)=0.642; C01 (191.38, 147.78)=0.596; C03 (186.60, 139.86)=0.563; C04 (184.80, 136.86)=0.546 | STOP | UNKNOWN / UNKNOWN |
| F06_end | C02 (204.87, 125.67)=0.673; C03 (207.97, 127.78)=0.529; C01 (202.39, 123.99)=0.363 | STOP | UNKNOWN / UNKNOWN |
| F07_start | C02 (206.00, 128.03)=0.630; C01 (208.17, 134.42)=0.560; C03 (204.23, 122.82)=0.329 | STOP | UNKNOWN / UNKNOWN |
| F07_end | C04 (230.79, 124.54)=0.505; C02 (228.16, 118.05)=0.470; C05 (232.19, 128.01)=0.454; C03 (230.13, 122.91)=0.436; C01 (227.04, 115.27)=0.412 | STOP | UNKNOWN / UNKNOWN |

For example, F03's end offers two almost tied routes: C01=0.518 and C02=0.514.
Neither is accepted. F02's end has a larger score margin but still fails curvature
consistency and the overall score gate. Thus STOP is not based only on ties.
The complete per-candidate rejection lists are preserved rather than reduced
to a single generic “low confidence” flag.

## Oracle availability and metric interpretation

The manually labeled cube identified in the scout has native ZYX bounds
`[2256:2512,2512:2768,4816:5072]`. Its nominal level-4 XY footprint is
**x=[301,317), y=[157,173)**. The present trace lies left of x=284, so it does
not overlap that cube. The orientation and P0.1 acquisition did not contain
any other registered manual layer segmentation of this region. This analysis
does not claim that no such annotation exists anywhere in the public archive.

No segmentation was used to construct the fragments or score candidates.
The analyst's construction masks and visual confidence labels are explicitly
not independent truth. Without a registered layer-identity oracle, report:

| Metric | Result / interpretation |
| --- | --- |
| Physically plausible fragments | 7 provisional candidates; 5 provisionally useful, 2 withheld uncertain; not oracle-certified |
| True inherited layer-jump count | **UNKNOWN**; 3 suspicious regions, not 3 confirmed jumps |
| Newly appended continuation links | **0** |
| False newly appended links | **0 by construction** because no links were appended |
| False continuation / layer-switch rate | **UNDEFINED (0 accepted links)**; inherited fragment errors also unknown |
| Correctly detected physical interruptions | **UNKNOWN**; 3 untraceable review regions are not confirmed tears |
| Missed true continuations | **UNKNOWN**; abstention may discard legitimate connections |
| Ambiguous physical-identity endpoints | **14/14** |
| Endpoints with multiple geometric candidates | **14/14** |
| Endpoint CONTINUE / STOP decisions | **0 / 14** |
| Oracle-resolved endpoints | **0** |

`metrics.json` uses null for oracle-dependent truth counts. `evaluation.json`
clarifies that zero *new false links* is a separate operational count from an
unknown truth rate. Do not interpret null as zero or report “0% layer switches.”

### Interior controls expose over-abstention

After the endpoint run, the same unchanged scorer was evaluated at the midpoint
of each retained fragment using only its preceding points as history. These
are plausible retained-path interiors, **not independent true-positive labels**.
All seven also returned STOP. Full measurements are in
`interior_control_probes.json`. This is negative evidence about selectivity:
an all-STOP policy does not establish correct interruption detection, low
missed-continuation rate, or nontrivial continuation accuracy. Thresholds were
not relaxed after this check, and no alternative scorer replaced the result.

## Reproduction, preservation, and artifacts

Run from `/home/wpotech/deployments/vesuvius-spiral-stack`:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p01b .venv/bin/python -B scripts/p0_1b_fragments.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p01b .venv/bin/python -B scripts/p0_1b_audit.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-p01b .venv/bin/python -B scripts/p0_1b_fragments.py --verify
```

The first command was executed once and refuses to overwrite completed metrics.
Use `--verify` for deterministic replay without replacing fragment evidence.
The audit script creates post-analysis support/control/disposition outputs;
it does not change the fragment mask, scores, or decisions. NumPy 2.5.3,
SciPy 1.18.1, Matplotlib 3.11.2, and Pillow 12.3.0 are recorded in evaluation.json.
No package installation or external model was used.

Key files under `outputs/p0_1b/`:

- `analysis_plan.json`, `review_intervals.csv`: disclosed masks and parameters.
- `trace_index_audit.png`, `fragments_overview.png`: source review and split trace.
- `F01_coordinates.csv` through `F07_coordinates.csv`: exact retained fragments.
- `Fxx_start_candidates.png`, `Fxx_end_candidates.png`: all endpoint proposals
  beside an unobstructed source view. Plotted curves are diagnostic, not accepted.
- `fragments.json`, `fragment_support.json`: coordinates, lengths, termination
  reasons, confidence evidence and useful/withheld disposition.
- `endpoints.json`, `candidate_scores.csv`, `endpoint_table.md`: every candidate,
  seven score components, raw evidence, rejection reasons and decisions.
- `metrics.json`, `evaluation.json`, `interior_control_probes.json`: counts,
  truth limitations and the negative selectivity control.
- `preservation_before.json`, `verification.json`, `SHA256SUMS`: preservation
  inventory, numerical replay checks and output hashes.

Verification exactly reproduced all fragment, endpoint and metric JSONs;
checked nonoverlapping CSV subsets of the original 556-point trace, all bounded
score components, and preserved hashes of 64 prior files. README.md, AGENTS.md,
orientation evidence, P0.1 evidence and prior reports remain unchanged.

## Conclusion and next recommended experiment

**Qualified outcome:** useful local fragments can be retained without imposing
global continuity, and uncertainty is now explicit. Five moderate-confidence
fragments are available for later evaluation; no inter-fragment link is accepted.
The scientific low-false-switch claim remains **NOT EVALUABLE**, and continuation
selectivity is **NOT DEMONSTRATED**. Neither is rescued by a visually plausible
overlay or by the zero-new-link count.

Next recommended experiment, not run: obtain independent same-plane layer
identity annotations for selected fragment interiors and endpoint candidates
(or a registered finer-resolution same-plane reference) before evaluating
continuation precision/recall. Include known safe continuations, true tears,
ambiguous contacts and wrong-layer alternatives; keep this retrospective
scorer's failures intact. Neighboring-Z correspondence remains reserved for
P0.2. **Stop after P0.1b.**
