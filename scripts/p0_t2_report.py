"""Publish checked T2 outputs and measured report without altering either evidence run."""
import csv, json, shutil, sys, importlib.metadata
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import p0_t2_recovery as m

ROOT=m.ROOT;OUT=ROOT/'outputs/p0_t2';RUN=OUT/'run_02'
metrics=json.loads((RUN/'metrics.json').read_text());verification=json.loads((OUT/'verification.json').read_text())
checks={r['z']:r for r in verification['slices']}
for z in m.ZS:
    for suffix in ['candidate_paths.png','best_continuous_thread.png','provenance.png','unwrapped_geometry.png','unwrapped_signal.png','samples.npz']:
        name=f'slice_{z}_{suffix}';assert not (OUT/name).exists();shutil.copyfile(RUN/name,OUT/name)
        assert m.sha(OUT/name)==m.sha(RUN/name)

# Summary and precise geometric-crossing overlays are diagnostics, not retracing or downstream masking.
fig,axs=plt.subplots(2,3,figsize=(18,12));volume=np.load(ROOT/m.PLAN['source'],mmap_mode='r')
for r,ax in zip(metrics,axs.flat):
    z=r['z'];a=np.load(OUT/f'slice_{z}_samples.npz');verts=a['vertex_xy'];offset=a['segment_offsets'];labels=a['segment_provenance']
    seg=[verts[offset[i]:offset[i+1]] for i in range(len(labels))]
    ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255)
    ax.add_collection(LineCollection(seg,colors=[m.COLORS[l] for l in labels],linewidths=.8))
    ax.scatter(*r['seed_xy'],marker='*',c='lime',s=50);ax.add_patch(Rectangle((430,150),30,60,fill=False,edgecolor='lime'))
    cross=checks[z]['candidate_checks'][0]['crossing_xy']
    for p in cross:ax.scatter(*p,s=90,facecolors='none',edgecolors='red',lw=1.5)
    ax.set(xlim=(150,475),ylim=(320,25),title=f"Z{z}: {r['observed_graph_arc_px']:.0f}px observed, {r['inferred_arc_px']:.0f}px inferred\n{r['observed_gain_vs_T1_raw']:.1f}× T1; {len(cross)} proper crossings")
axs.flat[-1].axis('off');axs.flat[-1].text(.02,.85,'Cyan: observed ridge geometry\nMagenta: inferred, signal masked\nOrange: local low confidence\nGreen: seed / exterior target\nRed circle: post-run crossing flag\n\nNo cross-Z tracking or stacking.\nNo downstream masking at uncertain links.',fontsize=13,va='top')
fig.tight_layout();fig.savefig(OUT/'maximum_recovery_overview.png',dpi=170);plt.close(fig)

summaries=[]
for r in metrics:
    c=checks[r['z']];s={k:v for k,v in r.items() if k!='transitions'}
    s.update(proper_self_crossings=c['candidate_checks'][0]['proper_self_crossings'],T1_raw_valid_samples=c['T1_raw_valid_samples'],
             sample_utilization_gain_vs_T1=c['sample_utilization_gain_vs_T1'],observed_arc_mm=r['observed_graph_arc_px']*.12656,
             inferred_arc_mm=r['inferred_arc_px']*.12656,usable_sample_fraction=r['usable_ct_sample_count']/sum(r['sample_counts_by_label']))
    summaries.append(s)
m.dump(OUT/'summary.json',summaries)
with (OUT/'measurements.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(summaries[0]));w.writeheader();w.writerows(summaries)
m.dump(OUT/'source_manifest.json',{'source':m.PLAN['source'],'source_sha256':m.sha(ROOT/m.PLAN['source']),
 'scripts':{str(p.relative_to(ROOT)):m.sha(p) for p in ROOT.glob('scripts/p0_t2*.py')},
 'python':sys.version,'versions':{n:importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib','pillow']},
 'selected_run':'run_02','initial_run':'run_01','publication':'root slice outputs byte-identical to run_02'})

table='\n'.join(f"| {r['z']} | {r['total_arc_px']:.2f} / {r['total_arc_mm']:.2f} | {r['observed_graph_arc_px']:.2f} | {r['inferred_arc_px']:.2f} | {r['observed_graph_percent']:.2f}% | {r['fragments']} | {r['bridges']} | {r['longest_bridge_px']:.2f} / {r['median_bridge_px']:.2f} | {r['low_confidence_transitions']} |" for r in summaries)
comparison='\n'.join(f"| {r['z']} | {r['T1_raw_observed_arc_px']:.2f} | {r['T1_audited_observed_arc_px']:.2f} | {r['observed_gain_vs_T1_raw']:.2f}× | {r['T1_raw_valid_samples']} → {r['usable_ct_sample_count']} ({r['sample_utilization_gain_vs_T1']:.2f}×) | {r['T1_radial_progress_px']:.2f} → {r['radial_progress_px']:.2f} | {r['complete_count']} | {r['proper_self_crossings']} |" for r in summaries)
positions='\n'.join(f"| {r['z']} | {r['seed_xy']} | {r['end_xy']} | {r['target_distance_reduction_px']:.2f} | YES | {sum(r['sample_counts_by_label'])} | {r['sample_counts_by_label']} |" for r in summaries)
alts=[]
for z in m.ZS:
    top=json.loads((RUN/f'slice_{z}_top_paths.json').read_text())
    for p,c in zip(top,checks[z]['candidate_checks']):alts.append(f"| {z} | {p['rank']} | {p['score']:.3f} | {len(p['path'])} | {c['proper_self_crossings']} |")
logical=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
report=f'''# P0-T2 — maximum-recovery continuous papyrus thread

## Verdict and scope

**Exploratory P0-T2 SUCCESS for maximum-recovery geometric construction, with explicit limitations.**
All five independent slices obtain a continuous inferred line from the exact T1 seed to a designated right peripheral papyrus region.
Observed graph-contour length grows **10.89–70.15×** relative to the original, unmasked T1 traces; usable CT sample counts grow **10.50–66.47×**.
Z144 and Z170 provide the clearest criterion examples: exterior reached, much greater real ridge-supported recovery,
and no proper geometric self-crossing in any of their three reported candidate paths.
Selected Z128/Z146/Z160 paths retain respectively1/2/1 proper crossing defects, clearly recorded rather than silently hidden.

This verdict is about the requested **continuous inferred maximum-recovery representation**, not proof that all fragments belong
to one original physical sheet. The CT-backed portions visibly follow irregular layered bands. Inferred links remain local and bounded,
not a direct arbitrary interior-to-exterior jump. Many are uncertain or appear to connect neighboring bands; this uncertainty is localized,
not used to discard subsequent usable CT. No minimum observed percentage was imposed. Measured graph support is85.83–88.18%.
These percentages are detector-based support proxies, not manual papyrus-label accuracy or confidence probabilities.

## Inputs, coordinates and preservation

- PHercParis4 / Scroll1, `54keV_7.91um_Scroll1A.zarr/4`.
- Source reference: https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4
- Local source `outputs/orientation/ct_level4_zyx.npy`, uint8, shape(128,493,506), ZYX; global Z origin128.
- Same five slices as T1: global128/144/146/160/170, local0/16/18/32/42. X right, Y down.
- Isotropic spacing126.56 µm; physical arc sampling Δs=0.5 px=63.28 µm.
- Reused unchanged T1 graph chains, DoG ridge maps, skeletons, exact seeds, original traces, scores and visual flags.
  T1 score thresholds are NOT continuation gates here; their decisions remain comparison/negative evidence.
  The two T1 suspect-zone locations flag nearby T2 transitions low-confidence without masking later segments.
- No downloads, new labels, cross-Z evidence, learned model, mesh, ideal spiral, OCR, ink inference or language prior.
  `spiral.png` remains a conceptual irregular-contour illustration, never registered truth.

All prior artifacts, README.md and AGENTS.md passed hash preservation. Both T2 runs and all rejected alternatives are retained.

## Interior/exterior convention

Start at the exact T1 seed, exploring both directions of its original contour chain. These are interior-region seeds,
not a claim to locate the anatomical innermost sheet; some sit well above the central cavity.
Exterior target is the fixed XY box **X430–460, Y150–210**, on the right peripheral layered papyrus.
This is an exterior REGION, not the exact outermost material boundary or a full cross-section segmentation.
Target choice was fixed before run01 using the already inspected CT context.
Endpoint distance to this box supplies a soft objective term; individual steps may move inward, leftward or away from the target.
There is no polar coordinate model, turn count or angular-coverage objective.

## Disclosed correction and negative evidence

`run_01/` used all T1 chains in its broad rectangular ROI. Z128 visibly followed the surrounding holder,
because long non-papyrus ridges also collect length reward. Thus its larger recovery number is contaminated, not credited as papyrus recovery.
All five run01 results were completed and preserved before moving them into that archive.

`run_02/` applies a single manually specified coarse XY papyrus-envelope polygon to all five slices:
`{json.loads((RUN/'plan.json').read_text())['domain_polygon_xy']}`.
Only whole chains and bridges inside it are admitted. This is an explicit user-task implementation judgment based on existing same-slice CT,
not trusted/manual sheet segmentation or cross-Z tracking. It can exclude real peripheral papyrus and does not remove every internal void.
The objective, bridge bound, weights, seeds, target, beam size and depth are unchanged. This correction is exploratory, made after observing
the holder failure; it is not a preregistered blind validation. Root-level required outputs are exact copies from run02.

## Graph and virtual-bridge evidence

T1 used Gaussian(0.6)−Gaussian(2) response>3 CT units, skeletonization and maximal degree-two chains retaining junctions.
T2 retains chains of arc>=3 px plus the seed chain, in both directions. Tiny chains remain in T1 evidence but are not nodes in this search.
Enumerate endpoint pairs within12 px, excluding transitions to the same original chain.
Use chords over up to five endpoint graph edges for local tangent/ridge orientation.
Layer-spacing proxy: normal DoG profile±10 px at0.5 px spacing, peaks prominence2 and minimum1 px apart;
median peak separation clipped2–10 px, fallback4 px. It is not sheet thickness.

For each bridge, record distance d, tangent cosine, alignment of the straight bridge with both endpoint tangents,
normal miss, layer spacing, endpoint turn, required angular bend per pixel, local gap ridge evidence and crossing proxy.
Gap evidence is mean clipped DoG/12 at five stations between20–80% of the connector.
Crossing proxy is mean skeleton occupancy there times min(normal_miss/spacing,3), zero for d<=2.
It penalizes suspected transversal ridge encounters but cannot prove physical sheet crossings.

Bridge penalty:

    2*d + 3*acos(tangent_cos)
    + 2*(acos(outgoing_chord_cos) + acos(incoming_chord_cos))
    + 4*normal_miss/spacing + 3*crossing_proxy
    + 2*(1-gap_ridge_support) [last term only for d>0]

The stored curvature proxy is summed required connector bend divided by max(d,1); the cost uses its integrated bend, not a fitted smooth curve.
Virtual bridges are straight polylines and can have endpoint corners. No fictitious CT is sampled on them, even if some gap stations look bright.
Coincident endpoints have no virtual length and get chord cosine1; tangent mismatch is still penalized.
Reject only strong immediate reversals (tangent cosine<−0.5 or either chord cosine<−0.35), excessive distance, domain violations,
and path-history conflicts. There is no T1-style certainty/margin STOP gate.

## Bounded global optimization

Path score:

    total observed graph arc − sum(bridge penalties)
    + 1.5*(initial target distance − endpoint target distance)

This rewards material recovery rather than shortest path. Target progress telescopes and is soft; long detours can be favored.
Search retains256 states per depth, with up to192 first-per-cell/orientation states for coverage diversity
(12 px endpoint cells, four direction quadrants), then fills by score.
Expand the12 lowest-penalty outgoing candidates per endpoint, up to180 transitions after the seed.
Forbid original-chain reuse and shared rasterized history pixels except the immediate joining endpoint;
test bridge interiors against previously occupied pixels. These constraints track history WITHIN one XY slice only.
The raster conflict test is approximate: it missed the proper intersections found by the independent geometric audit.
A whole chain can also contain a loop/touch; physical nonintersection is not guaranteed by the graph.

Complete candidates have final endpoints inside the target box and stop there. The algorithm tests chain endpoints,
not first intersection of the continuous curve with the target box; it may pass through and leave that box earlier.
All complete paths encountered are scored and retained, including close variants. The reported three alternatives are score-ordered,
with Jaccard chain-set overlap<0.9 preferred for diversity; if unavailable, distinct paths fill the slots.
Thus rank1 is the highest-scoring encountered complete path; ranks2/3 are diverse high-scoring alternatives,
not necessarily the second/third entries in the full unfiltered score list. The entire score-sorted list is preserved.

This is a bounded approximate prize-collecting search, NOT an exact global maximum.
Beam pruning, local outgoing cap, node filter, seed, target and domain affect the result.
Z144's selected path uses181 chains and hits the depth budget. Length reward can encourage unnecessary excursions;
the experiment discovers achievable recovery under this budget, not a canonical unwind coordinate.

## Provenance, continuity and CT sampling

Every stored vertex/segment and every uniform arc sample carries a provenance code:

| Code | Meaning | Signal handling |
| --- | --- | --- |
| 0 | no usable ridge-supported data at this sample | NaN |
| 1 | observed detector-supported contour | original CT, bilinear |
| 2 | inferred bridge | always NaN, continuous XY retained |
| 3 | low-confidence observed/transition neighborhood | CT only when ridge-supported |

Flag a transition low-confidence if tangent cosine<0.7, normal miss>0.45spacing, crossing proxy>0,
the best two local eligible penalties differ by<1, or its start lies within5 px of a preserved T1 suspect zone.
This margin is a local ambiguity descriptor, not a probability; global selection need not choose the locally best bridge.
Only the first5 px of the subsequent observed chain (or its entire shorter length) gets code3.
The rest returns to code1 if supported. There is **no downstream masking** caused by uncertain identity.

Construct one continuous polyline by exact shared endpoints and explicit straight bridges.
Parameterize by cumulative Euclidean pixel-center arc. Evaluate XY by linear interpolation at Δs0.5 px;
evaluate ORIGINAL uint8 CT by bilinear interpolation into float64 only for code1/3 where interpolated DoG>3.
Other samples remain NaN. Thus a short diagonal interval in an observed skeleton chain can yield code0 if ridge evidence dips.
Vertex/segment provenance describes the graph geometry; uniform-sample provenance additionally applies that evidence check.
No row resizing or signal enhancement. PNG signal is height1 with fixed0–255 rounded values, missing samples black;
the separate color provenance and lossless numeric labels distinguish missing black from observed black.
The geometry diagnostic plots XY against s and shows a provenance strip and signal breaks.
Physical arc length along an inferred bridge is an assumed chord, not measured torn-paper length; grid polylines also have discretization bias.

## Per-slice measurements — selected run02 paths

Lengths are level4 pixels except the explicit total-mm column. Observed graph arc includes low-confidence observed portions,
but excludes inferred bridges. Graph fragments are incorporated source chains, not independently certified physical pieces.

| Z | Total px /mm | Observed graph px | Inferred px | Graph-supported ratio | Fragments | Bridges | Longest /median bridge px | Low-confidence transitions |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
{table}

The largest bridge is{max(r['longest_bridge_px'] for r in summaries):.3f} px
({max(r['longest_bridge_px'] for r in summaries)*.12656:.3f} mm).
Total paths are long because they include extended observed contours and repeated inward/outward excursions, not a straight radial shortcut.
Neither bridge presence nor imperfect shape is itself a failure.

| Z | Seed XY | Final XY | Target-distance improvement px | Exterior reached | Row samples | Counts [0,1,2,3] |
| --- | --- | --- | --- | --- | --- | --- |
{positions}

## Direct observed-data utilization comparison with T1

The fair primary baseline is T1's original automatic trace, NOT its shorter post-audit masked prefix.
Both are reported. At the same63.28 µm interval, usable CT sample count is a second direct recovery metric.
T2 applies an additional interpolated-ridge-support check; T1 did not, so these counts do not artificially favor T2 by sampling its bridges.
Sample count×0.5 px is an approximate sampled-length equivalent, not an exact polyline integral (endpoint quantization applies).
Exact raw signal samples and separate graph arc make this distinction inspectable.

| Z | T1 original arc px | T1 audited prefix px | T2/T1 original graph arc | Usable samples T1→T2 (gain) | Radial progress T1→T2 px | Complete candidates evaluated | Selected proper crossings |
| --- | --- | --- | --- | --- | --- | --- | --- |
{comparison}

Radial progress is end-minus-start Euclidean distance from the old approximate center(240,140), reported descriptively only.
The exterior-distance improvement and explicit target intersection are stronger measurements of this experiment's stated endpoint convention.
No observation/bridge ratio was used as a pass/fail gate. Percentages do not establish physical sheet identity.

## Top candidate alternatives

All listed alternatives reach the exterior and have exact continuous geometry. The complete archives contain many highly correlated variants,
not independent experiments or independent evidence of physical correctness.

| Z | Diverse rank | Recomputed score | Incorporated chains | Proper crossings |
| --- | --- | --- | --- | --- |
{chr(10).join(alts)}

## Qualitative geometric assessment and remaining defects

All five best-path overlays and the high-penalty bridge panels for Z128/144/160/170 were inspected.
The main observed runs follow the visible bent/flattened layered bands; the paths can double back without being forced into an ideal spiral.
Z144 and Z170 connect many such arcs to the right exterior, with limited local connectors rather than long unsupported traversals.
Their worst-penalty bridge panels still show questionable cross-band links; they are explicit inferred geometry, not hidden observations.
This does not prove a single physical sheet, but it meets the current maximum-recovery purpose without making certainty its objective.

Some inferred links cross visible dark spaces or neighboring-layer directions, and length reward encourages detours and near-parallel returns.
These are distortions/identity risks, not reasons to erase subsequent CT. They are shown in magenta/orange and preserved with scores.
The independent exact line-segment test found selected proper crossings at:
Z128(361.6,111.6); Z146(291,184.25) and(383.3333,151.3333); Z160(363.5385,176.4615).
They are red-circled in `maximum_recovery_overview.png` and stored for every alternative in `verification.json`.
No path is secretly repaired or truncated. These three selected paths have geometric defects despite their increased recovery.
The exact test excludes endpoint touches and collinear overlap, so zero proper crossings is not a comprehensive non-self-contact proof.

Run01 holder contamination is a genuine negative result. Run02's coarse envelope removes the conspicuous support excursion,
but detector-supported CT can still contain noise, merged sheets or non-papyrus structures inside the polygon.
There is no registered oracle: true wrong-layer counts and physical correctness probabilities remain UNKNOWN.
The successful result is a provisional data-recovery geometry, not a validated document coordinate system or a P0 coherent stacked-image PASS.

## Artifacts, reproduction and verification

Required root-level outputs, for each128/144/146/160/170:

- `slice_<z>_candidate_paths.png`: three complete alternatives.
- `slice_<z>_best_continuous_thread.png`: selected line, seed, target, observed/inferred/uncertain portions and unused local alternatives.
- `slice_<z>_provenance.png`: one-row color provenance (black0, cyan1, magenta2, orange3).
- `slice_<z>_unwrapped_geometry.png`: arc-length coordinate/provenance/signal diagnostic.
- `slice_<z>_unwrapped_signal.png`: unenhanced one-row CT with inferred/missing data black.
- `slice_<z>_samples.npz`: lossless arc, XY, CT, sample labels and complete vertex/segment geometry with labels and offsets.
- `slice_<z>_bridge_audit.png`: six highest-penalty inferred links, not cherry-picked easy bridges.

`run_01/` and `run_02/` retain plans, preservation manifests, all virtual edges/features/admission flags,
oriented-chain mappings to T1, complete score-ranked paths, selected alternatives, search-depth logs, metrics and numeric outputs.
`run_02/slice_<z>_domain.npz` records the correction mask and admitted chain orientations.
Root `measurements.csv`, `summary.json`, `verification.json`, `source_manifest.json` and `SHA256SUMS` make the selected result inspectable.
Root required outputs are byte-identical publications of run02, not new tuned outputs.
At report generation, retained artifact logical size was{logical:,} bytes, excluding this report and final checksum file;
most space is explicit alternative-path archives. Network acquisition was0 bytes.

Executed in the repository root:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t2 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t2_recovery.py
mkdir outputs/p0_t2/run_01
mv outputs/p0_t2/*.json outputs/p0_t2/*.npz outputs/p0_t2/*.png outputs/p0_t2/*.csv outputs/p0_t2/run_01/
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t2 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t2_papyrus_domain.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t2 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t2_validate.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t2 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t2_report.py
```

Generation uses exclusive evidence directories; do not run creation commands over this completed result.
Validation checked all15 published candidate geometries for continuous joins, distinct original chains, target reach,
edge admission and independently recomputed objective. All five selected rows passed provenance/NaN rules,
fixed spacing and an independent four-corner bilinear calculation to1e-10; PNGs decoded.
An exact cached-graph Z128 search replay reproduced selected paths, scores and complete-candidate count.
Other searches were not fully rerun; their saved paths/scores/samples were checked independently.
All prior input/output hashes passed, including protected README/AGENTS. Software versions and source SHA-256 are recorded.

## Conclusion and next bounded question (not executed)

Maximum recovery is feasible here in the requested exploratory sense: continuous interior-to-exterior inferred paths can carry
substantially more actual CT-supported data while preserving blank signal on bridges and local uncertainty.
The useful positive evidence is not a claim that every connection is correct, and the defects are not hidden by downstream masking.
The next narrow question would be whether exact geometric crossing penalties and a tighter papyrus/holder discriminator preserve
this utilization gain without the identified excursions. No cross-Z refinement, stacking, learned correction or text inference was performed.
**Stop after P0-T2.**
'''
(ROOT/'docs/p0_t2_maximum_recovery_thread.md').write_text(report)
files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS')
(OUT/'SHA256SUMS').write_text(''.join(f'{m.sha(p)}  {p.relative_to(ROOT)}\n' for p in files))
print(json.dumps({'slices':len(summaries),'required_root_images':25,'gain_range':[min(r['observed_gain_vs_T1_raw'] for r in summaries),max(r['observed_gain_vs_T1_raw'] for r in summaries)],'artifact_files':len(files),'logical_bytes':logical,'report':'docs/p0_t2_maximum_recovery_thread.md'},indent=2))
