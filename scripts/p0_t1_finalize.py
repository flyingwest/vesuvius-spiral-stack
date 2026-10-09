"""Explicit post-run visual assessment; preserve automatic evidence byte-for-byte."""
import csv, json, hashlib
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t1'
metrics=json.loads((OUT/'metrics.json').read_text())
# Visual findings, not registered manual ground truth. Entire suffix becomes unknown.
FLAGS={146:(0,'Possible upward layer transition at (269,109); weak ridge contact does not establish identity.'),
       160:(2,'Possible layer transition in junction cluster (299,109)-(303,115); local loop/normal displacement in merged ridges.')}
NOTES={128:'Locally plausible bent ridge and short masked bridge; terminal alternatives not established.',
       144:'Short contour bordering central dark cavity; bridge plausible but coarse CT cannot prove physical identity.',
       146:'Seed lies farther from reference than Z144; accepted graph junction is suspect, not a verified bend.',
       160:'Longest automatic candidate; substantial outward displacement but junction cluster prevents certifying its entire identity.',
       170:'Short interior fragment ends in crowded weak ridge network; detection/scoring limitation, not proven physical tear.'}
volume=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')
review=[];diagnostics=[]
for m in metrics:
    z=m['z'];g=json.loads((OUT/f'slice_{z}_graph.json').read_text());a=np.load(OUT/f'slice_{z}_samples.npz')
    pieces=list(map(np.array,g['pieces']));lengths=[float(np.linalg.norm(np.diff(p,axis=0),axis=1).sum()) for p in pieces]
    row=a['ct'].copy();mask=a['missing'].copy();prefix=sum(lengths);flag=None
    if z in FLAGS:
        i,flag=FLAGS[z];cut=sum(lengths[:i+1])+sum(g['gap_lengths_px'][:i])
        mask|=a['s_px']>cut;row[mask]=np.nan;prefix=sum(lengths[:i+1])
    np.savez_compressed(OUT/f'slice_{z}_audited_samples.npz',s_um=a['s_um'],ct=row,missing=mask,
                        automatic_ct=a['ct'],visual_flag=np.array(flag or 'none identified; not oracle-confirmed'))
    Image.fromarray(np.nan_to_num(row,nan=0).round().astype('uint8')[None,:]).save(OUT/f'slice_{z}_audited_row.png')
    Image.fromarray(mask.astype('uint8')[None,:]*255).save(OUT/f'slice_{z}_audited_mask.png')
    r=dict(m);r.update(suspected_wrong_layer_jumps=int(z in FLAGS),confidence='low' if z in FLAGS or z==170 else 'moderate local geometry; identity unverified',
                      qualitative_notes=NOTES[z],visual_flag=flag,pre_flag_visible_arc_px=prefix,pre_flag_visible_arc_mm=prefix*.12656,
                      audited_mask_fraction=float(mask.mean()),true_wrong_layer_jumps=None)
    review.append(r)
    lens=np.array([np.linalg.norm(np.diff(np.array(p),axis=0),axis=1).sum() for p in g['chains']])
    terminal=g['decisions'][-1];reasons={}
    for c in terminal['candidates']:
        for reason in c['rejected_by']:reasons[reason]=reasons.get(reason,0)+1
    pts=np.concatenate(pieces);_,counts=np.unique(pts,axis=0,return_counts=True)
    diagnostics.append({'z':z,'chain_length_median_px':float(np.median(lens)),'chains_shorter_than_3px':int((lens<3).sum()),
                        'total_chains':len(lens),'terminal_reject_reasons_nonexclusive':reasons,
                        'terminal_viable':sum(not c['rejected_by'] for c in terminal['candidates']),
                        'vertices_repeated_including_shared_junctions':int((counts>1).sum())})
    # Independently assert no automatic sample lies strictly inside any gap.
    start=0.
    for i,L in enumerate(lengths[:-1]):
        start+=L;gap=g['gap_lengths_px'][i]
        inside=(a['s_px']>start+1e-9)&(a['s_px']<start+gap-1e-9)
        assert np.all(a['missing'][inside]);start+=gap
    assert Image.open(OUT/f'slice_{z}_ct.png').size==(506,493)
    assert np.array_equal(np.asarray(Image.open(OUT/f'slice_{z}_ct.png')),volume[z-128])
    assert not np.isfinite(row[mask]).any()

(OUT/'visual_assessment.json').write_text(json.dumps(review,indent=2)+'\n')
(OUT/'detection_diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
with (OUT/'audited_measurements.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(review[0]));w.writeheader();w.writerows(review)
fig,axs=plt.subplots(5,2,figsize=(13,12))
for r,(left,right) in zip(review,axs):
    z=r['z'];g=json.loads((OUT/f'slice_{z}_graph.json').read_text());p=np.concatenate(list(map(np.array,g['pieces'])))
    left.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255)
    for q in g['pieces']:q=np.array(q);left.plot(*q.T,color='cyan',lw=.8)
    left.scatter(*p[0],c='lime',marker='*',s=40)
    if z in FLAGS:
        ep=g['decisions'][FLAGS[z][0]]['endpoint'];left.scatter(*ep,facecolors='none',edgecolors='red',s=160,lw=2)
    left.set(xlim=(p[:,0].min()-10,p[:,0].max()+10),ylim=(p[:,1].max()+10,p[:,1].min()-10),title=f'Z{z}; red circle = suspected identity transition')
    a=np.load(OUT/f'slice_{z}_audited_samples.npz');cmap=plt.cm.gray.copy();cmap.set_bad('magenta')
    right.imshow(np.ma.masked_invalid(a['ct'][None,:]),cmap=cmap,vmin=0,vmax=255,aspect='auto',interpolation='nearest')
    right.set(title=f'Audited row; magenta = missing/unknown; Δs=63.28 µm',xlabel='sample (display height enlarged only)',yticks=[])
fig.tight_layout();fig.savefig(OUT/'audited_overview.png',dpi=150);plt.close(fig)
table='\n'.join(f"| {r['z']} | ({r['seed_x']:g},{r['seed_y']:g}) | {r['visible_arc_px']:.3f} / {r['visible_arc_mm']:.3f} | {r['gap_arc_lower_bound_px']:.3f} | {r['extent_px']:.3f} | {r['visible_fragments']} / {r['interruptions']} / {r['accepted_reconnections']} | {r['unresolved_endpoints']} | {r['suspected_wrong_layer_jumps']} | {r['radial_progress_px']:.3f} | {r['row_width']} | {r['masked_fraction']:.2%} |" for r in review)
audit_table='\n'.join(f"| {r['z']} | {r['pre_flag_visible_arc_px']:.3f} / {r['pre_flag_visible_arc_mm']:.3f} | {r['audited_mask_fraction']:.2%} | {r['confidence']} | {r['qualitative_notes']} |" for r in review)
report=f'''# P0-T1 — independent 2-D papyrus thread reconstruction

## Conclusion

**PARTIAL local-fragment recovery; the primary long interior-to-exterior thread criterion is NOT established.**
All five slices produce independent ridge graphs, provisional thread paths and fixed-spacing 1-D CT rows.
The longest automatic candidate is Z160: {review[3]['visible_arc_px']:.3f} px ({review[3]['visible_arc_mm']:.3f} mm), with three masked connectors.
However, visual review flags one possible identity transition in Z160 and one in Z146.
The audited rows mask the entire suffix after these flags; the automatic results remain unchanged as negative evidence.
This is not a claim of low physical layer-switch rate: no overlapping oracle establishes that rate.
No complete interior-to-exterior path has been demonstrated, and no stacked image is produced in T1.

## Hypothesis and conceptual reference

An irregular papyrus contour can be followed using local XY ridge geometry, allowing bends, kinks, reversals and explicit gaps,
then sampled by arc length. No ideal spiral, turns, angular coverage, polar ordering, or cross-Z persistence is required.
`spiral.png` was inspected: its individual colored contours are irregular layered paths with sharp bends and local reversals.
It is a conceptual illustration, not registered CT ground truth; none of its coordinates or colors constructs this trace.
The present result is a sparse subset of the desired contour reconstruction, not a recreation of that whole image.

## Input and selection

Exact source: `outputs/orientation/ct_level4_zyx.npy`, uint8, shape `(128,493,506)`, axes ZYX, global Z origin128.
Dataset: Scroll1 / PHercParis4, `54keV_7.91um_Scroll1A.zarr/4`.
Source URL: https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/4
Physical spacing is126.56 µm isotropic (published native7.91 µm ×16); X right, Y down.
Global Z128,144,146,160,170 correspond to local indices0,16,18,32,42.
No new data were downloaded. Full grayscale images retain every original pixel, including support material outside the detection ROI.

Selection was fixed before T1 execution: Z146 relatively traceable in the earlier pilot, Z170 difficult/no eligible polar seed,
Z128/144/160 representative cases spanning the existing region and including the historical central slice144.
These are convenience strata, not objectively ranked median/easy slices. Earlier polar difficulty does not predict contour difficulty:
Z160, not Z146, is the longest T1 candidate. No slices were replaced after results.

## Exact method and parameters

The pure numerical entry point is `process_slice(image)` in `scripts/p0_t1_threads.py`.
It receives one XY array only: no Z number, prior trace, F02/F03, labels, neighboring slice or track.

1. Convert CT to float; difference of Gaussian filters sigma0.6 and2 px, threshold response>3 CT units.
   Restrict detection to half-open X[25,475), Y[25,325), excluding much but not all support. This rectangle is not a papyrus label.
2. Skeletonize the thresholded mask. Form an eight-neighbor graph, removing diagonal shortcuts when an orthogonal intermediary exists.
   Decompose into maximal degree-two chains with shared junction endpoints; preserve every chain and cycle, including tiny/noisy ones.
   Unlike earlier branch-pruned extraction, junction pixels are not discarded. No closing, learned detection or CT interpolation fills gaps.
3. Select the nearest point to fixed interior-region reference(246,112) on a chain of arc>=20 px.
   Orient this chain toward the endpoint farther from approximate interior reference(240,140).
   Both references are conventions, not anatomical measurements; seed movement between slices is independent, not registration.
   Thereafter no radial/angle monotonicity or global spiral objective is imposed. Reversals inside a chain are allowed.
4. At each endpoint enumerate both orientations of unvisited chains with endpoint distance<=8 px and arc>=3 px.
   Tangents are chords over up to five graph edges. Measure tangent agreement, connector alignment to both tangents,
   endpoint normal displacement, local ridge strength, and local normal-profile layer-spacing proxy.
   Spacing uses DoG profile±10 px at0.25 px, peaks prominence2 and minimum1 px apart; median peak separation,
   clipped2–10 px, with4 px fallback. This is not a measured physical sheet thickness.
5. Bounded score:
   `0.30*max(0,tangent_cos) + 0.25*max(0,chord_cos) + 0.20*exp(-distance/8)
   + 0.15*exp(-normal_miss/spacing) + 0.10*clip(mean_first_six_ridge_response/12,0,1)`.
   Require score>=0.72, tangent cosine>=0.5, connector cosines>=0.65 for nonzero gaps,
   normal miss<=0.45spacing and no intervening skeleton pixel beyond1.5 px from both connector endpoints.
   Require winner margin>=0.08 among candidates passing gates. A stored margin1 means only one eligible candidate,
   NOT perfect confidence. Coincident junctions get chord cosine1 and zero normal miss, a significant weakness discussed below.
6. Accept a unique passing candidate; preserve any nonzero connector entirely as masked gap.
   Otherwise terminate UNKNOWN. Forbid chain reuse, cap200 chains (not reached).
   Connected branch transitions are not counted as missing gaps. Adjacent chain pieces sharing a vertex belong to one visible fragment.
   Colored overlay legend entries identify graph pieces, so there may be more colored pieces than visible fragments.

There is no separate curvature-field term, texture model, angular objective or cross-Z objective.
Short-range tangent history is the continuity proxy, not a physical identity proof.
Kinks inside chains are unrestricted; sharp junction turns can fail local tangent gates even if physically real.
This is a known sensitivity limitation, not a rule that papyrus must resemble a mathematical spiral.

## Arc length and sampling

Visible arc is cumulative Euclidean length of the pixel-center polyline, without smoothing/resizing.
The lattice can overestimate smooth physical length. Each inferred gap reserves its straight chord length as a lower bound,
not the unknown physical length of torn/folded papyrus. Total extent is visible length plus these placeholders.
Uniform Δs=0.5 px=63.28 µm. Linear interpolation in polyline arc length supplies XY; bilinear interpolation samples ORIGINAL CT into float64.
Samples strictly inside connectors are NaN and explicitly masked. One terminal masked sample marks open-ended UNKNOWN;
it does not estimate the unobserved remainder. Raw PNGs are one pixel high, rounded to uint8 using fixed0–255 display, gaps black.
Separate masks distinguish valid black pixels from missing samples. No contrast enhancement, text prior or ink probability is used.

## Automatic measurements (preserved, provisional)

Lengths/progress are level4 pixels unless marked mm. Interruptions count accepted nonzero gaps plus terminal UNKNOWN.
Suspected jumps are post-run visual flags, not oracle-confirmed errors. Each slice has one terminal unresolved endpoint;
flagged slices additionally have an internal identity ambiguity with its downstream suffix masked in audited outputs.

| Z | Seed XY | Visible arc px /mm | Missing lower-bound px | Total extent px | Fragments /interruptions /reconnections | Terminal UNKNOWN | Suspected jump zones | End-minus-start radius px | Row samples | Automatic masked |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
{table}

Radial progress is only a descriptive Euclidean distance change from(240,140), not a tracing objective or sheet topology metric.
It can grow while following an off-center layer and therefore does not prove interior-to-exterior unwinding.
Seeds range from near the central cavity to denser upper interior bands. No trace reaches the exterior.
Small automatic masked fractions refer only to short attempted domains; they emphatically do not indicate almost complete papyrus recovery.

## Visual audit and conservative diagnostic rows

The five full trace overlays were inspected. Raw CT versus graph decision panels were also inspected for Z144/146/160.
All decisions have preserved panels and exact scores in `decision_audit.json`.
Post-run assessment is explicitly exploratory and not a manual segmentation oracle; it does not retune or retrace the experiment.

| Z | Visible prefix before first flagged transition px /mm | Audited masked | Confidence | Notes |
| --- | --- | --- | --- | --- |
{audit_table}

Z146 decision0 at(269,109) accepts a branch heading toward(273,106), score0.823. The CT ridge contact is weak and the
path may move onto a neighboring band. This is a **suspected**, not confirmed, wrong-layer transition; a real kink is also possible.
Z160 decisions2–4 traverse the crowded junction cluster(299,109)–(303,115), with a small local loop and normal displacement
among tightly packed bands. Treat this as one suspected transition zone, not three confirmed jumps.
The apparently high scores there (roughly0.95–0.99) expose how coincident junctions make proximity/normal terms uninformative.
Z160's longer suffix therefore remains useful proposed geometry but cannot support a claim of a single validated physical thread.

At Z160's terminal(340,139), two passing candidates differ by only0.037, correctly triggering the algorithm's UNKNOWN rule.
The other four terminal decisions have no passing candidate. This establishes rule behavior, not oracle correctness of STOP.
No interruption is called a physical tear: detector fragmentation, merged layers, low contrast and real interruption remain confounded.
True jump count, false-link rate, oracle-correct continuation and missed-continuation count remain UNKNOWN.

## Failure diagnosis, not conservative-success rhetoric

The mask skeleton contains thousands of tiny chains and branch contacts. Many CT ridges visible to the eye are fragmented or
fused at this coarse126.56 µm sampling. Chains<3 px are preserved but not continuation candidates; the >=20 px seed filter
also moves seeds away from the fixed reference. `detection_diagnostics.json` quantifies fragmentation and terminal reject reasons.
Four slices stop after short arcs, so this is not a general long-thread success. Local scoring cannot certify identity at merged
nodes; spacing/normal gates are especially weak when endpoints coincide. Conversely tangent gates can reject a genuine kink.
The primary bottleneck is within-slice ridge graph fidelity and junction interpretation, not cross-Z persistence or lack of full turns.
Parameters were not adjusted after seeing these results; all automatic paths and adverse cases remain preserved.

## Artifacts and exact commands

Under `outputs/p0_t1/`, for each Z128/144/146/160/170:

- `slice_<z>_ct.png`: exact493×506 grayscale CT.
- `slice_<z>_ridge_candidates.png`: whole-ROI ridge graph visualization.
- `slice_<z>_thread_trace.png`: provisional automatic thread and seed.
- `slice_<z>_thread_trace_with_gaps.png`: fragments, dashed masked gaps, rejected alternatives, terminal endpoint.
- `slice_<z>_unwrapped_row.png` and `_row_mask.png`: original automatic row, height1; provisional, not corrected by visual audit.
- `slice_<z>_candidates.npz`: lossless response/skeleton; `_graph.json`: every chain, selected pieces, candidate scores and rejection reasons.
- `slice_<z>_samples.npz`: lossless s in pixels/µm, XY coordinates, CT samples, missing mask.
- `slice_<z>_decision_audit.png`: raw CT and candidate graph at every decision.
- `slice_<z>_audited_row.png`, `_audited_mask.png`, `_audited_samples.npz`: conservative reviewed diagnostics with suspect suffixes missing.

Global files: `plan.json`, `source_identity.json`, `preservation.json`, `metrics.json`, `measurements.csv`,
`visual_assessment.json`, `audited_measurements.csv`, `decision_audit.json`, `detection_diagnostics.json`,
`audited_overview.png`, `verification.json`, `SHA256SUMS`.

Executed from repository root using the existing .venv, with no dependency installation:

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t1_threads.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t1_audit.py
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t1_threads.py --verify
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t1_finalize.py
```

Generation refuses an existing evidence directory. Reverse-order replay reproduces all five paths/decisions and numeric samples exactly.
Independent four-corner bilinear calculation agrees within1e-10. PNG decoding, raw CT pixel identity, unsampled gap interiors,
and all prior artifact hashes are checked. README.md and AGENTS.md are unchanged.
The plan is a frozen exploratory configuration, not an independently reviewed or committed preregistration.

## Next recommendation, not executed

Evaluate same-slice ridge separation and junction discrimination on a small annotated set of bends/contacts, with current failing
decisions preserved as controls. If coarse sampling cannot separate layers, a bounded higher-resolution XY crop is a data-resolution
test, not a reason to introduce cross-Z tracking. Do not scale or stack these candidate rows as a validated single sheet yet.
No mesh, OCR, ink inference, learned reconstruction, language prior or neighboring-slice rescue was used. **Stop after P0-T1.**
'''
(ROOT/'docs/p0_t1_papyrus_thread_trace.md').write_text(report)
manifest=[]
for p in sorted(OUT.rglob('*')):
    if p.is_file() and p.name!='SHA256SUMS':manifest.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}')
(OUT/'SHA256SUMS').write_text('\n'.join(manifest)+'\n')
print(json.dumps({'visual_assessment':[{k:r[k] for k in ['z','pre_flag_visible_arc_px','suspected_wrong_layer_jumps','audited_mask_fraction']} for r in review],
                  'checks':'exact CT pixels, explicit gaps and audited masks PASS','report':'docs/p0_t1_papyrus_thread_trace.md'},indent=2))
