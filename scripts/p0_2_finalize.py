"""Assemble P0.2 reports and verify frozen evidence; no new tracking."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/p0_2'


def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2)+'\n')


def main():
    a=read('outputs/p0_2/run_01/metrics.json');b=read('outputs/p0_2/run_02/metrics.json')
    v=read('outputs/p0_2/validation/metrics.json');base=OUT/'validation'
    rows=['| Z | Run 01 fragments | Run 02 local proposals | Control-window / endpoint-window proposals |',
          '| --- | ---: | ---: | ---: |']
    for x,y in zip(a['fragment_counts'],b['fragment_counts']):
        assert x['z']==y['z']
        rows.append(f"| {x['z']} | {x['fragments']} | {y['all']} | {y['controls']} / {y['endpoint_windows']} |")
    counts='\n'.join(rows)
    endpoint_rows=['| Endpoint | Strongest support Z (ties with no votes are not winners) | Decision |','| --- | --- | --- |']
    for e in read('outputs/p0_2/run_02/endpoint_reconsideration.json'):
        best=e['candidates'][0]
        endpoint_rows.append(f"| {e['endpoint']} | {best['candidate']}: {best['z_support']} | {e['decision']} |")
    endpoints='\n'.join(endpoint_rows)
    report=f'''# P0.2 — bounded cross-Z fragment tracking

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

{counts}

| Run / candidate class | Accepted | Rejected | Unknown |
| --- | ---: | ---: | ---: |
| Run01 adjacent | 111 | 1727 | 419 |
| Run01 one missing slice | 2 | 1350 | 320 |
| Run02 adjacent | 524 | 818 | 832 |
| Run02 one missing slice | 27 | 242 | 446 |

Run01:1197 fragments,3929 candidates,1084 tracks.
Run02:1417 window-local fragments,2889 candidates,866 tracks.
Run02 track-length histogram (observed slices:count):
`{json.dumps(b['track_length_histogram'],sort_keys=True)}`.
There are25 run02 tracks containing a missing slice and27 accepted gap links.
Lengths count observations, not inclusive span.

Controls: F02/T0077 is observed137–150 (14 slices); F03/T0243 at
140–145,147–148 (8 observations); F04/T0467 at144–145 (2); F05/T0470
at144 (1); F06 has no unique center attachment (0). The ≥3/5 coverage gate
fails at2/5. F04 fails neighboring texture/displacement gates; F05 suffers
split detections and near-competing matches; F06 has competing near-anchor
attachments. Full follow-up diagnosis is in the validation report.

## Endpoint results

{endpoints}

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
'''
    path=ROOT/'docs/p0_2_cross_z_tracking.md'
    assert not path.exists(), 'Do not overwrite an earlier report'
    path.write_text(report)
    verdict={'READY_FOR_P0_3':'YES','scope':'bounded F02/F03 locally supported tracks with missing-data masks; not whole winding or full P0 PASS',
        'physical_consistency':'moderate CT-supported local identity, no registered oracle',
        'suspected_cross_z_layer_switches_F02_F03':0,'true_layer_switch_count':None,
        'accepted_links_reviewed':20,'gap_identity_recoveries':1,'recovered_missing_paths':0,
        'unresolved_short_control_cases':2,'unaccepted_outer_boundary_cases':3,
        'unchanged_P0_1b_unknown_endpoints':14,'original_P0_2_gate_met':False,
        'original_control_pass_count':2,'original_control_required':3,'stacking_performed':False,
        'F03_gap_cause':'independent detector normal-step splitting then minimum component length filter',
        'observations':'same local ridge before/after gap plausible; no missing geometry constructed'}
    save(base/'validation_verdict.json',verdict)
    with (base/'observations.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['control','track','z','node','radius_px','normal_px','arc_length_px','incoming_score','incoming_NCC','nearest_competitor_px'])
        writer.writeheader()
        for fid,m in v.items():
            for r in m['observations_detail']:
                row={k:value for k,value in r.items() if k!='competitors'}
                row.update(control=fid,track=m['track'],nearest_competitor_px=r['competitors'][0]['separation_px'] if r['competitors'] else '')
                writer.writerow(row)
    # Structural checks independent of the track-building routine.
    ns=read('outputs/p0_2/run_02/fragments.json');es=read('outputs/p0_2/run_02/associations.json');ts=read('outputs/p0_2/run_02/tracks.json')
    look={n['id']:n for n in ns};accepted=[e for e in es if e['status']=='ACCEPTED']
    assert len({e['source'] for e in accepted})==len(accepted)
    assert len({e['target'] for e in accepted})==len(accepted)
    for e in accepted:
        assert e['z_target']-e['z_source'] in [1,2]
        assert look[e['source']]['window']==look[e['target']]['window']
    memberships=[n for t in ts for n in t['nodes']]
    assert len(memberships)==len(set(memberships))==len(ns)
    assert len(ns)-len(ts)==len(accepted)
    assert len([e for e in accepted if e['window'] in ['control_F02','control_F03'] and any(e['source']==r['node'] for k in ['F02','F03'] for r in v[k]['observations_detail'])])==20
    preserved=read('outputs/p0_2/validation/preservation_before.json')
    for name,digest in preserved.items():assert sha(ROOT/name)==digest,name
    for run in ['run_01','run_02']:
        for name,digest in read(f'outputs/p0_2/{run}/preservation_before.json').items():assert sha(ROOT/name)==digest,name
    pngs=list(base.glob('*.png'))
    for p in pngs:
        with Image.open(p) as im:im.verify()
    assert len(pngs)==13
    save(base/'final_verification.json',{'prior_source_files_hash_preserved':len(preserved),
        'both_runs_exact_replay':'PASS (see run verification.json)',
        'node_partition_one_to_one_and_forward_z':'PASS','accepted_links':len(accepted),
        'F02_F03_links_reviewed':20,'validation_pngs_readable':len(pngs),
        'all_prior_run_preservation_inventories':'PASS','README_AGENTS':'unchanged',
        'report_sha256':sha(ROOT/'docs/p0_2_validation.md'),
        'scripts':{str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'scripts').glob('p0_2*.py')}})
    for folder in [base,OUT]:
        files=sorted(p for p in folder.rglob('*') if p.is_file() and p.name!='SHA256SUMS')
        (folder/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(folder)}\n' for p in files))
    print(json.dumps({'validation':verdict,'checks':read('outputs/p0_2/validation/final_verification.json')},indent=2))


if __name__=='__main__':main()
