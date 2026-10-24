"""Report and seal completed bounded extension; no scientific recomputation."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from p0_t13_control_extension import R,O,OLD,NEW,ALL,DATA,read,dump,sha,thread_dir

def main():
 assert not (O/'SHA256SUMS').exists();v=read(O/'verification.json');assert v['status']=='PASS';a=read(O/'trusted_overlap_audit.json');p=read(O/'extension_plan.json')
 for name,h in read(O/'preservation.json').items():assert sha(R/name)==h,name
 fig,axs=plt.subplots(4,2,figsize=(16,13))
 for ax,z in zip(axs.flat,NEW):
  image=np.load(O/f'ct_only_inputs/slice_{z}.npy')
  with np.load(thread_dir(z)/'thread.npz') as ar:xy=ar['vertices_xy'];labels=ar['edge_provenance']
  ax.imshow(image,cmap='gray',vmin=0,vmax=255);ax.add_collection(LineCollection(np.stack([xy[:-1],xy[1:]],1),colors=[{1:'cyan',2:'magenta',3:'orange'}[int(k)] for k in labels],linewidths=.7));ax.scatter(*xy[0],s=16,color='lime');ax.set_title(f'New Z{z}; cyan observed / magenta bridge / orange low confidence')
 fig.tight_layout();fig.savefig(O/'new_threads_overview.png',dpi=140);plt.close(fig)
 lines=['|Native Z|Source|Observed samples|Qualifying samples|Qualifying fraction|','|---|---|---:|---:|---:|']
 for row in a['per_slice']:lines.append(f"|{row['z']}|{row['source']}|{row['observed_sample_count']}|{row['qualifying_count']}|{100*row['qualifying_fraction']:.4f}%|")
 windows=['|Five-slice window|Per-slice counts|Minimum|Total shortfall|Pass|','|---|---|---:|---:|---|']
 for w in a['windows']:windows.append(f"|{w['start_z']}–{w['end_z']}|{[w[f'count_{j}'] for j in range(5)]}|{w['minimum_count']}|{w['total_deficit']}|{w['pass']}|")
 further=[]
 for side,b in a['boundary_continuation'].items():further.append(f"- {side}: boundary-to-interior counts {b['boundary_counts_outside_to_inside']}; qualifying boundary run={b['boundary_inward_qualifying_run']}; prospective ≥3-run continuation indicator={b['prospective_promising_rule_at_least3']}.")
 totalobs=sum(r['observed_sample_count'] for r in a['per_slice']);totalnear=sum(r['qualifying_count'] for r in a['per_slice']);best=a['best_window'];nearnew=sum(r['qualifying_count'] for r in a['per_slice'] if r['z'] in NEW)
 report=f'''# P0-T13 control prerequisite — bounded slice-range extension

```text
CONTROL_EXTENSION_COMPLETE = YES
CONTROL_EXTENSION_RESULT = {a['CONTROL_EXTENSION_RESULT']}
TOTAL_CONTROL_SLICES = 15
NEW_SLICES_RECONSTRUCTED = 8
BEST_5_SLICE_WINDOW = {best['start_z']}..{best['end_z']}
BEST_WINDOW_COUNTS = {a['BEST_WINDOW_COUNTS']}
TRUSTED_OVERLAP = {a['TRUSTED_OVERLAP']}
LEAKAGE_AUDIT = PASS
READY_TO_RESUME_P0_T13 = {a['READY_TO_RESUME_P0_T13']}
```

The bounded extension is complete. There are {a['passing_windows']} passing five-slice
windows under the unchanged per-slice gate. The original prerequisite and its
seven threads remain unchanged. No T13 correspondence, ribbon construction,
frame-error experiment, model loading or inference was performed. No farther
slice extension was acquired or reconstructed.

## Hypothesis, prospective limitation and exact inputs

The hypothesis is that nearby independent CT-only threads may extend useful
trusted-sheet overlap. The user specified exactly native Z44449..44463. Only
Z44449..44452 and44460..44463 are new; Z44453..44459 are sealed reuse.

Before acquisition or reconstruction, the extension plan recorded a mathematical
limitation: old endpoints44453 and44459 have zero qualifying points. Every newly
introduced five-slice window in this15-slice interval contains one of those
endpoints. The three wholly old windows already fail. Therefore this prescribed
range cannot yield PASS without changing old evidence or the gate. It remains
useful as the specifically requested bounded measurement of boundary continuation.
This was disclosed before execution, not inferred as an excuse after failure.

Use the exact same native scan as the sealed control prerequisite:

`{p['source']}`

PHercParis4-w02, trusted segment20231031143852, native scan20260411134726.
Native axes are ZYX, transverse XY. The unchanged rectangular native crop is
X[18944,21504), Y[15616,16640). Native Z spacing is2.4µm; exact4×4 XY mean with
half-even uint8 rounding gives9.6µm XY samples and4.8µm arc samples. The15-slice
center span is33.6µm. No Z pooling or row resampling occurs.

320 HTTP206 ranges acquired only the eight new native planes:
**20,971,520 bytes** under the frozen128MB acquisition ceiling. No existing CT
payload or old thread was replaced. The exact range receipts, Content-Range
headers, hashes and pooled input hashes are retained. The crop bounds are inherited
unchanged; no trusted geometry, labels or overlap scores were used to choose new
seeds, paths, branches, bridges or topology.

## Exact frozen reconstruction and independence

The mounted worker is the original
`scripts/p0_t13_control_reconstruct.py`, SHA256
`{p['worker_sha256']}`, verified byte-identical to commit576b562.
All414 prior host code/input/library hash entries also match. The pinned CPU
runtime image and deterministic thread settings are unchanged.

No reconstruction parameter was retuned: DoG(.6,2)>3; maximal skeleton chains;
3px minimum graph-chain arc; both orientations of eight longest eligible chains;
bounded beam128/depth64/outgoing12 maximum-recovery search; ≤12px bridges and
original tangent/chord/penalty rules; source-chain reuse and geometric intersection
guards; provenance1 OBSERVED,2 INFERRED_BRIDGE,3 LOW_CONFIDENCE;0.5px arc sampling;
NaN/mask on inferred bridges. This is a bounded search, not a global optimality or
physical sheet-identity proof.

Each new worker sees one CT input, unchanged worker source, libraries and a fresh
output directory. Docker has no network, read-only root, all capabilities dropped,
and no-new-privileges. Trusted geometry, labels, probability maps, overlap scores
and neighboring CT are not mounted. The entire set of eight new threads was
hash-frozen before post-reconstruction geometry access.

All eight new threads are independently reconstructed. The sealed T3 exhaustive
segment-pair auditor confirms zero proper crossings, non-adjacent contacts and
collinear overlaps in each. Coordinates, provenance, masks, search decisions,
metrics and all candidate links remain available. The original seven threads are
read-only references; they were never retraced or replaced.

## Post-freeze overlap and registration

Evaluation imports the unchanged sealed control helpers for trusted-section
intersection and nearest-point association. The trusted local grid and same-Z
section rule are identical to the first prerequisite. Qualifying samples are
observed/unmasked points at distance≤19.2µm with the inherited finite-frame
eligibility (scan-Z projected tangent norm≥0.25). Distance is to the nearest point
on the piecewise-linear trusted section at the **same native Z**; it is an upper
bound on unrestricted3D nearest-surface distance, not a newly fitted surface.
The CSV also reports distance-only counts so the inherited frame guard is visible.

All registered arrays for the seven old slices reproduce their sealed values
exactly. This is a downstream evaluation check, not a rerun of their reconstruction.
Native registration remains
`XYZ=(18944+4*x+1.5,15616+4*y+1.5,nativeZ)`, with no fitted transform or correction.
Maximum numerical forward/inverse residual is
**{a['coordinate_roundtrip_max_native_voxel']:.17g} native voxel**.
Physical registration uncertainty remains **uncalibrated**, stored as null.
Neither numerical roundtrip nor surface proximity estimates physical accuracy.
No already-available independent calibration was identified in the sealed source;
no calibration project was begun. The old pooling-cell, mesh-sampling and local
normal-variation limitations are retained.

{chr(10).join(lines)}

Across15 slices, {totalnear}/{totalobs} observed samples qualify
({100*totalnear/totalobs:.4f}%). The new slices contribute {nearnew} qualifying
samples. Aggregate counts are descriptive only and never substitute for the
per-slice gate. Proximity does not prove physical sheet identity or cross-slice
correspondence.

## Every consecutive five-slice window and frozen decision

All11 windows are evaluated. Each of the five counts must independently be≥25;
no averaging or total-window substitute is used. The inherited100µm central arc
span safeguard is reported too; it cannot change a count-gate PASS because25
unique4.8µm samples span at least115.2µm. No criterion was loosened.

The prospective selection rule is earliest passing window; if none, largest
minimum per-slice count, then smallest summed shortfall, then greatest total count,
then earliest start. A shortfall is max(25−count,0), reported without using it as
a substitute acceptance rule.

{chr(10).join(windows)}

Best window: **{best['start_z']}..{best['end_z']}**, counts
**{a['BEST_WINDOW_COUNTS']}**; {best['failing_slices']} failing slice(s),
minimum={best['minimum_count']}, total shortfall={best['total_deficit']} point(s).
The threshold near-miss does not justify changing25 to24. No accepted control
input manifest is emitted on FAIL; partial associations are retained as evaluation
evidence only. No T13 matching has been selected.

Visual inspection of `overlap_windows.png` shows isolated lower-boundary support
and a short upper run interrupted at44463. The minimum-count plot has one24-point
window and ten zero-minimum windows. This matches the numeric gate; appearance
is not used as an acceptance criterion or interpreted as ink.

## Boundary evidence and stopping decision

Before evidence, a descriptive continuation indicator was fixed: at least three
consecutive qualifying slices adjoining a new outer boundary would support a
separately authorized test of the remaining one or two needed slices. This is not
an acceptance gate or a prediction about unobserved data.

{chr(10).join(further)}

No boundary result authorizes automatic expansion. Failure here is not a proof
that no useful window exists elsewhere, and a positive continuation indicator
would not prove that the next unobserved slices qualify. The next recommended
externally controlled decision is to assess these boundary counts and the structural
limitation of this prescribed interval before authorizing any further data.
No such additional experiment is run in this milestone.

## Replay, leakage, preservation and reproducibility

All eight new graph/search/decision/thread results replay exactly in reverse order,
and again under syscall-level file tracing. Primary and both replay arrays and JSON
scientific records match exactly; ZIP-container timestamps are not the numeric
replay criterion. Each of24 successful workers has Python-open logs; all eight
instrumented replays also have raw filesystem syscall traces. Successful read paths,
code/input/library/runtime hashes and container mounts are audited. All scientific
input reads are the single assigned CT slice; leakage audit PASS.

Independent native range reassembly and4×4 pooling are exact; four-corner bilinear
sampling, arc/XY replay, provenance/masks, topology and frame arrays verify.
**{v['predecessor_entries_unchanged']} predecessor entries are byte-identical**,
including P0-T1 throughT12, blocked T13, all8fdf4cd evidence and the entire original
control prerequisite. All {len(v['prior_seals_verified'])} prior experiment root
SHA256SUMS manifests pass. The extension seal covers report/code, native payloads,
all numeric evidence, replays and audit logs. Verification PASS refers to integrity
and execution, separately from the scientific extension result.

Required artifacts are in `outputs/p0_t13_control_extension/`:
`extension_plan.json`, `new_independent_threads/`, `thread_metrics.csv`,
`per_slice_overlap.csv`, `five_slice_windows.csv`, `trusted_overlap_audit.json`,
`leakage_audit.json`, `verification.json`, and `SHA256SUMS`.
Additional evidence includes acquired CT inputs, range receipts, frozen manifests,
registered evaluation arrays, two replay directories, file-access inventories,
`overlap_windows.png` and `new_threads_overview.png`. The original control input
manifest remains unchanged. Large source ranges and redundant replay NPZs stay
local but sealed; small CT-only inputs, primary arrays and registered arrays are
published explicitly. `git_artifact_inventory.json` records this split.

Implementation entry points are `scripts/p0_t13_control_extension.py` (setup,
acquire,new_independent_threads,replay,audit_replay), `_evaluate.py`, `_verify.py`
and `_seal.py`. They write only the separate extension paths and refuse existing
stage directories. Do not overwrite sealed outputs when reproducing; use isolated
fresh destinations. Read-only checksum verification from repository root:

```sh
sha256sum -c outputs/p0_t13_control_extension/SHA256SUMS
```

The milestone is committed and pushed separately to origin/main without force.
Commit/remote identity is reported in the final execution response. Stop after
this extension; P0-T13 remains unexecuted.
'''
 reportpath=R/'docs/p0_t13_control_thread_extension.md';reportpath.write_text(report)
 files=[q for q in O.rglob('*') if q.is_file()]+list((R/'scripts').glob('p0_t13_control_extension*.py'))+[reportpath]+list(DATA.glob('*'))
 def publish(q):return not str(q.relative_to(R)).startswith('data/') and not(q.suffix=='.npz' and any(k in q.parts for k in ['replay','audit_replay']))
 dump(O/'git_artifact_inventory.json',[{'path':str(q.relative_to(R)),'bytes':q.stat().st_size,'planned_git':publish(q)} for q in sorted(files)]);files.append(O/'git_artifact_inventory.json');(O/'SHA256SUMS').write_text(''.join(f'{sha(q)}  {q.relative_to(R)}\n' for q in sorted(files)))
 for line in (O/'SHA256SUMS').read_text().splitlines():
  h,name=line.split('  ',1);assert sha(R/name)==h,name
 print('sealed',len(files),'files',flush=True)
if __name__=='__main__':main()
