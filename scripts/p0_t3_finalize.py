"""Seal existing P0-T3 run_04. No tracing, optimization or experiment edits."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import p0_t3_verify as verifier

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/p0_t3';RUN=BASE/'run_04';CHECK=BASE/'finalization_checks'
REPORT=ROOT/'docs/p0_t3_unwindable_threads.md'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    assert not REPORT.exists() and not (BASE/'SHA256SUMS').exists()
    metrics=json.loads((RUN/'metrics.json').read_text())
    authority=json.loads((BASE/'verification.json').read_text())
    frozen={str(p.relative_to(ROOT)):sha(p) for p in BASE.rglob('*') if p.is_file()}
    frozen.update({n:sha(ROOT/n) for n in ['README.md','AGENTS.md']})
    CHECK.mkdir(exist_ok=False)
    # Existing verifier only reads stored paths and source CT. Redirect its receipt;
    # do not overwrite the user's authoritative verification.json.
    verifier.RUN=RUN;verifier.OUT=CHECK;verifier.main()
    replay=[]
    for r in metrics:
        z=r['z'];p=np.load(RUN/f'slice_{z}/path_and_samples.npz')
        a=p['edge_start_xy'];b=p['edge_end_xy'];length=np.linalg.norm(b-a,axis=1)
        s=np.r_[0,np.cumsum(length)];ss=np.arange(0,s[-1]+1e-8,.5)
        i=np.minimum(np.searchsorted(s,ss,side='right')-1,len(a)-1)
        xy=a[i]+((ss-s[i])/length[i])[:,None]*(b[i]-a[i])
        assert np.array_equal(s,p['edge_s_px']) and np.array_equal(ss,p['s_px']) and np.array_equal(xy,p['xy'])
        assert np.array_equal(ss*126.56,p['s_um'])
        assert len(ss)==r['total_samples']
        assert np.isclose(length[p['edge_provenance']!=2].sum(),r['observed_arc_px'])
        assert np.isclose(length[p['edge_provenance']==2].sum(),r['bridge_arc_px'])
        assert np.array_equal(np.load(RUN/f'slice_{z}/provenance.npy'),p['provenance'])
        replay.append({'z':z,'arc_and_XY_regeneration':'bit-identical','stored_provenance_array':'identical','metrics':'consistent'})
    rerun=json.loads((CHECK/'verification.json').read_text())
    assert rerun['slices']==authority['slices']
    assert all(s[k]==0 for s in rerun['slices'] for k in ['proper_crossings','nonadjacent_point_contacts','collinear_overlaps'])
    for name,h in frozen.items():assert sha(ROOT/name)==h,name
    receipt={'mode':'saved-artifact replay only; no tracing or topology optimization rerun','run':'outputs/p0_t3/run_04',
             'checks':replay,'existing_validator':'PASS, reproduced authoritative per-slice results',
             'prior_evidence_and_protected_files':'byte-identical','full_experiment_replay':'NOT RUN; would rerun tracing/search/topology experiment'}
    (CHECK/'replay.json').write_text(json.dumps(receipt,indent=2)+'\n')
    table='\n'.join(f"| {r['z']} | {r['total_arc_mm']:.3f} | {r['observed_arc_mm']:.3f} | {r['bridge_arc_mm']:.3f} | {r['observed_graph_percent']:.3f}% | {r['observed_retained_vs_T2_percent']:.3f}% | {r['fragments_incorporated']} | YES |" for r in metrics)
    rows='\n'.join(f"| {r['z']} | {r['total_samples']} | {r['observed_samples']} | {r['masked_bridge_samples']} | {r['unresolved_samples']} | {r['low_confidence_samples']} |" for r in metrics)
    topology=[]
    for r in metrics:
        corrections=json.loads((RUN/f"slice_{r['z']}/corrections.json").read_text())
        proper=sum(c['crossing']['type']=='proper' for c in corrections)
        contacts=sum(c['crossing']['type']=='point_contact' for c in corrections)
        topology.append(f"| {r['z']} | {proper} → 0 | {contacts} → 0 | 0 | {r['corrections']} |")
    report=f'''# P0-T3 — continuous unwindable threads

## Authoritative result and conclusion

**P0_T3_COMPLETE = YES. READY_FOR_P0_T4 = YES.**

The authoritative final run is `outputs/p0_t3/run_04/`; the authoritative validation is
`outputs/p0_t3/verification.json`. All five selected slices reach the exterior with continuous,
physically spaced 1-D coordinates. Final exhaustive validation reports **zero proper crossings,
zero non-adjacent point contacts and zero collinear overlaps on every slice**. Z146's residual
point contact is resolved. Observed contour retention is **99.843–100% of P0-T2**.

This satisfies P0-T3's geometry-and-unwinding criterion without sacrificing maximum recovery.
Readiness means that the per-slice artifacts are available for a separately authorized next phase;
it is not a claim of cross-Z alignment, coherent stacked reconstruction, or physical sheet identity.
No P0-T4 work, stacking or new experiment was performed during finalization.

## Objective and input scope

Make each existing maximum-recovery PHercParis4 transverse XY path a continuous center-to-exterior
coordinate suitable for unwinding, removing obvious topology defects while keeping observed contour data.
The five global slice indices are **128, 144, 146, 160 and 170**. Processing is independent per slice.
The preserved input source is `outputs/orientation/ct_level4_zyx.npy`, referenced in the run's preservation
manifest. The finalization report uses only run_04 metrics, its preserved plan/corrections/arrays,
and the authoritative validation; no new tracing findings are introduced.

Physical pixel spacing is126.56 µm. Arc sampling is fixed at0.5 pixel = **63.28 µm**.
Coordinates retain X horizontally and Y vertically; each slice has its own seed and path order.
The unchanged inherited exterior-region endpoints are recorded in `metrics.json`; all five
`exterior_reached` results are true. This is a designated exterior region, not a newly measured
outermost material boundary.

## Method and topology cleanup

Start from the selected P0-T2 maximum-recovery path, retaining its observed/inferred/low-confidence
provenance. The preserved plan records that selection from existing archived crossing-free alternatives
would lose31–48% of observed contour in the affected slices. Removing entire crossing-bounded loops
would remove96–682 pixels of intervening path. Those approaches were not adopted as the final result.

Instead use a local **two-edge exchange**: replace intersecting edges `A→B` and `C→D` with
`A→C` and `B→D`, reversing the intervening path order. Its interior contour geometry is retained,
possibly traversed in reverse, rather than deleted. Both replacement links are explicit inferred
bridges, regardless of any bright CT beneath them. For candidate exchanges, require fewer
intersection events, then prefer less lost observed length and shorter new bridges; the preserved
plan caps corrections at20. No ideal spiral, global surface, cross-Z correspondence or learned
correction is used.

The final contact amendment extends the same operation to vertex-on-edge contacts. Observed geometry
within approximately5 pixels of a new bridge is locally relabeled low-confidence. This does not
mask or discard the downstream path. Z144 and Z170 retain their original selected T2 geometry;
Z128/Z146/Z160 receive respectively1/3/1 exchanges. Every correction's old/new endpoint coordinates,
provenance, reversed source-edge IDs and removed observed length is preserved in `corrections.json`.

### Before/after topology

Point-contact counts below are distinct contact locations, not incident segment-pair counts.
Final collinear overlaps are from the authoritative exhaustive validator.

| Z | Proper crossings before → after | Point-contact sites before → after | Final collinear overlaps | Exchanges |
| --- | --- | --- | --- | --- |
{chr(10).join(topology)}

**Z146 metadata clarification:** `run_04/metrics.json` records `crossings_before=4` because that
implementation counts two proper crossings plus two incident segment-pair contacts at one residual
point-contact location. It does not mean four proper crossings. The preserved contact amendment and
typed correction records distinguish these cases. The final counts are zero under both interpretations.

### Run identity and negative evidence

The phrase “all run03 pre-existing hashes PASS” in `verification.json` is stale descriptive text.
The inspected final validation wrapper selects **run_04**, and the validated slice results and saved
artifacts are run_04. The authoritative receipt is preserved byte-for-byte, rather than edited to
erase that historical wording. An isolated finalization recheck reproduced its per-slice results.
Earlier failed attempts and intermediate runs remain preserved as negative evidence and are not
the final result. In particular the run_04 plan documents why archive-only shortening and the
proper-crossing-only cleanup were insufficient.

## Recovery measurements

Observed length includes OBSERVED and LOW_CONFIDENCE contour edges, excluding inferred bridges.
These are detector-supported contours, not independently certified papyrus labels. Fragments are
incorporated source graph fragments, not distinct proven physical sheets.

| Z | Total path mm | Observed mm | Bridge mm | Observed / total | Observed retained vs T2 | Source fragments | Exterior reached |
| --- | --- | --- | --- | --- | --- | --- | --- |
{table}

The small observed losses are explicitly recorded in `metrics.json`: none for Z144/Z170,
1.414214 pixels each for Z128/Z160, and2.414214 pixels for Z146. There is no dramatic shortening
and no minimum observed-percentage gate. Approximately85.829–88.139% of final geometric length
is on detected observed contours; bridge length remains part of the coordinate.

## Path, provenance and continuous physical arc length

The final polyline is an ordered sequence of edges with exact shared endpoints. Its cumulative
Euclidean length includes **every inferred bridge**. Multiplying by126.56 converts pixels to µm.
Uniform samples start at s=0, advance63.28 µm, and stop at the last grid point within the path;
no row is resized to a convenient width and no endpoint sample is appended off that grid.

Each edge and each sampled point carries provenance:

| Code | Meaning | CT signal |
| --- | --- | --- |
| 0 | unresolved / no usable ridge-supported data | NaN, masked |
| 1 | observed CT-supported contour | bilinear original CT |
| 2 | inferred bridge | NaN, masked; continuous XY and s retained |
| 3 | low-confidence observed/transition neighborhood | bilinear CT where ridge-supported |

XY is interpolated linearly along the saved polyline at the fixed arc grid. Original grayscale
CT is sampled bilinearly only for codes1/3 with usable ridge evidence; the inherited interpolated
ridge-response check is>3. Inferred samples are never assigned invented CT intensity. A geometry
edge may remain observed while individual samples receive code0 when interpolated ridge support dips.
`mask.npy` is true for missing data (codes0/2); the uint8 mask PNG uses255 for masked and0 for valid.
The signal PNG is a height-one, fixed0–255 grayscale row, rounded to uint8; masked pixels display
black. Lossless arrays and masks distinguish missing black from genuine dark CT values.

| Z | Total samples / row width | Observed samples (1+3) | Masked bridge samples (2) | Other unresolved samples (0) | Low-confidence samples (3, included in observed) |
| --- | --- | --- | --- | --- | --- |
{rows}

Sample counts need not change exactly in proportion to retained arc: reversing/reconnecting a path
changes where the fixed arc grid lands on subsequent geometry. Slight increases in usable sample
count are sampling-phase effects, not invented observed material.

## Authoritative artifacts

All final per-slice files are under **`outputs/p0_t3/run_04/slice_<z>/`**, for each of the five Z values:

- `final_thread_overlay.png`: continuous final geometry, seed/exterior and corrected sites.
- `provenance_overlay.png`: observed/inferred/low-confidence portions distinguished.
- `unwrapped_row.png`, `unwrapped_mask.png`: actual one-row signal and explicit mask.
- `path_coordinates.csv`: ordered edge endpoints, cumulative physical s, provenance and original T2 edge IDs.
- `sample_provenance.csv`: each sample's physical s, XY, label, mask and CT value (blank when missing).
- `row.npy`, `mask.npy`, `provenance.npy`: lossless sample arrays.
- `path_and_samples.npz`: lossless complete edge geometry, source lineage, edge/sample provenance,
  cumulative coordinates, CT and mask.
- `corrections.json`: exact local topology changes; empty for unchanged paths.
- `topology_corrections.png`: before/after detail for slices with corrections.

Run-level `metrics.json`, `plan.json` and `preservation.json` retain measurements, method and prior hashes.
The authoritative validator receipt is `outputs/p0_t3/verification.json`.
The final checksum manifest is **`outputs/p0_t3/SHA256SUMS`** and covers this report, every file
under run_04, and the authoritative validation. It additionally covers finalization check receipts
and the finalization utility for reproducibility.

## Validation and non-tracing reproducibility

The authoritative validation and its final isolated recheck both establish:

- Five exterior-reaching, continuous paths; exact shared endpoints and continuous physical s.
- Zero proper crossings, non-adjacent point contacts and collinear overlaps, using exhaustive
  orientation/bounding-box segment-pair tests independent of the KD-tree cleanup detector.
- Retained observed edges are original T2 geometry, possibly reversed; only local1→3 confidence
  upgrades are allowed. New replacement edges are inferred, never relabeled observed.
- Bridge interiors retain their physical s extent and masked NaN samples.
- Independent four-corner bilinear CT sampling agrees within1e-10.
- Lossless arrays, CSV sample counts, PNG values/masks and image decoding are consistent.
- Preservation checks passed for earlier evidence, T2 documentation and protected README/AGENTS.

Finalization calls the existing verifier with `RUN=run_04` and redirects its new receipt to
`outputs/p0_t3/finalization_checks/verification.json`; it does not overwrite the authoritative receipt.
It also regenerates cumulative arc, uniform s and XY from saved final edges, obtaining bit-identical
arrays on all five slices, and checks metrics/provenance consistency. Results are recorded in
`finalization_checks/replay.json`. This is a **saved-artifact/coordinate and sampling replay**, not
a new tracing run or a full algorithm replay. Repeating the original search/topology construction
would rerun the experiment and was deliberately not done.

Finalization command (exclusive output creation; not intended to overwrite a completed seal):

```bash
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t3 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t3_finalize.py
```

The utility generates the manifest and immediately runs, from repository root:

```bash
sha256sum -c outputs/p0_t3/SHA256SUMS
```

This read-only checksum command can be repeated independently after sealing.

## Limitations and stopping point

Topology cleanliness is a geometric property of the saved inferred path, not proof that every
connection belongs to the same physical papyrus sheet. Two-edge exchanges can reverse long
intervening contour sequences; their ordering remains inferred even though almost all CT-supported
geometry is retained. Bridge lengths are geometric assumptions, not measured torn-material lengths.
Coarse CT ridge detection, polyline discretization and uncertain layer identity remain limitations.
Coordinates and rows are independent per slice; cross-Z correspondence or stack coherence has not
been demonstrated here. No OCR, ink inference, ideal spiral, mesh, learned correction or stacking
was introduced.

**Conclusion:** all five final paths are continuously parameterized and unwindable under the stated
2-D geometry criterion, with explicit gaps/provenance and essentially all maximum-recovery contour
retained. P0-T3 is complete and ready to hand off its sealed artifacts to a separately authorized
P0-T4. **Stop after P0-T3; no next-phase work executed.**
'''
    REPORT.write_text(report)
    files=[REPORT,BASE/'verification.json',Path(__file__).resolve()]
    files+=sorted(p for p in RUN.rglob('*') if p.is_file())
    files+=sorted(p for p in CHECK.rglob('*') if p.is_file())
    manifest=BASE/'SHA256SUMS'
    manifest.write_text(''.join(f'{sha(p)}  {p.relative_to(ROOT)}\n' for p in files))
    result=subprocess.run(['sha256sum','-c',str(manifest.relative_to(ROOT))],cwd=ROOT,capture_output=True,text=True)
    print(result.stdout,end='');assert result.returncode==0,result.stderr
    for name,h in frozen.items():assert sha(ROOT/name)==h,name
    print(json.dumps({'P0_T3_COMPLETE':'YES','READY_FOR_P0_T4':'YES','manifest_files':len(files),'manifest_verification':'PASS',
                     'replay':'saved-coordinate regeneration and existing non-tracing verification PASS','prior_evidence':'unchanged'}))

if __name__=='__main__':main()
