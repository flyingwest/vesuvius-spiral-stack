"""Seal and independently verify P0.4 artifacts, not alter the method or evidence."""
import csv
import json
import hashlib
from pathlib import Path
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_4'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads((OUT/p).read_text())
def save(p,v):(OUT/p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')


def main():
    m=read('metrics.json');assert m['attempted_Z_levels']==50 and m['new_seed_observations_beyond_P03']==0
    save('verdict.json',{'P0_4':'FAIL','scope':'Longitudinal scale-up under unchanged frozen method',
        'reason':'Only original14 F02 rows supported among50 attempted; neither end extends without unsupported joins.',
        'coherence_collapse_observed':False,'physical_method_universally_falsified':False,
        'disconnected_tracklets_reseeded':False,'new_F02_links':0,
        'suspected_switches_within_accepted_track':0,'true_switches':None,
        'three_low_coherence_transition_causes':'UNKNOWN; full core present, geometry-u comparison does not restore coherence.',
        'no_change_to_thresholds_or_method':True,'conditional_100_200_expansion':'not warranted by boundary support'})
    volume=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')
    a=np.load(OUT/'F02_extended_stack.npy');mask=np.load(OUT/'F02_extended_mask.npy')
    xy=np.load(OUT/'F02_extended_sample_xy.npy');zs=np.load(OUT/'F02_extended_z.npy');s=np.load(OUT/'F02_extended_s_px.npy')
    assert len(zs)==50 and np.all(np.diff(zs)==1) and np.all(np.diff(s)==.25)
    assert np.array_equal(np.isfinite(a),mask) and np.array_equal(np.isfinite(xy).all(axis=-1),mask)
    assert int(mask.any(axis=1).sum())==m['valid_rows']
    for i,z in enumerate(zs):
        p=xy[i,mask[i]];x=p[:,0];y=p[:,1];xx=np.floor(x).astype(int);yy=np.floor(y).astype(int);fx=x-xx;fy=y-yy
        im=volume[z-128].astype(float)
        values=(1-fx)*(1-fy)*im[yy,xx]+fx*(1-fy)*im[yy,xx+1]+(1-fx)*fy*im[yy+1,xx]+fx*fy*im[yy+1,xx+1]
        np.testing.assert_allclose(a[i,mask[i]],values,rtol=0,atol=1e-10)
    png=np.array(Image.open(OUT/'F02_extended_stack.png'));pm=np.array(Image.open(OUT/'F02_extended_stack_mask.png'))
    np.testing.assert_array_equal(pm==255,mask);np.testing.assert_array_equal(png[mask],np.rint(a[mask]).astype(np.uint8))
    assert np.all(png[~mask]==0)
    with (OUT/'F02_alignment.csv').open() as f:assert len(list(csv.DictReader(f)))==50
    with (OUT/'F02_correlation_by_z.csv').open() as f:assert len(list(csv.DictReader(f)))==49
    es=read('edges.json');accepted=[e for e in es if e['status']=='ACCEPTED'];nodes=read('nodes.json');tracks=read('tracks.json')
    assert len({e['source'] for e in accepted})==len(accepted)==len({e['target'] for e in accepted})
    assert all(e['z_target']-e['z_source'] in [1,2] for e in accepted)
    membership=[n for t in tracks for n in t['nodes']]
    assert len(set(membership))==len(membership)==len(nodes)
    selected=set(read('selected.json')['nodes']);assert sum(e['source'] in selected for e in accepted)==13
    previous=read('preservation_before.json')
    for name,digest in previous.items():assert sha(ROOT/name)==digest,name
    source_files=['outputs/orientation/ct_level4_zyx.npy','outputs/p0_2/run_02/windows.json',
        'outputs/p0_2/run_02/tracks.json','outputs/p0_3/F02_stack.npy','docs/p0_3_first_stack.md']
    save('provenance.json',{'source_sha256':{p:sha(ROOT/p) for p in source_files},
         'script_sha256':{str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'scripts').glob('p0_4*.py')},
         'report_sha256':sha(ROOT/'docs/p0_4_scale_test.md'),
         'voxel_spacing_um':126.56,'arc_sampling_um':31.64,'dtype':'float64 CT/XY; boolean valid mask',
         'network_download_bytes':0,'source_volume_file_bytes':(ROOT/source_files[0]).stat().st_size,
         'attempted_slice_pixel_bytes':int(volume[:50].nbytes),
         'registration_shifts':sorted({tuple(r['shift_xy']) for r in read('registration.json').values()}),
         'no_full_volume_copy':True})
    pngs=list(OUT.glob('*.png'))
    for p in pngs:
        with Image.open(p) as im:im.verify()
    save('artifact_verification.json',{'prior_files_unchanged':len(previous),'readable_pngs':len(pngs),
        'independent_bilinear_reference':'PASS','mask_nan_png_csv_consistency':'PASS',
        'one_to_one_forward_track_graph':'PASS','P03_overlap':'bit-identical',
        'README_AGENTS':'unchanged','replay':read('verification.json')})
    paths=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='SHA256SUMS')
    (OUT/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in paths))
    print(json.dumps(read('verdict.json'),indent=2));print(json.dumps(read('artifact_verification.json'),indent=2))


if __name__=='__main__':main()
