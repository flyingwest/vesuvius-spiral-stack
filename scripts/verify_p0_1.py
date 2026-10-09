"""Offline verification of preserved P0.1 numerical artifacts, not sheet identity."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import Image
from p0_1_single_slice import geometry, trace, resample, sample

root=Path(__file__).resolve().parents[1]
out=root/'outputs/p0_1'
metrics=json.loads((out/'metrics.json').read_text())
before=json.loads((out/'preservation_before.json').read_text())
for name,digest in before.items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
image=np.load(root/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')[16]
assert np.array_equal(image,np.asarray(Image.open(out/'source_slice.png')))
guide,normal=geometry()
path,_,_=trace(image,guide,normal)
assert np.array_equal(path,np.load(out/'optimizer_path_xy.npy'))
arc,xy=resample(path)
with (out/'trace_coordinates.csv').open() as f: rows=list(csv.DictReader(f))
assert len(rows)==metrics['trace_points']
for field,expected in [('s_pixels',arc),('s_mm',arc*.12656),('x_level4',xy[:,0]),('y_level4',xy[:,1])]:
    assert np.array_equal(np.array([float(r[field]) for r in rows]),expected),field
assert all(float(r['z_level4'])==144 and float(r['z_native_nominal'])==2304 for r in rows)
values=sample(image,xy)
assert np.array_equal(values,np.load(out/'sampling_values.npy'))
assert np.array_equal(np.asarray(Image.open(out/'unwrapped_row.png')),np.rint(values).astype(np.uint8)[None,:])
length=np.linalg.norm(np.diff(path,axis=0),axis=1).sum()
assert abs(length-metrics['arc_length_px'])<1e-10
assert np.all((xy[:,0]>=0)&(xy[:,0]<image.shape[1]-1)&(xy[:,1]>=0)&(xy[:,1]<image.shape[0]-1))
assert hashlib.sha256((root/'scripts/p0_1_single_slice.py').read_bytes()).hexdigest()==metrics['script_sha256']
result={'status':'PASS_NUMERICAL_AND_PRESERVATION_ONLY','prior_files_verified':len(before),
        'source_slice_exact':True,'deterministic_path_exact':True,'csv_matches_path':True,
        'png_matches_bilinear_sampling':True,'arc_length_recomputed':True,
        'physical_winding_identity':'NOT_ESTABLISHED; see report and trace_audit.png'}
(out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
(out/'SHA256SUMS').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n'
    for p in sorted(out.iterdir()) if p.is_file() and p.name!='SHA256SUMS'))
print(json.dumps(result,indent=2))
