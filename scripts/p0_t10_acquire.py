"""Acquire the frozen small official control, with an explicit HEAD budget."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import urllib.request
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t10'
DATA = ROOT / 'data/raw/p0_t10'

def main():
    selection = json.loads((OUT / 'crop_selection.json').read_text())['selected']
    candidates = json.loads((OUT / 'candidate_metadata.json').read_text())
    selected = next(row for row in candidates if row['segment'] == 'phercparis4-w02')
    y0, x0 = selection['y0'], selection['x0']
    plan = []
    for y in range(y0 // 128, y0 // 128 + 4):
        for x in range(x0 // 128, x0 // 128 + 4):
            plan.append({'url': f"{selected['surface_url']}/2/0/{y}/{x}", 'filename': f'data/raw/p0_t10/ct_0_{y}_{x}.bin', 'chunk': [0, y, x], 'kind': 'CT'})
    catalog = json.loads((OUT / 'segment_catalog.json').read_text())['mesh/']
    for item in catalog:
        if 'on-20260411134726-2.4um.tifxyz/' in item['key']:
            plan.append({'url': 'https://vesuvius-challenge-open-data.s3.amazonaws.com/' + item['key'], 'filename': 'data/raw/p0_t10/geometry_' + item['key'].rsplit('/', 1)[-1], 'kind': 'geometry'})
    def head(row):
        with urllib.request.urlopen(urllib.request.Request(row['url'], method='HEAD'), timeout=30) as response:
            row['expected_bytes'] = int(response.headers['Content-Length'])
            row['etag'] = response.headers.get('ETag')
        return row
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        plan = list(pool.map(head, plan))
    total = sum(row['expected_bytes'] for row in plan)
    assert total < 100000000
    (OUT / 'control_acquisition_plan.json').write_text(json.dumps({'cap_bytes': 100000000, 'total_bytes': total, 'objects': plan}, indent=2) + '\n')
    print('HEAD budget', total, 'bytes', flush=True)
    def acquire(row):
        destination = ROOT / row['filename']
        assert not destination.exists()
        with urllib.request.urlopen(row['url'], timeout=60) as response:
            payload = response.read(row['expected_bytes'] + 1)
        assert len(payload) == row['expected_bytes']
        destination.write_bytes(payload)
        row['sha256'] = hashlib.sha256(payload).hexdigest()
        return row
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        receipt = list(pool.map(acquire, plan))
    (OUT / 'control_acquisition.json').write_text(json.dumps(receipt, indent=2) + '\n')
    # A local, exact crop of the producer-rendered surface volume. No new geometry.
    zroot = OUT / 'source_crop.zarr'
    (zroot / '2').mkdir(parents=True, exist_ok=False)
    (zroot / '.zgroup').write_text('{"zarr_format":2}')
    (zroot / '.zattrs').write_text(json.dumps({'original_source': selected['surface_url'], 'original_level': 2, 'crop_origin_yx': [y0, x0]}))
    schema = dict(selected['/2/.zarray']); schema['shape'] = [109, 512, 512]
    (zroot / '2/.zarray').write_text(json.dumps(schema))
    for row in receipt:
        if row['kind'] != 'CT':
            continue
        _, y, x = row['chunk']
        target = zroot / f'2/0/{y-y0//128}/{x-x0//128}'
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = (ROOT / row['filename']).read_bytes()
        assert np.frombuffer(payload, np.uint8).size == 109 * 128 * 128
        target.write_bytes(payload)
    print('Acquisition complete; CT and trusted geometry preserved', flush=True)

if __name__ == '__main__':
    main()
