"""Bounded real-CT acquisition and orthogonal inspection; no reconstruction."""
import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import urllib.request

import numpy as np
import zarr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/orientation'
RAW = ROOT / 'data/raw/p0_orientation/scroll1.zarr'
BASE = ('https://dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/'
        'volumes_zarr_standardized/54keV_7.91um_Scroll1A.zarr/')
BOUNDS = [[128, 256], [0, 493], [0, 506]]
SPACING = 7.91 * 16
KEYS = ['.zattrs', '0/.zarray', '4/.zarray'] + [
    f'4/1/{y}/{x}' for y in range(4) for x in range(4)]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire():
    if RAW.exists() or (OUT / 'manifest.json').exists():
        raise RuntimeError('Acquisition already exists; use render/verify, do not overwrite.')
    OUT.mkdir(parents=True, exist_ok=True)
    heads = []
    for key in KEYS:
        with urllib.request.urlopen(urllib.request.Request(BASE + key, method='HEAD'), timeout=30) as r:
            heads.append({'key': key, 'url': BASE + key, 'bytes': int(r.headers['Content-Length'])})
    assert sum(h['bytes'] for h in heads) < 24 * 2**20, 'Download budget exceeded'
    (OUT / 'acquisition_plan.json').write_text(json.dumps({
        'hypothesis': 'Stored z is approximately longitudinal; XY shows rolled cross-section.',
        'scope': 'Orientation inspection, not P0 reconstruction or formal hypothesis test.',
        'bounds_zyx': BOUNDS, 'budget_bytes': 24 * 2**20, 'files': heads,
        'stop_rule': 'Stop after orientation is established; no tracing or unwrapping.',
    }, indent=2) + '\n')

    def fetch(h):
        path = RAW / h['key']
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(h['url'], timeout=60) as r:
            body = r.read(h['bytes'] + 1)
            assert len(body) == h['bytes'], 'HEAD/GET size mismatch'
            headers = dict(r.headers)
        path.write_bytes(body)
        return h | {'local_path': str(path.relative_to(ROOT)), 'sha256': sha(path),
                    'response_headers': headers}

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        files = list(pool.map(fetch, heads))
    attrs = json.loads((RAW / '.zattrs').read_text())['multiscales'][0]
    assert [a['name'] for a in attrs['axes']] == ['z', 'y', 'x']
    level = next(d for d in attrs['datasets'] if d['path'] == '4')
    assert level['coordinateTransformations'] == [{'scale': [16., 16., 16.], 'type': 'scale'}]
    manifest = {
        'dataset': 'Scroll1 / PHercParis4 / 54keV_7.91um_Scroll1A',
        'base_url': BASE, 'level': 4, 'axes': 'zyx', 'bounds_half_open_zyx': BOUNDS,
        'native_bounds_half_open_zyx': [[a*16, b*16] for a, b in BOUNDS],
        'shape_zyx': [128, 493, 506], 'dtype': 'uint8',
        'voxel_spacing_um_zyx': [SPACING]*3, 'native_spacing_um': 7.91,
        'spacing_source': 'external/villa/vesuvius/src/vesuvius/install/configs/scrolls.yaml:5-11',
        'local_store_warning': 'Partial cache only: read exactly the recorded bounds; absent chunks elsewhere are NOT acquired data.',
        'downloaded_files': files, 'download_payload_bytes': sum(f['bytes'] for f in files),
        'environment': {'python': platform.python_version(), **{
            p: importlib.metadata.version(p) for p in ['numpy', 'zarr', 'numcodecs', 'matplotlib', 'pillow']}},
        'script_sha256': sha(Path(__file__)),
        'protected_sha256': {p: sha(ROOT / p) for p in ['README.md', 'AGENTS.md']},
    }
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Acquired', manifest['download_payload_bytes'], 'bytes', flush=True)


def load():
    manifest = json.loads((OUT / 'manifest.json').read_text())
    for f in manifest['downloaded_files']:
        path = ROOT / f['local_path']
        assert path.stat().st_size == f['bytes'] and sha(path) == f['sha256']
    a = zarr.open_array(str(RAW / '4'), mode='r')[128:256, 0:493, 0:506]
    assert a.shape == (128, 493, 506) and a.dtype == np.uint8
    return a, manifest


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from PIL import Image
    a, manifest = load()
    np.save(OUT / 'ct_level4_zyx.npy', a)
    views = [('XY_z0144', a[16], 'x', 'y', 0, 0),
             ('XY_z0192', a[64], 'x', 'y', 0, 0),
             ('XY_z0240', a[112], 'x', 'y', 0, 0),
             ('XZ_y0246', a[:, 246, :], 'x', 'z', 0, 128),
             ('YZ_x0253', a[:, :, 253], 'y', 'z', 0, 128)]
    for name, v, horizontal, vertical, h0, v0 in views:
        Image.fromarray(v).save(OUT / (name + '_raw.png'))
        fig, ax = plt.subplots(figsize=(9, 9 if horizontal == 'x' and vertical == 'y' else 4))
        extent = [(h0-.5)*SPACING/1000, (h0+v.shape[1]-.5)*SPACING/1000,
                  (v0+v.shape[0]-.5)*SPACING/1000, (v0-.5)*SPACING/1000]
        ax.imshow(v, cmap='gray', vmin=0, vmax=255, origin='upper',
                  interpolation='nearest', extent=extent, aspect='equal')
        ax.set(xlabel=f'{horizontal} (mm, scan-coordinate origin)',
               ylabel=f'{vertical} (mm, scan-coordinate origin)',
               title=f'Scroll 1 | {name} | 126.56 um/voxel | uint8 0–255')
        fig.tight_layout()
        fig.savefig(OUT / (name + '.png'), dpi=150)
        plt.close(fig)
    # Descriptive continuity check in a fixed central cuboid, no tracing or fitted axis.
    roi = a[:, 80:420, 80:420].astype(np.float32)
    metrics = {'shape': list(a.shape), 'dtype': str(a.dtype), 'decoded_bytes': a.nbytes,
               'nonzero_fraction': float(np.mean(a != 0)),
               'intensity_quantiles': dict(zip(['0', '1', '50', '99', '100'],
                                              np.percentile(a, [0, 1, 50, 99, 100]).tolist())),
               'continuity_roi_global_zyx': [[128,256], [80,420], [80,420]],
               'adjacent_voxel_statistics': {}}
    for axis, name in enumerate('zyx'):
        left = np.take(roi, range(roi.shape[axis]-1), axis=axis)
        right = np.take(roi, range(1,roi.shape[axis]), axis=axis)
        metrics['adjacent_voxel_statistics'][name] = {
            'mean_absolute_difference': float(np.mean(np.abs(left-right))),
            'pearson_correlation': float(np.corrcoef(left.ravel(), right.ravel())[0,1]),
            'pairs': left.size}
    metrics['images'] = [name + suffix for name, *_ in views for suffix in ['.png', '_raw.png']]
    metrics['physical_fov_mm_zyx'] = (np.array(a.shape)*SPACING/1000).tolist()
    (OUT / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
    print(json.dumps(metrics, indent=2))


def verify():
    from PIL import Image
    a, manifest = load()
    assert np.array_equal(a, np.load(OUT / 'ct_level4_zyx.npy'))
    for name, expected in [('XY_z0144',a[16]), ('XY_z0192',a[64]), ('XY_z0240',a[112]),
                           ('XZ_y0246',a[:,246,:]), ('YZ_x0253',a[:,:,253])]:
        assert np.array_equal(np.asarray(Image.open(OUT / (name+'_raw.png'))), expected)
    for p, digest in manifest['protected_sha256'].items():
        assert sha(ROOT / p) == digest
    files = sorted(p for p in OUT.iterdir() if p.is_file() and p.name != 'SHA256SUMS')
    (OUT / 'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in files))
    print('PASS: source hashes/sizes, offline decoded array, all orthogonal raw pixels, protected files')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['acquire', 'render', 'verify'])
    args = parser.parse_args()
    {'acquire': acquire, 'render': render, 'verify': verify}[args.action]()
