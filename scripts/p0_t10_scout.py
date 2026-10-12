"""Bounded label-only crop selection, before CT/model evidence."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import time
import urllib.request
import urllib.error
import numpy as np
from numcodecs import Blosc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t10'
DATA = ROOT / 'data/raw/p0_t10'
BASE = 'https://huggingface.co/buckets/scrollprize/datasets/resolve/ink_9um/labels/aligned-scrollprizeorg-21slices/phercparis4-w02/phercparis4-w02_'

def fetch(item):
    kind, y, x = item
    url = f'{BASE}{kind}.zarr/0/0.{y}.{x}'
    destination = DATA / f'{kind}_{y}_{x}.bin'
    record = {'url': url, 'filename': str(destination.relative_to(ROOT)), 'chunk': [0, y, x]}
    try:
        with urllib.request.urlopen(urllib.request.Request(url + '?nonce=' + str(time.time_ns()), method='HEAD'), timeout=30) as response:
            size = response.headers.get('Content-Length')
        assert size is None or int(size) < 1000000
        with urllib.request.urlopen(url + '?nonce=' + str(time.time_ns()), timeout=30) as response:
            payload = response.read(1000001)
        assert len(payload) <= 1000000
        destination.write_bytes(payload)
        record.update(status=200, bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
        array = np.frombuffer(Blosc().decode(payload), np.uint8).reshape(21, 128, 128)[10].copy()
    except urllib.error.HTTPError as exc:
        if exc.code != 404:
            raise
        record.update(status=404, bytes=0, fill_value=0)
        array = np.zeros((128, 128), np.uint8)
    return item, array, record

def main():
    DATA.mkdir(parents=True, exist_ok=True)
    # Fixed center-out order; no CT or model response is used for selection.
    origins = [(4096, 5248), (3072, 4224), (3072, 6272), (5120, 4224), (5120, 6272)]
    receipts, candidates = [], []
    for y0, x0 in origins:
        arrays = {kind: np.zeros((512, 512), np.uint8) for kind in ['inklabels', 'supervision_mask']}
        items = [(kind, y, x) for kind in arrays for y in range(y0 // 128, y0 // 128 + 4) for x in range(x0 // 128, x0 // 128 + 4)]
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            for (kind, y, x), array, receipt in pool.map(fetch, items):
                arrays[kind][y * 128 - y0:y * 128 - y0 + 128, x * 128 - x0:x * 128 - x0 + 128] = array
                receipts.append(receipt)
        mask = arrays['supervision_mask'] > 0
        positive = arrays['inklabels'] > 0
        record = {'y0': y0, 'x0': x0, 'size': 512, 'supervised': int(mask.sum()), 'positive': int((mask & positive).sum()), 'negative': int((mask & ~positive).sum())}
        candidates.append(record)
        print(json.dumps(record), flush=True)
        if record['positive'] >= 1024 and record['negative'] >= 1024 and record['supervised'] >= 131072:
            np.savez_compressed(OUT / 'selected_labels.npz', **arrays)
            selection = {'selected': record, 'rule': 'first fixed center-out512crop with >=1024positives, >=1024negatives, >=50%supervision; labels only, no CT/predictions', 'candidate_order_yx': origins, 'candidates_inspected': candidates, 'label_values': {k: np.unique(v).tolist() for k, v in arrays.items()}}
            (OUT / 'crop_selection.json').write_text(json.dumps(selection, indent=2) + '\n')
            break
    (OUT / 'label_acquisition.json').write_text(json.dumps(receipts, indent=2) + '\n')
    assert (OUT / 'crop_selection.json').exists(), 'No qualifying crop in bounded selection'

if __name__ == '__main__':
    main()
