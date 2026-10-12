"""Diagnose shared NPZ reads; verify eager arrays against frozen archives."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t9/nemotron_recovery'

def main():
    failures = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for attempt in range(20):
            with np.load(ROOT / 'outputs/p0_t6/slice_128_profiles.npz') as archive:
                try:
                    list(pool.map(lambda j: (archive['xy'][j], archive['tangent'][j], archive['normal'][j]), range(64)))
                except Exception as exc:
                    failures.append({'attempt': attempt, 'type': type(exc).__name__, 'error': str(exc)})
        checks = []
        for z in range(128, 178):
            path = ROOT / f'outputs/p0_t6/slice_{z}_profiles.npz'
            with np.load(path) as archive:
                arrays = {key: archive[key] for key in archive.files}
            with np.load(path) as archive:
                for key, value in arrays.items():
                    assert np.array_equal(value, archive[key], equal_nan=True)
            indices = list(range(min(64, len(arrays['xy']))))
            actual = list(pool.map(lambda j: (arrays['xy'][j], arrays['tangent'][j], arrays['normal'][j]), indices))
            for j, triple in zip(indices, actual):
                for key, value in zip(['xy', 'tangent', 'normal'], triple):
                    assert np.array_equal(value, arrays[key][j])
            checks.append({'z': z, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    result = {'shared_archive_failures': failures, 'eager_arrays_serial_and_threaded_exact': True, 'archives_verified': checks}
    (OUT / 'archive_check.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'shared_archive_failures': len(failures), 'eager_archives_verified': len(checks)}))

if __name__ == '__main__':
    main()
