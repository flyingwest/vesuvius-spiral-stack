"""Seal completed T9 artifacts, including gitignored lossless arrays."""
from pathlib import Path
import hashlib
import json
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t9'

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def main():
    required = ['ink_probability_stack_raw.png', 'ink_probability_stack_display.png', 'mask_stack.png', 'provenance_stack.png', 'physical_aspect_stack.png', 'metrics.json', 'model_metadata.json', 'inference_parameters.json', 'verification.json', 'ink_probability.npy']
    assert all((OUT / name).is_file() for name in required)
    verification = json.loads((OUT / 'verification.json').read_text())
    assert verification['status'] == 'PASS' and verification['rows_verified'] == 50
    metrics = json.loads((OUT / 'metrics.json').read_text())
    assert metrics['scientific_result'] in ['PASS', 'PARTIAL', 'FAIL']
    probability = np.load(OUT / 'ink_probability.npy')
    mask = np.load(OUT / 'mask.npy')
    assert probability.shape == mask.shape and probability.shape[0] == 50
    assert np.array_equal(np.isnan(probability), mask)
    assert np.isfinite(probability[~mask]).all()
    assert ((probability[~mask] >= 0) & (probability[~mask] <= 1)).all()
    for name in required[:4]:
        with Image.open(OUT / name) as picture:
            assert picture.size == (probability.shape[1], 50)
    with Image.open(OUT / 'physical_aspect_stack.png') as picture:
        assert picture.size == (probability.shape[1], 100)
    prior = json.loads((OUT / 'preservation.json').read_text())
    assert all(sha(ROOT / path) == digest for path, digest in prior.items())
    original = json.loads((OUT / 'FROZEN_RUN_PLAN.json').read_text())
    amendment = json.loads((OUT / 'nemotron_recovery/infrastructure_amendment.json').read_text())
    for name, digest in original['code_hashes'].items():
        if name == 'scripts/p0_t9_infer.py':
            assert digest == amendment['original_sha256']
            assert sha(OUT / 'nemotron_recovery/p0_t9_infer_original.py') == digest
            assert sha(ROOT / name) == amendment['repaired_sha256']
        else:
            assert sha(ROOT / name) == digest, name
    assert sha(ROOT / 'data/raw/p0_t9/step-075000.pth') == original['checkpoint_sha256']
    acquisition_path = OUT / 'acquisition.json'
    recovery = json.loads((OUT / 'execution_recovery/fresh_validation.json').read_text())
    assert sha(acquisition_path) == recovery['acquisition_manifest_sha256']
    acquisition = json.loads(acquisition_path.read_text())['chunks']
    assert len(acquisition) == len({tuple(row['key']) for row in acquisition}) == 5523
    assert sum(row['expected_bytes'] for row in acquisition) == 4443682806
    for row in acquisition:
        path = ROOT / row['filename']
        assert path.stat().st_size == row['expected_bytes']
        assert sha(path) == row['sha256']
    restoration = json.loads((OUT / 'nemotron_recovery/restoration.json').read_text())
    assert restoration['operational_health'] == 'PASS'
    assert restoration['original_configuration_unchanged']
    files = sorted(path for path in OUT.rglob('*') if path.is_file() and path.name != 'SHA256SUMS')
    files += [ROOT / 'docs/p0_t9_pretrained_ink_model.md']
    files += sorted((ROOT / 'scripts').glob('p0_t9_*.py'))
    files += sorted((ROOT / 'scripts').glob('p0_t9_*.sh'))
    entries = [(path, sha(path)) for path in files]
    (OUT / 'SHA256SUMS').write_text(''.join(f'{digest}  {path.relative_to(ROOT)}\n' for path, digest in entries))
    assert all(sha(path) == digest for path, digest in entries)
    print(json.dumps({'status': 'PASS', 'files_sealed': len(entries), 'shape': list(probability.shape), 'valid_predictions': int((~mask).sum()), 'prior_files_preserved': len(prior), 'scientific_result': metrics['scientific_result'], 'nemotron_operational_health': restoration['operational_health'], 'requested_model_name_check': restoration['requested_model_name_check']}, indent=2))

if __name__ == '__main__':
    main()
