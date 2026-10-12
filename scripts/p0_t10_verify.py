"""Read-only evidence checks, with a T10-only verification receipt."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.stats import rankdata
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t10'

def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def read(name):
    return json.loads((OUT / name).read_text())

def main():
    preserved = read('preservation.json')
    for name, expected in preserved.items():
        assert sha(ROOT / name) == expected, name
    frozen = read('PRE_EVIDENCE_FREEZE.json')['source_sha256']
    for name, expected in frozen.items():
        assert sha(ROOT / name) == expected, name
    acquired = read('control_acquisition.json') + read('label_acquisition.json')
    for item in acquired:
        assert sha(ROOT / item['filename']) == item['sha256'], item['filename']
    p = np.load(OUT / 'official_control_probability.npy')
    mask = np.load(OUT / 'evaluation_mask.npy')
    with np.load(OUT / 'selected_labels.npz') as archive:
        y = archive['inklabels'][mask] > 0
    scores = p[mask]
    assert p.shape == (512, 512) and np.isfinite(p).all()
    assert p.min() >= 0 and p.max() <= 1
    npos, nneg = int(y.sum()), int((~y).sum())
    auc = float((rankdata(scores)[y].sum(dtype=np.float64) - npos*(npos+1)/2)/(npos*nneg))
    # Independent score-bin AP calculation, ascending bins then reverse sums.
    values, inverse = np.unique(scores, return_inverse=True)
    positives = np.bincount(inverse, weights=y, minlength=len(values))
    totals = np.bincount(inverse, minlength=len(values))
    tp = np.cumsum(positives[::-1])[::-1]
    count = np.cumsum(totals[::-1])[::-1]
    ap = float(np.sum(positives/npos * tp/count))
    metrics = read('official_control_metrics.json')
    assert abs(auc-metrics['ROC_AUC']) < 1e-12
    assert abs(ap-metrics['PR_AUC_average_precision']) < 1e-12
    assert metrics['status'] == 'PASS' and all(metrics['gates'].values())
    retry, replay, receipt = read('gpu_memory_retry.json'), read('replay.json'), read('inference_receipt.json')
    assert retry['status'] == replay['status'] == 'PASS'
    assert replay['all_raw_patches_exact'] and replay['official_pooling_recomputed_exact']
    assert replay['maximum_probability_error'] <= replay['tolerance']
    assert receipt['official_reference_max_uint8_error'] == 0
    assert receipt['official_dataset_vs_T9_normalization_exact']
    for name in ['official_control_probability.png', 'official_control_display.png']:
        assert Image.open(OUT/name).size == (512,512)
    for name in ['control_vs_t9_sampling.png','representative_patch_comparison.png','control_vs_t9_geometry.json']:
        assert (OUT/name).stat().st_size > 0
    result = {'status':'PASS', 'predecessor_files_byte_identical':len(preserved),
              'pre_evidence_files_unchanged':len(frozen), 'acquisition_files_hash_verified':len(acquired),
              'independent_rank_sum_ROC_AUC':auc, 'independent_score_bin_AP':ap,
              'finite_probability_fraction':float(np.isfinite(p).mean()),
              'probability_range':[float(p.min()),float(p.max())],
              'replay':replay, 'official_api_encoded_max_error':0,
              'T9_full_inference_rerun':False, 'local_llm_dependency':'NONE',
              'checksum_scope':'SHA256SUMS seals local raw/lossless evidence as well as tracked reports; run sha256sum -c from repository root'}
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
