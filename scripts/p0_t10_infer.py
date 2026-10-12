"""Fixed official-surface positive control and replay; no parameter search."""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import json
import sys
import numpy as np
import torch
import tifffile
import zarr
from vesuvius.ink_detection.inference import infer as official
from vesuvius.ink_detection.preprocessing.prepare_9um_isotropic_input import prepare_isotropic_input

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t10'
CHECKPOINT = ROOT / 'data/raw/p0_t9/step-075000.pth'

def main():
    replay = '--replay' in sys.argv
    assert hashlib.sha256(CHECKPOINT.read_bytes()).hexdigest() == 'e635558ae6a1a807a7e5ec1e83adfd45bc3c0ac53883ea43f1d4e085d62a9cab'
    torch.set_num_threads(4)
    torch.manual_seed(2910)
    torch.backends.cudnn.benchmark = False
    input_path = OUT / 'pooled.zarr'
    if not replay:
        prepare_isotropic_input(str(OUT / 'source_crop.zarr'), input_path, level='2', workers=1)
    source = zarr.open_array(str(OUT / 'source_crop.zarr/2'), mode='r')[:]
    independently_pooled = np.rint(source[13:97].astype(np.float32).reshape(21, 4, 512, 512).mean(1)).astype(np.uint8)
    pooled = zarr.open_array(str(input_path / '0'), mode='r')[:]
    assert np.array_equal(pooled, independently_pooled)
    configured = official.configure_model(SimpleNamespace(checkpoint=CHECKPOINT, amp_dtype='auto'))
    configured.model.cuda().eval()
    indices = official.select_layer_indices(21, layer_start=None, layer_end=None, output_depth=17, direction='forward')
    assert np.array_equal(indices, np.arange(2, 19))
    reader = official.FlatPatchReader(input_path=input_path, resolution='0', depth_axis_first=True, height=512, width=512, layer_indices=indices, output_depth=17, preprocessing=configured.preprocessing)
    blocks = official.iter_blocks((512, 512), 128, 64)
    dataset = official.FlatBlockDataset(reader=reader, blocks=blocks, patch_size=128, preprocessing=configured.preprocessing)
    weight = official.compute_importance_map_2d(patch_size=(128, 128), mode='hann').numpy()
    numerator = np.zeros((512, 512), np.float32)
    denominator = np.zeros_like(numerator)
    maps, patch_hashes = [], []
    for first in range(0, len(blocks), 8):
        items = [dataset[i] for i in range(first, min(first + 8, len(blocks)))]
        tensors = torch.stack([item[0] for item in items])
        with torch.inference_mode(), torch.autocast('cuda', dtype=torch.float16):
            probabilities = official.logits_to_probabilities(configured.model(tensors.cuda()), image_hw=(128, 128)).cpu().numpy()[:, 0]
        assert np.isfinite(probabilities).all() and probabilities.min() >= 0 and probabilities.max() <= 1
        for offset, (probability, item) in enumerate(zip(probabilities, items)):
            block = blocks[first + offset]
            raw = pooled[indices, block.y0:block.y0+128, block.x0:block.x0+128]
            # T9 preprocessing helper and official dataset must produce identical tensors.
            assert np.array_equal(item[0].numpy()[0], official.normalize_flat_patch(raw, configured.preprocessing))
            patch_hashes.append(hashlib.sha256(raw.tobytes()).hexdigest())
            maps.append(probability)
            ys, xs = slice(block.y0, block.y0+128), slice(block.x0, block.x0+128)
            numerator[ys, xs] += probability * weight
            denominator[ys, xs] += weight
    assert (denominator > 0).all()
    probability = numerator / denominator
    if replay:
        original = np.load(OUT / 'official_control_probability.npy')
        error = float(np.max(abs(original - probability)))
        recorded = json.loads((OUT / 'inference_receipt.json').read_text())
        assert patch_hashes == recorded['raw_patch_sha256']
        assert error <= .002
        (OUT / 'replay.json').write_text(json.dumps({'status': 'PASS', 'patches_replayed': len(blocks), 'all_raw_patches_exact': True, 'official_pooling_recomputed_exact': True, 'maximum_probability_error': error, 'tolerance': .002}, indent=2) + '\n')
        print('REPLAY PASS', error, flush=True)
        return
    np.save(OUT / 'official_control_probability.npy', probability)
    np.save(OUT / 'patch_probability_maps.npy', np.stack(maps))
    np.save(OUT / 'official_input.npy', pooled[indices])
    args = official.parse_args([str(input_path), str(CHECKPOINT), str(OUT / 'official_reference.tif'), '--batch-size', '8', '--num-workers', '0', '--no-compile'])
    official.infer_single_zarr(args=args, input_zarr=input_path, configured_model=configured, device=torch.device('cuda'), output_tiff=OUT / 'official_reference.tif', layer_direction='forward')
    reference = tifffile.imread(OUT / 'official_reference.tif')
    encoded = np.clip(probability * 255, 0, 255).astype(np.uint8)
    reference_error = int(np.max(abs(reference.astype(int)-encoded.astype(int))))
    assert reference_error <= 1
    receipt = {'checkpoint_sha256': hashlib.sha256(CHECKPOINT.read_bytes()).hexdigest(), 'source_z_half_open': [13, 97], 'pooling': 'official four-plane rounded mean', 'model_layer_indices': indices.tolist(), 'shape': list(probability.shape), 'patches': len(blocks), 'patch_shape': [17, 128, 128], 'stride': 64, 'batch_size': 8, 'amp': 'float16', 'preprocessing': configured.preprocessing, 'direction': 'forward', 'depth_offset': 0, 'tta': False, 'blend': 'official Hann', 'official_pooling_independent_exact': True, 'official_dataset_vs_T9_normalization_exact': True, 'official_reference_max_uint8_error': reference_error, 'raw_patch_sha256': patch_hashes}
    (OUT / 'inference_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('CONTROL INFERENCE COMPLETE', json.dumps({k:v for k,v in receipt.items() if k!='raw_patch_sha256'}), flush=True)

if __name__ == '__main__':
    main()
