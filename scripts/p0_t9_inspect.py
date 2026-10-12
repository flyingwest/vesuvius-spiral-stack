"""Fixed, evenly spaced visual windows and frozen sensitivity summaries."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / 'outputs/p0_t9'

def main():
    probability = np.load(OUT / 'ink_probability.npy')
    fig, axes = plt.subplots(4, 1, figsize=(15, 8))
    starts = np.linspace(0, probability.shape[1] - 512, 4, dtype=int)
    for ax, start in zip(axes, starts):
        ax.imshow(np.ma.masked_invalid((probability[:, start:start + 512] - .25) / .5), vmin=0, vmax=1, cmap='gray', interpolation='nearest', aspect=2, extent=[int(start), int(start) + 512, 177.5, 127.5])
        ax.set_title(f'Fixed evenly spaced window: samples {start}..{start + 511}')
        ax.set_ylabel('Z')
    fig.tight_layout()
    fig.savefig(OUT / 'fixed_window_inspection.png', dpi=140)
    plt.close(fig)
    probes = json.loads((OUT / 'orientation_depth_sensitivity.json').read_text())
    primary = np.array([row['primary'] for row in probes])
    result = {}
    for key in ['reverse', 'minus2', 'plus2']:
        alternative = np.array([row[key] for row in probes])
        result[key] = {'mean_absolute_difference': float(abs(primary - alternative).mean()), 'max_absolute_difference': float(abs(primary - alternative).max()), 'pearson_r': float(np.corrcoef(primary, alternative)[0, 1])}
    result['probes'] = len(probes)
    (OUT / 'sensitivity_summary.json').write_text(json.dumps(result, indent=2) + '\n')

if __name__ == '__main__':
    main()
