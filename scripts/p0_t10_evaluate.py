"""Frozen label discrimination and overlap metrics for the positive control."""
from pathlib import Path
import json
import numpy as np
from scipy.ndimage import distance_transform_edt
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t10'

def ranking(y, p, weights=None):
    order = np.argsort(-p, kind='stable')
    y, p = y[order], p[order]
    weights = np.ones(len(y)) if weights is None else weights[order]
    ends = np.r_[np.flatnonzero(np.diff(p)), len(p)-1]
    tp = np.cumsum(weights*y)[ends]
    fp = np.cumsum(weights*(1-y))[ends]
    if not tp[-1] or not fp[-1]:
        return np.nan, np.nan
    tpr, fpr = np.r_[0, tp/tp[-1]], np.r_[0, fp/fp[-1]]
    auc = float(np.sum((tpr[1:]+tpr[:-1])*np.diff(fpr)/2))
    ap = float(np.sum(np.diff(tpr)*tp/np.maximum(tp+fp, 1)))
    return auc, ap

def main():
    p = np.load(OUT/'official_control_probability.npy')
    with np.load(OUT/'selected_labels.npz') as archive:
        y = archive['inklabels'] > 0
        supervision = archive['supervision_mask'] > 0
    yy, xx = np.indices(p.shape)
    mask = supervision & (yy>=64) & (yy<448) & (xx>=64) & (xx<448)
    mask &= np.where(y, distance_transform_edt(y), distance_transform_edt(~y)) > 2
    labels, scores = y[mask].astype(int), p[mask]
    assert labels.sum() >= 1024 and (1-labels).sum() >= 1024
    auc, ap = ranking(labels, scores)
    blocks = ((yy//64)*8 + xx//64)[mask]
    ids = np.unique(blocks)
    rng = np.random.default_rng(2910)
    bootstrap = []
    for _ in range(1000):
        count = np.bincount(rng.choice(ids, len(ids), replace=True), minlength=64)
        bootstrap.append(ranking(labels, scores, count[blocks]))
    bootstrap = np.array(bootstrap)
    bootstrap_valid = np.isfinite(bootstrap).all(axis=1)
    np.save(OUT/'bootstrap_metrics.npy', bootstrap)
    positive, negative = scores[labels==1], scores[labels==0]
    pooled_sd = float(np.sqrt((positive.var()+negative.var())/2))
    effect = float((positive.mean()-negative.mean())/pooled_sd) if pooled_sd else 0.0
    predicted = scores > .5
    tp = int((predicted & (labels==1)).sum()); fp = int((predicted & (labels==0)).sum()); fn = int((~predicted & (labels==1)).sum()); tn = int((~predicted & (labels==0)).sum())
    prevalence = float(labels.mean())
    ci = np.percentile(bootstrap[bootstrap_valid,0], [2.5,97.5]).tolist() if bootstrap_valid.any() else [0.,1.]
    gates = {'auc_at_least_075': auc>=.75, 'AP_above_prevalence_by_015': ap>=prevalence+.15, 'effect_above_05': effect>.5, 'spatial_bootstrap_auc_lower_above_chance': ci[0]>.5}
    passed = all(gates.values())
    result = {'status': ('PASS' if passed else 'FAIL') if bootstrap_valid.all() else 'INCONCLUSIVE', 'invalid_bootstrap_draws':int((~bootstrap_valid).sum()), 'constant_probability_effect_convention': 'd=0 if pooled within-class SD is zero', 'shape': list(p.shape), 'evaluated_pixels': len(scores), 'ink_pixels': int(labels.sum()), 'nonink_pixels': int((1-labels).sum()), 'prevalence': prevalence, 'ROC_AUC': auc, 'PR_AUC_average_precision': ap, 'ROC_AUC_spatial_bootstrap95': ci, 'AP_spatial_bootstrap95': np.percentile(bootstrap[bootstrap_valid,1],[2.5,97.5]).tolist() if bootstrap_valid.any() else [0.,1.], 'effect_size_d': effect, 'positive_probability_quantiles': np.quantile(positive,[0,.05,.25,.5,.75,.95,1]).tolist(), 'negative_probability_quantiles': np.quantile(negative,[0,.05,.25,.5,.75,.95,1]).tolist(), 'positive_mean': float(positive.mean()), 'negative_mean': float(negative.mean()), 'threshold': .5, 'TP':tp,'FP':fp,'FN':fn,'TN':tn, 'IoU':tp/(tp+fp+fn), 'Dice':2*tp/(2*tp+fp+fn), 'gates':gates,'evaluation_mask':'supervision;64pixel crop border excluded; distance to opposite class>2px', 'bootstrap':'1000 draws of64x64spatial blocks, seed2910; same-crop uncertainty only', 'scope':'training-corpus positive control; not held-out generalization or causal proof about T9'}
    (OUT/'official_control_metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    np.save(OUT/'evaluation_mask.npy',mask)
    Image.fromarray(np.rint(np.clip(p,0,1)*65535).astype(np.uint16)).save(OUT/'official_control_probability.png')
    Image.fromarray(np.rint(np.clip((p-.25)/.5,0,1)*255).astype(np.uint8)).save(OUT/'official_control_display.png')
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    axs[0].imshow(np.ma.masked_where(~supervision,y),cmap='gray',vmin=0,vmax=1);axs[0].set_title('Producer ink labels; supervised only')
    axs[1].imshow(p,cmap='gray',vmin=.25,vmax=.75);axs[1].set_title('Fixed model display [.25,.75]')
    axs[2].hist(negative,bins=np.linspace(0,1,51),density=True,alpha=.6,label='non-ink');axs[2].hist(positive,bins=np.linspace(0,1,51),density=True,alpha=.6,label='ink');axs[2].legend();axs[2].set_title(f'AUC={auc:.3f}, AP={ap:.3f}')
    fig.tight_layout();fig.savefig(OUT/'registered_control_comparison.png',dpi=140);plt.close(fig)
    if passed:
        volume=np.load(OUT/'official_input.npy');fig,axs=plt.subplots(3,6,figsize=(18,8))
        for row,(z,start) in enumerate(zip([128,152,177],[128,192,256])):
            control=volume[:,start:start+128,start:start+128]
            with np.load(ROOT/f'outputs/p0_t9/smoke_{z}.npz') as archive: thread=archive['raw_patch'].copy()
            for group,(name,a) in enumerate([('Control',control),(f'T9 Z{z}',thread)]):
                for axis,plane in enumerate([a[8],a[:,64,:],a[:,:,64]]):
                    ax=axs[row,group*3+axis];ax.imshow(plane,cmap='gray',vmin=0,vmax=255,interpolation='nearest',aspect='auto');ax.set_title(f'{name}: '+['YX d=8','DX y=64','DY x=64'][axis]);ax.axhline(8,color='orange',lw=.5) if axis else None
        fig.tight_layout();fig.savefig(OUT/'representative_patch_comparison.png',dpi=140);plt.close(fig)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    main()
