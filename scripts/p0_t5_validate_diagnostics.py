"""Independent scalar replay of diagnostic lag correlations; no row writes."""
import csv,json
import numpy as np
from PIL import Image
from p0_t5_consecutive import OUT,dump
from p0_t5_evaluate import corr

def main():
    a=np.load(OUT/'stack.npy');pairs=json.loads((OUT/'adjacent_correlations.json').read_text())
    records=list(csv.DictReader((OUT/'local_lag_windows.csv').open()));max_error=0.
    for q in records:
        i=int(q['z_a'])-128;start=round(float(q['s_start_um'])/63.28);n=pairs[i]['shared_native_samples'];x=a[i,start:start+256];y=np.pad(a[i+1,:n],(64,64),constant_values=np.nan)
        rs=[]
        for lag in range(-64,65):
            r,count=corr(x,y[start+64+lag:start+64+lag+256],64);rs.append(r if r is not None else -2)
        best=min(range(129),key=lambda j:(-rs[j],abs(j-64),j));err=abs(rs[best]-float(q['best_r']));max_error=max(max_error,err)
        assert err<1e-10 and best-64==int(q['best_lag_samples'])
        assert abs(rs[64]-float(q['zero_r']))<1e-10
    metrics=json.loads((OUT/'lag_metrics.json').read_text());assert sum(metrics['raw_histogram'])==sum(metrics['best_histogram'])==49
    assert len(records)==metrics['local_windows']==1001
    display=json.loads((OUT/'display.json').read_text());mask=np.load(OUT/'mask.npy');p=np.load(OUT/'provenance.npy')
    contrast=np.zeros(a.shape,np.uint8);contrast[~mask]=np.rint(np.clip((a[~mask]-display['low'])/(display['high']-display['low']),0,1)*255).astype(np.uint8)
    assert np.array_equal(contrast,np.array(Image.open(OUT/'raw_stack_contrast.png')))
    palette=np.array([[0,0,0],[0,191,255],[255,0,255],[255,176,0]],np.uint8)
    assert np.array_equal(palette[p],np.array(Image.open(OUT/'provenance_stack.png')))
    assert np.array_equal(mask.astype(np.uint8)*255,np.array(Image.open(OUT/'mask_stack.png')))
    for name in ['lag_diagnostic.png','local_window_coherence.png','local_lag_examples.png']:
        with Image.open(OUT/name) as im:im.verify()
    dump(OUT/'diagnostic_validation.json',{'status':'PASS','local_windows_replayed':len(records),'lags_per_window':129,'max_scalar_vs_vectorized_error':max_error,'raw_histogram_count':49,'best_histogram_count':49,'contrast_provenance_mask_pngs':'exact','rows_modified':False})
    print('PASS:1001 local windows x129 lags scalar replay; histograms; contrast, provenance, mask PNGs',max_error)
if __name__=='__main__':main()
