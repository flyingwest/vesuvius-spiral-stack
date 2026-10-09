"""Stack sealed T3 rows without resampling, shifts, warps or cross-Z tracking."""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t4';SRC=ROOT/'outputs/p0_t3/run_04'
ZS=[128,144,146,160,170];DS=63.28;DZ=126.56
PLAN={'slices':ZS,'source':str(SRC.relative_to(ROOT)),'coordinate':'s=0 at each unchanged T3 seed; shared distance grid, not registered material coordinate',
      'ds_um':DS,'dz_um':DZ,'registration':'none; offsets0; no stretch/resizing/reversal',
      'z':'five observed rows in increasing Z, nonuniform spacing preserved in coordinate array; no missing Z interpolation',
      'mask':'retain source masks and labels; pad missing row tails with NaN, maskTrue, provenance0 and coverageFalse',
      'raw_display':'round/clip original CT to uint8 0..255; missing black',
      'contrast_display':'one pooled valid-data 1st/99th percentile window; visualization only',
      'diagnostics':'Pearson raw and contiguous first differences; 256-sample windows stride128 minimum64 jointly valid; diagnostic lag ±64 minimum128 pairs; never apply shifts',
      'null':'199 deterministic circular shifts of second row within shared native prefix, shifts>=128 and <=N-128; roll mask with data; descriptive autocorrelation-preserving comparison, not independent replicates',
      'scope':'representation and coherence assessment only, no changes to T3 paths; no OCR/ink/mesh/learned correction'}
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def corr(a,b,minimum=32):
    good=np.isfinite(a)&np.isfinite(b);n=int(good.sum())
    if n<minimum:return None,n
    x=a[good]-a[good].mean();y=b[good]-b[good].mean();den=np.linalg.norm(x)*np.linalg.norm(y)
    return (float(x@y/den) if den>0 else None),n
def difference(a):
    return np.diff(a) # NaNs prevent differences across gaps.
def csvwrite(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def build():
    native=[];extent=[];metrics={r['z']:r for r in json.loads((SRC/'metrics.json').read_text())}
    for z in ZS:
        p=SRC/f'slice_{z}';a=np.load(p/'row.npy');mask=np.load(p/'mask.npy');prov=np.load(p/'provenance.npy');coords=np.load(p/'path_and_samples.npz')
        assert np.array_equal(np.isnan(a),mask) and np.array_equal(mask,np.isin(prov,[0,2]))
        assert np.array_equal(coords['s_um'],np.arange(len(a))*DS)
        assert np.array_equal(a,coords['ct'],equal_nan=True)
        native.append((a,mask,prov))
        extent.append({'z':z,'relative_z_um':(z-ZS[0])*DZ,'length_um':metrics[z]['total_arc_mm']*1000,'last_sample_s_um':(len(a)-1)*DS,
                       'native_samples':len(a),'valid_samples':int((~mask).sum()),'bridge_samples':int((prov==2).sum()),
                       'unresolved_samples':int((prov==0).sum()),'low_confidence_samples':int((prov==3).sum()),'s_offset_um':0.})
    n=max(len(a) for a,_,_ in native);stack=np.full((len(ZS),n),np.nan);mask=np.ones((len(ZS),n),bool);prov=np.zeros((len(ZS),n),np.uint8);coverage=np.zeros((len(ZS),n),bool)
    for i,(a,b,c) in enumerate(native):
        stack[i,:len(a)]=a;mask[i,:len(a)]=b;prov[i,:len(a)]=c;coverage[i,:len(a)]=True
        extent[i].update(padding_samples=n-len(a),valid_fraction_native=float((~b).mean()),valid_fraction_rectangle=float((~mask[i]).mean()))
    return stack,mask,prov,coverage,extent

def diagnostics(stack,extent):
    pairs=[];windows=[];lagrows=[];nullrows=[]
    for i in range(4):
        n=min(extent[i]['native_samples'],extent[i+1]['native_samples']);a=stack[i,:n];b=stack[i+1,:n]
        r,count=corr(a,b);dr,dcount=corr(difference(a),difference(b));null=[];dnull=[]
        for shift in np.unique(np.linspace(128,n-128,199).astype(int)):
            br=np.roll(b,int(shift));nr,nn=corr(a,br);nd,ndn=corr(difference(a),difference(br))
            null.append(nr);dnull.append(nd);nullrows.append({'z_a':ZS[i],'z_b':ZS[i+1],'shift_samples':int(shift),'correlation':nr,'valid_overlap':nn,'difference_correlation':nd})
        lags=[]
        for lag in range(-64,65):
            aa=a[max(0,-lag):min(n,n-lag)];bb=b[max(0,lag):min(n,n+lag)];lr,ln=corr(aa,bb,128)
            lagrows.append({'z_a':ZS[i],'z_b':ZS[i+1],'lag_samples':lag,'lag_um':lag*DS,'correlation':lr,'overlap':ln});lags.append((lr if lr is not None else -2,lag,ln))
        best=max(lags,key=lambda t:(t[0],-abs(t[1])))
        local=[]
        for start in range(0,n-255,128):
            w,wn=corr(a[start:start+256],b[start:start+256],64)
            windows.append({'z_a':ZS[i],'z_b':ZS[i+1],'s_start_um':start*DS,'s_end_um':(start+255)*DS,'correlation':w,'valid_overlap':wn});local.append(w)
        good=[v for v in local if v is not None]
        pairs.append({'z_a':ZS[i],'z_b':ZS[i+1],'delta_z_um':(ZS[i+1]-ZS[i])*DZ,'shared_native_samples':n,'shared_s_last_um':(n-1)*DS,
                      'joint_valid_samples':count,'joint_valid_fraction_shared':count/n,'zero_lag_pearson':r,'difference_pearson':dr,'difference_pairs':dcount,
                      'shift_null_median':float(np.median(null)),'shift_null_p95':float(np.percentile(null,95)),
                      'null_fraction_at_least_observed':float(np.mean(np.array(null)>=r)),
                      'difference_null_p95':float(np.percentile(dnull,95)),'difference_null_fraction_at_least_observed':float(np.mean(np.array(dnull)>=dr)),
                      'best_diagnostic_lag_samples':best[1],'best_diagnostic_lag_um':best[1]*DS,'best_lag_correlation':best[0],
                      'window_median_correlation':float(np.median(good)) if good else None,'window_positive_fraction':float(np.mean(np.array(good)>0)) if good else None})
    return pairs,windows,lagrows,nullrows

def render(stack,mask,prov,coverage,extent,pairs,windows,lagrows):
    raw=np.nan_to_num(stack,nan=0).round().clip(0,255).astype(np.uint8)
    lo,hi=np.percentile(stack[~mask],[1,99]);contrast=np.zeros_like(raw);contrast[~mask]=np.rint(np.clip((stack[~mask]-lo)/(hi-lo),0,1)*255).astype(np.uint8)
    palette=np.array([[0,0,0],[0,191,255],[255,0,255],[255,176,0]],np.uint8)
    for name,a in [('raw_stack.png',raw),('raw_stack_contrast.png',contrast),('mask_stack.png',mask.astype(np.uint8)*255),('provenance_stack.png',palette[prov])]:Image.fromarray(a).save(OUT/name)
    fig,axs=plt.subplots(3,1,figsize=(16,8),sharex=True)
    cm=plt.cm.gray.copy();cm.set_bad('#bf42ba');smax=(stack.shape[1]-1)*DS/1000
    for ax,a,cmap,vmin,vmax,title in [(axs[0],np.ma.masked_array(stack,mask),cm,0,255,'Fixed CT mapping; magenta = masked (including padding)'),(axs[1],np.ma.masked_array(stack,mask),cm,lo,hi,'Global 1–99 percentile display; rows are NOT registered'),(axs[2],palette[prov],None,None,None,'Provenance: black0 / cyan1 / magenta2 / orange3')]:
        ax.imshow(a,aspect='auto',interpolation='nearest',extent=[0,smax,4.5,-.5],cmap=cmap,vmin=vmin,vmax=vmax)
        ax.set(yticks=range(5),yticklabels=ZS,ylabel='observed Z index',title=title)
    axs[-1].set_xlabel('Seed-relative physical arc s (mm)');fig.tight_layout();fig.savefig(OUT/'stack_inspection.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(15,5))
    for i,r in enumerate(extent):
        axs[0].plot([0,r['length_um']/1000],[r['z'],r['z']],lw=8);axs[0].text(r['length_um']/1000+1,r['z'],f"{r['valid_fraction_native']:.1%} valid",va='center',fontsize=8)
    axs[0].set(xlabel='Physical path extent s (mm)',ylabel='Actual Z index (nonuniform)',title='Independent origins; no lateral adjustment',xlim=(0,290))
    sparse=np.full((43,stack.shape[1]),np.nan);sparse[np.array(ZS)-128]=stack
    axs[1].imshow(np.ma.masked_invalid(sparse),cmap=cm,vmin=0,vmax=255,aspect='auto',interpolation='nearest',extent=[0,smax,(42.5)*DZ/1000,-.5*DZ/1000])
    axs[1].set(xlabel='s (mm)',ylabel='Z distance from Z128 (mm)',title='Actual Z spacing; unsampled levels remain blank')
    fig.tight_layout();fig.savefig(OUT/'row_extent_diagnostic.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(13,7))
    for p in pairs:
        w=[r for r in windows if r['z_a']==p['z_a']];l=[r for r in lagrows if r['z_a']==p['z_a']];label=f"{p['z_a']}–{p['z_b']}"
        axs[0].plot([(r['s_start_um']+r['s_end_um'])/2000 for r in w],[r['correlation'] for r in w],'.-',label=label)
        axs[1].plot([r['lag_um']/1000 for r in l],[r['correlation'] for r in l],label=label)
    axs[0].set(xlabel='s window center (mm)',ylabel='Pearson r',title='Unshifted local correlations; gaps excluded');axs[1].set(xlabel='Diagnostic relative lag (mm), NOT applied',ylabel='Pearson r',title='Offset sensitivity; no row transformations')
    for ax in axs:ax.axhline(0,c='gray',lw=.5);ax.legend()
    fig.tight_layout();fig.savefig(OUT/'coherence_diagnostic.png',dpi=160);plt.close(fig)
    return {'global_contrast_percentiles':[1,99],'low':float(lo),'high':float(hi),'raw_mapping':'uint8 rint clipped0..255; NaN→0'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    for line in (ROOT/'outputs/p0_t3/SHA256SUMS').read_text().splitlines():
        h,p=line.split('  ',1);assert sha(ROOT/p)==h,p
    stack,mask,prov,coverage,extent=build()
    if args.verify:
        for name,a in [('stack',stack),('mask',mask),('provenance',prov),('coverage',coverage)]:assert np.array_equal(np.load(OUT/f'{name}.npy'),a,equal_nan=True)
        assert np.array_equal(np.isnan(stack),mask) and np.array_equal(mask,np.isin(prov,[0,2]))
        assert np.array_equal(np.array(Image.open(OUT/'raw_stack.png')),np.nan_to_num(stack,nan=0).round().astype('uint8'))
        assert np.array_equal(np.array(Image.open(OUT/'mask_stack.png')),mask.astype('uint8')*255)
        assert diagnostics(stack,extent)[0]==json.loads((OUT/'neighbor_correlations.json').read_text())
        for name,h in json.loads((OUT/'preservation.json').read_text()).items():assert sha(ROOT/name)==h,name
        print('PASS: T3 seal; exact stack/mask/provenance/coverage replay; raw PNG; correlations; prior preservation');return
    OUT.mkdir(exist_ok=False)
    prior=[p for base in ['outputs','docs','scripts'] for p in (ROOT/base).rglob('*') if p.is_file() and OUT not in p.parents]+[ROOT/'README.md',ROOT/'AGENTS.md']
    dump(OUT/'preservation.json',{str(p.relative_to(ROOT)):sha(p) for p in prior});dump(OUT/'plan.json',PLAN)
    for name,a in [('stack',stack),('mask',mask),('provenance',prov),('coverage',coverage),('z_indices',np.array(ZS)),('z_relative_um',(np.array(ZS)-128)*DZ),('s_um',np.arange(stack.shape[1])*DS)]:np.save(OUT/f'{name}.npy',a)
    csvwrite(OUT/'row_extents.csv',extent)
    pairs,windows,lags,null=diagnostics(stack,extent)
    dump(OUT/'neighbor_correlations.json',pairs)
    for name,r in [('neighbor_correlations',pairs),('window_correlations',windows),('lag_correlations',lags),('shift_null',null)]:csvwrite(OUT/f'{name}.csv',r)
    display=render(stack,mask,prov,coverage,extent,pairs,windows,lags);dump(OUT/'display.json',display)
    metrics={'rows':5,'width':stack.shape[1],'z_indices':ZS,'z_span_um':(ZS[-1]-ZS[0])*DZ,'z_gaps_um':(np.diff(ZS)*DZ).tolist(),
             'common_grid_last_s_um':(stack.shape[1]-1)*DS,'common_grid_cell_extent_um':stack.shape[1]*DS,
             'all_rows_native_overlap_last_s_um':(min(r['native_samples'] for r in extent)-1)*DS,
             'valid_samples':int((~mask).sum()),'valid_fraction_rectangle':float((~mask).mean()),'masked_fraction_rectangle':float(mask.mean()),
             'bridge_samples':int((prov==2).sum()),'bridge_fraction_rectangle':float((prov==2).mean()),
             'padding_samples':int((~coverage).sum()),'unresolved_native_samples':int(((prov==0)&coverage).sum()),
             'low_confidence_samples':int((prov==3).sum()),'low_confidence_fraction_rectangle':float((prov==3).mean()),
             'low_confidence_fraction_valid':float((prov==3).sum()/(~mask).sum()),'valid_fraction_native_domains':float((~mask).sum()/coverage.sum()),
             'all_five_joint_valid_samples':int((~mask).all(axis=0).sum()),'rows':extent,'pairs':pairs}
    dump(OUT/'metrics.json',metrics);print(json.dumps(metrics,indent=2),flush=True)

if __name__=='__main__':main()
