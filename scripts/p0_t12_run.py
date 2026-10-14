"""Local fit, explicit support inventory, and fixed CT geometry comparison."""
from pathlib import Path
import hashlib,json,sys,itertools,csv,os
import numpy as np
from p0_t12_core import fit,coordinates,render,ct
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t12';T11=ROOT/'outputs/p0_t11/complete_support_attempt'
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def read(p):return json.loads(p.read_text())
def bases():
    samples=read(OUT/'samples.json');frames=read(T11/'frames/records.json');answer=[]
    assert len(samples)==72 and samples==read(ROOT/'outputs/p0_t11/samples.json')
    for i,p in enumerate(samples):
        f=next(r for r in frames if r['sample']==i and r['sigma']==2)
        for key in ['P_xyz','z','s_index','tangent_xyz','normal_xyz']:assert f[key]==p[key],(i,key)
        with np.load(T11/f'frames/patch_{i}_2.npz') as a:old=a['old'].copy();base=a['new'].copy() if f['accepted'] else old.copy()
        t=f['t_s'] if f['accepted'] else p['tangent_xyz'];y=f['image_y'] if f['accepted'] else [0,0,1];n=f['n'] if f['accepted'] else p['normal_xyz']
        answer.append((p,t,y,n,old,base,f['accepted']))
    return answer
def identity_check():
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:assert os.environ.get(key)=='1',key
    frozen=read(OUT/'PRE_EVIDENCE_FREEZE.json')['sha256']
    for name,h in frozen.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
def main():
    stage=sys.argv[1];identity_check();bs=bases()
    if stage in ['fit','fit_replay']:
        dest=OUT/stage;dest.mkdir(exist_ok=False);records=[]
        for i,(p,t,y,n,old,base,accepted) in enumerate(bs):
            result=fit(base);records.append({'sample':i,**p,'T11_orientation_accepted':accepted,'t':t,'y':y,'n':n,**result})
            print('fit',i+1,'/72',flush=True)
        volume=np.load(ROOT/'outputs/p0_t10/official_input.npy');control=[]
        for start in [128,192,256]:
            a=volume[:,start:start+128,start:start+128];control.append({'origin_yx':[start,start],'quality':ct.quality(a),'depth_position_proxy':float(np.argmax(a.mean((1,2)))-8),'fit':fit(a)})
        dump(dest/'records.json',records);dump(dest/'control.json',control)
        if stage=='fit_replay':
            assert records==read(OUT/'fit/records.json') and control==read(OUT/'fit/control.json');dump(OUT/'fit_replay.json',{'status':'PASS','target_fits':72,'control_fits':3,'all_records_exact':True})
        return
    records=read(OUT/'fit/records.json')
    if stage=='coverage':
        patches=[];missing={}
        for r in records:
            for method in ['plane','quadratic']:
                f=r[method]
                if not f['accepted']:continue
                c=coordinates(r['P_xyz'],r['t'],r['y'],r['n'],f['coefficients'],f['degree']);lo=np.floor(c.reshape(3,-1).min(1)).astype(int)[::-1];hi=(np.floor(c.reshape(3,-1).max(1)).astype(int)+2)[::-1]
                keys=list(itertools.product(*[range(a//128,(b-1)//128+1) for a,b in zip(lo,hi)]));absent=[]
                for key in keys:
                    name='_'.join(map(str,key))
                    if not any((ROOT/f'data/raw/{phase}/native_chunks/{name}').exists() for phase in ['p0_t9','p0_t11','p0_t12']):
                        absent.append(list(key));missing.setdefault(name,[]).append({'sample':r['sample'],'method':method})
                patches.append({'sample':r['sample'],'method':method,'bounds_zyx_half_open':[lo.tolist(),hi.tolist()],'keys':[list(k) for k in keys],'missing':absent})
        dump(OUT/'coverage.json',{'status':'PASS' if not missing else 'BLOCKED','fitted_patches':len(patches),'patches':patches,'missing_chunk_dependencies':missing});print('coverage',len(patches),'fit patches',len(missing),'missing chunks');return
    assert read(OUT/'coverage.json')['status']=='PASS','Missing CT support; do not pad'
    assert stage in ['geometry','geometry_replay'];dest=OUT/stage;dest.mkdir(exist_ok=False);metrics=[];arrays={};localrows=[];offsetrows=[]
    for i,(p,t,y,n,old,base,orientation) in enumerate(bs):
        r=records[i];qs={'T9':ct.quality(old),'T11':ct.quality(base)};actual={};accepted={}
        for policy in ['plane','plane_then_quadratic']:
            f=r['quadratic'] if policy=='plane_then_quadratic' and r['quadratic']['accepted'] else (r['plane'] if r['plane']['accepted'] else None)
            a=render(p['P_xyz'],t,y,n,f) if f else base.copy();actual[policy]=a;accepted[policy]=f is not None;qs[policy]=ct.quality(a)
            shift=float(f['coefficients'][0])*9.6 if f else 0.;point=np.asarray(p['P_xyz'])+np.asarray(n)*shift/7.91
            offsetrows.append({'sample':i,'z':p['z'],'s_index':p['s_index'],'policy':policy,'accepted':bool(f),'signed_offset_um':shift,'displacement_um':abs(shift),'original_xyz':json.dumps(p['P_xyz']),'fitted_center_xyz':json.dumps(point.tolist()),'normal_xyz':json.dumps(n),'confidence_proxy':float(1/(1+f['fit_residual_holdout'])) if f else 0.})
            localrows.append({'sample':i,'policy':policy,'accepted':bool(f),'degree':f['degree'] if f else 0,'ridge_fraction':f.get('ridge_fraction',0) if f else r['plane'].get('ridge_fraction',0),'center_residual':f['sheet_center_residual'] if f else None,'heldout_surface_residual':f['fit_residual_holdout'] if f else None,'conflict_fraction':f['competitor_conflict_fraction'] if f else None})
        hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in actual.items()}
        np.savez_compressed(dest/f'patch_{i}.npz',T9=old,T11=base,**actual)
        metrics.append({'sample':i,'quality':qs,'accepted':accepted,'hashes':hashes,'identity_reliable':r.get('identity_reliable',False)})
    gates=[]
    for policy in ['plane','plane_then_quadratic']:
        med={name:{k:float(np.median([a['quality'][name][k] for a in metrics])) for k in metrics[0]['quality']['T9']} for name in ['T9','T11',policy]}
        reliable=sum(r.get('identity_reliable',False) for r in records);accepted=sum(r['accepted'][policy] for r in metrics)
        checks={'at_least70pct_reliable_identity':reliable/72>=.7,'at_least70pct_accepted_fits':accepted/72>=.7,'depth_spread_improvement10pct':med[policy]['depth_peak_spread']<=.9*med['T9']['depth_peak_spread'],'no_increased_band_crossings':med[policy]['band_crossings']<=med['T9']['band_crossings'],'empty_volume_nondegradation':med[policy]['empty_fraction']<=med['T9']['empty_fraction']}
        localpass=[r['accepted'][policy] and r['quality'][policy]['depth_peak_spread']<=.9*r['quality']['T9']['depth_peak_spread'] and r['quality'][policy]['band_crossings']<=r['quality']['T9']['band_crossings'] and r['quality'][policy]['empty_fraction']<=r['quality']['T9']['empty_fraction'] for r in metrics]
        gates.append({'method':policy,'pass':bool(all(checks.values())),'checks':checks,'reliable_identity_count':reliable,'accepted_fit_count':accepted,'local_all_quality_pass_fraction':float(np.mean(localpass)),'medians':med})
    selected=next((r['method'] for r in gates if r['pass']),None);decision='D' if sum(r.get('identity_reliable',False) for r in records)/72<.7 else 'C'
    dump(dest/'metrics.json',metrics);dump(dest/'geometry_gate.json',{'status':'PASS' if selected else 'FAIL','selected_method':selected,'methods':gates,'stop_result':None if selected else decision,'jump_scope':'No ground-truth identity oracle; unique-nearest, bounded-displacement, residual, competitor and continuity safeguards. No systematic jump evidence is not proof of physical identity.'})
    for name,rows in [('local_fit_metrics.csv',localrows),('centering_offsets.csv',offsetrows)]:
        with (dest/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    if stage=='geometry_replay':
        assert metrics==read(OUT/'geometry/metrics.json');assert read(dest/'geometry_gate.json')==read(OUT/'geometry/geometry_gate.json');dump(OUT/'patch_replay.json',{'status':'PASS','samples':72,'policies':2,'all_patch_hashes_and_metrics_exact':True})
    print(json.dumps(read(dest/'geometry_gate.json')),flush=True)
if __name__=='__main__':main()
