"""Post-hoc CT validation of frozen local controls; no track modification."""
import csv
import json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_2_cross_z import ROOT, save, sha, sample


def main():
    out = ROOT / 'outputs/p0_2/validation'
    out.mkdir(exist_ok=False)
    base = ROOT / 'outputs/p0_2/run_02'
    read = lambda name: json.loads((base / name).read_text())
    nodes = read('fragments.json'); edges = read('associations.json')
    tracks = {t['id']: t for t in read('tracks.json')}
    controls = read('controls.json')[:4]
    ws = {w['id']: w for w in read('windows.json')}
    lookup = {n['id']: n for n in nodes}
    volume = np.load(ROOT / 'outputs/orientation/ct_level4_zyx.npy', mmap_mode='r')
    prior = [p for folder in ['outputs/orientation', 'outputs/p0_1', 'outputs/p0_1b',
                              'outputs/p0_2/run_01', 'outputs/p0_2/run_02', 'outputs/p0_2/audit', 'docs']
             for p in (ROOT / folder).rglob('*') if p.is_file()]
    prior += [ROOT / 'README.md', ROOT / 'AGENTS.md', Path(__file__)]
    save(out, 'preservation_before.json', {str(p.relative_to(ROOT)): sha(p) for p in prior})
    save(out, 'plan.json', {'type': 'retrospective validation, no parameter tuning or new accepted tracks',
        'radial_center_proxy_xy': [240, 140], 'spacing_um': 126.56,
        'radial_method': 'median radius on stations common to EVERY observation; never fragment centroid',
        'competitor_method': 'other retained fragments, >=3 common stations; median absolute normal separation',
        'gap_probe': 'interpolate endpoint geometry for diagnostic only; enumerate same detector peaks at each station',
        'physical_truth': 'visual CT consistency only; no registered layer oracle'})
    summary = {}; all_probes = []; all_edges = []
    shortcuts = json.loads((ROOT / 'outputs/p0_2/audit/two_step_consistency.json').read_text())
    orders = json.loads((ROOT / 'outputs/p0_2/audit/layer_order_checks.json').read_text())
    for c in controls:
        fid = c['id']; track = tracks[c['track']]; selected = [lookup[k] for k in track['nodes']]
        ids = set(track['nodes']); w = ws['control_' + fid]
        origin = np.array(w['origin']); tangent = np.array(w['tangent']); normal = np.array([-tangent[1], tangent[0]])
        common = sorted(set.intersection(*(set(n['stations']) for n in selected)))
        assert common
        records = []
        for n in selected:
            indices = [n['stations'].index(k) for k in common]
            xy = np.array(n['xy'])[indices]; vv = np.array(n['v'])[indices]
            alternatives = []
            for other in nodes:
                if other['z'] != n['z'] or other['window'] != n['window'] or other['id'] == n['id']: continue
                ks, ia, ib = np.intersect1d(n['stations'], other['stations'], return_indices=True)
                if len(ks) < 3: continue
                d = np.array(other['v'])[ib] - np.array(n['v'])[ia]
                alternatives.append({'id': other['id'], 'separation_px': float(np.median(abs(d))),
                                     'signed_separation_px': float(np.median(d)), 'common_stations': len(ks)})
            alternatives.sort(key=lambda x: x['separation_px'])
            incoming = next((e for e in edges if e['status']=='ACCEPTED' and e['target']==n['id']), None)
            records.append({'z': n['z'], 'node': n['id'], 'radius_px': float(np.median(np.linalg.norm(xy-[240,140], axis=1))),
                'normal_px': float(np.median(vv)), 'arc_length_px': n['arc_length_px'],
                'incoming_score': incoming['score'] if incoming else None,
                'incoming_NCC': incoming['NCC'] if incoming else None, 'competitors': alternatives})
        selected_edges = [e for e in edges if e['source'] in ids or e['target'] in ids]
        all_edges.extend(selected_edges)
        accepted = [e for e in selected_edges if e['status']=='ACCEPTED']
        check = [r for r in shortcuts if r['track']==track['id']]
        order = [r for r in orders if r['window']==w['id'] and any(k in ids for k in r['pair'])]
        radii = [r['radius_px'] for r in records]; vvals = [r['normal_px'] for r in records]
        near = [r['competitors'][0]['separation_px'] for r in records if r['competitors']]
        summary[fid] = {'track': track['id'], 'observed_z': track['z_observed'], 'missing_z': track['missing_z'],
            'observations': len(selected), 'span_slices': max(track['z_observed'])-min(track['z_observed'])+1,
            'common_u_px': [w['u'][k] for k in common], 'common_station_count': len(common),
            'radial_range_px': float(np.ptp(radii)), 'radial_net_px': float(radii[-1]-radii[0]),
            'normal_range_px': float(np.ptp(vvals)), 'normal_net_px': float(vvals[-1]-vvals[0]),
            'nearest_retained_competitor_min_px': min(near) if near else None,
            'nearest_retained_competitor_median_px': float(np.median(near)) if near else None,
            'max_accepted_normal_residual_px': max((e['normal_error_px'] for e in accepted),default=None),
            'accepted_gap_links': sum(e['z_target']-e['z_source']==2 for e in accepted),
            'touching_candidate_statuses': dict(Counter(e['status'] for e in selected_edges)),
            'two_step_shortcut_statuses': dict(Counter(r['shortcut_status'] for r in check)),
            'conflicting_accepted_shortcuts': [r for r in check if r['shortcut_competing_target']],
            'order_checks': len(order), 'order_reversals': sum(r['reversal'] for r in order), 'observations_detail': records}
        save(out, fid+'_measurements.json', summary[fid])
        zz = list(range(137,152)) if fid=='F02' else list(range(139,150)) if fid=='F03' else list(range(141,149))
        for overlay in [False, True]:
            fig, axes = plt.subplots(int(np.ceil(len(zz)/4)), 4, figsize=(16,4*int(np.ceil(len(zz)/4))))
            for ax,z in zip(axes.flat,zz):
                ax.imshow(volume[z-128], cmap='gray', vmin=0,vmax=255,interpolation='nearest')
                n = next((n for n in selected if n['z']==z),None)
                incoming = next((e for e in accepted if n and e['target']==n['id']),None)
                if overlay:
                    for other in nodes:
                        if other['z']==z and other['window']==w['id'] and other['id'] not in ids:
                            xy=np.array(other['xy']);ax.plot(*xy.T,color='cyan',lw=1,alpha=.85)
                            ax.text(*xy[len(xy)//2],other['id'].split('_')[-1],color='cyan',fontsize=7)
                    if n: ax.plot(*np.array(n['xy']).T,color='red',lw=1.5)
                label = 'unobserved / UNKNOWN' if n is None else 'provisional CT ridge'
                if incoming: label += f"; score={incoming['score']:.2f}"
                ax.set(xlim=(origin[0]-8,origin[0]+8),ylim=(origin[1]+8,origin[1]-8),
                    title=f'{fid} {track["id"]} Z={z}\n{label}',xlabel='X (level 4)',ylabel='Y')
            for ax in list(axes.flat)[len(zz):]:ax.axis('off')
            fig.suptitle(('Red: selected; cyan: competing retained proposals; scores are NOT probabilities' if overlay else 'Unobstructed CT, fixed 0–255 display; no interpolated missing geometry'),fontsize=13)
            fig.tight_layout(rect=(0,0,1,.97));fig.savefig(out/f'{fid}_{"overlays" if overlay else "raw"}.png',dpi=130);plt.close(fig)
        # Fixed geometry probes are diagnostic only, including detector's discarded peaks.
        reference = lookup[c['attachment']]
        probe_z = list(range(143,149)) if fid=='F03' else list(range(142,148))
        probe = np.array(reference['xy']); probe_stations = reference['stations']
        if fid=='F03':
            a=next(n for n in selected if n['z']==145);b=next(n for n in selected if n['z']==147)
            ks,ia,ib=np.intersect1d(a['stations'],b['stations'],return_indices=True)
            probe=(np.array(a['xy'])[ia]+np.array(b['xy'])[ib])/2;probe_stations=ks.tolist()
        vs = np.arange(-8,8.01,.25)
        fig,axes=plt.subplots(2,3,figsize=(15,8))
        for ax,z in zip(axes.flat,probe_z):
            im=np.asarray(volume[z-128],float);resp=gaussian_filter(im,.6)-gaussian_filter(im,2.)
            for k,p in zip(probe_stations,probe):
                profile=sample(resp,p[None,:]+vs[:,None]*normal)
                peaks,props=find_peaks(profile,prominence=1,distance=4)
                eligible=peaks[profile[peaks]>1]
                distance=float(min(abs(vs[eligible]))) if len(eligible) else None
                all_probes.append({'control':fid,'z':z,'station':int(k),'expected_xy':p.tolist(),
                    'response_at_expected':float(sample(resp,p[None,:])[0]),
                    'nearest_eligible_peak_distance_px':distance,
                    'eligible_peak_offsets_px':vs[eligible].tolist()})
                ax.plot(vs,profile,lw=.7,alpha=.6)
            ax.axvline(0,color='black',ls='--');ax.axhline(1,color='red',ls=':')
            ax.set(title=f'{fid} Z={z}: fixed-geometry probes',xlabel='Normal offset from probe (px)',ylabel='DoG CT units')
        fig.tight_layout();fig.savefig(out/f'{fid}_profile_diagnostics.png',dpi=130);plt.close(fig)
    save(out,'metrics.json',summary);save(out,'fixed_geometry_probes.json',all_probes)
    save(out,'control_candidate_associations.json',all_edges)
    for name,digest in json.loads((out/'preservation_before.json').read_text()).items():assert sha(ROOT/name)==digest,name
    save(out,'verification.json',{'prior_files_unchanged':len(prior),'tracks_modified':False,'new_associations':0,
        'fixed_station_metric_check':'nonempty intersection asserted for all four tracks'})
    print(json.dumps({f:{k:v for k,v in m.items() if k!='observations_detail'} for f,m in summary.items()},indent=2))


if __name__=='__main__': main()
