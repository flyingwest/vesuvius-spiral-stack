"""Post-run false-association indicators and visual audits; never changes tracks."""
from collections import Counter,defaultdict
import itertools
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_2_cross_z import ROOT,ZS,sha,save
from p0_2_local_ridges import associate,attach


def main():
    out=ROOT/'outputs/p0_2/audit';out.mkdir(exist_ok=True)
    base=ROOT/'outputs/p0_2/run_02'
    nodes=json.loads((base/'fragments.json').read_text());edges=json.loads((base/'associations.json').read_text())
    tracks=json.loads((base/'tracks.json').read_text());controls=json.loads((base/'controls.json').read_text())
    lookup={n['id']:n for n in nodes};membership={n:t['id'] for t in tracks for n in t['nodes']}
    accepted=[e for e in edges if e['status']=='ACCEPTED']
    byz={z:[n for n in nodes if n['z']==z] for z in ZS}
    regs=json.loads((ROOT/'outputs/p0_2/run_01/registration.json').read_text())
    # Audit identity consistency through independently rescored two-slice shortcuts.
    direct={}
    for z in ZS[:-2]:
        shift=np.array(regs[str(z)]['shift_xy'])+regs[str(z+1)]['shift_xy']
        for e in associate(byz[z],byz[z+2],shift):direct[(e['source'],e['target'])]=e
    direct_accepted={e['source']:e['target'] for e in direct.values() if e['status']=='ACCEPTED'}
    one_step={e['source']:e['target'] for e in accepted if e['z_target']==e['z_source']+1}
    twohop=[]
    for a,b in one_step.items():
        if b not in one_step:continue
        c=one_step[b];e=direct.get((a,c))
        competing=direct_accepted.get(a)
        twohop.append({'source':a,'via':b,'target':c,'track':membership[a],
                       'shortcut_status':e['status'] if e else 'NOT_ENUMERATED',
                       'shortcut_competing_target':competing if competing and competing!=c else None,
                       'shortcut_failures':e.get('decision_failures',e['hard_reject']) if e else ['not_enumerated']})
    save(out,'two_step_consistency.json',twohop)
    # Layer-order reversal audit at shared longitudinal stations within the same window.
    groups=defaultdict(list)
    for e in accepted:groups[(e['z_source'],e['z_target'],e['window'])].append(e)
    order_checks=[]
    for key,group in groups.items():
        for a,b in itertools.combinations(group,2):
            aa,ab,ba,bb=[lookup[n] for n in [a['source'],a['target'],b['source'],b['target']]]
            common=set(aa['stations'])&set(ab['stations'])&set(ba['stations'])&set(bb['stations'])
            if len(common)<3:continue
            maps=[dict(zip(n['stations'],n['v'])) for n in [aa,ab,ba,bb]]
            before=float(np.median([maps[0][k]-maps[2][k] for k in common]))
            after=float(np.median([maps[1][k]-maps[3][k] for k in common]))
            if min(abs(before),abs(after))<.5:continue
            order_checks.append({'pair':[a['source'],b['source']],'z_pair':key[:2],'window':key[2],
              'before_normal_order_px':before,'after_normal_order_px':after,'reversal':before*after<0})
    save(out,'layer_order_checks.json',order_checks)
    gap_examples=[]
    for e in accepted:
        if e['z_target']!=e['z_source']+2:continue
        a,b=lookup[e['source']],lookup[e['target']]
        common=np.intersect1d(a['stations'],b['stations']);ma=dict(zip(a['stations'],a['v']));mb=dict(zip(b['stations'],b['v']))
        mid_z=e['z_source']+1;matches=[]
        for n in byz[mid_z]:
            if n['window']!=e['window']:continue
            mn=dict(zip(n['stations'],n['v']));ks=[int(k) for k in common if int(k) in mn]
            if len(ks)<5:continue
            error=float(np.median([abs(mn[k]-(ma[k]+mb[k])/2) for k in ks]))
            matches.append({'node':n['id'],'median_expected_offset_error_px':error,'common_stations':len(ks)})
        matches.sort(key=lambda d:d['median_expected_offset_error_px'])
        detected=bool(matches and matches[0]['median_expected_offset_error_px']<=.75)
        gap_examples.append({'source':a['id'],'target':b['id'],'track':membership[a['id']],
          'window':e['window'],'z_missing':mid_z,'score':e['score'],'NCC':e['NCC'],
          'normal_error_px':e['normal_error_px'],
          'classification':'RETAINED_DETECTION_BUT_ASSOCIATION_GAP' if detected else 'NO_CLOSE_RETAINED_FRAGMENT',
          'middle_candidates':matches[:3],
          'physical_disappearance':'UNKNOWN; loss of detector/association evidence is not a proven tear'})
    save(out,'disappearance_reappearance.json',gap_examples)
    summary={'accepted_edges':len(accepted),'accepted_adjacent':sum(e['z_target']==e['z_source']+1 for e in accepted),
       'accepted_gap':len(gap_examples),'order_pairs_evaluable':len(order_checks),
       'order_reversals':sum(x['reversal'] for x in order_checks),
       'two_step_chains':len(twohop),'two_step_shortcut_statuses':dict(Counter(x['shortcut_status'] for x in twohop)),
       'two_step_conflicting_accepted_targets':sum(x['shortcut_competing_target'] is not None for x in twohop),
       'gap_classifications':dict(Counter(g['classification'] for g in gap_examples)),
       'accepted_NCC_min_median':[float(min(e['NCC'] for e in accepted)),float(np.median([e['NCC'] for e in accepted]))],
       'accepted_normal_error_max_px':float(max(e['normal_error_px'] for e in accepted)),
       'true_false_association_count':None,'true_false_association_rate':None,
       'limitation':'No registered physical layer oracle; gate compliance is not independent evidence of accuracy.'}
    save(out,'false_association_indicators.json',summary)
    volume=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')
    slice_manifest=[{'z_level4':z,'local_z':z-128,'native_nominal_z':z*16,'shape_yx':list(volume[z-128].shape),
                    'dtype':str(volume.dtype),'array_bytes':int(volume[z-128].nbytes),
                    'pixel_sha256':__import__('hashlib').sha256(volume[z-128].tobytes()).hexdigest()} for z in ZS]
    save(out,'slice_manifest.json',slice_manifest)
    chosen=[]
    for c in controls:
        if c['id'] in ['F02','F03']:
            track=next(t for t in tracks if t['id']==c['track'])
            zz=[137,140,144,148,150] if c['id']=='F02' else [143,144,145,146,147,148]
            fig,axes=plt.subplots(1,len(zz),figsize=(3.5*len(zz),4))
            point=np.array(c['anchor_xy'])
            for ax,z in zip(axes,zz):
                ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
                node=next((lookup[n] for n in track['nodes'] if lookup[n]['z']==z),None)
                if node:ax.plot(*np.array(node['xy']).T,color='red',lw=1.3)
                for n in byz[z]:
                    if n['window']=='control_'+c['id'] and (not node or n['id']!=node['id']):
                        ax.plot(*np.array(n['xy']).T,color='cyan',alpha=.7,lw=.6)
                ax.set(xlim=(point[0]-7,point[0]+7),ylim=(point[1]+7,point[1]-7),title=f"{c['id']} Z={z}"+(' UNKNOWN' if not node else ''))
            fig.tight_layout();fig.savefig(out/(c['id']+'_track_detail.png'),dpi=150);plt.close(fig)
    # A second disappearance example, preferring a retained middle detection.
    example=next(g for g in gap_examples if g['classification']=='RETAINED_DETECTION_BUT_ASSOCIATION_GAP' and not g['window'].startswith('control'))
    save(out,'selected_gap_example.json',example)
    point=np.array(lookup[example['source']]['xy']).mean(axis=0)
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,z in zip(axes,range(example['z_missing']-1,example['z_missing']+2)):
        ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
        for n in byz[z]:
            if n['window']==example['window']:
                color='red' if n['id'] in [example['source'],example['target']] else 'cyan'
                ax.plot(*np.array(n['xy']).T,color=color,lw=1)
        ax.set(xlim=(point[0]-8,point[0]+8),ylim=(point[1]+8,point[1]-8),title=f"{example['window']} Z={z}")
    fig.tight_layout();fig.savefig(out/'association_gap_example.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2))
    print('Selected gap',json.dumps(example,indent=2))
    for g in gap_examples:
        if g['window']=='control_F03':print('F03 gap',json.dumps(g,indent=2))


if __name__=='__main__':main()
