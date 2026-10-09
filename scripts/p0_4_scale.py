"""Scale the frozen F02 extraction/association and arc-length raster without tuning."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_2_cross_z import ROOT, sha, registration, tracks_from
from p0_2_local_ridges import extract_window, associate, PLAN as TRACK_PLAN
from p0_3_first_stack import arc_geometry, sample_row, correlation

OUT=ROOT/'outputs/p0_4'
PLAN={'initial_z_inclusive':[128,177],'initial_attempts':50,'available_z_inclusive':[128,255],
 'expansion':'Extend upper bound by50 (capped255) only if seed-connected track reaches within one plane of current upper bound. Lower128 is a cache boundary, not a physical ending; inspect if reached.',
 'source':'outputs/orientation/ct_level4_zyx.npy','z_offset':128,'spacing_um':126.56,
 'identity_seed':'z144_control_F02_r02','prior_identity':'T0077',
 'extraction_and_association':'Import unchanged P0.2 run02 extract_window/associate; original fixed window and competition gates; same registration. No reseeding disconnected tracks.',
 'gap':'Original maximum one absent slice, original stricter gap gates; no missing CT row interpolation.',
 'raster':'P0.3 fixed core u=[1,3.5], anchor2.25, ds0.25px=31.64um, piecewise-linear arc XY/bilinear raw CT. Partial new support clipped, never enlarged; anchor must be present.',
 'low_coherence':'Eligible adjacent NCC<0.2, same P0.3 correlation eligibility.',
 'classification':'PASS needs >=50 supported rows with P0.3 coherence gates and no suspected switching; PARTIAL needs >=2 seed-connected supported runs of>=3 rows with local coherence; otherwise FAIL longitudinal scale-up (not automatic falsification of physical persistence).',
 'transition_audit':[[139,140],[146,147],[147,148]],
 'diagnostic_alignment':'Compare frozen arc alignment versus homologous window-u samples within same core; no correction applied. Descriptive only.',
 'limits':'No frozen gates changed, no longer-gap joins, no full scroll, mesh, OCR, ink, language or learned correction.'}


def read(p):return json.loads((ROOT/p).read_text())
def save(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def csv_write(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def segments(zs):
    result=[]
    for z in zs:
        z=int(z)
        if not result or z!=result[-1][-1]+1:result.append([])
        result[-1].append(z)
    return result


def distribution(values):
    a=np.array(values,float)
    if not len(a):return {'count':0,'min':None,'median':None,'mean':None,'max':None,'sd':None,'q05_q25_q75_q95':None}
    return {'count':len(a),'min':float(min(a)),'median':float(np.median(a)),'mean':float(np.mean(a)),
            'max':float(max(a)),'sd':float(np.std(a)),'q05_q25_q75_q95':np.percentile(a,[5,25,75,95]).tolist()}


def tracking(zs,volume,w):
    images={z:np.asarray(volume[z-128],float) for z in zs}
    byz={z:extract_window(images[z],z,w) for z in zs}
    regs={z:registration(images[z],images[z+1]) for z in zs[:-1]}
    adjacent=[]
    for z in zs[:-1]:adjacent.extend(associate(byz[z],byz[z+1],np.array(regs[z]['shift_xy'])))
    outgoing={e['source'] for e in adjacent if e['status']=='ACCEPTED'}
    incoming={e['target'] for e in adjacent if e['status']=='ACCEPTED'}
    gaps=[]
    for z in zs[:-2]:
        shift=np.array(regs[z]['shift_xy'])+regs[z+1]['shift_xy']
        gaps.extend(associate([n for n in byz[z] if n['id'] not in outgoing],
                              [n for n in byz[z+2] if n['id'] not in incoming],shift,gap=True))
    nodes=[n for z in zs for n in byz[z]];edges=adjacent+gaps
    tracks,mem=tracks_from(nodes,edges);selected=next(t for t in tracks if t['id']==mem[PLAN['identity_seed']])
    return nodes,edges,tracks,selected,regs


def compute():
    volume=np.load(ROOT/PLAN['source'],mmap_mode='r')
    w=next(w for w in read('outputs/p0_2/run_02/windows.json') if w['id']=='control_F02')
    attempts=[];upper=177
    while True:
        zs=list(range(128,upper+1));nodes,edges,tracks,selected,regs=tracking(zs,volume,w)
        attempts.append({'z_inclusive':[128,upper],'seed_connected_z':selected['z_observed']})
        if max(selected['z_observed'])<upper-1 or upper==255:break
        upper=min(upper+50,255)
    look={n['id']:n for n in nodes};chosen={look[k]['z']:look[k] for k in selected['nodes']}
    old=next(t for t in read('outputs/p0_2/run_02/tracks.json') if t['id']=='T0077')
    assert set(old['nodes'])<=set(selected['nodes']), 'Original identity changed: stop for diagnosis'
    oldnodes={n['id']:n for n in read('outputs/p0_2/run_02/fragments.json')}
    for n in nodes:
        if 137<=n['z']<=151:assert n==oldnodes[n['id']]
    geometries={};rowmeta=[]
    for z,n in chosen.items():
        u=np.array(n['u']);keep=(u>=1)&(u<=3.5);u=u[keep];xy=np.array(n['xy'])[keep]
        if len(u)<2 or not u[0]<=2.25<=u[-1]:continue
        s,origin=arc_geometry(xy,u,2.25);geometries[z]={'u':u,'xy':xy,'s':s,'origin':origin,'node':n['id']}
    oldgrid=np.load(ROOT/'outputs/p0_3/F02_s_px.npy')
    low=min([oldgrid[0]]+[-g['origin'] for g in geometries.values()]);high=max([oldgrid[-1]]+[g['s'][-1]-g['origin'] for g in geometries.values()])
    grid=np.arange(int(np.floor(low/.25)),int(np.ceil(high/.25))+1)*.25
    stack=np.full((len(zs),len(grid)),np.nan);mask=np.zeros_like(stack,dtype=bool);coords=np.full(stack.shape+(2,),np.nan)
    for i,z in enumerate(zs):
        g=geometries.get(z);n=chosen.get(z)
        r={'z':z,'local_z':z-128,'identity_supported':n is not None,'node':n['id'] if n else '',
           'raster_supported':g is not None,'reason':'ACCEPTED_IDENTITY_AND_CORE' if g else 'CORE_OR_ANCHOR_MISSING' if n else 'IDENTITY_UNKNOWN',
           'arc_length_px':None,'arc_length_um':None,'alignment_offset_px':None,'adjacent_shift_px':None,
           's_start_px':None,'s_end_px':None,'start_x':None,'start_y':None,'end_x':None,'end_y':None,'valid_samples':0}
        if g:
            stack[i],mask[i],coords[i]=sample_row(volume[z-128],g['xy'],g['s'],g['origin'],grid)
            r.update(arc_length_px=float(g['s'][-1]),arc_length_um=float(g['s'][-1]*126.56),alignment_offset_px=-g['origin'],
              adjacent_shift_px=(-g['origin']+geometries[z-1]['origin']) if z-1 in geometries else None,
              s_start_px=-g['origin'],s_end_px=float(g['s'][-1]-g['origin']),
              start_x=float(g['xy'][0,0]),start_y=float(g['xy'][0,1]),end_x=float(g['xy'][-1,0]),end_y=float(g['xy'][-1,1]),valid_samples=int(mask[i].sum()))
        rowmeta.append(r)
    pairs=[]
    for i,z in enumerate(zs[:-1]):
        c=correlation(stack[i],stack[i+1],grid)
        pairs.append({'z_source':z,'z_target':z+1,'both_rows_supported':z in geometries and z+1 in geometries,
            'samples':c['samples'],'span_px':c['span_px'],'eligible':c['eligible'],'NCC':c['NCC'],
            'descriptive_NCC':c['descriptive_NCC'],'mean_absolute_difference_CT':c['mean_absolute_difference_CT'],
            'low_coherence':c['NCC']<.2 if c['NCC'] is not None else None})
    known=set(selected['nodes']);supportededges=[e for e in edges if e['status']=='ACCEPTED' and e['source'] in known]
    bounds=[]
    for side,z in [('lower',min(chosen)),('upper',max(chosen))]:
        end=chosen[z]['id'];related=[e for e in edges if (e['target']==end and side=='lower') or (e['source']==end and side=='upper')]
        bounds.append({'side':side,'z':z,'node':end,'at_acquired_boundary':z in [128,255],
                      'continuation_candidates':related,'decision':'UNKNOWN'})
    # Direct two-step and relative-order indicators; never changes accepted edges.
    direct={}
    for z in zs[:-2]:
        for e in associate([n for n in nodes if n['z']==z],[n for n in nodes if n['z']==z+2],np.array(regs[z]['shift_xy'])+regs[z+1]['shift_xy']):
            direct[e['source'],e['target']]=e
    one={e['source']:e['target'] for e in supportededges if e['z_target']==e['z_source']+1};shortcuts=[]
    for a,b in one.items():
        if b not in one:continue
        c=one[b];other=[e['target'] for e in direct.values() if e['source']==a and e['status']=='ACCEPTED' and e['target']!=c]
        shortcuts.append({'source':a,'via':b,'target':c,'direct_status':direct.get((a,c),{}).get('status','NOT_ENUMERATED'),'conflicting_accepted_targets':other})
    lengths=[r['arc_length_px'] for r in rowmeta if r['raster_supported']]
    shifts=[r['adjacent_shift_px'] for r in rowmeta if r['adjacent_shift_px'] is not None]
    correlations=[p['NCC'] for p in pairs if p['eligible']];runs=segments(sorted(geometries));runresults=[]
    for run in runs:
        pp=[p for p in pairs if p['z_source'] in run and p['z_target'] in run];cc=[p['NCC'] for p in pp if p['eligible']]
        coherent=bool(len(run)>=3 and len(cc)>=.75*len(pp) and np.mean(cc)>=.5 and np.mean(np.array(cc)>0)>=.75)
        runresults.append({'z_inclusive':[run[0],run[-1]],'rows':len(run),'mean_NCC':float(np.mean(cc)) if cc else None,'coherent_by_P03_gate':coherent})
    candidates_count=Counter(e['status'] for e in edges)
    summary={'attempted_Z_levels':len(zs),'attempted_z_inclusive':[zs[0],zs[-1]],'seed_connected_observations':len(chosen),
        'seed_connected_z':sorted(chosen),'valid_rows':len(geometries),'masked_rows':len(zs)-len(geometries),
        'longest_continuous_supported_run':max(map(len,runs),default=0),'supported_runs':runresults,
        'correlation_distribution':distribution(correlations),'positive_adjacent_fraction':float(np.mean(np.array(correlations)>0)) if correlations else None,
        'low_coherence_transition_count':sum(p['low_coherence'] is True for p in pairs),
        'unevaluable_adjacent_pairs':sum(not p['eligible'] for p in pairs),'alignment_shift_distribution_px':distribution(shifts),
        'row_length_distribution_px':distribution(lengths),'missing_data_fraction':float(1-mask.mean()),
        'missing_samples':int((~mask).sum()),'valid_samples':int(mask.sum()),'shape_z_s':list(stack.shape),
        'association_counts':dict(candidates_count),'fragments_total':len(nodes),'all_tracklets':len(tracks),
        'accepted_seed_links':len(supportededges),'accepted_seed_gap_links':sum(e['z_target']-e['z_source']==2 for e in supportededges),
        'two_step_conflicting_targets':sum(bool(x['conflicting_accepted_targets']) for x in shortcuts),
        'suspected_layer_switch_count':None,'true_layer_switch_count':None,
        'unconnected_tracklets_never_reseeded':len(tracks)-1,'lower_cache_boundary_reached':min(chosen)==128,
        'new_seed_observations_beyond_P03':len(set(chosen)-set(range(137,151))),
        'conditional_extension_executed':len(attempts)>1,'attempt_history':attempts}
    # Fixed-scale raster must reproduce every old CT sample exactly.
    oldstack=np.load(ROOT/'outputs/p0_3/F02_stack.npy');oldzs=np.load(ROOT/'outputs/p0_3/F02_z.npy')
    cols=[list(grid).index(x) for x in oldgrid]
    np.testing.assert_array_equal(stack[np.ix_([zs.index(int(z)) for z in oldzs],cols)],oldstack)
    assert np.array_equal(np.isfinite(stack),mask)
    for e in supportededges:assert e['target'] in known
    return {'nodes':nodes,'edges':edges,'tracks':tracks,'selected':selected,'registration':{str(k):v for k,v in regs.items()},
            'boundary_candidates':bounds,'shortcut_audit':shortcuts,'alignment':rowmeta,'pairs':pairs,'metrics':summary,
            'stack':stack,'mask':mask,'sample_xy':coords,'s_px':grid,'z':np.array(zs)}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if not args.verify:
        OUT.mkdir(exist_ok=False);save('plan.json',PLAN);save('frozen_tracking_parameters.json',TRACK_PLAN)
        prior=[p for folder in ['outputs/orientation','outputs/p0_1','outputs/p0_1b','outputs/p0_2','outputs/p0_3','docs','scripts'] for p in (ROOT/folder).rglob('*') if p.is_file()]+[ROOT/'README.md',ROOT/'AGENTS.md']
        save('preservation_before.json',{str(p.relative_to(ROOT)):sha(p) for p in prior})
    r=compute()
    for k in ['stack','mask','sample_xy','s_px','z']:
        path=OUT/f'F02_extended_{k}.npy'
        if args.verify:np.testing.assert_array_equal(np.load(path),r[k])
        else:np.save(path,r[k],allow_pickle=False)
    for k in ['nodes','edges','tracks','selected','registration','boundary_candidates','shortcut_audit','alignment','pairs','metrics']:
        if args.verify:assert read(f'outputs/p0_4/{k}.json')==r[k],k
        else:save(k+'.json',r[k])
    for name,digest in read('outputs/p0_4/preservation_before.json').items():assert sha(ROOT/name)==digest,name
    if args.verify:
        save('verification.json',{'exact_replay':'PASS','prior_hashes':'unchanged','P03_numeric_overlap':'bit-identical','independent_extraction_old_nodes':'exact','missing_mask':'NaN consistent'})
    else:
        csv_write('F02_alignment.csv',r['alignment']);csv_write('F02_correlation_by_z.csv',r['pairs'])
        a=r['stack'];m=r['mask'];p=np.zeros_like(a,dtype=np.uint8);p[m]=np.rint(a[m]).astype(np.uint8)
        Image.fromarray(p).save(OUT/'F02_extended_stack.png');Image.fromarray(m.astype(np.uint8)*255).save(OUT/'F02_extended_stack_mask.png')
        fig,axes=plt.subplots(1,4,figsize=(15,10),gridspec_kw={'width_ratios':[1.5,1,2,1.5]})
        s=r['s_px']*126.56;z=r['z'];ext=[s[0]-15.82,s[-1]+15.82,z[-1]+.5,z[0]-.5]
        cmap=plt.get_cmap('gray').copy();cmap.set_bad('magenta')
        axes[0].imshow(np.ma.masked_invalid(a),cmap=cmap,vmin=0,vmax=255,interpolation='nearest',aspect='auto',extent=ext)
        for row in r['alignment']:
            if row['raster_supported']:axes[0].plot(np.array([row['s_start_px'],row['s_end_px']])*126.56,[row['z']]*2,'c|',ms=5)
        axes[0].set(title='Raw CT; missing magenta\ncyan endpoints',xlabel='s (µm)',ylabel='Z level4')
        axes[1].imshow(m,cmap='gray',vmin=0,vmax=1,interpolation='nearest',aspect='auto',extent=ext);axes[1].set(title='Valid mask',xlabel='s (µm)')
        axes[2].plot([p['NCC'] if p['NCC'] is not None else np.nan for p in r['pairs']],[p['z_source']+.5 for p in r['pairs']],'o-')
        axes[2].set(xlim=(-1,1),ylim=(z[-1]+.5,z[0]-.5),title='Adjacent correlation vs Z\nblank = unevaluable',xlabel='NCC');axes[2].axvline(.2,color='red',ls=':');axes[2].grid(alpha=.3)
        axes[3].plot([row['alignment_offset_px']*126.56 if row['raster_supported'] else np.nan for row in r['alignment']],z,'o-')
        axes[3].set(ylim=(z[-1]+.5,z[0]-.5),title='Geometric row offsets',xlabel='Offset (µm)');axes[3].grid(alpha=.3)
        fig.suptitle('P0.4 F02: fixed method, no disconnected track reseeding; displayed panels enlarged')
        fig.tight_layout();fig.savefig(OUT/'F02_extended_diagnostic.png',dpi=150);plt.close(fig)
    print(json.dumps(r['metrics'],indent=2))


if __name__=='__main__':main()
