"""Retrospective single-XY fragment audit; conservative bounded continuation evidence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/p0_1b'
CENTER=np.array([240.,140.])
# Retrospective source-image audit intervals, inclusive original CSV indices.
# These are construction annotations, NOT independent manual segmentation truth.
EXCLUDED=[(0,12,'UNTRACEABLE','Open apex: merged ridges / insufficient unique continuation'),
          (44,63,'LIKELY_SWITCH','Abrupt cross-band movement on upper-right diagonal'),
          (89,104,'AMBIGUOUS','Right shoulder bend has competing ridge continuations'),
          (209,262,'LIKELY_SWITCH','Lower arc moves between near-parallel bands'),
          (304,357,'LIKELY_SWITCH','Lower-left diagonal crosses dark trough / competing bands'),
          (408,426,'UNTRACEABLE','Left-upper contact: ridge identity not resolved'),
          (456,479,'AMBIGUOUS','Upper-left zigzag: multiple ridge continuations'),
          (509,555,'UNTRACEABLE','Upper apex merges and strong guide-shift sensitivity')]
PLAN={'analysis':'Retrospective semi-manual fragmentation, exploratory scores, not calibrated probabilities',
      'z_level4':144,'source':'outputs/p0_1/source_slice.png','spacing_mm':.12656,
      'excluded_intervals':[{'start':a,'end':b,'category':c,'reason':r} for a,b,c,r in EXCLUDED],
      'minimum_fragment_points':15,'continuation_advance_px':6.,'normal_search_px':[-8,8,.25],
      'candidate_peak_prominence':1.5,'candidate_min_ridge_response':1.,
      'candidate_peak_separation_px':1.,'spacing_profile_px':[-8,8,.25],
      'spacing_peak_prominence':2.,'endpoint_fit_points':13,'ridge_sigma':[.6,2.],
      'gradient_sigma':1.,'structure_tensor_sigma':1.2,'center_proxy_xy':CENTER.tolist(),
      'score':'equal mean of seven [0,1] heuristic evidence components; not probability',
      'acceptance':{'total_score':.8,'winner_margin':.12,'coherence_min':.35,
                    'tangent_angle_max_deg':25,'curvature_delta_max':.08,
                    'endpoint_abs_curvature_max':.15,'lateral_spacing_fraction_max':.4,
                    'radial_spacing_fraction_max':.5,'texture_correlation_min':.6,
                    'mean_CT_difference_max':25,'bridge_positive_ridge_fraction_min':.9,
                    'bridge_abs_orientation_cosine_min':.85},
      'oracle':'Unavailable at candidate; do not equate construction annotations with ground truth',
      'scope':'Single PNG only; no neighboring z, no gap filling, no accepted link appended automatically'}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def sample(a,xy):
    xy=np.atleast_2d(xy)
    return map_coordinates(a,[xy[:,1],xy[:,0]],order=1,mode='nearest',prefilter=False)


def unit(x):return x/max(np.linalg.norm(x),1e-12)


def field(image):
    smooth=gaussian_filter(image,1.)
    gy,gx=np.gradient(smooth)
    xx=gaussian_filter(gx*gx,1.2);yy=gaussian_filter(gy*gy,1.2)
    xy=gaussian_filter(gx*gy,1.2)
    return xx,xy,yy


def orientation(fields,p,reference):
    xx,xy,yy=[float(sample(a,p)[0]) for a in fields]
    values,vectors=np.linalg.eigh([[xx,xy],[xy,yy]])
    t=vectors[:,0]
    if t@reference<0:t=-t
    coherence=(values[1]-values[0])/(values.sum()+1e-12)
    return t,float(coherence)


def endpoint_geometry(history):
    h=history[-13:]
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(h,axis=0),axis=1))];s-=s[-1]
    coef=np.polyfit(s,h,2)
    v=coef[1];acc=2*coef[0]
    tangent=unit(v)
    k=(v[0]*acc[1]-v[1]*acc[0])/max(np.linalg.norm(v)**3,1e-12)
    return h[-1],tangent,float(k)


def spacing(response,p,normal):
    offsets=np.arange(-8,8.01,.25)
    profile=sample(response,p+offsets[:,None]*normal)
    peaks,_=find_peaks(profile,prominence=2.,distance=4)
    peaks=peaks[profile[peaks]>1.]
    if len(peaks)<2:return None,offsets,profile,peaks
    return float(np.median(np.diff(offsets[peaks]))),offsets,profile,peaks


def texture(image,p,t):
    n=np.array([-t[1],t[0]])
    return sample(image,p+np.linspace(-3,3,13)[:,None]*n)


def score_endpoint(image,response,fields,history,eid):
    p,t,k=endpoint_geometry(history)
    n=np.array([-t[1],t[0]])
    predicted=p+6*t+.5*np.clip(k,-.15,.15)*36*n
    layer_spacing,offsets,profile,peaks=spacing(response,p,n)
    candidate_profile=sample(response,predicted+offsets[:,None]*n)
    qpeaks,_=find_peaks(candidate_profile,prominence=1.5,distance=4)
    qpeaks=qpeaks[candidate_profile[qpeaks]>1.]
    tsource,csource=orientation(fields,p,t)
    source_texture=texture(image,p,t)
    candidates=[]
    for j,idx in enumerate(qpeaks):
        offset=float(offsets[idx]);q=predicted+offset*n
        tq,coherence=orientation(fields,q,t)
        ta,_=orientation(fields,q-2*tq,tq);tb,_=orientation(fields,q+2*tq,tq)
        kq=float(np.arctan2(ta[0]*tb[1]-ta[1]*tb[0],ta@tb)/4)
        tangent_angle=float(np.degrees(np.arccos(np.clip(t@tq,-1,1))))
        curvature_delta=abs(k-kq)
        radial_delta=float(np.linalg.norm(q-CENTER)-np.linalg.norm(predicted-CENTER))
        tq_profile=texture(image,q,tq)
        corr=float(np.corrcoef(source_texture,tq_profile)[0,1]) if min(source_texture.std(),tq_profile.std())>1e-9 else 0.
        mean_difference=float(abs(source_texture.mean()-tq_profile.mean()))
        std_difference=float(abs(source_texture.std()-tq_profile.std()))
        distance=float(np.linalg.norm(q-p))
        # Candidate connector is a diagnostic cubic Hermite curve, not an accepted gap fill.
        u=np.linspace(0,1,33)[:,None]
        bridge=(2*u**3-3*u**2+1)*p+(u**3-2*u**2+u)*distance*t+(-2*u**3+3*u**2)*q+(u**3-u**2)*distance*tq
        bt=np.gradient(bridge,axis=0);bt/=np.linalg.norm(bt,axis=1)[:,None]
        ridge_fraction=float(np.mean(sample(response,bridge)>1.))
        orientation_agreement=float(np.mean([abs(orientation(fields,b,d)[0]@d) for b,d in zip(bridge,bt)]))
        local_spacing_q,_,_,_=spacing(response,q,np.array([-tq[1],tq[0]]))
        order_score=float(np.exp(-.5*(radial_delta/(.35*layer_spacing))**2)) if layer_spacing else 0.
        spacing_score=float(np.exp(-abs(np.log(local_spacing_q/layer_spacing)))) if layer_spacing and local_spacing_q else 0.
        components={'proximity':float(np.exp(-distance/12)),
                    'tangent':float((1+np.cos(np.radians(tangent_angle)))/2),
                    'curvature':float(np.exp(-curvature_delta/.08)),
                    'radial_order_proxy':order_score,'layer_spacing':spacing_score,
                    'gradient_ridge':float(orientation_agreement*np.sqrt(csource*coherence)*ridge_fraction),
                    'texture_intensity':float(max(0,corr)*np.exp(-mean_difference/25)*np.exp(-std_difference/20))}
        reasons=[]
        if not layer_spacing or not local_spacing_q:reasons.append('spacing_unknown')
        if tangent_angle>25:reasons.append('tangent_mismatch')
        if abs(k)>.15:reasons.append('endpoint_curvature_out_of_bound')
        if curvature_delta>.08:reasons.append('curvature_mismatch')
        if layer_spacing and abs(offset)>.4*layer_spacing:reasons.append('lateral_layer_order_risk')
        if layer_spacing and abs(radial_delta)>.5*layer_spacing:reasons.append('radial_order_risk')
        if min(csource,coherence)<.35:reasons.append('orientation_uncertain')
        if corr<.6 or mean_difference>25:reasons.append('texture_mismatch')
        if ridge_fraction<.9:reasons.append('bridge_untraceable')
        if orientation_agreement<.85:reasons.append('bridge_orientation_mismatch')
        total=float(np.mean(list(components.values())))
        if total<.8:reasons.append('score_below_threshold')
        candidates.append({'id':f'{eid}_C{j+1:02}','xy':q.tolist(),'normal_offset_px':offset,
             'score':total,'components':components,'distance_px':distance,
             'tangent_angle_deg':tangent_angle,'candidate_curvature':kq,'curvature_delta':curvature_delta,
             'radial_delta_px':radial_delta,'candidate_spacing_px':local_spacing_q,
             'coherence':coherence,'texture_correlation':corr,'mean_CT_difference':mean_difference,
             'std_CT_difference':std_difference,'bridge_ridge_fraction':ridge_fraction,
             'bridge_orientation_agreement':orientation_agreement,'rejection_reasons':reasons,
             'oracle_correct':None,'bridge_xy':bridge.tolist()})
    candidates.sort(key=lambda c:(-c['score'],c['id']))
    margin=candidates[0]['score']-candidates[1]['score'] if len(candidates)>1 else None
    accepted=bool(candidates and not candidates[0]['rejection_reasons'] and (margin is None or margin>=.12))
    return {'id':eid,'xy':p.tolist(),'outward_tangent':t.tolist(),'endpoint_curvature':k,
            'endpoint_coherence':csource,'local_spacing_px':layer_spacing,
            'spacing_peak_offsets_px':offsets[peaks].tolist(),'predicted_xy':predicted.tolist(),
            'candidate_count':len(candidates),'candidates':candidates,'winner_margin':margin,
            'decision':'CONTINUE' if accepted else 'STOP','identity_status':'UNKNOWN',
            'decision_correct_oracle':None,'oracle_correct_continuation':None,
            'stop_reasons':([] if accepted else (candidates[0]['rejection_reasons']+(['competing_candidate_scores'] if margin is not None and margin<.12 else []) if candidates else ['no_ridge_candidate_in_bounded_search']))}


def compute():
    image=np.asarray(Image.open(ROOT/PLAN['source'])).astype(float)
    table=np.genfromtxt(ROOT/'outputs/p0_1/trace_coordinates.csv',delimiter=',',names=True)
    path=np.c_[table['x_level4'],table['y_level4']]
    assert len(path)==556
    response=gaussian_filter(image,.6)-gaussian_filter(image,2.)
    fields=field(image)
    keep=np.ones(len(path),bool)
    for a,b,_,_ in EXCLUDED:keep[a:b+1]=False
    runs=np.flatnonzero(np.diff(np.r_[False,keep,False].astype(int))).reshape(-1,2)
    fragments=[];endpoints=[]
    for start,stop in runs:
        if stop-start<15:continue
        fid=f'F{len(fragments)+1:02}'
        frag=path[start:stop]
        length=float(np.linalg.norm(np.diff(frag,axis=0),axis=1).sum())
        reasons={}
        for side,index in [('start',start-1),('end',stop)]:
            reasons[side]=next((f'{cat}: {why}' for a,b,cat,why in EXCLUDED if a<=index<=b),'Input path boundary; no closure assumed')
        fragments.append({'id':fid,'start_index':int(start),'end_index':int(stop-1),
                          'start_xy':frag[0].tolist(),'end_xy':frag[-1].tolist(),
                          'points':len(frag),'arc_length_px':length,'arc_length_mm':length*.12656,
                          'confidence':'PROVISIONAL_VISUAL_PLAUSIBILITY; not calibrated or oracle-certified',
                          'termination_reasons':reasons})
        endpoints.extend([score_endpoint(image,response,fields,frag[::-1],fid+'_start'),
                          score_endpoint(image,response,fields,frag,fid+'_end')])
    summary={'plausible_fragments_provisional':len(fragments),'oracle_verified_fragments':None,
             'retained_points':sum(f['points'] for f in fragments),'source_points':len(path),
             'retained_arc_length_px':sum(f['arc_length_px'] for f in fragments),
             'excluded_review_intervals':len(EXCLUDED),
             'suspected_switch_regions':sum(c=='LIKELY_SWITCH' for _,_,c,_ in EXCLUDED),
             'true_layer_jump_count':None,'false_continuation_count':None,
             'false_layer_switch_rate':None,'correctly_detected_physical_interruptions':None,
             'missed_continuations':None,'oracle_resolved_endpoints':0,
             'ambiguous_identity_endpoints':len(endpoints),
             'multiple_geometric_candidate_endpoints':sum(e['candidate_count']>1 for e in endpoints),
             'CONTINUE_decisions':sum(e['decision']=='CONTINUE' for e in endpoints),
             'STOP_decisions':sum(e['decision']=='STOP' for e in endpoints),
             'candidate_count':sum(e['candidate_count'] for e in endpoints),
             'note':'Null denotes not measurable without independent registered layer identity; STOP correctness is not inferred from score.'}
    return image,path,fragments,endpoints,summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify:
        image,path,fragments,endpoints,summary=compute()
        for filename,data in [('fragments.json',fragments),('endpoints.json',endpoints),('metrics.json',summary)]:
            assert json.loads((OUT/filename).read_text())==data,filename
        before=json.loads((OUT/'preservation_before.json').read_text())
        for name,digest in before.items():assert sha(ROOT/name)==digest,name
        covered=[]
        for f in fragments:
            covered.extend(range(f['start_index'],f['end_index']+1))
            coords=np.genfromtxt(OUT/(f['id']+'_coordinates.csv'),delimiter=',',names=True)
            assert np.array_equal(np.c_[coords['x'],coords['y']],path[f['start_index']:f['end_index']+1])
        assert len(set(covered))==len(covered)==summary['retained_points']
        for e in endpoints:
            for c in e['candidates']:assert all(0<=v<=1.0000001 for v in c['components'].values())
            if e['decision']=='CONTINUE':assert not e['candidates'][0]['rejection_reasons']
        result={'deterministic_replay':'exact','prior_files_preserved':len(before),
                'fragment_partition':'nonoverlapping exact subsets of original CSV',
                'oracle_truth_claims':'none','single_input_plane_only':True,
                'script_sha256':sha(Path(__file__))}
        (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2));return
    assert not (OUT/'metrics.json').exists(),'Do not overwrite prior analysis'
    (OUT/'analysis_plan.json').write_text(json.dumps(PLAN,indent=2)+'\n')
    image,path,fragments,endpoints,summary=compute()
    for filename,data in [('fragments.json',fragments),('endpoints.json',endpoints),('metrics.json',summary)]:
        (OUT/filename).write_text(json.dumps(data,indent=2)+'\n')
    colors=plt.get_cmap('tab10')
    fig,ax=plt.subplots(figsize=(15,9));ax.imshow(image,cmap='gray',vmin=0,vmax=255,interpolation='nearest')
    ax.plot(*path.T,color='red',alpha=.3,lw=.7,label='original (not joined in new result)')
    for k,f in enumerate(fragments):
        q=path[f['start_index']:f['end_index']+1]
        ax.plot(*q.T,color=colors(k),lw=2,label=f['id'])
        ax.scatter(*q[[0,-1]].T,color=colors(k),s=14)
        mid=q[len(q)//2];ax.text(*mid,f['id'],color=colors(k),fontsize=11)
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(q,axis=0),axis=1))]
        with (OUT/(f['id']+'_coordinates.csv')).open('w') as file:
            writer=csv.writer(file);writer.writerow(['original_index','x','y','z','s_fragment_px','s_fragment_mm'])
            writer.writerows((idx,*point,144,s,s*.12656) for idx,point,s in zip(range(f['start_index'],f['end_index']+1),q,arc))
    ax.set(xlim=(165,295),ylim=(170,100),title='P0.1b: provisional fragments; no inferred bridge is appended')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'fragments_overview.png',dpi=160);plt.close(fig)
    # Small multiples at every endpoint show the source separately from proposed bridges.
    for e in endpoints:
        p=np.array(e['xy']);fig,axes=plt.subplots(1,2,figsize=(10,5))
        for ax in axes:
            ax.imshow(image,cmap='gray',vmin=0,vmax=255,interpolation='nearest')
            ax.set_xlim(p[0]-12,p[0]+12);ax.set_ylim(p[1]+12,p[1]-12)
            ax.scatter(*p,c='yellow',s=22)
        for c in e['candidates']:
            q=np.array(c['bridge_xy']);axes[1].plot(*q.T,lw=.7)
            axes[1].text(*q[-1],f"{c['id'].split('_')[-1]} {c['score']:.2f}",fontsize=7,color='orange')
        axes[0].set_title(e['id']+' source (no oracle)')
        axes[1].set_title(e['decision']+' / identity UNKNOWN')
        fig.tight_layout();fig.savefig(OUT/(e['id']+'_candidates.png'),dpi=130);plt.close(fig)
    with (OUT/'candidate_scores.csv').open('w') as file:
        columns=['endpoint','candidate','x','y','score','decision','proximity','tangent','curvature','radial_order_proxy','layer_spacing','gradient_ridge','texture_intensity','rejection_reasons']
        writer=csv.DictWriter(file,fieldnames=columns);writer.writeheader()
        for e in endpoints:
            for c in e['candidates']:writer.writerow({'endpoint':e['id'],'candidate':c['id'],
                'x':c['xy'][0],'y':c['xy'][1],'score':c['score'],'decision':e['decision'],
                **c['components'],'rejection_reasons':';'.join(c['rejection_reasons'])})
    print(json.dumps(summary,indent=2))
    print(json.dumps(fragments,indent=2))
    for e in endpoints:print(e['id'],e['decision'],e['candidate_count'],e['stop_reasons'])


if __name__=='__main__':main()
