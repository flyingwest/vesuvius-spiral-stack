"""Post-analysis support audit; never retune fragments or continuation thresholds."""
import csv
import importlib.metadata
import json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
from p0_1b_fragments import ROOT, OUT, EXCLUDED, field, orientation, sample, score_endpoint, sha


def main():
    image=np.asarray(Image.open(ROOT/'outputs/p0_1/source_slice.png')).astype(float)
    table=np.genfromtxt(ROOT/'outputs/p0_1/trace_coordinates.csv',delimiter=',',names=True)
    path=np.c_[table['x_level4'],table['y_level4']]
    base=np.load(ROOT/'outputs/p0_1/optimizer_path_xy.npy')
    base_s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(base,axis=0),axis=1))]
    variants=[]
    for name in ['sigma_0.4','sigma_0.8','guide_shift_-0.5','guide_shift_+0.5']:
        q=np.load(ROOT/'outputs/p0_1'/f'{name}_path_xy.npy')
        variants.append(np.column_stack([np.interp(table['s_pixels'],base_s,q[:,i]) for i in range(2)]))
    dev=np.max([np.linalg.norm(q-path,axis=1) for q in variants],axis=0)
    response=gaussian_filter(image,.6)-gaussian_filter(image,2.)
    fields=field(image)
    fragments=json.loads((OUT/'fragments.json').read_text())
    support=[];interior=[]
    for f in fragments:
        start,end=f['start_index'],f['end_index']
        q=path[start:end+1];tangent=np.gradient(q,axis=0);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
        agreement=np.array([abs(orientation(fields,p,t)[0]@t) for p,t in zip(q,tangent)])
        p95=float(np.percentile(dev[start:end+1],95))
        support.append({'fragment':f['id'],'positive_ridge_fraction':float(np.mean(sample(response,q)>1)),
            'mean_abs_ridge_orientation_cosine':float(agreement.mean()),
            'variant_max_displacement_px':float(dev[start:end+1].max()),
            'variant_p95_displacement_px':p95,
            'confidence':'LOW' if p95>2 else 'MODERATE',
            'use_status':'WITHHELD_UNCERTAIN_HYPOTHESIS' if p95>2 else 'PROVISIONALLY_USEFUL',
            'confidence_basis':'Retrospective visual plausibility plus >2 px p95 sensitivity warning; not calibrated or independently validated'})
        midpoint=len(q)//2
        probe=score_endpoint(image,response,fields,q[:midpoint+1],f['id']+'_interior_control')
        probe['control_note']='Known retained-path interior, NOT independent physical truth; diagnostic of blanket abstention.'
        interior.append(probe)
    (OUT/'fragment_support.json').write_text(json.dumps(support,indent=2)+'\n')
    (OUT/'interior_control_probes.json').write_text(json.dumps(interior,indent=2)+'\n')
    endpoints=json.loads((OUT/'endpoints.json').read_text())
    summary={'useful_provisional_fragments':sum(x['use_status']=='PROVISIONALLY_USEFUL' for x in support),
        'withheld_uncertain_fragments':sum(x['use_status']=='WITHHELD_UNCERTAIN_HYPOTHESIS' for x in support),
        'confidence_rule_status':'Post-analysis exploratory warning, not a retrospectively preregistered threshold',
        'new_links_appended':0,'false_new_links_appended':0,
        'false_link_rate':None,'rate_denominator_accepted_links':0,
        'true_inherited_layer_jump_count':None,'correctly_detected_physical_interruptions':None,
        'missed_true_continuations':None,'endpoint_decision_biological_correctness':'UNKNOWN at every endpoint',
        'interior_control_continue':sum(p['decision']=='CONTINUE' for p in interior),
        'interior_control_stop':sum(p['decision']=='STOP' for p in interior),
        'retained_arc_length_mm':sum(f['arc_length_mm'] for f in fragments),
        'useful_arc_length_mm':sum(f['arc_length_mm'] for f,s in zip(fragments,support) if s['use_status']=='PROVISIONALLY_USEFUL'),
        'source_reads':'source_slice.png, P0.1 CSV and preserved same-slice variant NPYs only',
        'environment':{p:importlib.metadata.version(p) for p in ['numpy','scipy','matplotlib','pillow']}}
    (OUT/'evaluation.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (OUT/'review_intervals.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['start_index','end_index','start_x','start_y','end_x','end_y','category','reason','physical_cause_oracle'])
        for a,b,cat,reason in EXCLUDED:writer.writerow([a,b,*path[a],*path[b],cat,reason,'UNKNOWN'])
    # Render all candidate evidence into compact report-ready tables.
    lines=['## All endpoint candidates\n',
           'Scores are heuristic evidence, not probabilities. Candidate IDs are stable within each endpoint.\n',
           '| Endpoint | Candidate XY and total score (descending) | Decision | Oracle-correct / biological correctness |',
           '| --- | --- | --- | --- |']
    for e in endpoints:
        labels='; '.join(f"{c['id'].split('_')[-1]} ({c['xy'][0]:.2f}, {c['xy'][1]:.2f})={c['score']:.3f}" for c in e['candidates'])
        lines.append(f"| {e['id']} | {labels} | {e['decision']} | UNKNOWN / UNKNOWN |")
    (OUT/'endpoint_table.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary,indent=2))
    print('Interior decisions',[(x['id'],x['decision']) for x in interior])


if __name__=='__main__':main()
