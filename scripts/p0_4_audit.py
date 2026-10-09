"""Diagnose scale-test boundaries and three existing low-coherence transitions."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_2_cross_z import ROOT, sample, ncc
from p0_3_first_stack import arc_geometry

OUT=ROOT/'outputs/p0_4'
def read(name):return json.loads((OUT/name).read_text())
def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    nodes=read('nodes.json');look={n['id']:n for n in nodes};selected=read('selected.json')
    chosen={look[k]['z']:look[k] for k in selected['nodes']};edges=read('edges.json');pairs=read('pairs.json')
    volume=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')
    summaries=[]
    for za,zb in [(139,140),(146,147),(147,148)]:
        a,b=chosen[za],chosen[zb];e=next(e for e in edges if e['source']==a['id'] and e['target']==b['id'])
        pp=next(p for p in pairs if p['z_source']==za)
        ks,ia,ib=np.intersect1d(a['u'],b['u'],return_indices=True);keep=(ks>=1)&(ks<=3.5);u=ks[keep];ia=ia[keep];ib=ib[keep]
        xy_a=np.array(a['xy'])[ia];xy_b=np.array(b['xy'])[ib]
        ca=sample(np.asarray(volume[za-128],float),xy_a);cb=sample(np.asarray(volume[zb-128],float),xy_b)
        sa,oa=arc_geometry(xy_a,u,2.25);sb,ob=arc_geometry(xy_b,u,2.25)
        result={'z_pair':[za,zb],'source':a['id'],'target':b['id'],'arc_aligned_NCC':pp['NCC'],
          'homologous_u_NCC_descriptive':ncc(ca,cb),'homologous_samples':len(u),
          'homologous_u_px':u.tolist(),'source_xy':xy_a.tolist(),'target_xy':xy_b.tolist(),
          'source_CT_at_u':ca.tolist(),'target_CT_at_u':cb.tolist(),
          'mean_CT_change_at_u':float(np.mean(cb-ca)),
          'std_CT_at_u':[float(np.std(ca)),float(np.std(cb))],
          'max_anchored_arc_residual_px':float(max(abs((sa-oa)-(sb-ob)))),
          'acceptance_evidence':e,'core_has_detector_gap':False,
          'material_change':'Raw CT changes observed; biological cause not established.',
          'alignment_error':'Pure origin translation cannot guarantee material correspondence; diagnostic u comparison is not ground truth.',
          'detector_fragmentation':'Not the direct cause: full frozen u=1..3.5 core exists in both observations.',
          'possible_layer_switch':'No positive switch evidence from accepted CT review; true identity unverified without oracle.',
          'insufficient_support':'Original bounded shift scan cannot retain required 2px support at all shifts.',
          'primary_cause':'UNKNOWN: CT/material variation versus geometric sampling/correspondence error cannot be separated.'}
        summaries.append(result)
        fig,axes=plt.subplots(1,3,figsize=(13,4))
        center=(xy_a.mean(axis=0)+xy_b.mean(axis=0))/2
        for ax,z,node in zip(axes[:2],[za,zb],[a,b]):
            ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
            for candidate in nodes:
                if candidate['z']!=z:continue
                q=np.array(candidate['xy']);ax.plot(*q.T,color='red' if candidate['id']==node['id'] else 'cyan',lw=1)
            keep=(np.array(node['u'])>=1)&(np.array(node['u'])<=3.5);q=np.array(node['xy'])[keep]
            ax.plot(*q.T,'yo',ms=3)
            ax.set(xlim=(center[0]-7,center[0]+7),ylim=(center[1]+7,center[1]-7),title=f'Z={z}: selected red / candidates cyan\nyellow = frozen raster core',xlabel='X',ylabel='Y')
        axes[2].plot(u,ca,'o-',label=f'Z={za}');axes[2].plot(u,cb,'o-',label=f'Z={zb}')
        axes[2].set(title=f'Geometry-u comparison only\nNCC={ncc(ca,cb):.3f}; no correction applied',xlabel='Frozen window u (px)',ylabel='Raw CT units');axes[2].legend()
        fig.tight_layout();fig.savefig(OUT/f'transition_{za}_{zb}.png',dpi=150);plt.close(fig)
    save('low_coherence_transition_audit.json',summaries)
    fig,axes=plt.subplots(2,3,figsize=(13,8))
    origin=np.array([270.91603599636494,132.74196951005848])
    for ax,z in zip(axes.flat,[135,136,137,150,151,152]):
        ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
        for n in nodes:
            if n['z']!=z:continue
            q=np.array(n['xy']);color='red' if z in chosen and n['id']==chosen[z]['id'] else 'cyan'
            ax.plot(*q.T,color=color,lw=1.2);ax.text(*q[len(q)//2],n['id'].split('_')[-1],color=color,fontsize=8)
        ax.set(xlim=(origin[0]-8,origin[0]+8),ylim=(origin[1]+8,origin[1]-8),
               title=f'Z={z}: '+('F02 accepted (red)' if z in chosen else 'identity UNKNOWN; candidates only'),xlabel='X',ylabel='Y')
    fig.tight_layout();fig.savefig(OUT/'boundary_CT_overlays.png',dpi=150);plt.close(fig)
    save('visual_review_scope.json',{'selected_z_reviewed':list(range(137,151)),
        'prior_full_CT_review':'outputs/p0_2/validation/F02_overlays.png and F02_raw.png, unchanged input track',
        'new_CT_review':[135,136,137,139,140,146,147,148,150,151,152],
        'new_accepted_seed_links':0,'suspected_switches_in_accepted_seed_track':0,
        'true_switch_count':None,'unresolved_new_boundary_directions':2,
        'interpretation':'No extrapolation of physical identity to disconnected tracklets.'})
    for r in summaries:
        print(r['z_pair'],'arc NCC',r['arc_aligned_NCC'],'u NCC',r['homologous_u_NCC_descriptive'],
              'CTmeanchange',r['mean_CT_change_at_u'],'arcresidual',r['max_anchored_arc_residual_px'],
              'texture NCC',r['acceptance_evidence']['NCC'],'normal residual',r['acceptance_evidence']['normal_error_px'])


if __name__=='__main__':main()
