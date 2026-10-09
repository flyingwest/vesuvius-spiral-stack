"""Post-run diagnostics only; does not revise paths or scores."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/p0_t1'
volume=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')
summary=[]
for z in [128,144,146,160,170]:
    d=json.loads((OUT/f'slice_{z}_graph.json').read_text());pieces=list(map(np.array,d['pieces']))
    selected=[]
    fig,axs=plt.subplots(len(d['decisions']),2,figsize=(10,3.5*len(d['decisions'])),squeeze=False)
    for i,dec in enumerate(d['decisions']):
        x,y=dec['endpoint'];c=next((c for c in dec['candidates'] if c['chain']==dec['selected'] and not c['rejected_by']),None)
        selected.append({'index':i,'endpoint':dec['endpoint'],'decision':dec['decision'],'selected':c,'margin':dec['margin'],
                         'candidate_count':len(dec['candidates']),'viable_count':sum(not q['rejected_by'] for q in dec['candidates']),
                         'top_candidates':dec['candidates'][:3]})
        for ax in axs[i]:
            ax.imshow(volume[z-128],cmap='gray',vmin=0,vmax=255,interpolation='nearest')
            ax.set(xlim=(x-12,x+12),ylim=(y+12,y-12))
        axs[i,0].set_title(f'Z{z} decision {i}: raw CT at ({x:g},{y:g})')
        ax=axs[i,1]
        p=pieces[i];ax.plot(*p.T,'c-',lw=1)
        if c:
            q=pieces[i+1];ax.plot(*q.T,'lime',lw=1);ax.plot([x,q[0,0]],[y,q[0,1]],'m--')
        for alt in dec['candidates'][:6]:
            q=np.array(d['chains'][alt['chain']]);ax.plot(*q.T,color='orange',lw=.6,alpha=.65)
        ax.scatter(x,y,c='red',s=10)
        ax.set_title(f"{dec['decision']}; score={c['score']:.3f}; margin={dec['margin']:.3f}" if c else f"UNKNOWN; {len(dec['candidates'])} candidates")
    fig.tight_layout();fig.savefig(OUT/f'slice_{z}_decision_audit.png',dpi=140);plt.close(fig)
    summary.append({'z':z,'decisions':selected})
(OUT/'decision_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
for row in summary:
    print('Z',row['z'])
    for d in row['decisions']:
        c=d['selected']
        print(d['index'],d['endpoint'],d['decision'], 'to',c['xy'] if c else None,'score',round(c['score'],3) if c else None,'margin',round(d['margin'],3), 'candidates',d['candidate_count'],'viable',d['viable_count'])
