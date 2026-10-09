"""Expose discarded independent ridge components in F03's frozen detector."""
import json
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p0_2_cross_z import ROOT, sample, save


def main():
    out=ROOT/'outputs/p0_2/validation'
    base=ROOT/'outputs/p0_2/run_02'
    w=next(w for w in json.loads((base/'windows.json').read_text()) if w['id']=='control_F03')
    ns=json.loads((base/'fragments.json').read_text())
    volume=np.load(ROOT/'outputs/orientation/ct_level4_zyx.npy',mmap_mode='r')
    origin=np.array(w['origin']);t=np.array(w['tangent']);normal=np.array([-t[1],t[0]])
    us=np.array(w['u']);vs=np.arange(-8,8.01,.25);grid=origin+us[:,None,None]*t+vs[None,:,None]*normal
    result=[];fig,axes=plt.subplots(1,3,figsize=(15,5))
    for ax,z in zip(axes,[145,146,147]):
        im=np.asarray(volume[z-128],float);response=gaussian_filter(im,.6)-gaussian_filter(im,2.)
        values=sample(response,grid);peaks=[]
        for i,profile in enumerate(values):
            idx,_=find_peaks(profile,prominence=1.,distance=4);idx=idx[profile[idx]>1.]
            row=[]
            for k in idx:
                denom=profile[k-1]-2*profile[k]+profile[k+1]
                correction=np.clip(.5*(profile[k-1]-profile[k+1])/denom,-.5,.5) if abs(denom)>1e-9 else 0.
                row.append((i,float(vs[k]+correction*.25)))
            peaks.append(row)
        graph={p:[] for row in peaks for p in row};breaks=[]
        for a,b in zip(peaks[:-1],peaks[1:]):
            if not a or not b:continue
            d=np.abs(np.array([p[1] for p in a])[:,None]-np.array([p[1] for p in b]))
            for i,p in enumerate(a):
                j=int(d[i].argmin());row=np.sort(d[i]);col=np.sort(d[:,j]);fail=[]
                if d[:,j].argmin()!=i:fail.append('not_mutual_nearest')
                if d[i,j]>.75:fail.append('normal_step_over_0.75')
                if (len(row)>1 and row[1]-row[0]<.125) or (len(col)>1 and col[1]-col[0]<.125):fail.append('competitor_margin_under_0.125')
                if fail:
                    if abs(p[1])<2 and abs(b[j][1])<2:
                        breaks.append({'u_from':float(us[p[0]]),'u_to':float(us[b[j][0]]),
                            'v_from':p[1],'v_to':b[j][1],'step':float(d[i,j]),'failures':fail})
                    continue
                graph[p].append(b[j]);graph[b[j]].append(p)
        seen=set();components=[]
        for seed in sorted(graph):
            if seed in seen:continue
            comp=set();todo=[seed]
            while todo:
                p=todo.pop()
                if p in comp:continue
                comp.add(p);todo.extend(q for q in graph[p] if q not in comp)
            seen.update(comp);ordered=sorted(comp);u=np.array([us[p[0]] for p in ordered]);v=np.array([p[1] for p in ordered])
            xy=origin+u[:,None]*t+v[:,None]*normal;arc=float(np.linalg.norm(np.diff(xy,axis=0),axis=1).sum())
            retained=len(comp)>=7 and arc>=3
            components.append({'u':u.tolist(),'v':v.tolist(),'xy':xy.tolist(),'points':len(comp),'arc_px':arc,
                'retained':retained,'near_expected_ridge':bool(np.median(abs(v))<2)})
        expected=[n for n in ns if n['z']==z and n['window']==w['id']]
        kept=[c for c in components if c['retained']]
        assert len(kept)==len(expected)
        for a,b in zip(kept,expected):assert np.allclose(a['xy'],b['xy'],rtol=0,atol=1e-12)
        result.append({'z':z,'components':components,'near_ridge_breaks':breaks})
        ax.imshow(im,cmap='gray',vmin=0,vmax=255,interpolation='nearest')
        for c in components:
            if not c['near_expected_ridge']:continue
            q=np.array(c['xy']);ax.plot(*q.T,'.-',color='red' if c['retained'] else 'magenta',lw=1,ms=4)
            ax.text(*q[len(q)//2],f"n={c['points']}",color='yellow',fontsize=9)
        ax.set(xlim=(264,277),ylim=(161,151),title=f'F03 Z={z}\nRed: retained; magenta: discarded peaks',xlabel='X',ylabel='Y')
    fig.tight_layout();fig.savefig(out/'F03_gap_discarded_components.png',dpi=150);plt.close(fig)
    save(out,'F03_gap_detector_diagnosis.json',result)
    for r in result:
        print('Z',r['z'],'near components',[(c['points'],c['arc_px'],c['retained']) for c in r['components'] if c['near_expected_ridge']])
        print('breaks',r['near_ridge_breaks'])


if __name__=='__main__':main()
