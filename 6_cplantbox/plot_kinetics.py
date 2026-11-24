"""Growth over time: measured (PimpiBig2 tracings, days 1-5, control; medians scaled to the panel day-5 size) vs CPlantBox (v2)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, rsa_model as m
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
ORDER=['C','B','dT']; NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}; COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}
P=json.load(open('params_v2.json')); K=json.load(open('kinetics_targets.json'))
fig,axs=plt.subplots(1,4,figsize=(14,3.6)); days=np.arange(1,6)
for t in ORDER:
    F=json.load(open(f'fit_v2_{t}.json')); t0=F['x']['t0']; tg=F['target']
    sims=[]
    for s in range(12):
        p=m.build(P[t]['3D'],100+s,stochastic=True); tt=0.; row=[]
        for d in days:
            while tt<t0+d-1e-9: st=min(0.1,t0+d-tt); p.simulate(st,False); tt+=st
            o=m.measure(p); L=o['lats']; L=L[L[:,1]>0.1]; row.append((o['MR'],len(L),L[:,1].sum()))
        sims.append(row)
    S=np.array(sims)
    for k,(ax,key) in enumerate(zip(axs[:3],['MR','nLR','tot'])):
        ax.plot(days,tg[key],'o-',color=COL[t],lw=2,label=f'{NAME[t]}: measured' if k==0 else None)
        md=np.median(S[:,:,k],0); ax.plot(days,md,'--',color=COL[t],lw=1.5,label=f'{NAME[t]}: CPlantBox' if k==0 else None)
        ax.fill_between(days,np.percentile(S[:,:,k],10,0),np.percentile(S[:,:,k],90,0),color=COL[t],alpha=0.12,lw=0)
    a=np.arange(1,5); axs[3].plot(a,[K[t]['rate_by_age'][str(i)] for i in a],'o-',color=COL[t],lw=2)
    x=F['x']; axs[3].plot(a,x['rL']+x['sL']*a,'--',color=COL[t],lw=1.5)
for ax,lab in zip(axs,['Main-root length (cm)','Number of laterals (> 1 mm)','Total lateral length (cm)','Lateral elongation rate (cm d⁻¹)']): ax.set_ylabel(lab)
for ax in axs[:3]: ax.set_xlabel('Imaging day'); ax.set_xticks(days)
axs[3].set_xlabel('Days since the lateral was first seen'); axs[3].set_xticks([1,2,3,4])
axs[0].legend(fontsize=7,frameon=False)
fig.suptitle('Growth over time: PimpiBig2 tracings of the same plants on days 1–5 (control; solid, medians scaled to the panel day-5 size) vs CPlantBox (dashed, median and 10–90% of 12 replicates)',fontsize=9)
fig.tight_layout(rect=(0,0,1,0.92))
for e in ('png','pdf'): fig.savefig(f'figs_v2/kinetics_sim_vs_measured.{e}',dpi=200)
