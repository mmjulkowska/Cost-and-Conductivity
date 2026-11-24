"""Simulated (equal total root length) vs measured (control, day 5): main root length, lateral number, lateral length by depth, angles."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, pandas as pd, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT=os.environ.get('OUT','sims'); FIG=os.environ.get('FIG','figs')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
ORDER=['C','B','dT']; NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}; COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}
R=pd.read_csv(f'{OUT}/replicate_metrics.csv')
d=pd.read_csv(CP.inp('01_panel_root_geometry.csv')); d=d[d.condition=='control']; d['TRL']=d.MR_length+d.total_lateral_root_length
Lm=pd.read_csv('rsml_laterals.csv'); Lm=Lm[Lm.cond=='control']
fig,axs=plt.subplots(2,4,figsize=(13,6.6))
mets=[('TRL','TRL','Total root length (cm)'),('MR','MR_length','Main root length (cm)'),('nLR','n_lateral_roots','Number of laterals'),('meanLR','mean_lateral_root_length','Mean lateral length (cm)'),
      ('LR_ins_angle','LR_insertion_angle_deg','Lateral insertion angle (°)'),('LR_angle_from_vertical','LR_angle_from_vertical_deg','Lateral angle from vertical (°)')]
for ax,(cs,cm,lab) in zip(axs.flat[:6],mets):
    for i,t in enumerate(ORDER):
        m=d.loc[d.architectural_type==t,cm].dropna(); b=ax.boxplot([m],positions=[i*3],widths=0.8,showfliers=False,patch_artist=True)
        b['boxes'][0].set(facecolor='#dddddd',edgecolor='#777'); [x.set(color='#777') for k in ('whiskers','caps') for x in b[k]]; b['medians'][0].set(color='k')
        for g,off,mk in (('plate',1,'o'),('3D',1.8,'^')):
            s=R.loc[(R.type==t)&(R.geometry==g),cs]; ax.scatter(i*3+off+np.random.default_rng(1).uniform(-.15,.15,len(s)),s,s=12,color=COL[t],marker=mk,alpha=.8,lw=0)
    ax.set_xticks([i*3+0.9 for i in range(3)]); ax.set_xticklabels([NAME[t].replace(' telephone','\ntelephone') for t in ORDER],fontsize=8); ax.set_ylabel(lab)
axs.flat[0].set_title('grey box: measured (control, day 5)\n● simulated 2D plate   ▲ simulated 3D',fontsize=8,loc='left')
# lateral length vs insertion depth
ax=axs.flat[6]; gs=axs[1,2].get_gridspec(); axs.flat[6].remove(); axs.flat[7].remove(); ax=fig.add_subplot(gs[1,2:])
bins=np.arange(0,9.01,0.5); mid=bins[:-1]+0.25
for t in ORDER:
    m=Lm[Lm.type==t]; n=m.arbor.nunique(); h,_=np.histogram(m.ins,bins,weights=m.len); ax.plot(mid,h/n,color=COL[t],lw=2,label=f'{NAME[t]}: measured (n = {n})')
    sims=[]
    for s in range(1,int(R.rep.max())+1):
        pl=np.load(f'{OUT}/plate_{t}_rep{s:02d}_polylines.npy',allow_pickle=True); mr=pl[0]-pl[0][0]; cum=np.r_[0,np.cumsum(np.linalg.norm(np.diff(mr,axis=0),axis=1))]
        ins=[cum[np.argmin(np.linalg.norm(mr-(q[0]-pl[0][0]),axis=1))] for q in pl[1:]]; ln=[np.linalg.norm(np.diff(q,axis=0),axis=1).sum() for q in pl[1:]]
        sims.append(np.histogram(ins,bins,weights=ln)[0])
    ax.plot(mid,np.mean(sims,0),color=COL[t],lw=1.5,ls='--',label=f'{NAME[t]}: simulated 2D')
ax.set_xlabel('Distance of lateral insertion from main-root base (cm)'); ax.set_ylabel('Lateral root length per plant\nin each 0.5 cm of main root (cm)'); ax.legend(fontsize=7,frameon=False,ncol=2)
fig.suptitle('Simulated root systems (all grown to 30 cm total root length) vs the measured panel',fontsize=11)
fig.tight_layout(rect=(0,0,1,0.96))
for e in ('png','pdf'): fig.savefig(f'{FIG}/validation_sim_vs_measured.{e}',dpi=200)
